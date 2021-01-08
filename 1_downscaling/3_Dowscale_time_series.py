#! /usr/local/bin/python3
# Coding: latin1
import os
import codecs
import pandas as pd
import numpy as np
import multiprocessing as mp
import calendar
import copy

###
# This code computes the actual downscaled time series baes on the output of
# the script '2_Mean_year_computation.py'. Historical data should be provided
# in the SMET format
# (see https://models.slf.ch/docserver/meteoio/SMET_specifications.pdf), or
# user can redefine the __init__ function of the ReadHistorical class.
# User can use DOY data from other source than the script
# '2_Mean_year_computation.py' by updating te ReadDOY class.
# The function set_units in the class DownscaledData performs some unit conversion.
# It can be modified to fit user needs.
# Variables in the first section of the 'main' should be modified to correspond
# to the user paths, periods and scenarios.
# Climate change scenario data should be provided in the same format than CH2018,
# or user can redefine the first 3 function of the ReadScenario class.
# Parallel execution can be enabled by uncommenting lines at the end of
# the 'main'.
#
# This code is related to the paper: Climate change scenarios at hourly time-step
# over Switzerland from an enhanced temporal downscaling approach,
# currently under review in the International Journal of Climatology.
#
# Copyright (C) 2021, Adrien Michel, EPFL
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
###

########## HARMONIC FITTING FUNCTIONS ###################

PI = 3.14159
def getA(data, q, length):
  a = []
  for j in range(1, q + 1):
    s = 0
    w = j / length
    for t in range(length):
      s = s + 2 / length * data[t] * np.cos(2 * PI * w * (t + 1))
    a.append(s)
  return(a)

def getB(data, q, length):
  b = []
  for j in range(1, q + 1):
    s = 0
    w = j / length
    for t in range(length):
      s = s + 2 / length * data[t] * np.sin(2 * PI * w * (t + 1))
    b.append(s)
  return(b)

def getX(t, a, b, a0, length, q):
  x = np.full(len(t), a0)
  for j in range(q):
    w = (j + 1) / length
    x = x + a[j] * np.cos(2 * PI * w * t) + b[j] * np.sin(2 * PI * w * t)
  return(x)

def getHarmonic(var, q):
  length = len(var)
  a = getA(var, q, length)
  b = getB(var, q, length)
  a0 = np.mean(var)
  return(getX(np.arange(1, 366) * length / 365, a, b, a0, length, q))

########## END OF HARMONIC FITTING FUNCTIONS ###################

# Object to read and store historical data set
class ReadHistorical:
  ##### TO BE REPLACED BY NEW FUNCTIONS DEPENDONG ON INPUT DATA PROVIDED ####
  ##### The function should store a panda data frame in self.data with dates as indices
  ##### and variables names as column names
  ##### Should also contain:
  ##### __station = station name
  ##### __header = header from input file to be copied in output files, might be empty
  ##### variables = list of variables in the data frame

  def __init__(self, station, meteo_path):

    self.__station = station

    # Read header and fields
    with open(meteo_path + station + ".smet", 'r') as file:
      line = next(file)
      fields = []
      pos_read = 1
      self.__header=[]
      while not line.startswith("[DATA]"):
        line = line.strip()
        if line.startswith("fields"):
          fields = line.split("=")[1].strip().split("\t")
        else:
          self.__header.append(line)
        line = next(file)
        pos_read = pos_read + 1
    self.__header.append("[DATA]")

    # Read data
    self.data = pd.read_csv(meteo_path + station + ".smet", delimiter='\t| |;|,', dtype=np.float64,
                            skiprows=pos_read, engine='python', names=fields[1:], index_col=0, parse_dates=[0])

    # Remove leap days, set na value, compute daily mean
    self.data.replace(-999, np.NaN, inplace=True)
    self.data = self.data[~(
        (self.data.index.month == 2) & (self.data.index.day == 29))]

    # Store variables to be compared to variables available in climate change scenarios
    self.variables = self.data.columns
  ##### END OF FUNCTIONS TO REPLACE ####

  def get_header(self):
    return(self.__header)

# Object containing the downscaled time series
class DownscaledData:

  def __init__(self,historical_data, scenario_doy, period,historical_period,dates,station,scenario):
    # Store vriables for file print out
    self.__station=station
    self.__scenario=scenario
    self.__header=historical_data.get_header()
    self.__historical_period=historical_period
    self.__period=period

    # Copy data from historical time series taht are both in historical data and scenarios data
    self.variables_to_downscale=np.intersect1d(historical_data.variables,scenario_doy.variables)
    self.data = historical_data.data[self.variables_to_downscale].copy()

    # Change timestep to the future period
    name=str(period[0])+"_"+str(period[1])
    index=dates[name]
    self.data.set_index(index,inplace=True)

    # Store the range of years (used later to add leap days)
    self.years=range(period[0],period[1]+1)

  ##### TO BE REPLACED BY NEW FUNCTIONS DEPENDONG ON INPUT DATA PROVIDED ####
  # Change units of some variables
  def set_units(self):
    for var in self.variables_to_downscale:
      if var == "TA":
        self.data[var] = self.data[var] + 273.15
      if var == "RH":
        self.data[var] = self.data[var]/100
  ##### END OF FUNCTIONS TO REPLACE ####


  def complete_leap(self):
    for year in range(self.__period[0],self.__period[1]+1):
      if calendar.isleap(year):
        index=self.data.index
        leap_day=self.data[(index.year==year) & (index.month==2) & (index.day==28)].copy()
        leap_day.rename( index=lambda ind: ind + pd.DateOffset(days=1), inplace=True)
        self.data=self.data.append(leap_day).sort_index()

  ##### TO BE REPLACED BY NEW FUNCTIONS DEPENDONG ON OUTPUT DATA DESIRED ####
  def print_file(self,output):
    os.makedirs(output, exist_ok=True)
    os.makedirs(output+"/"+self.__station, exist_ok=True)
    os.makedirs(output+"/"+self.__station+"/"+str(self.__period[0])+"_"+str(self.__period[1]), exist_ok=True)
    file_name=output+"/"+self.__station+"/"+str(self.__period[0])+"_"+str(self.__period[1])+"/"+self.__station \
              +"_"+self.__scenario+"_"+str(self.__period[0])+"_"+str(self.__period[1])+".smet"
    header_loc=copy.copy(self.__header)
    header_loc.insert(-1,"info = Generated with XXX. Data for scenario "+self.__scenario+" and period " \
                        +str(self.__period[0])+"_"+str(self.__period[1]) + " based on historical period " \
                        +str(self.__historical_period[0])+"_"+str(self.__historical_period[1]))

    header_loc.insert(-1,"fields = timestamp "+" ".join(self.data.columns))

    with open(file_name,"w") as f:
      f.write("\n".join(header_loc)+"\n")
    self.data.to_csv(file_name,header=False, index=True,mode="a", sep='\t',float_format="%.3f",na_rep="-999",date_format="%Y-%m-%dT%H:%M")
  ##### END OF FUNCTIONS TO REPLACE ####

# Object to read and store DOY means for historical data and climate change scenarios
class ReadDOY:
  ##### TO BE REPLACED BY NEW FUNCTIONS DEPENDONG ON INPUT DATA PROVIDED ####
  ##### The function should store a panda data frame in self.data containg DOY means
  ##### with variables names as column names, Should also strore list of variables in
  ##### self.variables
  def __init__(self, doy_path):
    # Read fields
    self.data = pd.read_csv(doy_path, delimiter='\t', dtype=np.float64)
    self.variables = self.data.columns
  ##### END OF FUNCTIONS TO REPLACE ####

  def __str__(self):
    return(self.data.to_string())


# Function retunring timestamps for periods where downscaling is applied
def get_dates(historical_data,periods,historical_period):
  dates=dict()
  for i in range(len(periods)):
    name=str(periods[i][0])+"_"+str(periods[i][1])
    time_delta=periods[i][0]-historical_period[0]
    dates[name]=historical_data.data.index+pd.DateOffset(years=time_delta)
  return(dates)


# Do the actual dowscaling for the provided sttion, scenario and periods
def do_downscale(station,historical_data,historical_doy,scenario,periods,hist_period,dates,h,input_DOY,output):
  for i in range(len(periods)):
    print(station,h,scenario,str(periods[i][0])+"_"+str(periods[i][1]))
    # Read DOY means
    path = input_DOY + station + "/" + scenario + "/" + station + "_" + \
       scenario + "_" + str(periods[i][0]) + "_" + str(periods[i][1])
    scenario_doy = ReadDOY(path)
    # Create object to store dowscaled data
    dowscaled_data=DownscaledData(historical_data, scenario_doy, periods[i],hist_period,dates,station, scenario)

    for var in dowscaled_data.variables_to_downscale:
     # Do the downscaling
     if var == "TA":
       deltas = getHarmonic(
           scenario_doy.data[var], h) - getHarmonic(historical_doy.data[var], h)
       deltas=np.repeat(deltas,24)
       deltas=np.tile(deltas,periods[i][1]-periods[i][0]+1)
       dowscaled_data.data[var] = dowscaled_data.data[var] + deltas
     else:
       deltas = getHarmonic(
           scenario_doy.data[var], h)/getHarmonic(historical_doy.data[var], h)
       deltas=np.repeat(deltas,24)
       deltas=np.tile(deltas,periods[i][1]-periods[i][0]+1)
       dowscaled_data.data[var] = dowscaled_data.data[var] * deltas
    # Put variables in the right units
    dowscaled_data.set_units()
    # Complete leam years
    dowscaled_data.complete_leap()
    # Print output
    dowscaled_data.print_file(output)


#### MAIN ####

def main():
  currDir = os.path.dirname(os.path.realpath(__file__))

  ### VARIABLES TO BE SET ###
  input_meteo = currDir + '/SMET_MCH_COMPLETE_REDUCED_30/'
  input_DOY = currDir + '/DOY_30/'
  h_table = [3,5,7,9,11,13,15]

  hist_period=[1985, 2015]

  periods = [[1980,2010],[2010,2040],[2040,2070],[2070,2100]]

  files=[]
  for file in os.listdir(input_meteo):
    if (os.path.isfile(input_meteo+"/"+file) and file.endswith(".smet")):
     files.append(file)

  stations=[x.replace('.smet','') for x in files]

  for h in h_table:
    output = currDir + '/output_30_'+str(h)+'/'

    scenarios=["CLMCOM-CCLM4_ECEARTH_EUR11_RCP45", \
          "CLMCOM-CCLM4_ECEARTH_EUR11_RCP85", \
          "CLMCOM-CCLM4_HADGEM_EUR11_RCP45", \
          "CLMCOM-CCLM4_HADGEM_EUR11_RCP85", \
          "CLMCOM-CCLM4_HADGEM_EUR44_RCP85", \
          "CLMCOM-CCLM4_MPIESM_EUR11_RCP45", \
          "CLMCOM-CCLM4_MPIESM_EUR11_RCP85", \
          "CLMCOM-CCLM4_MPIESM_EUR44_RCP45", \
          "CLMCOM-CCLM4_MPIESM_EUR44_RCP85", \
          "CLMCOM-CCLM5_ECEARTH_EUR44_RCP85", \
          "CLMCOM-CCLM5_HADGEM_EUR44_RCP85", \
          "CLMCOM-CCLM5_MIROC_EUR44_RCP85", \
          "CLMCOM-CCLM5_MPIESM_EUR44_RCP85", \
          "DMI-HIRHAM_ECEARTH_EUR11_RCP26", \
          "DMI-HIRHAM_ECEARTH_EUR11_RCP45", \
          "DMI-HIRHAM_ECEARTH_EUR11_RCP85", \
          "DMI-HIRHAM_ECEARTH_EUR44_RCP45", \
          "DMI-HIRHAM_ECEARTH_EUR44_RCP85", \
          "ICTP-REGCM_HADGEM_EUR44_RCP85", \
          "KNMI-RACMO_ECEARTH_EUR44_RCP45", \
          "KNMI-RACMO_ECEARTH_EUR44_RCP85", \
          "KNMI-RACMO_HADGEM_EUR44_RCP26", \
          "KNMI-RACMO_HADGEM_EUR44_RCP45", \
          "KNMI-RACMO_HADGEM_EUR44_RCP85", \
          "MPICSC-REMO1_MPIESM_EUR11_RCP26", \
          "MPICSC-REMO1_MPIESM_EUR11_RCP45", \
          "MPICSC-REMO1_MPIESM_EUR11_RCP85", \
          "MPICSC-REMO1_MPIESM_EUR44_RCP26", \
          "MPICSC-REMO1_MPIESM_EUR44_RCP45", \
          "MPICSC-REMO1_MPIESM_EUR44_RCP85", \
          "MPICSC-REMO2_MPIESM_EUR11_RCP26", \
          "MPICSC-REMO2_MPIESM_EUR11_RCP45", \
          "MPICSC-REMO2_MPIESM_EUR11_RCP85", \
          "MPICSC-REMO2_MPIESM_EUR44_RCP26", \
          "MPICSC-REMO2_MPIESM_EUR44_RCP45", \
          "MPICSC-REMO2_MPIESM_EUR44_RCP85", \
          "SMHI-RCA_CCCMA_EUR44_RCP45", \
          "SMHI-RCA_CCCMA_EUR44_RCP85", \
          "SMHI-RCA_CSIRO_EUR44_RCP45", \
          "SMHI-RCA_CSIRO_EUR44_RCP85", \
          "SMHI-RCA_ECEARTH_EUR11_RCP26", \
          "SMHI-RCA_ECEARTH_EUR11_RCP45", \
          "SMHI-RCA_ECEARTH_EUR11_RCP85", \
          "SMHI-RCA_ECEARTH_EUR44_RCP26", \
          "SMHI-RCA_ECEARTH_EUR44_RCP45", \
          "SMHI-RCA_ECEARTH_EUR44_RCP85", \
          "SMHI-RCA_GFDL_EUR44_RCP45", \
          "SMHI-RCA_GFDL_EUR44_RCP85", \
          "SMHI-RCA_HADGEM_EUR11_RCP45", \
          "SMHI-RCA_HADGEM_EUR11_RCP85", \
          "SMHI-RCA_HADGEM_EUR44_RCP26", \
          "SMHI-RCA_HADGEM_EUR44_RCP45", \
          "SMHI-RCA_HADGEM_EUR44_RCP85", \
          "SMHI-RCA_IPSL_EUR11_RCP45", \
          "SMHI-RCA_IPSL_EUR11_RCP85", \
          "SMHI-RCA_IPSL_EUR44_RCP45", \
          "SMHI-RCA_IPSL_EUR44_RCP85", \
          "SMHI-RCA_MIROC_EUR44_RCP26", \
          "SMHI-RCA_MIROC_EUR44_RCP45", \
          "SMHI-RCA_MIROC_EUR44_RCP85", \
          "SMHI-RCA_MPIESM_EUR11_RCP45", \
          "SMHI-RCA_MPIESM_EUR11_RCP85", \
          "SMHI-RCA_MPIESM_EUR44_RCP26", \
          "SMHI-RCA_MPIESM_EUR44_RCP45", \
          "SMHI-RCA_MPIESM_EUR44_RCP85", \
          "SMHI-RCA_NORESM_EUR44_RCP26", \
          "SMHI-RCA_NORESM_EUR44_RCP45", \
          "SMHI-RCA_NORESM_EUR44_RCP85"]

### END OF VARIABLES TO BE SET ###

    for station in stations:
      # Read historical data and DOY
      historical_data = ReadHistorical(station, input_meteo)
      historical_doy = ReadDOY(input_DOY + station + "/" + station +
                             "_historical_" + str(hist_period[0]) + "_" + str(hist_period[1]))
      # Correct precipitations
      if("PSUM" in historical_doy.variables):
        historical_doy.data["PSUM"] = historical_doy.data["PSUM"] * 24

      # Compute only once all the new timestamps fo all the futur periods, to be assigned to generated df
      dates=get_dates(historical_data,periods,hist_period)

      # Loop over scenarios and do downscaling
      ### UNCOMMENT LINE BELOW FOR PARALLEL EXECUTION ###
      # pool = mp.Pool(10)
      for scenario in scenarios:
        ### COMMENT LINE FOR PARALLEL EXECUTION ###
        do_downscale(station,historical_data,historical_doy,scenario,periods,hist_period,dates,h,input_DOY,output)
        ### UNCOMMENT LINES BELOW FOR PARALLEL EXECUTION ###
        #pool.apply_async(do_downscale,[station,historical_data,historical_doy,scenario,periods,hist_period,dates,h,input_DOY,output])
      #pool.close()
      #pool.join()

if __name__ == "__main__":
  main()
