#! /usr/local/bin/python3
# Coding: latin1
import os
import codecs
import pandas as pd
import numpy as np
import multiprocessing as mp


###
# This code computes DOY means from given historical and climate change scenarios
# output variables. Historical data should be provided in the SMET format
# (see https://models.slf.ch/docserver/meteoio/SMET_specifications.pdf), or
# user can redefine the __init__ function of the ReadHistorical class.
# Climate change scenario data should be provided in the same fromat than CH2018,
# or user can redefine the first 3 function of the ReadScenario class.
# Variables in the first section of the 'main' should be modified to correspond
# to the user paths, periods and scenarios.
# Parallel execution can be enabled by uncommenting lines at the end of
# the 'main'.
#
# This code is related to the paper: Climate change scenarios at hourly time-step
# over Switzerland from an enhanced temporal downscaling approach,
# published in the International Journal of Climatology, doi:10.1002/joc.7032
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
# along with this program. If not, see <https://www.gnu.org/licenses/>.
###

class ReadHistorical:

  ##### TO BE REPLACED BY NEW FUNCTIONS DEPENDING ON INPUT DATA PROVIDED ####

  def __init__(self,station,meteo_path):
    self.__station=station
    # Read fields
    with open(meteo_path+station+".smet",'r') as file:
      line = next(file)
      fields=[]
      pos_read=1
      while not line.startswith("[DATA]") :
        line=line.strip()
        if line.startswith("fields"):
          fields=line.split("=")[1].strip().split("\t")
        line=next(file)
        pos_read=pos_read+1
    # Read data
    data_tmp=pd.read_csv(meteo_path+station+".smet", delimiter = '\t| |;|,',dtype=np.float64, \
                         skiprows=pos_read, engine='python', names=fields[1:], index_col=0, parse_dates=[0])
    # Romove leap days, reduce to period, set na value, compute daily mean
    data_tmp = data_tmp[~((data_tmp.index.month == 2) & (data_tmp.index.day == 29))]
    data_tmp = data_tmp.replace(-999,np.NaN)
    # Srote final data (daily mean)
    self.data = data_tmp.resample('D').mean()

  ##### END OF FUNCTIONS TO BE REPLACED #####

  def __cut_data(self,period):
    self.data = self.data[((self.data.index.year >= period[0]) & (self.data.index.year <= period[1]))]

  def __compute_doy_mean(self):
    self.data=self.data.groupby([self.data.index.month, self.data.index.day]).mean().dropna()

  def __print_to_file(self,output_path,period):
    path=output_path+self.__station+"/"
    file_name=self.__station+"_historical"+"_"+str(period[0])+"_"+str(period[1])
    self.data.to_csv(path+file_name,sep="\t",index=False,float_format="%.3f")


  def process_data(self,period,output_path):
    self.__cut_data(period)
    self.__compute_doy_mean()
    self.__print_to_file(output_path,period)

########## INPUT READING FUNCTIONS #####################
class ReadScenario:

    ##### TO BE REPLACED BY NEW FUNCTIONS DEPENDING ON INPUT DATA PROVIDED ####

  __vars={"TA":"tas","PSUM":"pr","ISWR":"rsds","RH":"hurs","VW":"sfcWind"}

  def __init__(self,station,scenario_path,variables,scenario):
    self.__station=station
    self.__scenario_path=scenario_path
    self.__variables_path=variables
    self.__scenario=scenario

    # Read data for all variables and store in a dictionary
    to_concatenate=dict()
    for var in variables:
      data_tmp=self.__read_scenario_data_file(var)
      if data_tmp is not None:
        to_concatenate[var]=data_tmp
    if(len(to_concatenate)==0):
      print("[E] No CH2018 data found for station "+station+" and scenario "+scenario+". Probably files do not exist.")
      self.data=pd.DataFrame()
      return

    # Transform the dictionary into a PD dataframe
    self.data=self.__concatenate_scenario(to_concatenate)

  def __read_scenario_data_file(self,variable):
    var=ReadScenario.__vars[variable]
    name="CH2018_"+var+"_"+self.__scenario+"_QMstations_1981-2099"
    file_name=self.__scenario_path+var+"/"+name+"_csv/"+name+"_"+self.__station+".csv"
    # Count number of lines in header
    if(os.path.isfile(file_name)):
      with codecs.open(file_name, 'r', encoding='utf-8',
      errors='ignore') as file:
        lines = file.readlines()
        fields=[]
        pos_read=0
        while(True):
          line=lines[pos_read].strip()
          pos_read=pos_read+1
          if line.startswith("DATE;VALUE"):
            break
      # Read data and remove leap days
      data=pd.read_csv(file_name, delimiter = '[ \t]*;[ \t]*', dtype=np.float64,na_values='NA', \
                       skiprows=pos_read, engine='python', names=[variable], index_col=0, parse_dates=[0])
      data = data[~((data.index.month == 2) & (data.index.day == 29))]

      return data
    else:
      return None

  def __concatenate_scenario(self,to_concatenate):
    variables=list(to_concatenate.keys())
    data=to_concatenate[variables[0]]
    if len(variables)>1:
      for v in variables[1:]:
        data[v]=to_concatenate[v][v]
    return data

  ##### END OF FUNCTIONS TO BE REPLACED #####

  def __cut_data(self,period):
    return self.data[((self.data.index.year >= period[0]) & (self.data.index.year <= period[1]))]

  def __compute_doy_mean(self,data_tmp):
    return data_tmp.groupby([data_tmp.index.month, data_tmp.index.day]).mean().dropna()

  def __print_to_file(self,data,output_path,period):
    path=output_path+self.__station+"/"+self.__scenario+"/"
    file_name=self.__station+"_"+self.__scenario+"_"+str(period[0])+"_"+str(period[1])
    data.to_csv(path+file_name,sep="\t",index=False,float_format="%.3f")

  def process_data(self,periods,output_path):
    for period in periods:
      data_tmp = self.__cut_data(period)
      data_tmp = self.__compute_doy_mean(data_tmp)
      self.__print_to_file(data_tmp,output_path,period)



def compute_scenario(station,input_CH2018,variables,periods,scenario,output):
  print("\tProcessing:",station,scenario)
  if not os.path.exists(output+station+"/"+scenario):
    os.makedirs(output+station+"/"+scenario)
  scenario_data=ReadScenario(station,input_CH2018,variables,scenario)
  if scenario_data.data.empty:
    return
  scenario_data.process_data(periods,output)


def main():
  currDir = os.path.dirname(os.path.realpath(__file__))

  # Variables to be edited by user

  input_meteo = currDir+'/SMET_MCH_COMPLETE_REDUCED_10/'
  input_CH2018 = currDir+'/../CH2018/QMstations/'

  output = currDir+'/DOY_10/'

  hist_period=[2005,2015]

  periods=[[1980,1990],[1990,2000],[2000,2010],[2010,2020],[2020,2030],[2030,2040],\
           [2040,2050],[2050,2060],[2060,2070],[2070,2080],[2080,2090],[2090,2100]]

  files=[]
  for file in os.listdir(input_meteo):
    if (os.path.isfile(input_meteo+"/"+file) and file.endswith(".smet")):
     files.append(file)
  stations=[x.replace('.smet','') for x in files]

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

  # End of edition section

  for station in stations:
    print("Processing:",station)
    if not os.path.exists(output+station):
      os.makedirs(output+station)
    historical_data = ReadHistorical(station,input_meteo)
    historical_data.process_data(hist_period,output)
    variables=historical_data.data.columns

    ### UNCOMMENT LINE BELOW FOR PARALLEL EXECUTION ###
    # Set number of precessors to use
    #pool = mp.Pool(16)
    for scenario in scenarios:
      ### COMMENT LINE FOR PARALLEL EXECUTION ###
      compute_scenario(station,input_CH2018,variables,periods,scenario,output)
      ### UNCOMMENT LINES BELOW FOR PARALLEL EXECUTION ###
      #pool.apply_async(compute_scenario,[station,input_CH2018,variables,periods,scenario,output])
    #pool.close()
    #pool.join()

if __name__ == "__main__":
    main()
