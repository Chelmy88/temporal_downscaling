#! /usr/local/bin/python3
# Coding: latin1
import os
import pandas as pd
import numpy as np
import multiprocessing as mp
import codecs

###
# This code allows to extract seasonnal means from raw CH2018 data.
# The seasonnal mean computed with this script should be further used in the
# script plot_assessment.R
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

LOOKUP_SEASON = {
    11: '4Autumn',
    12: '1Winter',
    1: '1Winter',
    2: '1Winter',
    3: '2Spring',
    4: '2Spring',
    5: '2Spring',
    6: '3Summer',
    7: '3Summer',
    8: '3Summer',
    9: '4Autumn',
    10: '4Autumn'}

def get_season(x):
    return(LOOKUP_SEASON[x.month])

# Object to read and store the dataset and compute seasonnal means
class ReadScenario:

    ##### TO BE REPLACED BY NEW FUNCTIONS DEPENDING ON INPUT DATA PROVIDED ####

  __vars={"TA":"tas","PSUM":"pr","ISWR":"rsds","RH":"hurs","VW":"sfcWind"}

  def __init__(self,station,scenario_path,scenario,periods):
    self.__station=station
    self.__scenario_path=scenario_path
    self.__scenario=scenario
    variables=["TA","PSUM","ISWR","RH","VW"]
    # Read data for all variables and store in a dictionary
    self.output_all=""
    to_concatenate=dict()
    print("Reading",station,scenario)
    for var in variables:
      data_tmp=self.__read_scenario_data_file(var)
      if data_tmp is not None:
        to_concatenate[var]=data_tmp
    if len(to_concatenate)==0:
      return

    # Transform the dictionnary into a PD dataframe
    data=self.__concatenate_scenario(to_concatenate)
    cols=data.columns

    data=data.astype(np.float64)
    for p in periods:
      # Get start and end years from period name
      start=int(p[0:4])
      end=int(p[5:9])
      # Cut the data to the desired period
      sub_data=data.loc[(data.index.year >= start) & (data.index.year <= end)]
      # Compute the seasonnal mean
      means=sub_data.groupby(get_season).mean()

      self.output=[]
      if 'TA' in cols:
        self.output = self.output+ [str(x) for x in means['TA'].values]
      else:
        self.output.append(['NaN','NaN','NaN','NaN'])
      if 'PSUM' in cols:
        self.output = self.output+ [str(x) for x in means['PSUM'].values]
      else:
        self.output = self.output+['NaN','NaN','NaN','NaN']
      if 'RH' in cols:
        self.output = self.output+ [str(x) for x in means['RH'].values]
      else:
        self.output = self.output+['NaN','NaN','NaN','NaN']
      if 'ISWR' in cols:
        self.output = self.output+ [str(x) for x in means['ISWR'].values]
      else:
        self.output = self.output+['NaN','NaN','NaN','NaN']
      if 'VW' in cols:
        self.output = self.output+ [str(x) for x in means['VW'].values]
      else:
        self.output = self.output+['NaN','NaN','NaN','NaN']

      self.output=self.output+[scenario,p]
      self.output_all=self.output_all+"\t".join(self.output)+'\n'

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
      #Read data and remove leap days
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



def do_work(station,periods,scenarios,length):
  currDir = os.path.dirname(os.path.realpath(__file__))
  # Path of the raw CH2018 data
  scenario_path=currDir+'/../CH2018/QMstations/'
  output='seas_mean_'+str(length)+'_ds_raw'
  os.makedirs(output, exist_ok=True)
  with open(output+"/"+station+".txt","w") as file:
    file.write("\t".join(["TA_DJF","TA_MAM","TA_JJA","TA_SON", \
                          "PSUM_DJF","PSUM_MAM","PSUM_JJA","PSUM_SON", \
                          "RH_DJF","RH_MAM","RH_JJA","RH_SON", \
                          "ISWR_DJF","ISWR_MAM","ISWR_JJA","ISWR_SON", \
                          "VW_DJF","VW_MAM","VW_JJA","VW_SON", \
                          "scenario","period"])+"\n")
    for s in scenarios:
      ch2018=ReadScenario(station,scenario_path,s,periods)
      if len(ch2018.output_all)>0:
        file.write(ch2018.output_all)

def main():
  # Sample corresponding to a set of 20 MCH stations for 10 years periods
  stations=["ABO","BAS","CDF","CHU","COV","DIS","GVE","KLO","LUG","LUZ","OTL","PAY","ROB","SAM","SBE","SCU","SHA","SIO","WFJ","ZER"]
  periods=["1980_2010","2010_2040","2040_2070","2070_2100"]
  length=30

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

  # Loop over scenarios and do the computation
  ### UNCOMMENT LINE BELOW FOR PARALLEL EXECUTION ###
  # pool = mp.Pool(10)
  for station in stations:
    ### COMMENT LINE FOR PARALLEL EXECUTION ###
    do_work(station,periods,scenarios,length)
    ### UNCOMMENT LINES BELOW FOR PARALLEL EXECUTION ###
    #pool.apply_async(do_work,[station,periods,scenarios,length])
  #pool.close()
  #pool.join()

if __name__ == "__main__":
  main()
