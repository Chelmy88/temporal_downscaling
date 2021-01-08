#! /usr/local/bin/python3
# Coding: latin1
import os
import pandas as pd
import numpy as np
import multiprocessing as mp

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
# Object to read and store historical data set
class ReadSmet:
  ##### TO BE REPLACED BY NEW FUNCTIONS DEPENDONG ON INPUT DATA PROVIDED ####
  ##### The function should store a panda data frame in self.data with dates as indices
  ##### and variables names as column names
  ##### Should also contain:
  ##### __station = station name
  ##### __header = header from input file to be copied in output files, might be empty
  ##### variables = list of variables in the data frame

  def __init__(self,currDir,station,scenario,period,h,length):

    data_file=currDir+"/output_"+str(length)+"_"+str(h)+"/"+station+"/"+period+"/"+station+"_"+scenario+"_"+period+".smet"
    # Read header and fields
    print("Reading ",data_file)
    with open(data_file, 'r') as f:
      line = next(f)
      fields = []
      pos_read = 1
      self.__header=[]
      while not line.startswith("[DATA]"):
        line = line.strip()
        if line.startswith("fields"):
          fields = line.split("=")[1].strip().split(" ")
        else:
          self.__header.append(line)
        line = next(f)
        pos_read = pos_read + 1
    self.__header.append("[DATA]")
    data = pd.read_csv(data_file, delimiter='\t', dtype=np.float64,
                       skiprows=pos_read, names=fields[1:], index_col=0, parse_dates=[0])

    cols=data.columns

    data=data.astype(np.float64)
    data.replace(-999, np.NaN, inplace=True)
    means=data.groupby(get_season).mean()

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

    self.output=self.output+[scenario,period,h]


def do_work(station,periods,scenarios,h_table,length):
  currDir = os.path.dirname(os.path.realpath(__file__))
  output='seas_mean_10_ds'
  with open(output+"/"+station+".txt","w") as file:
    file.write("\t".join(["TA_DJF","TA_MAM","TA_JJA","TA_SON", \
                          "PSUM_DJF","PSUM_MAM","PSUM_JJA","PSUM_SON", \
                          "RH_DJF","RH_MAM","RH_JJA","RH_SON", \
                          "ISWR_DJF","ISWR_MAM","ISWR_JJA","ISWR_SON", \
                          "VW_DJF","VW_MAM","VW_JJA","VW_SON", \
                          "scenario","period","h"])+"\n")



  #os.makedirs(output, exist_ok=True)
    for h in h_table:
      print(h_table)
      for s in scenarios:
        for p in periods:
          smet=ReadSmet(currDir,station,s,p,h,length)
          file.write("\t".join(smet.output)+'\n')

def main():
  stations=["ABO","BAS","CDF","CHU","COV","DIS","GVE","KLO","LUG","LUZ","OTL","PAY","ROB","SAM","SBE","SCU","SHA","SIO","WFJ","ZER"]
  h_table=['3','5',"7",'9','11','13','15']
  #periods=["1980_2010","2010_2040","2040_2070","2070_2100"]
  #length=30
  periods=["1980_1990","1990_2000","2000_2010","2010_2020","2020_2030","2030_2040","2040_2050","2050_2060","2060_2070","2070_2080","2080_2090","2090_2100"]
  length=10

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

  pool = mp.Pool(10)
  for station in stations:
    #do_work(station,periods,scenarios,h_table,length)
    pool.apply_async(do_work,[station,periods,scenarios,h_table,length])
  pool.close()
  pool.join()

if __name__ == "__main__":
  main()
