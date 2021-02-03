setwd("~")
# Set the path above the the directory containing the directories with the
# output of the scripts: Get_seas_error_raw.py and Get_seas_error.py

rm(list=ls())
library(data.table)
library(lubridate)
library(fields)
library(raster)

###
# This script allows to reproduce the seasonnal and natural variability
# analysis. The paths in the script should be adapted.
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

scenarios=c("CLMCOM-CCLM4_ECEARTH_EUR11_RCP45",
            "CLMCOM-CCLM4_ECEARTH_EUR11_RCP85",
            "CLMCOM-CCLM4_HADGEM_EUR11_RCP45",
            "CLMCOM-CCLM4_HADGEM_EUR11_RCP85",
            "CLMCOM-CCLM4_HADGEM_EUR44_RCP85",
            "CLMCOM-CCLM4_MPIESM_EUR11_RCP45",
            "CLMCOM-CCLM4_MPIESM_EUR11_RCP85",
            "CLMCOM-CCLM4_MPIESM_EUR44_RCP45",
            "CLMCOM-CCLM4_MPIESM_EUR44_RCP85",
            "CLMCOM-CCLM5_ECEARTH_EUR44_RCP85",
            "CLMCOM-CCLM5_HADGEM_EUR44_RCP85",
            "CLMCOM-CCLM5_MIROC_EUR44_RCP85",
            "CLMCOM-CCLM5_MPIESM_EUR44_RCP85",
            "DMI-HIRHAM_ECEARTH_EUR11_RCP26",
            "DMI-HIRHAM_ECEARTH_EUR11_RCP45",
            "DMI-HIRHAM_ECEARTH_EUR11_RCP85",
            "DMI-HIRHAM_ECEARTH_EUR44_RCP45",
            "DMI-HIRHAM_ECEARTH_EUR44_RCP85",
            "ICTP-REGCM_HADGEM_EUR44_RCP85",
            "KNMI-RACMO_ECEARTH_EUR44_RCP45",
            "KNMI-RACMO_ECEARTH_EUR44_RCP85",
            "KNMI-RACMO_HADGEM_EUR44_RCP26",
            "KNMI-RACMO_HADGEM_EUR44_RCP45",
            "KNMI-RACMO_HADGEM_EUR44_RCP85",
            "MPICSC-REMO1_MPIESM_EUR11_RCP26",
            "MPICSC-REMO1_MPIESM_EUR11_RCP45",
            "MPICSC-REMO1_MPIESM_EUR11_RCP85",
            "MPICSC-REMO1_MPIESM_EUR44_RCP26",
            "MPICSC-REMO1_MPIESM_EUR44_RCP45",
            "MPICSC-REMO1_MPIESM_EUR44_RCP85",
            "MPICSC-REMO2_MPIESM_EUR11_RCP26",
            "MPICSC-REMO2_MPIESM_EUR11_RCP45",
            "MPICSC-REMO2_MPIESM_EUR11_RCP85",
            "MPICSC-REMO2_MPIESM_EUR44_RCP26",
            "MPICSC-REMO2_MPIESM_EUR44_RCP45",
            "MPICSC-REMO2_MPIESM_EUR44_RCP85",
            "SMHI-RCA_CCCMA_EUR44_RCP45",
            "SMHI-RCA_CCCMA_EUR44_RCP85",
            "SMHI-RCA_CSIRO_EUR44_RCP45",
            "SMHI-RCA_CSIRO_EUR44_RCP85",
            "SMHI-RCA_ECEARTH_EUR11_RCP26",
            "SMHI-RCA_ECEARTH_EUR11_RCP45",
            "SMHI-RCA_ECEARTH_EUR11_RCP85",
            "SMHI-RCA_ECEARTH_EUR44_RCP26",
            "SMHI-RCA_ECEARTH_EUR44_RCP45",
            "SMHI-RCA_ECEARTH_EUR44_RCP85",
            "SMHI-RCA_GFDL_EUR44_RCP45",
            "SMHI-RCA_GFDL_EUR44_RCP85",
            "SMHI-RCA_HADGEM_EUR11_RCP45",
            "SMHI-RCA_HADGEM_EUR11_RCP85",
            "SMHI-RCA_HADGEM_EUR44_RCP26",
            "SMHI-RCA_HADGEM_EUR44_RCP45",
            "SMHI-RCA_HADGEM_EUR44_RCP85",
            "SMHI-RCA_IPSL_EUR11_RCP45",
            "SMHI-RCA_IPSL_EUR11_RCP85",
            "SMHI-RCA_IPSL_EUR44_RCP45",
            "SMHI-RCA_IPSL_EUR44_RCP85",
            "SMHI-RCA_MIROC_EUR44_RCP26",
            "SMHI-RCA_MIROC_EUR44_RCP45",
            "SMHI-RCA_MIROC_EUR44_RCP85",
            "SMHI-RCA_MPIESM_EUR11_RCP45",
            "SMHI-RCA_MPIESM_EUR11_RCP85",
            "SMHI-RCA_MPIESM_EUR44_RCP26",
            "SMHI-RCA_MPIESM_EUR44_RCP45",
            "SMHI-RCA_MPIESM_EUR44_RCP85",
            "SMHI-RCA_NORESM_EUR44_RCP26",
            "SMHI-RCA_NORESM_EUR44_RCP45",
            "SMHI-RCA_NORESM_EUR44_RCP85")


periods_30=c("1980_2010","2010_2040","2040_2070","2070_2100")
data_30=list()
stations_30=c("ABO","BAS","CDF","CHU","COV","DIS","GVE","KLO","LUG","LUZ","OTL","PAY","ROB","SAM","SBE","SCU","SHA","SIO","WFJ","ZER")
h_table=c("3","5","7","9","11","13","15")

for (station in stations_30)
{
  print(paste("Reading",station))
  seas_values_all=as.data.frame(fread(paste0("seas_mean_30_ds/",station,".txt")))
  seas_values_all_raw=as.data.frame(fread(paste0("seas_mean_30_ds_raw/",station,".txt")))
  for(h in h_table)
  {
    seas_values_all_h=seas_values_all[which(seas_values_all$h==as.numeric(h)),]
    for (scenario in scenarios)
    {
      seas_values_all_s=seas_values_all_h[which(seas_values_all_h$scenario==scenario),]
      seas_values_all_s_raw=seas_values_all_raw[which(seas_values_all_raw$scenario==scenario),]
      seas_values_all_s
      for(p in periods_30){
        seas_values_all_p=seas_values_all_s[which(seas_values_all_s$period==p),]
        seas_values_all_p_raw=seas_values_all_s_raw[which(seas_values_all_s_raw$period==p),]
        seas_values_all_p
        seas_values_all_p_raw
        seas=data.frame(TA=unlist(seas_values_all_p[1:4]),
                        PSUM=unlist(seas_values_all_p[5:8]),
                        RH=unlist(seas_values_all_p[9:12]),
                        ISWR=unlist(seas_values_all_p[13:16]),
                        VW=unlist(seas_values_all_p[17:20]))
        rownames(seas)=c("DJF","MAM","JJA","SON")
        data_30[[station]][[h]][[scenario]][[p]][["downscaled"]]=seas
        seas=data.frame(TA=unlist(seas_values_all_p_raw[1:4]),
                        PSUM=unlist(seas_values_all_p_raw[5:8]),
                        RH=unlist(seas_values_all_p_raw[9:12]),
                        ISWR=unlist(seas_values_all_p_raw[13:16]),
                        VW=unlist(seas_values_all_p_raw[17:20]))
        rownames(seas)=c("DJF","MAM","JJA","SON")
        data_30[[station]][[h]][[scenario]][[p]][["RAW"]]=seas
      }
    }
  }
}

all_raw_error_30=list()

vars=c("TA","PSUM","ISWR","RH","VW")
cols_seas=list("DJF"="blue","MAM"="green","JJA"="orange","SON"="brown")
rcp_pch=list("RCP26"=1,"RCP45"=2,"RCP85"=3)

for (station in stations_30)
{
  all_raw_error_30[[station]]=list()
  for (h in h_table)
  {
    all_raw_error_30[[station]][[h]]=list()

    for(var in vars)
    {
      all_raw_error_30[[station]][[h]][[var]]=c(NA)
    }
  }
}


for (station in stations_30)
{
  for (h in h_table)
  {
    for (scenario in scenarios)
    {
      vars=colnames(data_30[[station]][[h]][[scenario]][[1]]$downscaled)
      for(p in periods_30)
      {
        for(var in vars)
        {

          raw=data_30[[station]][[h]][[scenario]][[p]]$RAW[[var]]
          downscaled=data_30[[station]][[h]][[scenario]][[p]]$downscaled[[var]]
          if(var=="TA")
          {
            downscaled=downscaled-273.15
          }else if(var=="RH")
          {
            downscaled=downscaled*100
          }
          if(var=="PSUM")
          {
            downscaled=downscaled*24
          }
          raw_error=raw-downscaled
          all_raw_error_30[[station]][[h]][[var]]=c(all_raw_error_30[[station]][[h]][[var]],raw_error)
        }
      }
    }
  }
}


#############################################################################
#############################################################################
### COMPUTPE NATURAL VARIABILITY FOR ALL SCENARION, STATIONS, AND PERIODS ###
#############################################################################
#############################################################################

# Computes the MSE for 1 sample with the approximation obtained from the 10
# other samples, as in Bosshard

all_mse=list()
mse=list()

satation="PAY"
raw_cc_path="<PATH TO RAW CH2018>"
var_conv=list("RH"="hurs","PSUM"="pr","ISWR"="rsds","VW"="sfcWind","TA"="tas")

var="hurs"
scenario="DMI-HIRHAM_ECEARTH_EUR11_RCP45"
cc=as.data.frame(fread(paste0(raw_cc_path,"/CH2018_",var,"_",scenario,"_QMstations_1981-2099_",station,".csv"),skip=17))
leap=which(strftime(cc$DATE,"%d%m")!="2902")
cc=cc[leap,]
p1=which(strftime(cc$DATE,"%Y")>=1980 & strftime(cc$DATE,"%Y")<2010 )
p2=which(strftime(cc$DATE,"%Y")>=2010 & strftime(cc$DATE,"%Y")<2040 )
p3=which(strftime(cc$DATE,"%Y")>=2040 & strftime(cc$DATE,"%Y")<2070 )
p4=which(strftime(cc$DATE,"%Y")>=2070 & strftime(cc$DATE,"%Y")<2100 )

h_table=c("3","5","7","9","11","13","15")
stations=c("ABO","BAS","CDF","CHU","COV","DIS","GVE","KLO","LUG","LUZ","OTL","PAY","ROB","SAM","SBE","SCU","SHA","SIO","WFJ","ZER")
for (station in stations)
{
  print(station)
  mse[[station]]=list()
  all_mse[[station]]=list()

  for (var in c("RH","TA","PSUM","ISWR","VW"))
  {
    loc_var=var_conv[[var]]
    print(var)
    mse[[station]][[var]]=list()
    all_mse[[station]][[var]]=list()

    for(h in h_table){
      all_mse[[station]][[var]][[h]]=c()
    }

    for(scenario in scenarios){
      mse[[station]][[var]][[scenario]]=list()

      file=paste0(raw_cc_path,"/CH2018_",loc_var,"_",scenario,"_QMstations_1981-2099_",station,".csv")
      if(file.exists(file)){
        cc=as.data.frame(fread(file,skip=17))
      }
      else{
        next
      }
      cc=cc[leap,]

      cc1_tmp=cc[p1,]
      cc2_tmp=cc[p2,]
      cc3_tmp=cc[p3,]
      cc4_tmp=cc[p4,]

      cc1=matrix(0,nrow=365,ncol=10)
      cc2=matrix(0,nrow=365,ncol=10)
      cc3=matrix(0,nrow=365,ncol=10)
      cc4=matrix(0,nrow=365,ncol=10)

      for (i in c (1:10))
      {
        tmp=cbind(cc1_tmp$VALUE[c(1:365)+(3*i-3)*365],cc1_tmp$VALUE[c(1:365)+(3*i-2)*365],cc1_tmp$VALUE[c(1:365)+(3*i-1)*365])
        cc1[,i]=rowMeans(tmp, na.rm=TRUE)
        tmp=cbind(cc2_tmp$VALUE[c(1:365)+(3*i-3)*365],cc2_tmp$VALUE[c(1:365)+(3*i-2)*365],cc2_tmp$VALUE[c(1:365)+(3*i-1)*365])
        cc2[,i]=rowMeans(tmp, na.rm=TRUE)
        tmp=cbind(cc3_tmp$VALUE[c(1:365)+(3*i-3)*365],cc3_tmp$VALUE[c(1:365)+(3*i-2)*365],cc3_tmp$VALUE[c(1:365)+(3*i-1)*365])
        cc3[,i]=rowMeans(tmp, na.rm=TRUE)
        tmp=cbind(cc4_tmp$VALUE[c(1:365)+(3*i-3)*365],cc4_tmp$VALUE[c(1:365)+(3*i-2)*365],cc4_tmp$VALUE[c(1:365)+(3*i-1)*365])
        cc4[,i]=rowMeans(tmp, na.rm=TRUE)
      }

      for(h in h_table){
        mse[[station]][[var]][[scenario]][["cc1"]][[h]]=0
        mse[[station]][[var]][[scenario]][["cc2"]][[h]]=0
        mse[[station]][[var]][[scenario]][["cc3"]][[h]]=0
        mse[[station]][[var]][[scenario]][["cc4"]][[h]]=0

        for (i in c (1:10))
        {
          sample=rowMeans(cc1[,-i])
          approx=getHarmonic(sample,strtoi(h))
          mse[[station]][[var]][[scenario]][["cc1"]][[h]]=
            mse[[station]][[var]][[scenario]][["cc1"]][[h]]+sqrt(sum((approx-cc1[,i])*(approx-cc1[,i]))/365)
          sample=rowMeans(cc2[,-i])
          approx=getHarmonic(sample,strtoi(h))
          mse[[station]][[var]][[scenario]][["cc2"]][[h]]=
            mse[[station]][[var]][[scenario]][["cc2"]][[h]]+sqrt(sum((approx-cc2[,i])*(approx-cc2[,i]))/365)
          sample=rowMeans(cc3[,-i])
          approx=getHarmonic(sample,strtoi(h))
          mse[[station]][[var]][[scenario]][["cc3"]][[h]]=
            mse[[station]][[var]][[scenario]][["cc3"]][[h]]+sqrt(sum((approx-cc3[,i])*(approx-cc3[,i]))/365)
          sample=rowMeans(cc4[,-i])
          approx=getHarmonic(sample,strtoi(h))
          mse[[station]][[var]][[scenario]][["cc4"]][[h]]=
            mse[[station]][[var]][[scenario]][["cc4"]][[h]]+sqrt(sum((approx-cc4[,i])*(approx-cc4[,i]))/365)

        }
        mse[[station]][[var]][[scenario]][["cc1"]][[h]]=(mse[[station]][[var]][[scenario]][["cc1"]][[h]])/10
        mse[[station]][[var]][[scenario]][["cc2"]][[h]]=(mse[[station]][[var]][[scenario]][["cc2"]][[h]])/10
        mse[[station]][[var]][[scenario]][["cc3"]][[h]]=(mse[[station]][[var]][[scenario]][["cc3"]][[h]])/10
        mse[[station]][[var]][[scenario]][["cc4"]][[h]]=(mse[[station]][[var]][[scenario]][["cc4"]][[h]])/10
        all_mse[[station]][[var]][[h]]=c(all_mse[[station]][[var]][[h]],
                                         mse[[station]][[var]][[scenario]][["cc1"]][[h]],
                                         mse[[station]][[var]][[scenario]][["cc2"]][[h]],
                                         mse[[station]][[var]][[scenario]][["cc3"]][[h]],
                                         mse[[station]][[var]][[scenario]][["cc4"]][[h]])
      }
    }
  }
}


# Substract mean for each scenario
diff_mse=list()
for (station in names(all_mse))
{
  for (var in names(all_mse[[station]])){
    means=rowMeans(as.data.frame(all_mse[[station]][[var]]))
    for (h in h_table)
    {
      diff_mse[[station]][[var]][[h]]=all_mse[[station]][[var]][[h]]-means
    }
  }
}


##############################
## DO THE ASSESSMENT PLOTS ###
##############################

y_units=list()
y_units[["TA"]]=expression("Centered mean RMSE ("*degree*"C)")
y_units[["PSUM"]]="Centered mean RMSE (mm/h)"
y_units[["RH"]]="Centered mean RMSE (%)"
y_units[["ISWR"]]=expression("Centered mean RMSE (W/m"^2*")")
y_units[["VW"]]="Centered mean RMSE (m/s)"

x_units=list()
x_units[["TA"]]=expression("Mean abs seas. raw error ("*degree*"C)")
x_units[["PSUM"]]="Mean abs seas. raw error  (mm/h)"
x_units[["RH"]]="Mean abs seas. raw error  (%)"
x_units[["ISWR"]]=expression("Mean abs seas. raw error  (W/m"^2*")")
x_units[["VW"]]="Mean abs seas. raw error  (m/s)"


pdf(paste0("plots/ALL_ASSES_all_stns_30.pdf"),width=10,height=6)
par(mfrow=c(5,1),mar=c(3,3,2,1),mgp=c(1.7,0.5,0),oma=c(0,0,2,0))
layout(matrix(c(1,1,2,2,3,3,4,4,0,5,5,0),byrow = TRUE,ncol = 4))
for (station in stations_30)
{
  ord=match(c("TA","PSUM","RH","ISWR","VW"),names(all_mse[[station]]))
  variables=names(all_mse[[station]])[ord]

  for (var in variables)
  {
    mean_raw=c()
    mean_mse=c()
    var_raw=c()
    var_mse=c()
    for(h in h_table){
      mean_raw=c(mean_raw, mean(abs(all_raw_error_30[[station]][[h]][[var]]),na.rm=T))
      mean_mse=c(mean_mse,mean(diff_mse[[station]][[var]][[h]],na.rm=T))
      var_raw=c(var_raw, sqrt(var(abs(all_raw_error_30[[station]][[h]][[var]]),na.rm=T)))
      var_mse=c(var_mse,sqrt(var(diff_mse[[station]][[var]][[h]],na.rm=T)))
    }
    x_lim=c(min(mean_raw-var_raw),max(mean_raw+var_raw))
    y_lim=c(min(mean_mse-var_mse),max(mean_mse+var_mse))
    plot(mean_raw,mean_mse,pch=1,ylim=y_lim,xlim=x_lim,main=var,ylab=y_units[[var]],xlab=x_units[[var]],type="n")
    arrows(x0=mean_raw, y0=mean_mse-var_mse, x1=mean_raw, y1=mean_mse+var_mse, angle=90,code=3,
           length=0.04,lwd=0.4,col="grey")
    arrows(x0=mean_raw-var_raw, y0=mean_mse, x1=mean_raw+var_raw,y1=mean_mse, angle=90, code=3,
           length=0.04,lwd=0.4,col="grey")
    text(mean_raw,mean_mse,h_table,pos=3)
    points(mean_raw,mean_mse,pch=1)
  }
  mtext(station,side=3,outer=T,font=2)
}
dev.off()
