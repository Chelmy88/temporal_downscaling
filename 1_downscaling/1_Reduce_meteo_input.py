#! /usr/local/bin/python3
import os
import glob
import datetime
import sys
import math
import multiprocessing as mp

import datetime
import re
import warnings

###
# This code allows to cut historical meteorological time series between the
# desires dates. It will also remove variables where no data count exceeds the
# threshold of 1000 missing data. Data should be provided in the SMET format
# (see https://models.slf.ch/docserver/meteoio/SMET_specifications.pdf)
# Variables in the first section of the code should be modified.
#
# This code is related to the paper: Climate change scenarios at hourly time-step
# over Switzerland from an enhanced temporal downscaling approach,
# currently under review in the International Journal of Climatology.
#
# Copyright (C) 2020, Adrien Michel, EPFL
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.
###


#### DEFINE GLOBAL VARIABLES ###
#### Variables to be edited by user ###
currDir = os.path.dirname(os.path.realpath(__file__))
input = currDir + '/<input directory>'
output = currDir + '/<output directory>/'
loaded_data = {}

os.makedirs(output, exist_ok=True)

date_format = '<date fromat string of input dates>, e.g. %Y-%m-%dT%H:%M'
starting_date = '<start date in format %Y-%m-%dT%H:%M> e.g. 2005-01-01T00:00'
ending_date = '<end date in format %Y-%m-%dT%H:%M> e.g. 2005-01-01T00:00'


########## INPUT READING FUNCTIONS #####################

def loadData(filename):
  global date_format
  with open(filename, 'r') as file:
    pos_read = 0
    header = []
    fields = []
    fields_to_keep = []

    isheader = True
    lines = file.readlines()
    while(isheader):
      line = lines[pos_read].strip()
      if line.find("[DATA]") != 0:
        pos_read = pos_read + 1
        if line.find("fields") == 0:
          fields = str(line.split("=")[1]).split()
          for v in ["TA", "PSUM", "VW", "ISWR", "ILWR", "RH"]:
            if v in fields:
              fields_to_keep.append(fields.index(v))
          continue
        header.append(line)
        continue
      else:
        isheader = False
        header.append("[DATA]")
        pos_read = pos_read + 1

    if len(fields_to_keep) == 0:
      return None, None, None

    dates = []
    vars = dict()
    for f in [fields[x] for x in fields_to_keep]:
      vars[f] = []

    while(pos_read < len(lines)):
      line = lines[pos_read]
      line = line.strip().split()
      dates.append(datetime.datetime.strptime(line[0], date_format))
      for x in fields_to_keep:
        vars[fields[x]].append(line[x])
      pos_read += 1

    return header, vars, dates


def cutData(header, vars, dates):
  global starting_date
  global ending_date
  try:
    start_date = dates.index(
        datetime.datetime.strptime(starting_date, '%Y-%m-%dT%H:%M'))
  except:
    return None, None
  try:
    end_date = dates.index(
        datetime.datetime.strptime(ending_date, '%Y-%m-%dT%H:%M'))
  except:
    return None, None

  dates = dates[start_date:end_date]
  for var in vars:
    vars[var] = vars[var][start_date:end_date]

  vars_to_remove = []
  for var in vars:
    i = 0
    count_nodata = vars[var].count('-999')
    if(count_nodata >= 1000):
      print("\t[W]", var, "removed because of to much inner no data")
      vars_to_remove.append(var)

  for var in vars_to_remove:
    del vars[var]
  if(len(vars) == 0):
    return None, None
  writer = []
  for i in range(len(dates)):
    line = datetime.datetime.strftime(dates[i], '%Y-%m-%dT%H:%M')
    for var in vars:
      line = line + "\t" + vars[var][i]
    writer.append(line)

  fields = "fields = timestamp\t" + "\t".join(vars.keys())
  header.insert(-2, fields)

  return header, writer


def writeData(file, header, writer):
  with open(file, 'w') as write_file:
    for line in header:
      write_file.write(line + "\n")
    for line in writer:
      write_file.write(line + "\n")

##########################


#### MAIN ####

def main():
  global input
  global loaded_data
  print("[I] Checking input")

  # Gather available data files and legend files
  files = []
  for file in os.listdir(input):
    if (os.path.isfile(input + "/" + file) and file.endswith(".smet")):
      files.append(file)

  # Load data
  print("[I] Loading Data")
  for file in files:
    print("[I] Reading file " + file)
    [header, vars, dates] = loadData(input + "/" + file)
    if header is None and vars is None and dates is None:
      print("\t[W] Nothing to keep for file " + file + ", station ignored")
    else:
      [header, writer] = cutData(header, vars, dates)
      if header is None and writer is None:
        print("\t[W] Variables too short for file " +
              file + ", station ignored")
      else:
        writeData(output + "/" + file, header, writer)


if __name__ == "__main__":
  main()
