# README

Source code for:
Climate change scenarios at hourly time-step over Switzerland
from an enhanced temporal downscaling approach

International Journal of Climatology, 2021, doi:10.1002/joc.7032

Adrien Michel, Varun Sharma, Michael Lehning, and Hendrik Huwald

School of Architecture, Civil and Environmental Engineering, Ecole Polytechnique
Fédérale de Lausanne (EPFL), Lausanne, Switzerland

WSL Institute for Snow and Avalanche Research SLF, Davos, Switzerland

## Introduction
This directory contains the source code to perform temporal donswnscaling of
climate change daily time series to hourly time series by using a delta change
method as described in the paper cited above. A careful reading of the paper is required
in order to understand the different steps performed by the scripts provided.

The first script is used to cut and filter input historical meteorological
data for missing values. The second script computes the DOY means.
The third script performs the actual downscaling.

## Usage
Please see the instruction at the beginning of each script. The scripts are
adapted to the datasets used in the paper. Instructions about where and how to
modify the scripts in order to use different data source are given in the scripts.

## Citation
This repository should be cited as:
Michel, A., 2021. Source code for: Climate change scenarios at hourly time-step
over Switzerland from an enhanced temporal downscaling approach. EnviDat.
doi:10.16904/envidat.203.

## License
Copyright (C) 2021, Adrien Michel, EPFL

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <https://www.gnu.org/licenses/>.
