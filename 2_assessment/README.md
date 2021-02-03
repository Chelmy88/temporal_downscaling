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
This directory contains the source code to perform assessment of the downscaled
time series. A careful reading of the paper cited above is required in order to understand the
different steps performed by the scripts provided.

The two python scripts are used to extract the seasonal mean from the raw CH2018
data (Get_seas_error_raw.py) and from the downscaled time series (Get_seas_error.py).
These values are then used in the first part of the R script (plot_Assesments.R)
in order to perform the seasonal analysis described in Section 4.1 in the paper.
The second part of the R script describes the added variability analysis
presented in Section 4.2 of the paper. The last part of the code produces the
assessment plots presented in Section 4.3 and S4 of the paper and SI.

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
