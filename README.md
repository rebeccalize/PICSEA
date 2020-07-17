# PICSEA

PICSEA - Predicting the Impacts of Cyclones in South-East Africa - is a research project undertaken at the National Centre for Atmospheric Science at the University of Reading, funded by the Natural Environment Research Council and Department for International Development under the SHEAR consortium. 

This repository contains code used for data analysis as part of the PICSEA project, for example to assess the skill of the UK Met Office and ECMWF foreacsts of tropical cyclone track and intensity, and tropical cyclone rainfall. It makes use of data such as forecast and observed tropical cyclone tracks, tracked using Kevin Hodges' TRACK software, rainfall observations, satellite rainfall data and forecast rainfall data, typically in text or netcdf format.

The repository is split into folders containing code used to complete a specific task, such as calculating the track statistics for a certain forecasting system, or producing precipitation composites or plotting data. It's also split into JASMIN (UKMO) and ECMWF, and many of the scripts are the same/very similar between the two, because I needed copies on both machines where I was running the UKMO tracking and analysis on JASMIN, and the ECMWF tracking and analysis at ECMWF. 
