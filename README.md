# CodeDemo-For-TerraQuanta

## Description
This collection of codes is used to determine whether there are any mesoscale convective systems (MCSs) 
using FY2E TBB data.

- Programs
  - *select_MCS_FY2E_Tianshan.f*
    - main code
    - read in FY2E data
    - invoke R and NCL programs
  - *cluster_MCS.R*
    - cluster grids that meet conditions
  - *write_MCS_to_ncdata.ncl*
    - write MCS cluster results into NetCDF files
    - plot MCSs

- DATA
  - FY2E_TBB_IR1_OTG_20130601_1800.AWX
    - FY2E example data
  - GTOPO30_5MIN.CDF
    - topographic data fo plotting
  - MCS_grid_FY2E_Xinjiang_2013060118.nc
    - output data
- PIC
  - example picture

## Notice
For ifort, use '-assume byterecl' to compile 'select_MCS_FY2E_Xinjiang.f'.
