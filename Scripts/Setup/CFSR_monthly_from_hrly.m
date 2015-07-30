%====================== CFSR_monthly_from_hourly.m =======================%
%=                                                                       =%
%=  Written by Austin C Toddd (NCSU, 2014) for personal use.             =%
%=                                                                       =%
%=  Program is designed to read CFSR hourly (or 6-hrly average) data     =%
%=  from file and compute monthly averages of the net heat flux and      =%
%=  net momentum fluxes.                                                 =%
%=                                                                       =%
%=========================================================================%

addpath(genpath('/Users/actodd/MYMATLAB/'));

%=========================================================================%
%=                     Set directories and filenames                     =%
%=========================================================================%
datadir='/Volumes/Black_box/Data/CFSR/';

varnames={'dlswfc';'dswsfc';'lhtfl';'shtfl';'ulwsfc';'uswsfc'};
filevars={'DLWRF_L1_Avg_1';'DSWRF_L1_Avg_1';...
          'LHTFL_L1_Avg_1';'sHTFL_L1_Avg_1';...
          'ULWRF_L1_Avg_1';'USWRF_L1_Avg_1'};


%=========================================================================%
%=                     Set directories and filenames                     =%
%=========================================================================%


