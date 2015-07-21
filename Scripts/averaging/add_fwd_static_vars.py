#===============================================================================
#
# SCRIPT :  add_fwd_vars.py
#
# PURPOSE : Ingest data from US East water age forward model and average over
#           the five years of forward model state.  Several variables already
#           have been averaged over the 5 year period for use in the LTRANS
#           simulations.  Several new netCDF dimensions are also added to the
#           existing file.
#
# METHOD :
#
# HISTORY : Created by Austin Todd on 11 February 2015.
#
#===============================================================================

#-------------------------------------------------------------------------------
# Define all required libraries, routines, and modules
#-------------------------------------------------------------------------------
import numpy as np
import netCDF4
import sys
import os

#-------------------------------------------------------------------------------
# Define all required sub-libraries, sub-routines, and sub-modules
#-------------------------------------------------------------------------------
from netCDF4 import Dataset

#-------------------------------------------------------------------------------
# Print message to screen
#-------------------------------------------------------------------------------
print "----------------------------------------------------------------------"
print " RUNNING: Script add_fwd_vars.py"
print "----------------------------------------------------------------------"
print ""

################################################################################
#                                                                              #
#                            User Defined Variables                            #
#                                                                              #
infile1 = '/gpfs_common/he_share/actodd/USeast-age/output/clim/fwdrun/part1/'+\
          'useast_his.nc'
#infile2 = '/gpfs_common/he_share/actodd/USeast-age/output/clim/fwdrun/part2/'+\
#          'useast_his.nc'
outfile = '/gpfs_common/he_share/actodd/USeast-age/output/clim/fwdrun/temp.nc'#+\
#           'daily_avg.nc'
#                                                                              #
################################################################################

#-------------------------------------------------------------------------------
# Read in additional grid variables from forward file
#-------------------------------------------------------------------------------
indata1=Dataset(inpfile1,mode='r')
#indata2=Dataset(inpfile2,mode='r')

#-------------------------------------------------------------------------------
# Open output file and create dimensions and variables
#-------------------------------------------------------------------------------
print " "
print "Creating output file: %s" % output_file
print " "

outdata=Dataset(outfile,mode='w')
outdata.createDimension('xi_psi'  ,len(indata1.dimensions['xi_psi'  ]))
outdata.createDimension('eta_psi' ,len(indata1.dimensions['eta_psi' ]))
outdata.createDimension('tracer'  ,len(indata1.dimensions['tracer'  ]))
outdata.createDimension('boundary',len(indata1.dimensions['boundary']))

#-------------------------------------------------------------------------------
# Create several of the static variables
#-------------------------------------------------------------------------------
print 'ntimes'
var= outdata.createVariable('ntimes',indata1.variables['ntimes'].dtype)
setattr(var,'long_name','number of long time-steps')
var = indata1.variables['ntimes']

print 'ndtfast'
var= outdata.createVariable('ndtfast',indata1.variables['ndtfast'].dtype)
setattr(var,'long_name','number of short time-steps')
var = indata1.variables['ndtfast']

print 'dt'
var= outdata.createVariable('dt',indata1.variables['dt'].dtype)
setattr(var,'long_name','size of long time-steps')
setattr(var,'units','second')
var = indata1.variables['dt']

print 'dtfast'
var= outdata.createVariable('dtfast',indata1.variables['dtfast'].dtype)
setattr(var,'long_name','size of short time-steps')
setattr(var,'units','second')
var = indata1.variables['dtfast']

print 'dstart'
var= outdata.createVariable('dstart',indata1.variables['dstart'].dtype)
setattr(var,'long_name','time stamp assigned to model initilization')
setattr(var,'units','days since 2004-01-01 00:00:00')
var = indata1.variables['dstart']

print 'nHIS'
var= outdata.createVariable('nHIS',indata1.variables['nHIS'].dtype)
setattr(var,'long_name','number of time-steps between history records')
var = indata1.variables['nHIS']

print 'nRST'
var= outdata.createVariable('nRST',indata1.variables['nRST'].dtype)
setattr(var,'long_name','number of time-steps between restart records')
var = indata1.variables['nRST']

print 'Falpha'
var= outdata.createVariable('Falpha',indata1.variables['Falpha'].dtype)
setattr(var,'long_name','Power-law shape barotropic filter parameter')
var = indata1.variables['Falpha']

print 'Fbeta'
var= outdata.createVariable('Fbeta',indata1.variables['Fbeta'].dtype)
setattr(var,'long_name','Power-law shape barotropic filter parameter')
var = indata1.variables['Fbeta']

print 'Fgamma'
var= outdata.createVariable('Fgamma',indata1.variables['Fgamma'].dtype)
setattr(var,'long_name','Power-law shape barotropic filter parameter')
var = indata1.variables['Fgamma']

print 'tnu2'
var= outdata.createVariable('tnu2',indata1.variables['tnu2'].dtype,('tracer'))
setattr(var,'long_name','Laplacian mixing coefficient for tracers')
setattr(var,'units','meter2 second-1')
var = indata1.variables['tnu2']

print 'visc2'
var= outdata.createVariable('visc2',indata1.variables['visc2'].dtype)
setattr(var,'long_name','Laplacian mixing coefficient for momentum')
setattr(var,'units','meter2 second-1')
var = indata1.variables['visc2']

print 'Akt_bak'
var= outdata.createVariable('Akt_bak',indata1.variables['Akt_bak'].dtype,('tracer'))
setattr(var,'long_name','background vertical mixing coefficient for tracers')
setattr(var,'units','meter2 second-1')
var[:] = indata1.variables['Akt_bak'][:]

print 'Akv_bak'
var= outdata.createVariable('Akv_bak',indata1.variables['Akv_bak'].dtype)
setattr(var,'long_name','background vertical mixing coefficient for momentum')
setattr(var,'units','meter2 second-1')
var = indata1.variables['Akv_bak']

print 'Akk_bak'
var= outdata.createVariable('Akk_bak',indata1.variables['Akk_bak'].dtype)
setattr(var,'long_name','background vertical mixing coefficient for turbulent energy')
setattr(var,'units','meter2 second-1')
var = indata1.variables['Akk_bak']

print 'Akp_bak'
var= outdata.createVariable('Akp_bak',indata1.variables['Akp_bak'].dtype)
setattr(var,'long_name','background vertical mixing coefficient for length scale')
setattr(var,'units','meter2 second-1')
var = indata1.variables['Akp_bak']

print 'rdrg'
var= outdata.createVariable('rdrg',indata1.variables['rdrg'].dtype)
setattr(var,'long_name','linear drag coefficient')
setattr(var,'units','meter second-1')
var = indata1.variables['rdrg']

print 'rdrg2'
var= outdata.createVariable('rdrg2',indata1.variables['rdrg2'].dtype)
setattr(var,'long_name','quadratic drag coefficient')
var = indata1.variables['rdrg2']

print 'Zob'
var= outdata.createVariable('Zob',indata1.variables['Zob'].dtype)
setattr(var,'long_name','bottom roughness')
setattr(var,'units','meter')
var = indata1.variables['Zob']

print 'Zos'
var= outdata.createVariable('Zos',indata1.variables['Zos'].dtype)
setattr(var,'long_name','surface roughness')
setattr(var,'units','meter')
var = indata1.variables['Zos']

print 'Znudg'
var= outdata.createVariable('Znudg',indata1.variables['Znudg'].dtype)
setattr(var,'long_name','free-surface nudging/relaxation inverse time scale')
setattr(var,'units','day-1')
var = indata1.variables['Znudg']

print 'M2nudg'
var= outdata.createVariable('M2nudg',indata1.variables['M2nudg'].dtype)
setattr(var,'long_name','2D momentum nudging/relaxation inverse time scale')
setattr(var,'units','day-1')
var = indata1.variables['M2nudg']

print 'M3nudg'
var= outdata.createVariable('M3nudg',indata1.variables['M3nudg'].dtype)
setattr(var,'long_name','3D momentum nudging/relaxation inverse time scale')
setattr(var,'units','day-1')
var = indata1.variables['M3nudg']

print 'Tnudg'
var= outdata.createVariable('Tnudg',indata1.variables['Tnudg'].dtype)
setattr(var,'long_name','Tracers nudging/relaxation inverse time scale')
setattr(var,'units','day-1')
var = indata1.variables['Tnudg']

print 'FSobc_in'
var= outdata.createVariable('FSobc_in',indata1.variables['FSobc_in'].dtype,('boundary'))
setattr(var,'long_name','free-surface inflow, nudging inverse time scale')
setattr(var,'units','second-1')
var[:] = indata1.variables['FSobc_in'][:]

print 'FSobc_out'
var= outdata.createVariable('FSobc_out',indata1.variables['FSobc_out'].dtype,('boundary'))
setattr(var,'long_name','free-surface outflow, nudging inverse time scale')
setattr(var,'units','second-1')
var[:] = indata1.variables['FSobc_out'][:]

print 'M2obc_in'
var= outdata.createVariable('M2obc_in',indata1.variables['M2obc_in'].dtype,('boundary'))
setattr(var,'long_name','2D momentum inflow, nudging inverse time scale')
setattr(var,'units','second-1')
var[:] = indata1.variables['M2obc_in'][:]

print 'M2obc_out'
var= outdata.createVariable('M2obc_out',indata1.variables['M2obc_out'].dtype,('boundary'))
setattr(var,'long_name','2D momentum outflow, nudging inverse time scale')
setattr(var,'units','second-1')
var[:] = indata1.variables['M2obc_out'][:]

print 'Tobc_in'
var= outdata.createVariable('Tobc_in',indata1.variables['Tobc_in'].dtype,('boundary','tracer'))
setattr(var,'long_name','Tracers inflow, nudging inverse time scale')
setattr(var,'units','second-1')
var[:] = indata1.variables['Tobc_in'][:]

print 'Tobc_out'
var= outdata.createVariable('Tobc_out',indata1.variables['Tobc_out'].dtype,('boundary','tracer'))
setattr(var,'long_name','Tracers inflow, nudging inverse time scale')
setattr(var,'units','second-1')
var[:] = indata1.variables['Tobc_out'][:]

print 'M3obc_in'
var= outdata.createVariable('M3obc_in',indata1.variables['M3obc_in'].dtype,('boundary'))
setattr(var,'long_name','3D momentum inflow, nudging inverse time scale')
setattr(var,'units','second-1')
var[:] = indata1.variables['M3obc_in'][:]

print 'M3obc_out'
var= outdata.createVariable('M3obc_out',indata1.variables['M3obc_out'].dtype,('boundary'))
setattr(var,'long_name','3D momentum outflow, nudging inverse time scale')
setattr(var,'units','second-1')
var[:] = indata1.variables['M3obc_out'][:]

print 'rho0'
var= outdata.createVariable('rho0',indata1.variables['rho0'].dtype)
setattr(var,'long_name','mean density used in Boussinesq approximation')
setattr(var,'units','kilogram meter-3')
var = indata1.variables['rho0']

print 'gamma2'
var= outdata.createVariable('gamma2',indata1.variables['gamma2'].dtype)
setattr(var,'long_name','slipperiness parameter')
var = indata1.variables['gamma2']

print 'spherical'
var= outdata.createVariable('spherical',indata1.variables['spherical'].dtype)
setattr(var,'long_name','grid type logical switch')
setattr(var,'option_T','spherical')
setattr(var,'option_F','Cartesian')
var = indata1.variables['spherical']

print 'xl'
var= outdata.createVariable('xl',indata1.variables['xl'].dtype)
setattr(var,'long_name','domain length in the XI-direction')
setattr(var,'units','meter')
var = indata1.variables['xl']

print 'el'
var= outdata.createVariable('el',indata1.variables['el'].dtype)
setattr(var,'long_name','domain length in the ETA-direction')
setattr(var,'units','meter')
var = indata1.variables['el']

print 'Vtransform'
var= outdata.createVariable('Vtransform',indata1.variables['Vtransform'].dtype)
setattr(var,'long_name','vertical terrain-following transformation equation')
var = indata1.variables['Vtransform']

print 'Vstretching'
var= outdata.createVariable('Vstretching',indata1.variables['Vstretching'].dtype)
setattr(var,'long_name','vertical terrain-following stretching function')
var = indata1.variables['Vstretching']

print 'theta_s'
var= outdata.createVariable('theta_s',indata1.variables['theta_s'].dtype)
setattr(var,'long_name','S-coordinate surface control parameter')
var = indata1.variables['theta_s']

print 'theta_b'
var= outdata.createVariable('theta_b',indata1.variables['theta_b'].dtype)
setattr(var,'long_name','S-coordinate bottom control parameter')
var = indata1.variables['theta_b']

print 'Tcline'
var= outdata.createVariable('Tcline',indata1.variables['Tcline'].dtype)
setattr(var,'long_name','S-coordinate surface/bottom layer width')
setattr(var,'units','meter')
var = indata1.variables['Tcline']

print 'hc'
var= outdata.createVariable('hc',indata1.variables['hc'].dtype)
setattr(var,'long_name','S-coordinate parameter, critical depth')
setattr(var,'units','meter')
var = indata1.variables['hc']

#-------------------------------------------------------------------------------
# Create several of the grid variables
#-------------------------------------------------------------------------------
print 'f'
var = outdata.createVariable('f','float32',('eta_rho','xi_rho'))
setattr(var,'long_name','Coriolis parameter at RHO-points')
setattr(var,'units','second-1')
setattr(var,'coordinates','lon_rho lat_rho')
var[:] = indata1.variables['f'][:]

print 'pm'
var = outdata.createVariable('pm','float32',('eta_rho','xi_rho'))
setattr(var,'long_name','curvilinear coordinate metric in XI')
setattr(var,'units','meter-1')
setattr(var,'coordinates','lon_rho lat_rho')
var[:] = indata1.variables['pm'][:]

print 'pn'
var = outdata.createVariable('pn','float32',('eta_rho','xi_rho'))
setattr(var,'long_name','curvilinear coordinate metric in ETA')
setattr(var,'units','meter-1')
setattr(var,'coordinates','lon_rho lat_rho')
var[:] = indata1.variables['pn'][:]

print 'lon_u'
var = outdata.createVariable('lon_u','float32',('eta_u','xi_u'))
setattr(var,'long_name','longitude of U-points')
setattr(var,'units','degree_east')
var[:] = indata1.variables['lon_u'][:]

print 'lat_u'
var = outdata.createVariable('lat_u','float32',('eta_u','xi_u'))
setattr(var,'long_name','latitude of U-points')
setattr(var,'units','degree_north')
var[:] = indata1.variables['lat_u'][:]

print 'lon_v'
var = outdata.createVariable('lon_v','float32',('eta_v','xi_v'))
setattr(var,'long_name','longitude of V-points')
setattr(var,'units','degree_east')
var[:] = indata1.variables['lon_v'][:]

print 'lat_v'
var = outdata.createVariable('lat_v','float32',('eta_v','xi_v'))
setattr(var,'long_name','latitude of V-points')
setattr(var,'units','degree_north')
var[:] = indata1.variables['lat_v'][:]

print 'lon_psi'
var = outdata.createVariable('lon_psi','float32',('eta_psi','xi_psi'))
setattr(var,'long_name','longitude of PSI-points')
setattr(var,'units','degree_east')
var[:] = indata1.variables['lon_psi'][:]

print 'lat_psi'
var = outdata.createVariable('lat_psi','float32',('eta_psi','xi_psi'))
setattr(var,'long_name','latitude of PSI-points')
setattr(var,'units','degree_north')
var[:] = indata1.variables['lat_psi'][:]

print 'angle'
var = outdata.createVariable('angle','float32',('eta_rho','xi_rho'))
setattr(var,'long_name','angle between XI-axis and EAST')
setattr(var,'units','radians')
setattr(var,'coordinates','lon_rho lat_rho')
var[:] = indata1.variables['angle'][:]

print 'mask_u'
var = outdata.createVariable('mask_u','float32',('eta_u','xi_u'))
setattr(var,'long_name','mask on U-points')
setattr(var,'option_0','land')
setattr(var,'option_1','water')
setattr(var,'coordinates','lon_u lat_u')
var[:] = indata1.variables['mask_u'][:]

print 'mask_v'
var = outdata.createVariable('mask_v','float32',('eta_v','xi_v'))
setattr(var,'long_name','mask on V-points')
setattr(var,'option_0','land')
setattr(var,'option_1','water')
setattr(var,'coordinates','lon_v lat_v')
var[:] = indata1.variables['mask_v'][:]

print 'mask_psi'
var = outdata.createVariable('mask_psi','float32',('eta_psi','xi_psi'))
setattr(var,'long_name','mask on PSI-points')
setattr(var,'option_0','land')
setattr(var,'option_1','water')
setattr(var,'coordinates','lon_psi lat_psi')
var[:] = indata1.variables['mask_psi'][:]

indata1.close()
#indata2.close()
outdata.close()




