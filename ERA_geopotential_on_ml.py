from netCDF4 import Dataset
import numpy as np

 
def ERA_geopotential_on_ml(temp,spec_hum,press_half_level,press_mid_level,geopotential_sfc):
  # calculates the geopotential on model levels for ERA-Interim based on 
  # Eq. 2.20, 2.21, 2.22, 2.23, in Chapter 2 Basic equations and discretisation, of Part III Dynamics and numerical procedures, at:
  # https://www.ecmwf.int/sites/default/files/elibrary/2007/9220-part-iii-dynamics-and-numerical-procedures.pdf
  # be aware that indexing in the document starts at the lowest level but indexing in the routine starts at the uppermost level (this is what you get when reading the data with Python netCDF4) 
  # due to the table at https://www.ecmwf.int/en/forecasts/documentation-and-support/60-model-levels geopotential for the uppermost level is undefined (only Major Tom knows...)
  # input is 
  # temp: temperature[ntime,nlevel,n_lat,nlon] on model levels [K]
  # spec_hum: specific_humidity[ntime,nlevel,n_lat,nlon] on model levels [kg/kg]
  # press_half_level: pressure on half model levels [ntime,nlevel+1,n_lat,nlon] [Pa]
  # press_mid_level: pressure on mid/full model levels [ntime,nlevel,n_lat,nlon] [Pa]
  # geopotential_sfc: surface geopotential [ntime,1,n_lat,nlon] [m**2 s**-2]
  # output is:
  # geopotential on full model levels [ntime,nlevel,n_lat,nlon] [m**2 s**-2]
  # geopotential on half model levels [ntime,nlevel+1,n_lat,nlon] [m**2 s**-2]   
  # written by Klemens Barfus (klemens.barfus at gmx.de)

  import numpy as np  

  Rd = 287.06
  g = 9.80665 
  ndim = temp.ndim
  if(ndim == 1):
    temp = temp[np.newaxis,:,np.newaxis,np.newaxis]
    spec_hum = spec_hum[np.newaxis,:,np.newaxis,np.newaxis]  
    press_half_level = press_half_level[np.newaxis,:,np.newaxis,np.newaxis]
    press_mid_level = press_mid_level[np.newaxis,:,np.newaxis,np.newaxis]
    if(np.isscalar(geopotential_sfc)):
      geopotential_sfc = np.asarray([geopotential_sfc])
    geopotential_sfc = geopotential_sfc[np.newaxis,:,np.newaxis,np.newaxis]
  if(ndim == 3):
    temp = temp[np.newaxis,:,:,:]
    spec_hum = spec_hum[np.newaxis,:,:,:]
    press_half_level = press_half_level[np.newaxis,:,:,:]
    press_mid_level = press_mid_level[np.newaxis,:,:,:]
    geopotential_sfc = geopotential_sfc[np.newaxis,:,:,:]
  # check order of arrays
  if(press_mid_level[0,0,0,0] > press_mid_level[0,10,0,0]): # low -> up order
    order = "ltu"
  else:
    order = "utl"
    temp = temp[:,::-1,:,:]
    spec_hum = spec_hum[:,::-1,:,:]
    press_half_level = press_half_level[:,::-1,:,:]
    press_mid_level = press_mid_level[:,::-1,:,:]
    geopotential_sfc = geopotential_sfc[:,::-1,:,:]
  nxyz = temp.shape
  n_time = nxyz[0]
  n_level = nxyz[1]+1
  n_lat = nxyz[2]
  n_lon = nxyz[3]
  geopot_half_level = np.zeros((n_time,n_level,n_lat,n_lon))
  geopot_half_level.fill(np.nan)
  geopot_full_level = np.zeros((n_time,n_level-1,n_lat,n_lon))
  geopot_full_level.fill(np.nan)
  # compute on half levels
  for i_time in range(0, n_time):
    #print(i_time, n_time)  
    for i_lat in range(0, n_lat):
      for i_lon in range(0, n_lon):
        i_level = 0
        while(i_level <= n_level-2):
           if(i_level == 0):
             geopot_half_level[i_time,i_level,i_lat,i_lon] = geopotential_sfc[i_time,0,i_lat,i_lon]
           else:
             p1 = press_half_level[i_time,i_level-1,i_lat,i_lon]
             p2 = press_half_level[i_time,i_level,i_lat,i_lon]
             t_level = temp[i_time,i_level-1,i_lat,i_lon]
             q_level = spec_hum[i_time,i_level-1,i_lat,i_lon]
             t_level = t_level * (1.+0.609133*q_level)
             if((p1 != 0.0) and (p2 != 0.0)):
               delta_z = Rd * t_level * np.log(p1/p2)
             else:
               delta_z = np.nan
             geopot_half_level[i_time,i_level,i_lat,i_lon] = geopot_half_level[i_time,i_level-1,i_lat,i_lon] + delta_z
           i_level = i_level + 1
  # compute on full levels
  for i_time in range(0, n_time):
    #print(i_time, n_time)
    for i_lat in range(0, n_lat):
      for i_lon in range(0, n_lon):
        for i_level in range(0, n_level-1):
          t_level = temp[i_time,i_level,i_lat,i_lon]
          q_level = spec_hum[i_time,i_level,i_lat,i_lon]
          t_level = t_level * (1.+0.609133*q_level)
          if(i_level == n_level-2):
            a = np.log(2.0)
          else:
            p1 = press_half_level[i_time,i_level+1,i_lat,i_lon]
            p2 = press_half_level[i_time,i_level,i_lat,i_lon]
            delta_p = p2 - p1 # Eq 2.13
            a = 1.0 - (p1/delta_p) * np.log(p2/p1)
          geopot_full_level[i_time,i_level,i_lat,i_lon] =  geopot_half_level[i_time,i_level,i_lat,i_lon] + a * Rd * t_level
          # print(geopot_half_level[i_time,i_level,i_lat,i_lon],geopot_full_level[i_time,i_level,i_lat,i_lon],geopot_half_level[i_time,i_level+1,i_lat,i_lon])
  if(order == "utl"):
    geopot_full_level = geopot_full_level[:,::-1,:,:]
    geopot_half_level = geopot_half_level[:,::-1,:,:]
  if(ndim == 1):
    geopot_full_level = geopot_full_level[0,:,0,0]
    geopot_half_level = geopot_half_level[0,:,0,0]
  if(ndim == 3):
    geopot_full_level = geopot_full_level[0,:,:,:]
    geopot_half_level = geopot_half_level[0,:,:,:]

  return geopot_full_level, geopot_half_level
  
