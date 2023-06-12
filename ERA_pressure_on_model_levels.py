# get pressure at model levels
def ERA_pressure_on_model_levels(ps, hyam, hybm, hyai, hybi):
  ## calculates the pressure [Pa] on models levels for ERA-files ###
  # input is var152 [time,lev,lat,lon] or [lev,lat,lon] or [lev] : logarithm of surface pressure [Pa]
  # or var134 [time,lev,lat,lon] or [lev,lat,lon] or [lev] : surface pressure [Pa] 
  import numpy as np
  ndim = ps.ndim
  if(ndim == 3):
    ps = ps[np.newaxis,:] # add time axis
  else:
    if(ndim == 1):
      ps = ps[np.newaxis,:,np.newaxis,np.newaxis]
  nxy = ps.shape
  n_hym = len(hyam) # half level coefficients
  n_hyi = len(hyai) # full level coefficients
  n_time = nxy[0]
  n_lat = nxy[2]
  n_lon = nxy[3] 
  press_mid_level = np.zeros((n_time,n_hym,n_lat,n_lon))
  press_half_level = np.zeros((n_time,n_hyi,n_lat,n_lon))
  for i_time in range(0, n_time):
    if(ps[i_time,:,:,:] < 50.): # data are log_pressure(Pa)
      for i_level in range(0, n_hym):
        press_mid_level[i_time,i_level,:,:] = np.exp(ps[i_time,:,:,:]) * hybm[i_level] + hyam[i_level]
      for i_level in range(0, n_hyi):
        press_half_level[i_time,i_level,:,:] = np.exp(ps[i_time,:,:,:]) * hybi[i_level] + hyai[i_level]
    else: # data is pressure [Pa]
      for i_level in range(0, n_hym):
        press_mid_level[i_time,i_level,:,:] = ps[i_time,:,:,:] * hybm[i_level] + hyam[i_level]
      for i_level in range(0, n_hyi):
        press_half_level[i_time,i_level,:,:] = ps[i_time,:,:,:] * hybi[i_level] + hyai[i_level]
  if(ndim < 4):
    press_mid_level = np.squeeze(press_mid_level)
    press_half_level = np.squeeze(press_half_level) 

  return press_mid_level, press_half_level
