def AtmModel_get_times(time):
  import datetime
  # gets a dictionary with times for ERA-Interim dataset
  # input is the variable as one get's it by 'time = f.variables['time']'
  # written by K.Barfus 2/2019
  
  time_units = time.units  # hours since 2007-3-1- 00:00:00
  temp_str = time_units.split()
  time_var = time[:]
  date = temp_str[2] # date
  clock_time = temp_str[3] # time
  date_str = date.split("-")
  ref_year = int(date_str[0])
  ref_month = int(date_str[1])
  ref_day = int(date_str[2])
  clock_time_str = clock_time.split(":")
  ref_hour = int(clock_time_str[0])
  ref_minute = int(clock_time_str[1])
  ref_second = int(clock_time_str[2])
  ref_date = datetime.datetime(ref_year,ref_month,ref_day,ref_hour,ref_minute,ref_second)
  n_times = len(time_var)
  dates = []
  for i_times in range(0, n_times):
    d_time = datetime.timedelta(hours=time_var[i_times])
    rec_date = ref_date + d_time
    dates.append(rec_date)
  return dates
