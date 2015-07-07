import numpy as np
import math
import coolmath

# ==============================================================================
def median_(lst):
  return np.median(np.array(lst))
# ==============================================================================

# ==============================================================================
def sort_array(x):
  y = sorted(x)
  return y
# ==============================================================================

# ==============================================================================
def mid_point(bin):
  y = int(math.floor( float(bin)*0.5 ) + 1)
  return y
# ==============================================================================

# ==============================================================================
def median_deviation(x):
 xmed = median_(x)
 res=[]
 for i in range(len(x)):
   res.append(abs(x[i]-xmed))
 y = median_(res)
 return y
# ==============================================================================

# ==============================================================================
def median_model(needs_sort,bin,x):
  # Check if sorting needed
  if needs_sort == 1:
    sort_array(x)
  # Compute median model
  n = len(x)
  jmax = n + 1 - bin
  midpoint = mid_point(bin)
  t_temp=[0 for i in range(bin)]
  f_temp=[0 for i in range(bin)]
  tmed=[]
  fmed=[]
  for j in range(jmax):
     jstart = j
     jend = j + bin - 1
     # Define temp lists
     for i in range(bin):
       t_temp[i] = x[jstart+i-1][0]
       f_temp[i] = x[jstart+i-1][1]
     tmed.append(t_temp[midpoint])
     fmed.append(median_(f_temp))
  # Define output array
  y=[[0 for j in range(2)] for i in range(n)]
  # First-segment
  for i in range(midpoint-1):
    y[i][0] = x[i][0]
    y[i][1] = x[i][1]
  # Mid-segment
  for i in range(n+2-2*midpoint):
    y[i+midpoint-1][0] = tmed[i]
    y[i+midpoint-1][1] = fmed[i]
  # Last-segment
  for i in range(n-midpoint+1,n):
    y[i][0] = x[i][0]
    y[i][1] = x[i][1]
  return y
# ==============================================================================

# ==============================================================================
def median_filter(needs_sort,bin,x):
  # Check if sorting needed
  if needs_sort == 1:
    sort_array(x)
  # Compute median model
  y = median_model(0,bin,x)
  # Define z array
  midpoint = mid_point(bin)
  zlen = len(x) + 2 - 2*midpoint
  z=[[0 for j in range(3)] for i in range(zlen)]
  for i in range(zlen):
    z[i][0] = x[i+midpoint-1][0]
    z[i][1] = x[i+midpoint-1][1]/y[i+midpoint-1][1]
    z[i][2] = x[i+midpoint-1][2]/y[i+midpoint-1][1]
  return z
# ==============================================================================

# ==============================================================================
def identify_outliers(scale_errors,needs_sort,bin,x):
  # Useful constants
  roottwo = 1.414213562373095
  madfac = 1.482602218505602
  n = len(x)
  # Get the residuals
  y = median_model(needs_sort,bin,x)
  res = []
  for i in range(n):
    res.append(y[i][1] - x[i][1])
  # Scale the errors
  if scale_errors == 1:
    maderror = madfac*median_deviation(res)
    errors = [row[2] for row in x] 
    mederror = median_(errors)
    sigmafac = maderror/mederror
    print 'Scaling errors by ',sigmafac
  else:
    sigmafac = 1.0
  # Number of sigmas to clip
  sigmas = coolmath.inverf(1.0-(1.0/n))
  # Identify outliers
  outlier=[0 for i in range(n)]
  for i in range(n):
    if abs(res[i]) >= sigmas*sigmafac*x[i][2]:
      outlier[i] = 1
    else:
      outlier[i] = 0
  print 'Number of outliers = ',sum(outlier),' [',100.0*sum(outlier)/n,'%]'
  return outlier
# ==============================================================================
  
# ==============================================================================
def detrend(needs_sort,scale_errors,bin_in,x):
  n = len(x)
  # Do not allow even bin sizes
  if (bin_in % 2 == 0): #even
    bin = bin_in + 1
    print 'Increasing bin size to ',bin
  else:
    bin = bin_in
  # Sort the array, if needed
  if needs_sort == 1:
    sort_array(x)
  # Find outliers
  outlier = identify_outliers(scale_errors,0,bin,x)
  # Define xclean
  nclean = n - sum(outlier)
  xclean=[[0 for j in range(3)] for i in range(nclean)]
  m = 0
  for i in range(n):
    if outlier[i] == 0:
      m = m + 1
      xclean[m-1][0] = x[i][0]
      xclean[m-1][1] = x[i][1]
      xclean[m-1][2] = x[i][2]
  # Execute median filter
  xnorm = median_filter(0,bin,xclean)
  # Check for any remaining outliers
  xnorm_outlier = identify_outliers(scale_errors,0,bin,xnorm)
  # Define final xout array
  nout = len(xnorm) - sum(xnorm_outlier)
  xout=[[0 for j in range(3)] for i in range(nout)]
  m = 0
  for i in range(len(xnorm)):
    if xnorm_outlier[i] == 0:
      m = m + 1
      xout[m-1][0] = xnorm[i][0]
      xout[m-1][1] = xnorm[i][1]
      xout[m-1][2] = xnorm[i][2]
  return xout
# ==============================================================================