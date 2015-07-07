import math

def inverf(x):
  # Constants
  a_winitzki = 0.1400122886866666
  b_winitzki = 4.546884979448284
  # Compute factor
  if x<0.0:
    factor = -1.0
  else:
    factor = 1.0
  # Compute LDCs
  logxsq = math.log(1.0-x*x)
  y = b_winitzki + 0.5*logxsq
  y = math.sqrt( y*y - (logxsq/a_winitzki) ) - y
  y = factor*math.sqrt(y)
  return y
