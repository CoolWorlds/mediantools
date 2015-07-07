def simple_read(name,ncols):
  f = open(name, 'r')
  flen=len(f.readlines())
  f = open(name, 'r')
  line = f.readlines()
  x=[[0 for j in range(ncols)] for i in range(flen)]
  for i in range(flen):
    line[i] = line[i].strip()
    columns = line[i].split()
    for j in range(ncols):
      x[i][j] = float(columns[j])
  return x
