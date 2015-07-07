import mediantools
import simpleread

x = simpleread.simple_read("raw.dat",3)

bin = 51
needs_sort = 1
scale_errors = 1

#y = mediantools.median_model(needs_sort,bin,x)
#y = mediantools.median_filter(needs_sort,bin,x)
#y = mediantools.identify_outliers(scale_errors,needs_sort,bin,x)
y = mediantools.detrend(needs_sort,scale_errors,bin,x)

with open("raw.med.dat", 'w') as f:
    f.writelines(' '.join(str(j) for j in i) + '\n' for i in y)
