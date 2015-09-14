import numpy
import matplotlib
from matplotlib import pyplot as plt

f = open('kkpl.dat', 'r')

nlist = []

for line in f.readlines():
      row = []
    for valstring in line.split('\t'):
                row.append(float(valstring))
    nlist.append(row)

                arr = numpy.array(nlist)

                  # color map determines the colors used on the plot
             mycmap = plt.get_cmap(name='rainbow')    # See matplotlib gallery

                  # normalisation determines how the color used depends on the value
                  #   - eg logarithmic, linear, etc
            mynorm = matplotlib.colors.LogNorm()     # See matplotlib docs

    plt.imshow(arr, cmap=mycmap, norm=mynorm)

