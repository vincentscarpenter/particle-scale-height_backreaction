import pencil
import numpy
from matplotlib import pylab
import sys

ts = pencil.read_ts()
ts_dictionary = ts.__dict__

t_start = ts.t[0]
t_stop  = ts.t[-1]

if(len(sys.argv) == 1):
    print("Need timeseries value to plot. Options are: ")
    print(ts_dictionary.keys())
    quantity_name = raw_input("Please make a selection, or type Q to quit: ")
    if(quantity_name == "Q"):
        sys.exit("Quitting...")
elif(len(sys.argv) == 2):
    quantity_name = sys.argv[1]
else:
    quantity_name = sys.argv[1]
    if(int(sys.argv[2]) < t_start or int(sys.argv[3]) > t_stop):
        sys.exit("Time out of bounds.")
    else:
        t_start = int(sys.argv[2])
        t_stop  = int(sys.argv[3])

try:
    quantity = ts_dictionary[quantity_name]
except KeyError:
    sys.exit(quantity_name + " is not the name of a quantity in the time series.")

i_start = numpy.where(numpy.abs(t_start - ts.t) == numpy.min(numpy.abs(t_start - ts.t)))[0][0]
i_stop  = numpy.where(numpy.abs(t_stop  - ts.t) == numpy.min(numpy.abs(t_stop  - ts.t)))[0][0]

figure  = pylab.figure()
subplot = figure.add_subplot(111)
subplot.plot(ts.t[i_start:i_stop],quantity[i_start:i_stop])

subplot.set_title(quantity_name + " vs t")
subplot.set_ylabel(quantity_name)
subplot.set_xlabel("t (code units)")

pylab.show()
