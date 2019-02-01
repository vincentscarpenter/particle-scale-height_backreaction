import pencil
import numpy
from matplotlib import pylab
import sys

ts          = pencil.read_ts()
ts_single   = pencil.read_ts("../../single-particle_2048x2048/data/time_series.dat")

tmax        = max(ts.t[-1],ts_single.t[-1])

start_params = pencil.read_param()
grav         = start_params.gravx
tausp        = start_params.tausp

time_single = ts_single.t/tausp
time        = ts.t/tausp

vpxm        = ts.vpxm*2000
vpxmin      = ts.vpxmin*2000
vpxmax      = ts.vpxmax*2000
vpx_single  = ts_single.vpxm*2000

figure  = pylab.figure()
subplot = figure.add_subplot(111)
subplot.set_xlim([0.0,tmax/tausp])

subplot.plot(numpy.linspace(0.0,tmax/tausp,1000),numpy.repeat(2.45*tausp*2000,1000),linestyle="-",color="black",label=r"$\tau g$")
#subplot.plot(time_single,vpx_single,linestyle="--",color="gray",label="single particle run")
#subplot.plot(time,vpxmin,linestyle="none",marker=".",markeredgecolor="red",markerfacecolor="none",markersize=0.5,label="min")
#subplot.plot(time,vpxmax,linestyle="none",marker=".",markeredgecolor="green",markerfacecolor="none",markersize=0.5,label="max")
#subplot.plot(time,vpxm,linestyle="none",marker=".",markeredgecolor="black",markerfacecolor="none",markersize=0.5,label="mean")

subplot.plot(time,vpxmin,linestyle="-",color="red",label="min")
subplot.plot(time,vpxmax,linestyle="-",color="green",label="max")
subplot.plot(time,vpxm,linestyle="-",color="black",label="mean")

subplot.set_ylabel(r"$v_{rel}$ (mm/s)",fontsize=35)
subplot.set_xlabel("t (friction times)",fontsize=30)
subplot.legend(fontsize=30)

for tick in subplot.xaxis.get_major_ticks():
    tick.label.set_fontsize(25)
for tick in subplot.yaxis.get_major_ticks():
    tick.label.set_fontsize(25)

pylab.show()
