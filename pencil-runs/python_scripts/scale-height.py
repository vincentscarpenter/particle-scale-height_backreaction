import pencil
import numpy
from matplotlib import pyplot
from matplotlib import cm
import sys
import os
from scipy.optimize import curve_fit

def gauss(x, amplitude, mean, sigma):
    return amplitude*numpy.exp(-((x - mean)**2)/(2.*sigma**2))

base_path,script_name   = os.path.split(sys.argv[0])
scratch,simulation_name = os.path.split(base_path)
script_name             = script_name[:-3]

data        = pencil.read_var(trimall=True)
pdata       = pencil.read_pvar()
time_series = pencil.read_ts()
parameters  = pencil.read_param()

npar   = len(pdata.ipars)
xgrid  = data.x
ygrid  = data.y
zgrid  = data.z
dxgrid = data.dx
dygrid = data.dy
dzgrid = data.dz
x0 = xgrid[0]
y0 = ygrid[0]
z0 = zgrid[0]
x1 = xgrid[-1]
y1 = ygrid[-1]
z1 = zgrid[-1]
last_snap = numpy.floor(time_series.t[-1]/parameters.tausp) + 1

ivar_lower = 0
ivar_upper = last_snap
if(len(sys.argv) == 3):
    ivar_lower = int(sys.argv[1])
    ivar_upper = int(sys.argv[2])
    nvar       = ivar_upper - ivar_lower + 1
if(ivar_lower < 0 or ivar_upper > last_snap or nvar < 1):
    print("ivar range is not valid. Exiting...")
    sys.exit()
else:
    print("Tracking scale heights from ivar = " + str(ivar_lower) + " to " + str(ivar_upper))

resolutions   = [100,250,500,750,1000]
nres          = len(resolutions)
scale_height  = numpy.zeros((nvar,nres))
axes          = numpy.empty(nres)
fits_fig,axes = pyplot.subplots(1,nres)
fits_fig.set_size_inches(10.8,10.8/nres)
fits_fig.set_dpi(100)
for ivar in range(ivar_lower,ivar_upper+1):
    print("----------- ivar = " + str(ivar) + " -----------")
    var   = "PVAR" + str(ivar)
    pdata = pencil.read_pvar(varfile=var)
    xp = pdata.xp
    yp = pdata.yp
    zp = pdata.zp
    for ires in range(nres):
        print("       **** ires = " + str(ires) + " ****")
        nzsteps     = resolutions[ires]
        z_coords,dz = numpy.linspace(start=z0,stop=z1,num=nzsteps,endpoint=False,retstep=True)
        densities   = numpy.zeros(nzsteps)
        for it in range(nzsteps):
            z_loop        = z_coords[it]
            densities[it] = len(numpy.intersect1d(numpy.where(zp >= z_loop)[0],numpy.where(zp < z_loop + dz)[0]))
            progress = numpy.round(it/(nzsteps-1.0),2)
            sys.stdout.flush()
            sys.stdout.write("density loop progress: " + str(progress) + "\r")
#
# Cut out the very low density parts
#
        density_threshold = 0.1
        density_max       = densities.max()
        above_threshold   = numpy.where(densities > density_threshold*density_max)
        densities_trunc   = densities[above_threshold]
        z_coords_trunc    = z_coords[above_threshold] + dz/2.
#
# Take reasonable guesses for the fit
#
        amplitude_guess                 = density_max
        mean_guess                      = z_coords_trunc[numpy.where(densities_trunc == density_max)][0]
        width_guess                     = 2**(-0.5)*z_coords_trunc[numpy.where(densities_trunc >= density_max/numpy.e)[0][-1]]
        fit_coefficients_guess          = [amplitude_guess, mean_guess, width_guess]
#
# Passing guesses and data into the fitting thing
#
        print("Starting fit with " + str(nzsteps) + " steps......")
        fit_coefficients,fit_covariance = curve_fit(gauss,z_coords_trunc,densities_trunc,p0=fit_coefficients_guess)
        print("Done!")
        amplitude    = fit_coefficients[0]
        mean         = fit_coefficients[1]
        width        = fit_coefficients[2]
        fit          = gauss(z_coords_trunc,amplitude,mean,width)
#
# width = 2*sigma
# sigma**2 = 2*H**2
# so width = 2**1.5*H
#
        ivar_scaled = ivar - ivar_lower
        scale_height[ivar_scaled,ires] = 2.**(-1.5)*width
#
# Plot the fit
#
        ax = axes[ires]
        ax.plot(z_coords_trunc,densities_trunc,linestyle="none",marker="*",markersize=10.0,markerfacecolor="green",markeredgecolor="green")
        ax.plot(z_coords_trunc,numpy.repeat(density_max/numpy.e,len(z_coords_trunc)),linestyle="--",color="gray")
        ax.axvline(mean + scale_height[ivar_scaled,ires],linestyle="--",color="gray")
        ax.axvline(mean - scale_height[ivar_scaled,ires],linestyle="--",color="gray")
        ax.plot(z_coords_trunc,fit,color="black")
        ax.set_xlabel(r"$z$")
        ax.set_ylabel(r"N$_{p}$")
        ax.set_title(str(ivar) + " orbits, " + str(nzsteps) + " points")
    fits_file_title = simulation_name + "_" + script_name + "_fits_t-" + str(ivar) + ".png"
    fits_fig.savefig(fits_file_title,bbox_inches="tight")
    fits_fig.clear()
pyplot.close(fits_fig)
#
# closeness vs time
#

hvt_fig, ((ax1)) = pyplot.subplots(1,1)
hvt_fig.set_size_inches(19.2,10.8)
hvt_fig.set_dpi(100)

if(nvar == 1):
    ax1.set_xlim(ivar_lower-0.5,ivar_lower+0.5)
else:
    ax1.set_xlim(ivar_lower,ivar_upper)
markerstyle = "*"
markersize   = 30.0

ymin = 1.0e10
ymax = 0
time       = numpy.array(range(ivar_lower,ivar_upper+1))
colors     = cm.get_cmap('viridis')(numpy.linspace(0.0,1.0,nres))
for ires in range(nres):
    scale_height_array = scale_height[:,ires]
    npoints            = resolutions[ires]
    markercolor        = colors[ires]
    markerlabel        = str(npoints) + " points"
    if(scale_height_array.min() < ymin):
        ymin = scale_height_array.min()
    if(scale_height_array.max() > ymax):
        ymax = scale_height_array.max()
    ax1.plot(time,scale_height_array,linestyle="none",marker=markerstyle,markersize=markersize,markeredgecolor=markercolor,markerfacecolor="none",label=markerlabel)

ax1.set_xlim([ivar_lower,ivar_upper])
if(ymax < 0.0):
    ymax = ymax + 0.1*(ymax-ymin)
ymin = ymin - 0.1*(ymax-ymin)
ax1.set_ylim(ymin,ymax)
ax1.set_yscale("log")

ax1.set_xlabel(r"t (orbits)",fontsize=24)
ax1.set_ylabel(r"h/H",fontsize=24)
ax1.legend(fontsize=24)

for tick in ax1.xaxis.get_major_ticks():
    tick.label.set_fontsize(19)
for tick in ax1.yaxis.get_major_ticks():
    tick.label.set_fontsize(19)

hvt_file_title = simulation_name + "_" + script_name + "_" + str(ivar_lower) + "-" + str(ivar_upper) + "_scale-height-vs-time.png"
hvt_fig.savefig(hvt_file_title,bbox_inches="tight")

pyplot.show()
