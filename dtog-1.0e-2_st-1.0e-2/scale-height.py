import pencil
import numpy
from matplotlib import pyplot
import sys
from scipy.optimize import curve_fit

def gauss(x, amplitude, mean, sigma):
    return amplitude*numpy.exp(-((x - mean)**2)/(2.*sigma**2))

var = "pvar.dat"
data = pencil.read_var()
t = data.t/(2*numpy.pi)
if len(sys.argv) > 1:
    var = sys.argv[1]
    if(var != "pvar.dat"):
        t    = int(var[4:])
        data = pencil.read_var(ivar=t,trimall=True)

parameters = pencil.read_param()
pdata      = pencil.read_pvar(varfile=var)

xp = pdata.xp
yp = pdata.yp
zp = pdata.zp

x0 = parameters.xyz0[2] ; y0 = parameters.xyz0[1] ; z0 = parameters.xyz0[0]
x1 = parameters.xyz1[2] ; y1 = parameters.xyz1[1] ; z1 = parameters.xyz1[0]
nx = len(data.x)        ; ny = len(data.y)        ; nz = len(data.z)

nzsteps = 100
if len(sys.argv) == 3:
    nzsteps = int(sys.argv[2])
    if(nzsteps < 2):
        print("Please choose a number of steps > 1.")
        sys.exit("Exiting...")
z_coords,dz = numpy.linspace(start=z0,stop=z1,num=nzsteps,endpoint=False,retstep=True)
print("Dividing z = {" + str(z0) + "," + str(z1) + "} into " + str(nzsteps) + " steps, giving dz = " + str(dz))

densities                 = numpy.zeros(nzsteps)
for it in range(nzsteps):
    z_loop        = z_coords[it]
    densities[it] = len(numpy.intersect1d(numpy.where(zp >= z_loop)[0],numpy.where(zp < z_loop + dz)[0]))

density_threshold = 0.1
density_max       = densities.max()
above_threshold   = numpy.where(densities > density_threshold*density_max)
densities_trunc   = densities[above_threshold]
z_coords_trunc    = z_coords[above_threshold] + dz/2.

amplitude_guess                 = density_max
mean_guess                      = z_coords_trunc[numpy.where(densities_trunc == density_max)][0]
width_guess                     = 2**(-0.5)*z_coords_trunc[numpy.where(densities_trunc >= density_max/numpy.e)[0][-1]]
fit_coefficients_guess          = [amplitude_guess, mean_guess, width_guess]
print("Guessing Gaussian fit parameters: ")
print("  amplitude =  " + str(fit_coefficients_guess[0]))
print("  mean      =  " + str(fit_coefficients_guess[1]))
print("  width     =  " + str(fit_coefficients_guess[2]))

fit_coefficients,fit_covariance = curve_fit(gauss,z_coords_trunc,densities_trunc,p0=fit_coefficients_guess)
amplitude    = fit_coefficients[0]
mean         = fit_coefficients[1]
width        = fit_coefficients[2]
fit          = gauss(z_coords_trunc,amplitude,mean,width)
scale_height = 2.**(0.5)*width
print("Gaussian fit parameters: ")
print("  amplitude =  " + str(fit_coefficients[0]))
print("  mean      =  " + str(fit_coefficients[1]))
print("  width     =  " + str(fit_coefficients[2]))
print("In principle, H = 1.41*width = " + str(scale_height))

fig, ((ax1)) = pyplot.subplots(1,1)

ax1.plot(z_coords_trunc,densities_trunc,linestyle="none",marker="*",markerfacecolor="green",markeredgecolor="green")
ax1.plot(z_coords_trunc,numpy.repeat(density_max/numpy.e,len(z_coords_trunc)),linestyle="--",color="gray")
ax1.axvline(mean + scale_height,linestyle="--",color="gray")
ax1.axvline(mean - scale_height,linestyle="--",color="gray")
ax1.plot(z_coords_trunc,fit,color="black")

ax1.set_xlabel(r"$z$")
ax1.set_ylabel(r"$\rho_{p}$")
ax1.set_title("t = " + str(t) + " orbits")

pyplot.show()
