import pencil
import numpy
from matplotlib import pyplot
import sys

var = "pvar.dat"
data = pencil.read_var()
t = data.t/(2*numpy.pi)
if len(sys.argv) > 1:
    var = sys.argv[1]
    if(var != "pvar.dat"):
        t  = int(var[4:])

parameters = pencil.read_param()
pdata      = pencil.read_pvar(varfile=var)

xp = pdata.xp
yp = pdata.yp
zp = pdata.zp

x0 = parameters.xyz0[0]
y0 = parameters.xyz0[1]
z0 = parameters.xyz0[2]
x1 = parameters.xyz1[0]
y1 = parameters.xyz1[1]
z1 = parameters.xyz1[2]

y_lower_frac = 0.01
y_upper_frac = 0.01
y_center     = 0.00
if len(sys.argv) == 4:
    y_lower_frac = float(sys.argv[2])
    y_upper_frac = float(sys.argv[3])
    if(y_lower_frac > 1.0 or y_lower_frac < 0.0 or y_upper_frac > 1.0 or y_upper_frac < 0.0):
        print("Slice size arguments are fractions of the length from the slice center")
        print("to the domain edge, and should be between 0.0 and 1.0.")
        sys.exit("Exiting...")
elif len(sys.argv) == 5:
    y_lower_frac = float(sys.argv[2])
    y_upper_frac = float(sys.argv[3])
    y_center     = float(sys.argv[4])
    if(y_lower_frac > 1.0 or y_lower_frac < 0.0 or y_upper_frac > 1.0 or y_upper_frac < 0.0):
        print("Slice size arguments are fractions of the length from the slice center")
        print("to the domain edge, and should be between 0.0 and 1.0.")
        sys.exit("Exiting...")
    if(y_center >= y1 or y_center <= y0):
        print("Slice is centered at or outside of the domain. Please pick a slice")
        print("center between " + str(y0) + " and " + str(y1) + ".")
        sys.exit("Exiting...")

fig, ((ax1)) = pyplot.subplots(1,1)

yslice_lower = y_center - y_lower_frac*(y_center-y0)
yslice_upper = y_center + y_upper_frac*(y1-y_center)
yslice       = numpy.intersect1d(numpy.where(yp >= yslice_lower),numpy.where(yp <= yslice_upper))

markersize = .01

ax1.set_xlim([x0,x1])
ax1.set_ylim([z0,z1])

#ax1.scatter(xp[yslice],zp[yslice],s=markersize)
ax1.contourf(data.x,data.z,data.rho[:,0,:])

ax1.set_xlabel(r"$x$")
ax1.set_ylabel(r"$z$")
ax1.set_title("Particle positions, y={" + str(yslice_lower) + "," + str(yslice_upper) + "}, t = " + str(t) + " orbits")

pyplot.show()
