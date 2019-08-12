import pencil
import numpy
from matplotlib import pyplot
import sys

var          = "pvar.dat"
if len(sys.argv) > 1:
    varnum_str = sys.argv[1]
    var = "PVAR" + varnum_str

parameters = pencil.read_param()
pdata = pencil.read_pvar(varfile=var)

t  = int(varnum_str)
xp  = pdata.xp
yp  = pdata.yp
zp  = pdata.zp
vpx = pdata.vpx
vpy = pdata.vpy
vpz = pdata.vpz

x0 = parameters.xyz0[2]
y0 = parameters.xyz0[1]
z0 = parameters.xyz0[0]
x1 = parameters.xyz1[2]
y1 = parameters.xyz1[1]
z1 = parameters.xyz1[0]

x_lower_frac = 1.0
x_upper_frac = 1.0
x_center     = 0.0
y_lower_frac = 1.0
y_upper_frac = 1.0
y_center     = 0.0
z_lower_frac = 1.0
z_upper_frac = 1.0
z_center     = 0.0
if len(sys.argv) == 4:
    z_lower_frac = float(sys.argv[2])
    z_upper_frac = float(sys.argv[3])
    if(z_lower_frac > 1.0 or z_lower_frac < 0.0 or z_upper_frac > 1.0 or z_upper_frac < 0.0):
        print("Box size arguments are fractions of the length from the box center")
        print("to the domain edge, and should be between 0.0 and 1.0.")
        sys.exit("Exiting...")
elif len(sys.argv) == 5:
    z_lower_frac = float(sys.argv[2])
    z_upper_frac = float(sys.argv[3])
    z_center     = float(sys.argv[4])
    if(z_lower_frac > 1.0 or z_lower_frac < 0.0 or z_upper_frac > 1.0 or z_upper_frac < 0.0):
        print("Box size arguments are fractions of the length from the box center")
        print("to the domain edge, and should be between 0.0 and 1.0.")
        sys.exit("Exiting...")
    if(z_center >= z1 or z_center <= z0):
        print("Box is centered at or outside of the domain. Please pick a box")
        print("center between " + str(y0) + " and " + str(y1) + ".")
        sys.exit("Exiting...")

xbox_lower = x_center - x_lower_frac*(x_center-x0)
xbox_upper = x_center + x_upper_frac*(x1-x_center)
xbox       = numpy.intersect1d(numpy.where(xp >= xbox_lower),numpy.where(xp <= xbox_upper))
ybox_lower = y_center - y_lower_frac*(y_center-y0)
ybox_upper = y_center + y_upper_frac*(y1-y_center)
ybox       = numpy.intersect1d(numpy.where(yp >= ybox_lower),numpy.where(yp <= ybox_upper))
zbox_lower = z_center - z_lower_frac*(z_center-z0)
zbox_upper = z_center + z_upper_frac*(z1-z_center)
zbox       = numpy.intersect1d(numpy.where(zp >= zbox_lower),numpy.where(zp <= zbox_upper))
xybox     = numpy.intersect1d(xbox,ybox)
xyzbox     = numpy.intersect1d(xybox,zbox)

vpz_avg = numpy.mean(vpz[xyzbox])
print("Mean z velocity in x = {" + str(xbox_lower) + "," + str(xbox_upper) + "}, y = {" + str(ybox_lower) + "," + str(ybox_upper) + "}, z = {" + str(zbox_lower) + "," + str(zbox_upper) + "} is:")
print(vpz_avg)

#fig, ((ax1)) = pyplot.subplots(1,1)
#
#ax1.set_xlim([x0,x1])
#ax1.set_ylim([z0,z1])
#
#ax1.scatter(xp[yslice],zp[yslice],s=markersize)
#
#ax1.set_xlabel(r"$x$")
#ax1.set_ylabel(r"$z$")
#ax1.set_title("Particle positions, y={" + str(yslice_lower) + "," + str(yslice_upper) + "}, t = " + str(t) + " orbits")
#
#pyplot.show()
