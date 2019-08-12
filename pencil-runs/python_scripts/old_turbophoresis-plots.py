import numpy
import matplotlib.pyplot as pyplot
import sys
import pencil
import pickle

parameters = pencil.read_param()

snap = -1
var  = "pvar.dat"
if len(sys.argv) == 2:
    snap = int(sys.argv[1])
    var  = "PVAR" + sys.argv[1]
data  = pencil.read_var(ivar=snap,magic="vorticity")
pdata = pencil.read_pvar(varfile=var)

x   = data.x[3:-3]
y   = data.y[3:-3]
z   = data.z[3:-3]
ux  = data.ux[3:-3,3:-3,3:-3]
uy  = data.uy[3:-3,3:-3,3:-3]
uz  = data.uz[3:-3,3:-3,3:-3]
rho = data.rho[3:-3,3:-3,3:-3]

xp  = pdata.xp
yp  = pdata.yp
zp  = pdata.zp
vpx = pdata.vpx
vpy = pdata.vpy
vpz = pdata.vpz

contours_figure, ((ax1,ax2,ax3)) = pyplot.subplots(1,3)
plots_figure, ((ax4,ax5,ax6))    = pyplot.subplots(1,3)

iyslice = 0
markersize = 0.1

# rms velocity (mean subtracted?) contour
ux_ref  = numpy.zeros([len(x),len(y),len(z)])
uy_ref  = numpy.transpose(numpy.broadcast_to(-0.5*x,(len(x),len(y),len(z))))
#uy_ref  = numpy.zeros([len(x),len(y),len(z)])
#ux_shear_profile = -0.5*x
#for i in range(len(y)):
#    for j in range(len(z)):
#        uy[:,i,j] = ux_shear_profile
uz_ref  = numpy.zeros([len(x),len(y),len(z)])
sumthin = ((ux - ux_ref)**2 + (uy - uy_ref)**2 + (uz - uz_ref)**2)**(0.5)
ax1.pcolormesh(x,z,sumthin[:,iyslice,:])
ax1.set_title("gas velocity residuals at y=0")
ax1.set_xlabel("x")
ax1.set_ylabel("z")

# vorticity contour
#velocity_array = numpy.array([ux,uy,uz])
#vorticity_raw  = pencil.curl(velocity_array,data.dx,data.dy,data.dz)
#vorticity_x    = vorticity_raw[0][3:-3,3:-3,3:-3]
#vorticity_y    = vorticity_raw[1][3:-3,3:-3,3:-3]
#vorticity_z    = vorticity_raw[2][3:-3,3:-3,3:-3]
vorticity      = data.vort[:,3:-3,3:-3,3:-3]
vorticity_mag  = (vorticity[0]**2 + vorticity[1]**2 + vorticity[2]**2)**(0.5)
ax2.pcolormesh(x,z,vorticity_mag[:,iyslice,:])
ax2.set_title("magnitude of vorticity at y=0")
ax2.set_xlabel("x")
ax2.set_ylabel("z")

# dust to gas ratio contour
filename = "rhop_orbit-" + str(snap) + ".pickle"
try:
    rhop_file = open(filename,"r")
    rhop = pickle.load(rhop_file)
except IOError:
    rhop_file = open(filename,"w")
    rhop = pencil.particles_to_density(xp,yp,zp,data.x,data.y,data.z)[3:-3,3:-3,3:-3]
    pickle.dump(rhop,rhop_file)
    rhop_file.close()
epsilon_d = parameters.eps_dtog*rhop/rho
ax3.pcolormesh(x,z,epsilon_d[:,iyslice,:])
ax3.set_title("dust to gas ratio at y=0")
ax3.set_xlabel("x")
ax3.set_ylabel("z")

# typical particle velocity vs epsilon for each grid cell
filename = "vp_rms_orbit-" + str(snap) + ".pickle"
try:
    vp_rms_file = open(filename,"r")
    vp_rms = pickle.load(vp_rms_file)
    where_particles = numpy.where(vp_rms != 0.0)
except IOError:
    vp_rms_file = open(filename,"w")
    npar   = numpy.zeros([len(x),len(y),len(z)])
    vp_rms = numpy.zeros([len(x),len(y),len(z)])
    for ip in range(len(pdata.ipars)):
        xpar  = xp[ip]
        ypar  = yp[ip]
        zpar  = zp[ip]
        ix = numpy.where(numpy.abs(xpar - x) == numpy.abs(xpar - x).min())[0][0]
        iy = numpy.where(numpy.abs(ypar - y) == numpy.abs(ypar - y).min())[0][0]
        iz = numpy.where(numpy.abs(zpar - z) == numpy.abs(zpar - z).min())[0][0]
        vp_rms[ix,iy,iz] += (vpx[ip]**2 + vpy[ip]**2 + vpz[ip]**2)**(0.5)
        npar[ix,iy,iz] += 1
        prog = float(ip)/float(len(pdata.ipars))
        sys.stdout.write("\r" + str(prog))
        sys.stdout.flush()
    where_particles = numpy.where(npar != 0)
    vp_rms[where_particles] = vp_rms[where_particles]/npar[where_particles]
    pickle.dump(vp_rms,vp_rms_file)
    vp_rms_file.close()
vp_rms_pars    = vp_rms[where_particles]
epsilon_d_pars = epsilon_d[where_particles]
order = numpy.unravel_index(numpy.argsort(epsilon_d_pars),epsilon_d_pars.shape)
ax4.scatter(epsilon_d_pars[order],vp_rms_pars[order]**2,s=markersize,color="black")
ax4.set_title("rms particle velocity vs dust to gas ratio")
ax4.set_xlabel(r"$\epsilon_{d}$")
ax4.set_ylabel(r"$v_{p,rms}$")

# voticity vs epsilon for each grid cell

vorticity_mag_pars = vorticity_mag[where_particles]
ax5.scatter(epsilon_d_pars[order],vorticity_mag_pars[order],s=markersize,color="black")
ax5.set_title("magnitude of vorticity vs dust to gas ratio")
ax5.set_xlabel(r"$\epsilon_{d}$")
ax5.set_ylabel(r"|$\omega$|")

# ratio of typical velocites for particles and gas in each grid cell

umag = (data.ux[3:-3,3:-3,3:-3]**2 + data.uy[3:-3,3:-3,3:-3]**2 + data.uz[3:-3,3:-3,3:-3]**2)**(0.5)
velocity_ratio = vp_rms_pars/umag[where_particles]
ax6.scatter(epsilon_d_pars[order],velocity_ratio[order],s=markersize,color="black")
ax6.axhline(1.0)
ax6.set_title("(rms particle velocity)/(magnitude of gas velocity) vs dust to gas ratio")
ax6.set_xlabel(r"$\epsilon_{d}$")
ax6.set_ylabel(r"|$\omega$|")

pyplot.show()
