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
if(snap == -1):
    orbit = data.t/(2*numpy.pi)
else:
    orbit = snap

x   = data.x[3:-3]
y   = data.y[3:-3]
z   = data.z[3:-3]
ux  = numpy.transpose(data.ux[3:-3,3:-3,3:-3])
uy  = numpy.transpose(data.uy[3:-3,3:-3,3:-3])
uz  = numpy.transpose(data.uz[3:-3,3:-3,3:-3])
rho = numpy.transpose(data.rho[3:-3,3:-3,3:-3])

xp  = pdata.xp
yp  = pdata.yp
zp  = pdata.zp
vpx = pdata.vpx
vpy = pdata.vpy
vpz = pdata.vpz

m_d = parameters.mp_swarm*len(x)*len(y)*len(z)

epsilon_d_figure, ((ax1)) = pyplot.subplots(1,1)
turbulence_figure, ((ax2,ax3)) = pyplot.subplots(1,2)

iyslice = 0
markersize = 0.1

# dust to gas ratio contour
filename = "rhop_orbit-" + str(orbit) + ".pickle"
try:
    rhop_file = open(filename,"r")
    rhop = pickle.load(rhop_file)
    rhop_file.close()
except IOError:
    rhop_file = open(filename,"w")
    rhop = numpy.transpose(pencil.particles_to_density(xp,yp,zp,data.x,data.y,data.z)[3:-3,3:-3,3:-3])
    pickle.dump(rhop,rhop_file)
    rhop_file.close()
epsilon_d = m_d*rhop/rho
epsilon_d_contour = ax1.pcolormesh(x,z,numpy.transpose(epsilon_d[:,iyslice,:]))
ax1.set_title("dust to gas ratio at y=0")
ax1.set_xlabel("x")
ax1.set_ylabel("z")
epsilon_d_figure.colorbar(epsilon_d_contour,ax=ax1)

# typical particle velocity vs epsilon for each grid cell
vp_rms_file_found = True
npar_file_found = True
filename = "vp_rms_orbit-" + str(orbit) + ".pickle"
try:
    vp_rms_file = open(filename,"r")
    vp_rms = pickle.load(vp_rms_file)
    vp_rms_file.close()
except IOError:
    vp_rms_file = open(filename,"w")
    vpx_mean = numpy.zeros([len(x),len(y),len(z)])
    vpy_mean = numpy.zeros([len(x),len(y),len(z)])
    vpz_mean = numpy.zeros([len(x),len(y),len(z)])
    ngpx     = numpy.zeros([len(pdata.ipars)],dtype=numpy.int)
    ngpy     = numpy.zeros([len(pdata.ipars)],dtype=numpy.int)
    ngpz     = numpy.zeros([len(pdata.ipars)],dtype=numpy.int)
    vp_rms   = numpy.zeros([len(x),len(y),len(z)])
    vp_rms_file_found = False
filename = "npar_orbit-" + str(orbit) + ".pickle"
try:
    npar_file = open(filename,"r")
    npar = pickle.load(npar_file)
    npar_file.close()
except IOError:
    npar_file = open(filename,"w")
    npar      = numpy.zeros([len(x),len(y),len(z)])
    npar_file_found = False
if((not npar_file_found) and (vp_rms_file_found)):
    for ip in range(len(pdata.ipars)):
        xpar  = xp[ip]
        ypar  = yp[ip]
        zpar  = zp[ip]
        ix = numpy.where(numpy.abs(xpar - x) == numpy.abs(xpar - x).min())[0][0]
        iy = numpy.where(numpy.abs(ypar - y) == numpy.abs(ypar - y).min())[0][0]
        iz = numpy.where(numpy.abs(zpar - z) == numpy.abs(zpar - z).min())[0][0]
        npar[ix,iy,iz] += 1
    pickle.dump(npar,npar_file)
    npar_file.close()
elif((npar_file_found) and (not vp_rms_file_found)):
    for ip in range(len(pdata.ipars)):
        xpar  = xp[ip]
        ypar  = yp[ip]
        zpar  = zp[ip]
        ix = numpy.where(numpy.abs(xpar - x) == numpy.abs(xpar - x).min())[0][0]
        iy = numpy.where(numpy.abs(ypar - y) == numpy.abs(ypar - y).min())[0][0]
        iz = numpy.where(numpy.abs(zpar - z) == numpy.abs(zpar - z).min())[0][0]
        ngpx[ip] = ix
        ngpy[ip] = iy
        ngpz[ip] = iz
        vpx_mean[ix,iy,iz] += vpx[ip]
        vpy_mean[ix,iy,iz] += vpy[ip]
        vpz_mean[ix,iy,iz] += vpz[ip]
    where_particles = numpy.where(npar != 0)
    vpx_mean[where_particles] = vpx_mean[where_particles]/npar[where_particles]
    vpy_mean[where_particles] = vpy_mean[where_particles]/npar[where_particles]
    vpz_mean[where_particles] = vpz_mean[where_particles]/npar[where_particles]
    for ip in range(len(pdata.ipars)):
        ix  = ngpx[ip]
        iy  = ngpy[ip]
        iz  = ngpz[ip]
        vp_rms[ix,iy,iz] += (vpx[ip] - vpx_mean[ix,iy,iz])**2 + (vpy[ip] - vpy_mean[ix,iy,iz])**2 + (vpz[ip] - vpz_mean[ix,iy,iz])**2
    vp_rms[where_particles] = vp_rms[where_particles]/npar[where_particles]
    pickle.dump(vp_rms,vp_rms_file)
    vp_rms_file.close()
elif((not npar_file_found) and (not vp_rms_file_found)):
    for ip in range(len(pdata.ipars)):
        xpar  = xp[ip]
        ypar  = yp[ip]
        zpar  = zp[ip]
        ix = numpy.where(numpy.abs(xpar - x) == numpy.abs(xpar - x).min())[0][0]
        iy = numpy.where(numpy.abs(ypar - y) == numpy.abs(ypar - y).min())[0][0]
        iz = numpy.where(numpy.abs(zpar - z) == numpy.abs(zpar - z).min())[0][0]
        ngpx[ip] = ix
        ngpy[ip] = iy
        ngpz[ip] = iz
        vpx_mean[ix,iy,iz] += vpx[ip]
        vpy_mean[ix,iy,iz] += vpy[ip]
        vpz_mean[ix,iy,iz] += vpz[ip]
        npar[ix,iy,iz] += 1
    where_particles = numpy.where(npar != 0)
    vpx_mean[where_particles] = vpx_mean[where_particles]/npar[where_particles]
    vpy_mean[where_particles] = vpy_mean[where_particles]/npar[where_particles]
    vpz_mean[where_particles] = vpz_mean[where_particles]/npar[where_particles]
    for ip in range(len(pdata.ipars)):
        ix  = ngpx[ip]
        iy  = ngpy[ip]
        iz  = ngpz[ip]
        vp_rms[ix,iy,iz] += (vpx[ip] - vpx_mean[ix,iy,iz])**2 + (vpy[ip] - vpy_mean[ix,iy,iz])**2 + (vpz[ip] - vpz_mean[ix,iy,iz])**2
    vp_rms[where_particles] = vp_rms[where_particles]/npar[where_particles]
    pickle.dump(npar,npar_file)
    pickle.dump(vp_rms,vp_rms_file)
    npar_file.close()
    vp_rms_file.close()
num_particles_cut = numpy.where((npar >= 3)) #& (vp_rms > 1e-15) & (epsilon_d > 1e-15))
vp_rms_pars    = vp_rms[num_particles_cut]
epsilon_d_pars = epsilon_d[num_particles_cut]
order = numpy.unravel_index(numpy.argsort(epsilon_d_pars),epsilon_d_pars.shape)
ax2.set_xscale("log")
ax2.set_yscale("log")
ax2.scatter(epsilon_d_pars[order],vp_rms_pars[order],s=markersize,color="black")
ax2.set_title("$v_{p,rms}$ vs $\epsilon_{d}$")
ax2.set_xlabel(r"$\epsilon_{d}$")
ax2.set_ylabel(r"$v_{p,rms}$")

# voticity vs epsilon for each grid cell
vorticity      = numpy.transpose(data.vort[:,3:-3,3:-3,3:-3])
vorticity_mag  = (vorticity[:,:,:,0]**2 + vorticity[:,:,:,1]**2 + vorticity[:,:,:,2]**2)**(0.5)
vorticity_mag_pars = vorticity_mag[num_particles_cut]
ax3.set_xscale("log")
ax3.set_yscale("log")
ax3.scatter(epsilon_d_pars[order],vorticity_mag_pars[order],s=markersize,color="black")
ax3.set_title("$|\omega|$ vs $\epsilon_{d}$")
ax3.set_xlabel(r"$\epsilon_{d}$")
ax3.set_ylabel(r"|$\omega$|")

pyplot.show()
