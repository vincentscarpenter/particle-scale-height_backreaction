import numpy
import matplotlib.pyplot as pyplot
import sys
import pencil
import pickle

pi = numpy.pi

parameters  = pencil.read_param()
time_series = pencil.read_ts()

time    = time_series.t/(2*pi)
rhopmax = time_series.rhopmax

figure, ((ax1)) = pyplot.subplots(1,1)
ax1.plot(time,rhopmax,linewidth=0.1,color="black")
ax1.set_title("time evolution of particle density, dust to gas ratio, vorticity, particle velocity")
ax1.set_xlabel(r"t")
ax1.set_ylabel(r"rhopmax")

ax2 = ax1.twinx()
#ax2.set_ylim([0.0,10.0])
markersize = 10.0

snap_start = 0
snap_stop  = numpy.floor(time[-1]/20*pi)
if len(sys.argv) == 3:
    snap_start = int(sys.argv[1])
    snap_stop  = int(sys.argv[2])

for snap in range(snap_start,snap_stop+1):
    data  = pencil.read_var(ivar=snap,magic="vorticity")
    x     = data.x[3:-3]
    y     = data.y[3:-3]
    z     = data.z[3:-3]
    ux    = data.ux[3:-3,3:-3,3:-3]
    uy    = data.uy[3:-3,3:-3,3:-3]
    uz    = data.uz[3:-3,3:-3,3:-3]
    rho   = data.rho[3:-3,3:-3,3:-3]
    m_d = parameters.mp_swarm*len(x)*len(y)*len(z)
    compute_rhop   = False
    compute_vp_rms = False
    filename = "rhop_orbit-" + str(snap) + ".pickle"
    try:
        rhop_file = open(filename,"r")
        rhop = pickle.load(rhop_file)
        rhop_file.close()
    except IOError:
        rhop_file = open(filename,"w")
        compute_rhop = True
    filename = "vp_rms_orbit-" + str(snap) + ".pickle"
    try:
        vp_rms_file = open(filename,"r")
        vp_rms = pickle.load(vp_rms_file)
        vp_rms_file.close()
        where_particles = numpy.where(vp_rms != 0.0)
    except IOError:
        vp_rms_file = open(filename,"w")
        compute_vp_rms = True
    if(compute_rhop or compute_vp_rms):
        var   = "PVAR" + str(snap)
        pdata = pencil.read_pvar(varfile=var)
        xp  = pdata.xp
        yp  = pdata.yp
        zp  = pdata.zp
        vpx = pdata.vpx
        vpy = pdata.vpy
        vpz = pdata.vpz
    if(compute_rhop):
        rhop = pencil.particles_to_density(xp,yp,zp,data.x,data.y,data.z)[3:-3,3:-3,3:-3]
        pickle.dump(rhop,rhop_file)
        rhop_file.close()
    if(compute_vp_rms):
        npar   = numpy.zeros([len(x),len(y),len(z)])
        vpx_mean = numpy.zeros([len(x),len(y),len(z)])
        vpy_mean = numpy.zeros([len(x),len(y),len(z)])
        vpz_mean = numpy.zeros([len(x),len(y),len(z)])
        ngpx     = numpy.zeros([len(pdata.ipars)])
        ngpy     = numpy.zeros([len(pdata.ipars)])
        ngpz     = numpy.zeros([len(pdata.ipars)])
        vp_rms   = numpy.zeros([len(x),len(y),len(z)])
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
            npar[ix,iy,iz]     += 1
        where_particles = numpy.where(npar != 0)
        vpx_mean[where_particles] = vpx_mean[where_particles]/npar[where_particles]
        vpy_mean[where_particles] = vpy_mean[where_particles]/npar[where_particles]
        vpz_mean[where_particles] = vpz_mean[where_particles]/npar[where_particles]
        for ip in range(len(pdata.ipars)):
            ix  = int(ngpx[ip])
            iy  = int(ngpy[ip])
            iz  = int(ngpz[ip])
            vp_rms[ix,iy,iz] += (vpx[ip] - vpx_mean[ix,iy,iz])**2 + (vpy[ip] - vpy_mean[ix,iy,iz])**2 + (vpz[ip] - vpz_mean[ix,iy,iz])**2
        pickle.dump(vp_rms,vp_rms_file)
        vp_rms_file.close()
    where_rhopmax  = numpy.where(rhop == rhop.max())
    epsilon_d      = m_d*rhop/rho
#    vorticity      = data.vort[:,3:-3,3:-3,3:-3]
#    vorticity_mag  = (vorticity[0]**2 + vorticity[1]**2 + vorticity[2]**2)**(0.5)
#    ax2.plot(snap,epsilon_d[where_rhopmax],color="red",linestyle=None,markersize=markersize,marker="^",markerfacecolor="none")
    ax2.plot(snap,vp_rms[where_rhopmax],color="blue",linestyle=None,markersize=markersize,marker="^",markerfacecolor="none")
#    ax2.plot(snap,vorticity_mag[where_rhopmax],color="green",linestyle=None,markersize=markersize,marker="*")

pyplot.show()
