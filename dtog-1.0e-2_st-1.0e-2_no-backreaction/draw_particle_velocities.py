import pencil
import numpy
import matplotlib
from matplotlib import pyplot
from matplotlib import colors
from matplotlib import cm
import sys

if(len(sys.argv) == 1):
    pdata = pencil.read_pvar(varfile="pvar.dat")
elif(len(sys.argv) == 2):
    var_file = "PVAR" + sys.argv[1]
    pdata    = pencil.read_pvar(varfile=var_file)
else:
    print("Too many args brah")
    sys.exit()

xp  = pdata.xp  ; yp  = pdata.yp
vpx = pdata.vpx ; vpy = pdata.vpy

vmag = (vpx**2 + vpy**2)**(0.5)

#fig, ((ax0)) = pyplot.subplots(1,1)
fig = pyplot.figure()
ax0 = fig.add_axes([0.1,0.1,0.8,0.8])

ax0.scatter(xp,yp,c="black",edgecolor="black")

scale_factor = 0.02

lenx = scale_factor*vpx
leny = scale_factor*vpy

colormap           = pyplot.cm.viridis
colornormalization = colors.Normalize(vmin=numpy.min(vmag),vmax=numpy.max(vmag))

scalarmap = cm.ScalarMappable(norm=colornormalization,cmap=colormap)

for i in range(len(xp)):
    xi = xp[i]   ; yi = yp[i]
    Lx = lenx[i] ; Ly = leny[i]
    icolor = scalarmap.to_rgba(vmag[i])
    ax0.arrow(xi,yi,Lx,Ly,width=0.0001,head_width=0.0002,color=icolor)

cbar_ax = fig.add_axes([0.91,0.1,0.05,0.8])
matplotlib.colorbar.ColorbarBase(cbar_ax,cmap=colormap,norm=colornormalization,orientation="vertical")

pyplot.show()
