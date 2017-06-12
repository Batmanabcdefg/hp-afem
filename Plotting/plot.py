import sys
import time
from inspect import isfunction
import random
from sympy import sympify, lambdify, symbols
from math import ceil, floor
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import matplotlib as mpl
from matplotlib.collections import PatchCollection
import matplotlib.patches as mpatches
import numpy as np

only_do_error_quotient = True
draw_3d = False
draw_arrows = False
draw_dofs = False
draw_trisid = False
fix_max_degree = (False, 20)
zoomlocation = (1,1)

cmap_errors = 'magma'
cmap_dofs = 'YlGnBu'

def G( x1, y1, x2, y2, x3, y3, x, y):
    return ((x2-x1)*x + (x3-x1)*y + x1, (y2-y1)*x + (y3-y1)*y + y1)

def trinum(deg):
    return (deg+2)*(deg+1)/2

def readsol(fp):
    x, y = symbols('x y')
    solline = fp.readline()
    if len(solline):
        sol = lambdify((x,y), sympify(solline).expand())
    else:
        sol = []
    return sol

def genpoints( p, tri):
    fp = float(p)
    points = [()]*trinum(p)
    innie = 0.1
    # vertex DOFs
    for i in range(3):
        points[i] = tuple(tri[i].tolist() + ['green'])

    if p >= 2:
        for d in range(1, p+1):
            tn = trinum(d)
            l = d/fp
            if d >= 3:
                for r1 in range(d-2):
                    #find x, y on ref triangle
                    xr = (r1 + 1.0)/fp
                    yr = (fp - d + 1.0)/fp
                    x,y = G(tri[0][0], tri[0][1], tri[1][0], tri[1][1], tri[2][0], tri[2][1], xr, yr)
                    points[tn - (d-2) + r1] = tuple([x, y, 'red'])
            if d < p:
                for e in range(3):
                    x = (tri[e][0] * l + tri[(1 + e)%3][0] *(1-l))*(1 - innie/2) + innie/2 * tri[(e - 1)%3][0]
                    y = (tri[e][1] * l + tri[(1 + e)%3][1] *(1-l))*(1 - innie/2) + innie/2 * tri[(e - 1)%3][1]
                    points[tn + e] = tuple([x, y, 'blue'])
    return points

def degree(dim):
    p = 0
    while trinum(p) <= dim:
        p = p+1
    return p-1

def readbases(dirname):
    funcs = [[]]*8
    for i in range(8):
        fn = "%s/basis_%i.mat" % (dirname, i)
        f = open(fn)
        lines = f.readlines()
        f.close()
        numfuncs = int(lines[0])
        funcs[i] = [0]*numfuncs
        for idx, line in enumerate(lines[1:numfuncs+1]):
            func = line.split(" ")
            funcdeg = int(func[0])
            numcoeffs = (funcdeg+1)**2
            funcs[i][idx] = np.array(map( float, func[1:numcoeffs+1])).reshape((funcdeg+1, funcdeg+1))

    return funcs

def basisfunc(x, y, func):
    n, m = func.shape
    return sum(map(sum, [[x**a * y**b * func[a,b] for a in range(n)] for b in range(m)]))

def readverts(fp, settings):
    if len(settings) == 1:
        try:
            nverts = int(settings[0])
        except:
            nverts = int(fp.readline())
    else:
        nverts = int(fp.readline())
    sol = np.empty(nverts, dtype=np.double)
    verts = np.empty([nverts, 2], dtype=np.double)
    dofs = np.empty(nverts, dtype=np.int)
    for i in range(nverts):
        line = fp.readline().split(" ")
        px = line[0]
        py = line[1]
        if "linsol" in settings:
            sol[i] = float(line[2])
        elif "dof" in settings:
            dofs[i] = int(line[2])
        verts[i,:] = [float(px), float(py)]

    return verts, sol, dofs

def readtris(fp, settings, vdofs):
    dims = []
    sol = []
    types = []
    dofs = []
    terror = []
    ntris = int(fp.readline())
    if "tridim" in settings:
        dims = np.empty(ntris, dtype=np.int)
        if "sol" in settings:
            assert("tritype" in settings)
            sol = [[]]*ntris
        if "dof" in settings:
            dofs = [[]]*ntris
    if "error" in settings:
        terror = [0]*ntris
    if "tritype" in settings:
        types = np.empty(ntris, dtype=np.int)
    tris = np.empty([ntris, 3], dtype=np.int)
    for i in range(ntris):
        line = fp.readline().split(" ")
        if "tritype" in settings:
            types[i] = int(line[3 + int("tridim" in settings)])
        if "tridim" in settings:
            dim, v0, v1, v2 = line[0:4]
            dims[i] = int(dim)
            if "sol" in settings:
                sol[i] = [float(x) for x in line[5:]]
            if "dof" in settings:
                dofs[i] = [int(x) for x in line[5:]]
                assert(dofs[i][0] == vdofs[int(v0)] and 
                       dofs[i][1] == vdofs[int(v1)] and 
                       dofs[i][2] == vdofs[int(v2)])
        if "error" in settings:
            v0, v1, v2 = line[int("tridim" in settings):(3 + int("tridim" in settings))]
            terror[i] = float(line[3+int("tridim" in settings)])
        if not "error" in settings and not "tridim" in settings:
            v0, v1, v2 = line[0:3]
        tris[i,:] = [int(v0), int(v1), int(v2)]
    return tris, dims, types, sol, dofs, terror

def createarrows(verts, tris):
    x = verts[:,0]
    y = verts[:,1]

    patches = []
    for t in tris:
        midpoint_edge = ((x[t[1]] + x[t[2]])/2, (y[t[1]] + y[t[2]])/2)
        arrow = mpatches.Arrow( x[t[0]], y[t[0]], -0.5 * (x[t[0]] - midpoint_edge[0]), -0.5*(y[t[0]] - midpoint_edge[1]), width=0.1)
        patches.append( arrow)
    collection = PatchCollection( patches)

    return collection

def ploterror(verts, tris, terror, settings, dims, fig, ax, colorbar=True, zoom=1.0, location=(0,0)):
    print sum(dims)
    assert("error" in settings)
    x = verts[:,0]
    y = verts[:,1]
    ax.set_aspect('equal')
    ax.set_ylim((min(y)-0.2*max(abs(y)))/zoom + location[1], (max(y)+0.2*max(abs(y)))/zoom + location[1])
    ax.set_xlim((min(x)-0.2*max(abs(x)))/zoom + location[0], (max(x)+0.2*max(abs(x)))/zoom + location[0])

    img = ax.tripcolor(x, y, tris, facecolors=np.log10(np.array(np.sqrt(terror))),
        edgecolors='k', linewidth=0.1, cmap=cmap_errors)
    if colorbar:
        fig.colorbar(img, ax=ax)

def plotdofs(verts, vdofs, tris, tdofs, settings, dims, fig, ax):
    assert("tridim" in settings)
    x = verts[:,0]
    y = verts[:,1]
    ax.set_aspect('equal')
    ax.set_ylim(min(y)-0.2*max(abs(y)), max(y)+0.2*max(abs(y)))
    ax.set_xlim(min(x)-0.2*max(abs(x)), max(x)+0.2*max(abs(x)))

    img = ax.tripcolor(x, y, tris, facecolors=np.array(map(degree, dims)), edgecolors='b', cmap=cmap_dofs)
    #fig.colorbar(img, ax=ax)
    fig.colorbar(img, ax=ax, ticks=np.linspace(1, 16, 16, endpoint=True))

    if draw_arrows:
        ax.add_collection(arrows)

    pointtexts = {}
    for i, t in enumerate(tris):
        dofs = tdofs[i]
        pts = genpoints(degree(dims[i]), [verts[t[0]], verts[t[1]], verts[t[2]]])
        for j in range(dims[i]):
            pointtexts[pts[j]] = dofs[j]
    for key, val in pointtexts.iteritems():
        ax.text(key[0], key[1], "%i" % val, color=key[2], ha='right', va='center')

def plotmesh(verts, tris, settings, dims, fig, ax, colorbar=True, zoom=1.0, location=(0,0)):
    cmap = cmap_dofs
    x = verts[:,0]
    y = verts[:,1]
    ax.set_aspect('equal')
    ax.set_ylim((min(y)-0.2*max(abs(y)))/zoom + location[1], (max(y)+0.2*max(abs(y)))/zoom + location[1])
    ax.set_xlim((min(x)-0.2*max(abs(x)))/zoom + location[0], (max(x)+0.2*max(abs(x)))/zoom + location[0])
    if "tridim" in settings:
        if fix_max_degree[0]:
            maxx = fix_max_degree[1]
            bounds = np.linspace(1, maxx, maxx)
            norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

            img = ax.tripcolor(x, y, tris, norm=norm, facecolors=np.array(map(degree, dims)), edgecolors='k', linewidth=0.1, cmap=cmap)
            if colorbar:
                fig.colorbar(img, ax=ax, cmap=cmap, norm=norm, ticks=np.linspace(1, maxx, maxx))
        else:
            img = ax.tripcolor(x, y, tris, facecolors=np.array(map(degree, dims))**(1.0/1.0), edgecolors='k', linewidth=0.1, cmap=cmap)
            if colorbar:
                fig.colorbar(img, ax=ax, cmap=cmap)
    else:
        ax.triplot(x, y, tris)
    if draw_arrows:
        ax.add_collection(arrows)
    if draw_trisid:
        for i, t in enumerate(tris):
            xf,yf = G(x[t[0]], y[t[0]], x[t[1]], y[t[1]], x[t[2]], y[t[2]], x=1/3.0, y=1/3.0)
            ax.text(xf, yf, "%i" % i, color='green', size=20)

    if draw_dofs:
        for i, t in enumerate(tris):
            xf,yf = G(x[t[0]], y[t[0]], x[t[1]], y[t[1]], x[t[2]], y[t[2]], x=1/3.0, y=1/3.0)
            for j in range(3):
                ax.text(verts[t[j]][0], verts[t[j]][1], "%i" % t[j], color='red')
            ax.text(xf, yf, "%i" % i, color='green')

def plotlinsol(verts, linsol, tris, settings, dims, fig, ax):
    x = verts[:,0]
    y = verts[:,1]
    triv = tri.Triangulation(x, y, tris)
    ax.plot_trisurf(triv, linsol, cmap='viridis', edgecolor='none')
    ax.set_aspect('equal')

def plotsol( verts, tris, types, realsol, sol, settings, dims, fix, ax):
    def Ginv( vol, x1, y1, x2, y2, x3, y3, x, y):
        return (1/(2.0 * vol) * ((y3-y1)*(x-x1) + (x1-x3)*(y-y1)), 1/(2.0*vol)*((y1-y2)*(x-x1) + (x2-x1)*(y-y1)));

    def getvols( x, y, tris):
        vols = []
        for i, t in enumerate(tris):
            vols.append( 1/2.0 * abs(-x[t[1]]*y[t[0]] + x[t[2]]*y[t[0]] + x[t[0]]*y[t[1]] - x[t[2]]*y[t[1]] - x[t[0]]*y[t[2]] + x[t[1]]*y[t[2]]));
        return vols

    maxdeg = max([degree(d) for d in dims])
    f = readbases("../Bases2/degree%i" % maxdeg)
    x = verts[:,0]
    y = verts[:,1]
    v = getvols(x, y, tris)
    NUM = 1000
    xx = []
    yy = []
    zz = []
    isrealsol = False
    if isfunction(realsol):
        isrealsol = True

    for i, tri in enumerate( tris):
        for _ in range(int(ceil(NUM*v[i]))):
            X = random.random()
            Y = random.random()
            while X + Y > 1:
                X = random.random()
                Y = random.random()
            xf,yf = G(x[tri[0]], y[tri[0]], x[tri[1]], y[tri[1]], x[tri[2]], y[tri[2]], x=X, y=Y)
            xx.append( xf)
            yy.append( yf)
            zf = 0
            #realsol = (xf*yf*(1-xf-yf))**2
            #realsol = xf*yf*(1-xf)*(1-yf)
            #realsol = xf*yf*(1-xf*xf)*(1-yf*yf)
            for j in range( dims[i]):
                zf = zf + sol[i][j]*basisfunc(X, Y, f[types[i]][j])
            if isrealsol:
                zz.append( abs(zf - realsol(xf, yf)))
            else:
                zz.append(zf)

    if draw_3d:
        ax.scatter( xx, yy, zz, c=zz)
    else:
        ax.scatter( xx, yy, c=zz)
    ax.set_aspect('equal')

def doplot(fp):
    settings = fp.readline()[:-1].split(" ")
    print "read settings ", settings

    verts, linsol, vdofs = readverts(fp, settings)
    tris, dims, types, sol, tdofs, terror = readtris(fp, settings, vdofs)
    realsol = readsol(fp)

    if draw_arrows:
        arrows = createarrows(verts, tris)

    fig = plt.figure()
    if "linsol" in settings:
        ax1 = fig.add_subplot(121, projection='3d')
        ax2 = fig.add_subplot(122)
        plotmesh(verts, tris, settings, dims, fig, ax2)
        plotlinsol(verts, linsol, tris, settings, dims, fig, ax1)
    elif "sol" in settings:
        if draw_3d:
            ax1 = fig.add_subplot(121, projection='3d')
            ax2 = fig.add_subplot(122)
            plotsol(verts, tris, types, realsol, sol, settings, dims, fig, ax1)
        else:
            ax2 = plt.subplot2grid((3,4), (0,0), colspan=3, rowspan=3)
            ax2.set_title("Polynomial degrees in triangulation")
            ax1 = plt.subplot2grid((3,4), (0,3))
            ax3 = plt.subplot2grid((3,4), (1,3))
            #ax4 = plt.subplot2grid((3,4), (2,3))
            plotmesh(verts, tris, settings, dims, fig, ax1, colorbar=False, zoom=10**1, location=zoomlocation)
            plotmesh(verts, tris, settings, dims, fig, ax3, colorbar=False, zoom=10**2, location=zoomlocation)
            #plotmesh(verts, tris, settings, dims, fig, ax4, colorbar=False, zoom=4*10**5, location=zoomlocation)
            ax1.set_title(r"zoom $10^1$")
            ax3.set_title(r"zoom $10^2$")
            #ax4.set_title(r"zoom $10^6$")
            ax1.get_xaxis().set_ticks([])
            ax1.get_yaxis().set_ticks([])
            ax3.get_xaxis().set_ticks([])
            ax3.get_yaxis().set_ticks([])
            #ax4.get_xaxis().set_ticks([])
            #ax4.get_yaxis().set_ticks([])
        plotmesh(verts, tris, settings, dims, fig, ax2)
    elif "dof" in settings:
        ax = fig.add_subplot(111)
        if draw_dofs:
            plotdofs(verts, vdofs, tris, tdofs, settings, dims, fig, ax)
    elif "error" in settings:
        print np.sqrt(max(terror)/min(terror))
        ax2 = plt.subplot2grid((3,4), (0,0), colspan=3, rowspan=3)
        ax2.set_title("10-log of estimated errors in triangulation")
        ax1 = plt.subplot2grid((3,4), (0,3))
        ax3 = plt.subplot2grid((3,4), (1,3))
        #ax4 = plt.subplot2grid((3,4), (2,3))
        ploterror(verts, tris, terror, settings, dims, fig, ax1, colorbar=False, zoom=10**1, location=zoomlocation)
        ploterror(verts, tris, terror, settings, dims, fig, ax3, colorbar=False, zoom=10**2, location=zoomlocation)
        #ploterror(verts, tris, terror, settings, dims, fig, ax4, colorbar=False, zoom=10**3, location=zoomlocation)
        ax1.set_title(r"zoom $10^1$")
        ax3.set_title(r"zoom $10^2$")
        #ax4.set_title(r"zoom $10^3$")
        ax1.get_xaxis().set_ticks([])
        ax1.get_yaxis().set_ticks([])
        ax3.get_xaxis().set_ticks([])
        ax3.get_yaxis().set_ticks([])
        #ax4.get_xaxis().set_ticks([])
        #ax4.get_yaxis().set_ticks([])
        ploterror(verts, tris, terror, settings, dims, fig, ax2)
    else:
        ax = fig.add_subplot(111)
        plotmesh(verts, tris, settings, dims, fig, ax)

    return fig
    #fig.savefig('lshaped_%d.pdf' % sum(dims))

if __name__ == '__main__':
    fig = doplot(sys.stdin)
    plt.tight_layout()
    plt.show()
