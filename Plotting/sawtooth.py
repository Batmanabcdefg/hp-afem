from __future__ import division
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors
import sys
import numpy as np

show_legend = True
plot_subits = True
do_progression = False
progression_filter_oits = (False, [1, 5, 14, 16, 25, 33])
oits_start = 16
override_legend = (False, [
    r"$\mu = \frac{1}{2}, \omega = 2$", 
    r"$\mu = \frac{1}{2}, \omega = 4$", 
    r"$\mu = \frac{1}{2}, \omega = 8$",
    r"$\mu = \frac{1}{4}, \omega = 4$", 
    r"$\mu = \sqrt{\frac{1}{2}}, \omega = 4$", 
    ])

def trinum(deg):
    return (deg+2)*(deg+1)/2

def degree(dim):
    p = 0
    while trinum(p) <= dim:
        p = p+1
    return p-1

def parse_progression(obj):
    output = []
    for n in range(len(obj['oits'])):
        if obj['iits'][n] == 0:
            if progression_filter_oits[0] and obj['oits'][n] not in progression_filter_oits[1]:
                continue
            verts, linsol, vdofs, tris, dims, types, sol, tdofs, terror = read_meshfile(
                    "testLshaped_%d.%d_%d.sol" % (
                        obj['oits'][n], obj['iits'][n], obj['dofs3'][n]))
            oit = obj['oits'][n]            # k
            compl = sum(dims)               # \#\C(\D_k^*))
            ndof = obj['dofs3'][n]          # \#DoF(\C(\D_k^*))
            val = obj['vals'][n]            # \hsno{u - u_{\C\D_k^*}}
            broken = obj['brokens'][n]      # E_{\D_k^*}(u_{k-1})^.5
            pinf = obj['pinfs'][n]          # \norm{p_{\D_k^*}}_{\infty}
            preval = obj['prevals'][n]      # \hsno{u - u_{k-1}}
            nbcompl = obj['nbcompls'][n]    # \# \D_k^*
            nreduce = obj['riits'][n]       # number of iterations in reduce
            fracerror = 10 #max(terror)/min(terror) TODO
            fractime = (obj['times'][n] - obj['times'][n-1] + 0.0)/(obj['times'][n] - obj['times'][n-obj['iits'][n-1]])
            log_fraction = 1.0/(((1 + np.log(pinf))**(1.5))/((val - preval)/broken))
            comparable = (pinf + 0.0)/min(map(degree, dims))
            conf_ref_frac = float(sum([trinum(degree(dim) + 1)-1 for dim in dims]))/nbcompl #this is >= \#C(D)/\# D
            output.append([oit, pinf, log_fraction, comparable, conf_ref_frac, nreduce, fracerror, fractime, ndof])
    res = np.array(output).transpose()
    return {'oits': res[0], 'pinfs': res[1], 'log_frac': res[2], 'comp': res[3],
            'conf_ref_frac': res[4], 'nred': res[5], 'frac_error': res[6], 'frac_time': res[7], 'ndofs': res[8]}

def read_meshfile(fn):
    f = open(fn)
    settings = f.readline()[:-1].split(" ")

    verts, linsol, vdofs = readverts(f, settings)
    tris, dims, types, sol, tdofs, terror = readtris(f, settings, vdofs)
    f.close()
    return verts, linsol, vdofs, tris, dims, types, sol, tdofs, terror

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

def find_time_partition(obj):
    total = max(obj['times'])
    n = 0
    sumtimes = 0
    tm1 = 0
    for j in range(1, max(obj['oits'])+1):
        t0 = 0
        while n < len(obj['oits']) and obj['oits'][n] == j:
            if obj['iits'][n] == 0:
                t0 = obj['times'][n]
            if obj['riits'][n] == 0:
                sumtimes += obj['times'][n] - t0
                print "nearbest time", (t0-tm1)/(obj['times'][n]-tm1)
                tm1 = obj['times'][n]
            n = n + 1
    sys.exit(0)

def parse_file(lines):
    def find_in_lines(search):
        for line in lines:
            if line.startswith(search):
                return line.split(": ")[1][:-1]
        return False

    # first, parse the header
    meshfile = find_in_lines("mesh file")
    rhsfile = find_in_lines("rhs file")
    initdeg = int(find_in_lines("initial degree"))
    inittrinum = int(find_in_lines("number of triangles"))
    esthref = int(find_in_lines("h-refines"))
    estpref = int(find_in_lines("p-refines"))
    initerror = float(find_in_lines("initial error"))
    theta = float(find_in_lines("theta"))
    omega = float(find_in_lines("omega"))
    mu = float(find_in_lines("mu"))

    header_end = -1
    for i, line in enumerate(lines):
        if line.startswith("---------"):
            header_end = i+1
            break

    lines = lines[header_end:]

    # then, parse the rest of the lines
    oits = np.zeros(len(lines), dtype=np.int64)
    iits = np.zeros(len(lines), dtype=np.int64)
    dofs = np.zeros(len(lines), dtype=np.int64)
    vals = np.zeros(len(lines), dtype=np.float64)
    times = np.zeros(len(lines), dtype=np.int64)
    brokens = np.zeros(len(lines), dtype=np.float64)
    nbcompls = np.zeros(len(lines), dtype=np.int64)
    pinfs = np.zeros(len(lines), dtype=np.int64)
    prevals = np.zeros(len(lines), dtype=np.float64)

    i = 0
    for linestr in lines:
        line = linestr.split(" ")
        if int(line[0]) == 666:
            break
        if int(line[0]) == 3:
            nbcompls[i] = int(line[1][2:])
            brokens[i] = float(line[6])
        if int(line[0]) == 1:
            prevals[i] = float(line[1])
        if int(line[0]) == 2:
            pinfs[i] = int(line[1])
        if int(line[0]) == 0:
            if line[1][-1] == ':':
                times[i] = int(line[1][:-1])
                line = line[2:]
            if int(line[0]) < oits_start:
                continue
            oits[i] = int(line[0])
            iits[i] = int(line[1])
            dofs[i] = int(line[2])
            vals[i] = float(line[3])
            i = i + 1
        else:
            continue
    oits = oits[:i]
    iits = iits[:i]
    dofs = dofs[:i]
    vals = vals[:i]
    times = times[:i]
    times = ((times - times[0])/3600)**(1.0/3)
    brokens = brokens[:i]
    pinfs = pinfs[:i]
    nbcompls = nbcompls[:i]
    prevals = prevals[:i]
    epsilons = initerror * mu**oits
    broken_epsilons = (1 + omega)* initerror * mu**(oits - 1)

    riits = np.array(iits)
    n = 0
    for j in range(1,max(oits)+1):
        while n < len(oits) and oits[n] == j:
            riits[n] = max(iits[oits == j]) - iits[n]
            n = n + 1

    return {'params': {'theta': theta, 'omega': omega, 'mu': mu, 'epsilon0': initerror},
            'oits': oits, 'iits': iits, 'dofs': dofs**(1/3), 'dofs3': dofs,
            'vals': vals, 'times': times, 'brokens': brokens, 'riits': riits,
            'pinfs': pinfs, 'prevals': prevals, 'epsilons': epsilons,
            'nbcompls': nbcompls, 'broken_epsilons': broken_epsilons}

fps = []

if len(sys.argv) > 1:
    for fn in sys.argv[1:]:
        fps.append(open(fn))
else:
    fps.append(sys.stdin)

lines = [[] for fp in fps]
objs = [0 for fp in fps]

for i, fp in enumerate(fps):
    for line in fp:
        lines[i].append(line)
    objs[i] = parse_file(lines[i])

cmap = cm.get_cmap("hot")
sm = cm.ScalarMappable(norm=colors.Normalize(vmin=0, vmax=1.2*len(fps)), cmap=cmap)
#colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
legends = {
    'oits': r"$hp$-AFEM iteration number",
    'iits': "Reduce iteration number",
    'dofs': r"#DoF$^{1/3}$",
    'dofs3': r"#DoF",
    'vals': "Estimated error",
    'times': "Time (hours)",
    'brokens': "Broken error",
    'epsilons': r"$\epsilon_k$",
    'riits': r"Reverse Reduce iteration number",
    'broken_epsilons': r"$(\omega + C_L) \epsilon_{k-1}$",
    'pinfs': r"$||p_D||_\infty$",
    'log_frac': r"$\chi(D)$", #higher is better
    'comp': r"Max nbr degree quotient", #lower is better
    'conf_ref_frac': r"$\# C(D)/\# D$", #lower is better
    'nred': r"# iterations in Reduce", #lower is better
    'frac_time': "Fraction of time inside NearBest"
}

default_params = {'theta': 0.8, 'omega': 4, 'mu': 0.5}
pretty_params = {'theta': r"$\theta$", 'omega': r"$\omega$", 'mu': r"$\mu$"}
different_params = []

def iterplot(i, obj, key, ax):
    if override_legend[0]:
        label = override_legend[1][i]
    else:
        label = ""
        for param, val in obj['params'].iteritems():
            if param in different_params:
                label = label + ", %s=%g" % (pretty_params[param], val)
        label = label[2:]

    prog = obj['progression']

    handle, = ax.plot(prog['oits'], prog[key], c=sm.to_rgba(i), linewidth=0.5, marker='.', markersize=1.5, label=label)#, marker='-')

    ax.set_xlabel(legends['oits'])
    ax.set_ylabel(legends[key])
    return handle

def doplot(i, obj, a, b, ax, plot_subit, reverse=False, marker='x'):
    if override_legend[0]:
        label = override_legend[1][i]
    else:
        label = ""
        for param, val in obj['params'].iteritems():
            if param in different_params:
                label = label + ", %s=%g" % (pretty_params[param], val)
        label = label[2:]

    check = 'riits' if reverse else 'iits'
    """
    xs = np.vstack([obj[a][obj[check] == 0], np.ones(sum(obj[check]==0))]).T
    ys = np.log10(obj[b][obj[check] == 0])
    print xs, ys
    sol, _, _, _ = np.linalg.lstsq(xs,ys)
    print sol
    ax.semilogy(obj[a][obj[check] == 0], 10**(sol[0] * obj[a][obj[check] == 0] + sol[1]))
    """
    if plot_subit:
        handle, = ax.semilogy(obj[a][obj[check] == 0], obj[b][obj[check] == 0],
                              c=sm.to_rgba(i), linestyle='None', marker=marker,
                              label=label)
        ax.semilogy(obj[a][obj[check] != 0], obj[b][obj[check] != 0], c=sm.to_rgba(i),
                    linestyle='None', marker='.', markersize=1.5)
        ax.semilogy(obj[a], obj[b], c=sm.to_rgba(i), linewidth=0.5)#, marker='-')
    else:
        handle, = ax.semilogy(obj[a][obj[check] == 0], obj[b][obj['iits'] == 0],
                              c=sm.to_rgba(i), linestyle='None', marker=marker,
                              label=label)
        ax.semilogy(obj[a][obj[check] == 0], obj[b][obj['iits'] == 0],
                    c=sm.to_rgba(i), linewidth=0.5)#, marker='--')

    ax.set_xlabel(legends[a])
    ax.set_ylabel(legends[b])
    return handle

handles = []

for param in default_params:
    for i, obj in enumerate(objs):
        if obj['params'][param] != default_params[param]:
            different_params.append(param)

xs = np.linspace(1000, 110000, 1000)
ys = 500 * xs**(-2)

if do_progression:
    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex=True)
    progression = [0 for fp in fps]

    for i, obj in reversed(list(enumerate(objs))):
        obj['progression'] = parse_progression(objs[i])
        print map(int, obj['progression']['ndofs'])
    for i, obj in reversed(list(enumerate(objs))):
        #handles.append(doplot(i, obj, 'dofs', 'vals', ax1, plot_subits))
        #doplot(i, obj, 'times', 'vals', ax2, plot_subits)
        handles.append(iterplot(i, obj, 'log_frac', ax1))
        iterplot(i, obj, 'log_frac', ax2)
        iterplot(i, obj, 'frac_time', ax3)
        iterplot(i, obj, 'comp', ax4)
    ax1.set_xlabel("")
    ax2.set_xlabel("")
else:
    f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
    for i, obj in reversed(list(enumerate(objs))):
        handles.append(doplot(i, obj, 'dofs', 'vals', ax1, plot_subits))
        doplot(i, obj, 'times', 'vals', ax2, plot_subits)
    ax2.set_ylabel("")
if show_legend:
    plt.legend(handles=handles)
#h = ax1.semilogy(xs**(1/3), ys, 'r--')
plt.tight_layout()
plt.show()
