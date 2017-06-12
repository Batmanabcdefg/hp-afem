from __future__ import division
import random
import time
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors
import sys
import os
import numpy as np
from multiprocessing import Pool
from subprocess import call

oits_start = 0
settings_start = 0

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
    times = ((times - times[0])/3600)
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

    return {'meshfile': meshfile, 'rhsfile': rhsfile,
            'params': {'theta': theta, 'omega': omega, 'mu': mu, 'epsilon0': initerror},
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

estimator_settings = [
        {'a': 1},
]

all_to_run = []

def parameters_to_thing(runline):
    print runline[1][0], runline[1][1], runline[0]
    return "-m %s -r %s -i %s -a %d" % (runline[1][0], runline[1][1], "output/testLshaped_%d.%d_%d.sol" % (runline[1][2], runline[1][3], runline[1][4]), runline[0]['a'])

for obj in objs:
    for setting in estimator_settings[settings_start:]:
        n = 0
        for j in range(1, max(obj['oits'])+1):
            while n < len(obj['oits']) and obj['oits'][n] == j:
                if obj['iits'][n] == 0:
                    all_to_run.append([setting, [obj['meshfile'], obj['rhsfile'], obj['oits'][n], obj['iits'][n], obj['dofs3'][n]]])
                n = n + 1

def run_process(innn):
    FNULL = open(os.devnull, 'w')
    call(("./oneshot " + innn).split(), stdout=FNULL)

commands_to_run = map(parameters_to_thing, all_to_run)
pool = Pool(processes=1)              # start 4 worker processes
pool.map(run_process, commands_to_run)
