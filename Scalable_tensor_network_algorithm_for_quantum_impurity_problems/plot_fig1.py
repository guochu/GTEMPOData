
import json
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import rc
from matplotlib.colors import LogNorm
from numpy import linspace, asarray, loadtxt, linspace, sqrt
from numpy.linalg import norm
import math
import numpy
import matplotlib.transforms as mtransforms

import matplotlib.ticker as ticker 
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import ScalarFormatter

rc('text', usetex=True)
plt.rc('font',family='Times New Roman')


def get_analytic():
    data_path = "toulouse1/thouless_analytic_beta10_mu1.0_N2000.json"
    with open(data_path, 'r') as f:
        data = f.read()
        data = json.loads(data)
    return asarray(data["gf"])
gf_ana = get_analytic()

def get_gtempo(dtau, chi):
    N = int(10 / dtau)
    data_path = "toulouse1/thouless_tempo_beta10_mu1.0_N%s_chi%s.json"%(N, chi)
    with open(data_path, 'r') as f:
        data = f.read()
        data = json.loads(data)
    return asarray(data["gf"])

def get_error(gf):
    n = int((len(gf_ana)-1) / (len(gf)-1))
    return sqrt(norm(gf_ana[::n] - gf)**2 / len(gf))

# dtaus = (0.4, 0.2, 0.1, 0.05, 0.01, 0.005)
# chis = (20, 30, 40, 50, 60, 70, 80, 90)
dtaus = (0.2, 0.1, 0.05, 0.01, 0.005)
Ddtaus =  (0.4, 0.2, 0.1, 0.02, 0.01)
chis = (40, 50, 60, 70, 80, 90)
# dtaus = (0.4, 0.2, 0.1, 0.05, 0.01)
# chis = (40, 50, 60, 70, 80, 90)

n1, n2 = len(dtaus), len(chis)
errs = numpy.zeros((n1, n2))
try:
    for i in range(n1):
        for j in range(n2):
            gf = get_gtempo(dtaus[i], chis[j])
            errs[i, j] = get_error(gf)
except:
    print()

print(errs)

plt.rcParams.update({'font.size': 14})
# Greater and Lesser, it and ed


fig, ax = plt.subplots(figsize=(6,5))
# im = ax.imshow(errs, cmap='viridis_r')
im = ax.imshow(errs, cmap='Blues')
ax.set_xticks(range(n2), labels=chis, fontsize=16)
ax.set_yticks(range(n1), labels=Ddtaus, fontsize=16)
ax.set_xlabel(r"$\chi$", fontsize=16)
ax.set_ylabel(r"$D \delta\tau$", fontsize=16)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1)  # size 是宽度, pad 是间距
cbar = plt.colorbar(im, cax=cax)

# cbar = fig.colorbar(im, ax=ax, shrink=0.8)
cbar.ax.tick_params(labelsize=14)



cbar.formatter.set_powerlimits((0, 0))

plt.tight_layout()

plt.savefig("pics/toulouse1.pdf")
plt.show()
