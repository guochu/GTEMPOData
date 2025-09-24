
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

rc('text', usetex=True)
plt.rc('font',family='Times New Roman')


def get_analytic(beta):
    N = int(beta / 0.05)
    data_path = "toulouse2/thouless_analytic_beta%s_mu1.0_N%s.json"%(beta,N)
    with open(data_path, 'r') as f:
        data = f.read()
        data = json.loads(data)
    return asarray(data["gf"])


def get_gtempo(beta, chi):
    N = int(beta / 0.05)
    data_path = "toulouse2/thouless_tempo_beta%s_mu1.0_N%s_chi%s.json"%(beta, N, chi)
    with open(data_path, 'r') as f:
        data = f.read()
        data = json.loads(data)
    return asarray(data["gf"])

def get_error(gf1, gf2):
    return sqrt(norm(gf1 - gf2)**2 / len(gf1))



plt.rcParams.update({'font.size': 14})
# Greater and Lesser, it and ed

betas = (10, 50)
chis = (20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300)

fontsize=20
fontsize2=14
fontsize3=22
fig, ax = plt.subplots(1,2, figsize=(8,3))

n1, n2 = len(betas), len(chis)
errs = numpy.zeros((n1, n2))
for i in range(n1):
    for j in range(n2):
        try:
            gf1 = get_analytic(betas[i])
            gf2 = get_gtempo(betas[i], chis[j])
            errs[i, j] = get_error(gf1, gf2)
        except:
            print()
print(errs)

xticks = ([0,4,8,12,16,20], [0, 20, 40, 60, 80, 100])
inyticks1 = [[0, 0.01, 0.02], [0, 0.01, 0.02, 0.03]]
inyticks2 = [[0, 0.003, 0.006], [0, 0.02, 0.04]]
ns = [11, 25]
ns_gf = [5, 19]
titles = (r"$\beta=10$", r"$\beta=50$")
for i in range(2):
    chi = chis[ns_gf[i]-1]
    print("chi of", chi)
    gf1 = get_analytic(betas[i])
    gf2 = get_gtempo(betas[i], chi)
    taus = numpy.arange(0,len(gf1)) * 0.05 * 2
    ax[i].plot(taus, -gf1, color="darkgray", label="analytic")
    ax[i].plot(taus, -gf2, ls=(0, (5, 5)), color="k", label="GTEMPO")
    
    ax[i].set_xlabel(r"$D\tau$", labelpad=0.1, fontsize=fontsize)
    ax[i].set_ylabel(r"$G(\tau)$", labelpad=0.0, fontsize=fontsize)
    ax[i].set_title(titles[i], fontsize=fontsize)

    ax[i].minorticks_on()
    ax[i].tick_params("both", which='both', direction='in', labelsize=fontsize)
    ax[i].set_xlim((-0.1,betas[i]*2+0.1))
    tops = [0.0, 0.05]
    ax[i].set_ylim(top=tops[i])
    ax[i].set_xticks(xticks[i])

    left, bottom, width, height = 0.18,0.2,0.3,0.3
    ax_ins = ax[i].inset_axes([left,bottom,width,height])
    ax_ins.plot(chis[:ns[i]], errs[i][:ns[i]], color="dodgerblue")#, marker="o", markerfacecolor='none', markersize=4)
    ax_ins.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
    ax_ins.set_xlabel(r"$\chi$", labelpad=-3, fontsize=fontsize2)
    ax_ins.set_ylabel("Mean Error", labelpad=0.1, fontsize=fontsize2)
    ax_ins.minorticks_on()
    ax_ins.tick_params("both", which='both', direction='in', labelsize=fontsize2)
    ax_ins.set_yticks(inyticks1[i])

    left, bottom, width, height = 0.5,0.2,0.3,0.3
    ax_ins = ax[i].inset_axes([left,bottom,width,height])
    ax_ins.plot(taus, abs(gf1-gf2), color="dodgerblue")
    ax_ins.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
    ax_ins.set_xlabel(r"$D\tau$", labelpad=-3, fontsize=fontsize2)
    ax_ins.set_ylabel(r"Abs. Error", labelpad=2, fontsize=fontsize2)
    ax_ins.yaxis.tick_right()
    ax_ins.yaxis.set_label_position("right")
    ax_ins.minorticks_on()
    ax_ins.tick_params("both", which='both', direction='in', labelsize=fontsize2)
    ax_ins.set_yticks(inyticks2[i])

    ax_ins.set_xlim((-0.1,betas[i]*2+0.1))
    ax_ins.set_ylim(bottom=0)

trans = mtransforms.ScaledTranslation(-20/72, 7/72, fig.dpi_scale_trans)
xpoint, ypoint = -0.2, 0.85
labels = ["(a)", "(b)"]
for i in range(2):
    ax[i].text(xpoint, ypoint, labels[i], transform=ax[i].transAxes + trans, fontsize=fontsize3, va='bottom')

ax[0].legend(loc=(0.2,0.6), fontsize=fontsize2)

plt.tight_layout(h_pad=0.2, w_pad=-0.5, pad=0.0)
plt.savefig("pics/toulouse2.pdf")
plt.show()
