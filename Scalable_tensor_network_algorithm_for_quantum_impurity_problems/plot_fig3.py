
import json
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import rc
from matplotlib.colors import LogNorm
from numpy import linspace, asarray, loadtxt, linspace, sqrt
from numpy.linalg import norm
import numpy
import matplotlib.transforms as mtransforms

from triqs.gf import *
from h5 import *
import os

from matplotlib.ticker import MultipleLocator


rc('text', usetex=True)
plt.rc('font',family='Times New Roman')



beta_chi = {10:60, 50:200}
def read_gtempo(beta, chi2, U):
	chi = beta_chi[beta]
	filename = "orb1/beta%s_δτ0.05_chi%s_chi2%s_U%sGtau.dat"%(beta, chi, chi2, U)
	if not os.path.exists(filename):
		return None
	print("Read from", filename)
	gtau = loadtxt(filename)
	return numpy.array(gtau)


def read_triqs(beta, U):
	filename = "cthybn/bands1/beta%s-U%.2f.h5"%(beta,U)
	with HDFArchive(filename, 'r') as A:
		Gtau_qmc = A['G_tau-0']
	n = int(500 / beta)
	return Gtau_qmc['up'].data[:,0,0][::n]

def get_error(gf1, gf2):
    return sqrt(norm(gf1 - gf2)**2 / len(gf1))



plt.rcParams.update({'font.size': 14})
# Greater and Lesser, it and ed

fontsize=20
fontsize2=14
fontsize3=22

Us = (1, 10)
betas = (10, 50)
chi2s = (6, 12, 18, 24, 30, 36, 42, 48, 54, 60)
fig, ax = plt.subplots(2,2, figsize=(8,6))

errs_Ubeta = []
for i in range(2):
	for j in range(2):
		U, beta = Us[i], betas[j]

		g1 = read_triqs(float(beta), float(U))
		errs = []
		for k in range(len(chi2s)):
			chi2 = chi2s[k]
			g2 = read_gtempo(beta, chi2, U)
			errs.append(get_error(g1, g2))
		errs_Ubeta.append(errs)

titles = [[r"$\beta=10,U=1$", r"$\beta=50,U=1$"], [r"$\beta=10,U=10$", r"$\beta=50,U=10$"]]
xticks = ([[0,4,8,12,16,20], [0, 20, 40, 60, 80, 100]], [[0,4,8,12,16,20], [0, 20, 40, 60, 80, 100]])
for i in range(2):
	for j in range(2):
		U, beta = Us[i], betas[j]

		g1 = read_triqs(float(beta), float(U))
		g2 = read_gtempo(beta, chi2s[-1], U)
		taus = numpy.arange(0,len(g1)) * 0.05 * 2
		ax[i,j].plot(taus, g1, color="darkgray", label="CTQMC")
		ax[i,j].plot(taus, g2, ls=(0, (5, 5)), color="k", label="GTEMPO")

		ax[i,j].set_xlabel(r"$D\tau$", labelpad=0.1, fontsize=fontsize)
		ax[i,j].set_ylabel(r"$G(\tau)$", labelpad=0.1, fontsize=fontsize)
		ax[i,j].set_title(titles[i][j], fontsize=fontsize)

		ax[i,j].minorticks_on()
		ax[i,j].tick_params("both", which='both', direction='in', labelsize=fontsize)
		ax[i,j].set_xlim((-0.1,beta*2+0.1))
		ax[i,j].set_ylim(top=0)
		ax[i,j].set_xticks(xticks[i][j])

		errs = []
		for k in range(len(chi2s)):
			chi2 = chi2s[k]
			g2 = read_gtempo(beta, chi2, U)
			errs.append(get_error(g1, g2))

		left, bottom, width, height = 0.23,0.2,0.3,0.3
		ax_ins = ax[i,j].inset_axes([left,bottom,width,height])
		ax_ins.plot(chi2s, errs, color="dodgerblue")
		ax_ins.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
		ax_ins.set_xlabel(r"$\chi_2$", labelpad=-2, fontsize=fontsize2)
		ax_ins.set_ylabel("Mean Error", labelpad=0.1, fontsize=fontsize2)
		ax_ins.minorticks_on()
		ax_ins.tick_params("both", which='both', direction='in', labelsize=fontsize2)
		if i==0 and j==0:
			ax_ins.set_yticks([0.0035, 0.0033, 0.0031])


		left, bottom, width, height = 0.55,0.2,0.3,0.3
		ax_ins = ax[i,j].inset_axes([left,bottom,width,height])
		if i==1 and j==0:
			ax_ins.plot(taus, abs(g1-g2), color="dodgerblue", lw=0.5)
		else:
			ax_ins.plot(taus, abs(g1-g2), color="dodgerblue")
		ax_ins.set_xlim((-0.1,beta*2+0.1))
		ax_ins.set_ylim(bottom=0)
		ax_ins.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
		ax_ins.set_xlabel(r"$D\tau$", labelpad=-2, fontsize=fontsize2)
		ax_ins.set_ylabel(r"Abs. Error", labelpad=2, fontsize=fontsize2)
		ax_ins.minorticks_on()
		ax_ins.tick_params("both", which='both', direction='in', labelsize=fontsize2)

		ax_ins.yaxis.tick_right()
		ax_ins.yaxis.set_label_position("right")

trans = mtransforms.ScaledTranslation(-20/72, 7/72, fig.dpi_scale_trans)
xpoint, ypoint = -0.2, 0.85
labels = [["(a)", "(b)"], ["(c)", "(d)"]]
for i in range(2):
	for j in range(2):
		ax[i,j].text(xpoint, ypoint, labels[i][j], transform=ax[i,j].transAxes + trans, fontsize=fontsize3, va='bottom')

ax[0,1].legend(loc=(0.3,0.6), fontsize=fontsize2)

plt.tight_layout(h_pad=0.5, w_pad=0, pad=0.2)
plt.savefig("pics/SIAM.pdf")
plt.show()


