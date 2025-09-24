
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




import matplotlib.pyplot as plt
import numpy as np

import matplotlib.collections as mcol
from matplotlib.legend_handler import HandlerLineCollection, HandlerTuple
from matplotlib.lines import Line2D

class HandlerDashedLines(HandlerLineCollection):
    """
    Custom Handler for LineCollection instances.
    """
    def create_artists(self, legend, orig_handle,
                       xdescent, ydescent, width, height, fontsize, trans):
        # figure out how many lines there are
        numlines = len(orig_handle.get_segments())
        xdata, xdata_marker = self.get_xdata(legend, xdescent, ydescent,
                                             width, height, fontsize)
        leglines = []
        # divide the vertical space where the lines will go
        # into equal parts based on the number of lines
        ydata = np.full_like(xdata, height / (numlines + 1))
        # for each line, create the line at the proper location
        # and set the dash pattern
        for i in range(numlines):
            legline = Line2D(xdata, ydata * (numlines - i) - ydescent)
            self.update_prop(legline, orig_handle, legend)
            # set color, dash pattern, and linewidth to that
            # of the lines in linecollection
            try:
                color = orig_handle.get_colors()[i]
            except IndexError:
                color = orig_handle.get_colors()[0]
            try:
                dashes = orig_handle.get_dashes()[i]
            except IndexError:
                dashes = orig_handle.get_dashes()[0]
            try:
                lw = orig_handle.get_linewidths()[i]
            except IndexError:
                lw = orig_handle.get_linewidths()[0]
            if dashes[1] is not None:
                legline.set_dashes(dashes[1])
            legline.set_color(color)
            legline.set_transform(trans)
            legline.set_linewidth(lw)
            leglines.append(legline)
        return leglines



def get_error(gf1, gf2):
    return sqrt(norm(gf1 - gf2)**2 / len(gf1))


def read_gtempo_orb2(chi2, idx):
	path = "orb2/norb2_beta10_δτ0.05_chi60_chi2%s/"%chi2
	filename = path+"Gtau-%s.dat"%idx
	if not os.path.exists(filename):
		return None
	print("Load from:", filename)	
	gtau = loadtxt(filename)
	return numpy.array(gtau)[::500]

def read_bethe_triqs_orb2(i):
	filename = "cthybn/bands2/beta10.0-U4.00-J1.00.h5"
	with HDFArchive(filename, 'r') as A:
		Gtau_qmc = A['Gtau-%s'%i]
	return Gtau_qmc['up_0'].data[:,0,0][::50]




def read_gtempo_orb3(chi2):
	filename = "orb3/anderson_tempo1_norb3_beta10.0_U4.0_J1.0_mu5.0_N200_chi60_chi2%s.json"%chi2

	print("Load from:", filename)
	if not os.path.exists(filename):
		return None
	
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gtau = data['gf']
	return -numpy.array(gtau)

def read_triqs_orb3():
	filename = "cthybn/bands3/beta10.0-U4.00-J1.00.h5"
	with HDFArchive(filename, 'r') as A:
		Gtau_qmc = A['Gtau-0']
	return Gtau_qmc['up_0'].data[:,0,0][::50]




def read_bethe_gtempo_orb3(i):
	filename = 'orb3/beta10_δτ0.05_chi60_chi2700/Gtau-%s.dat'%(i+1)
	print("Load from:", filename)
	data = loadtxt(filename)
	return data

def read_bethe_triqs_orb3(i):
	filename = "cthybn/bands3/beta10.0-U4.00-J1.00.h5"
	with HDFArchive(filename, 'r') as A:
		Gtau_qmc = A['Gtau-%s'%i]
	return Gtau_qmc['up_0'].data[:,0,0]

def get_bethe_orb3(i):
	beta = 10.0
	Gtau_gtempo = read_bethe_gtempo_orb3(i)
	taus_gtempo = numpy.linspace(0, beta, num=len(Gtau_gtempo))

	Gtau_qmc = read_bethe_triqs_orb3(i)
	taus_qmc = numpy.linspace(0, beta, num=len(Gtau_qmc))

	return taus_qmc, Gtau_qmc, Gtau_gtempo[::10]


plt.rcParams.update({'font.size': 14})

if True:
	fontsize=20
	fontsize2=14
	fontsize3=22

	beta = 10.0
	fig, ax = plt.subplots(1,2, figsize=(8,3))

	colors = ["r", "g", "b"]

	Gtau_qmc = read_bethe_triqs_orb2(0)
	taus_qmc = numpy.linspace(0, beta, num=len(Gtau_qmc)) * 2
	ax[0].plot(taus_qmc, Gtau_qmc, label="CTQMC", color="darkgray")
	chis = (100, 200, 300)
	for i in range(len(chis)):
		chi2 = chis[i]
		Gtau_gtempo = read_gtempo_orb2(chi2, 1)
		taus_gtempo = numpy.linspace(0, beta, num=len(Gtau_gtempo)) * 2
		ax[0].plot(taus_gtempo, Gtau_gtempo, ls = (0, (5, 5)), color=colors[i], label=r"$\chi_2$=%s"%chi2)
	ax[0].set_ylim(top=0)


	axin = ax[0].inset_axes([0.15,0.1,0.4,0.4])
	for i in range(len(chis)):
		chi2 = chis[i]
		Gtau_gtempo = read_gtempo_orb2(chi2, 1)
		taus_gtempo = numpy.linspace(0, beta, num=len(Gtau_gtempo)) * 2
		axin.plot(taus_gtempo, numpy.abs(Gtau_qmc-Gtau_gtempo), color=colors[i])
	axin.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
	axin.set_ylabel('Abs. Error', fontsize=fontsize2, labelpad=-0.1)
	axin.minorticks_on()
	axin.tick_params("both", which='both', direction='in', labelsize=fontsize2)
	axin.set_xlim((-0.1,20+0.1))
	axin.set_ylim(bottom=0)



	Gtau_qmc = read_triqs_orb3()
	chis = [300, 500, 700, 900]
	chis = [300, 500, 700]
	ax[1].plot(taus_qmc, Gtau_qmc, label="CTQMC", color="darkgray")
	for i in range(len(chis)):
		chi2 = chis[i]
		Gtau_gtempo = read_gtempo_orb3(chi2)
		taus_gtempo = numpy.linspace(0, beta, num=len(Gtau_gtempo)) * 2
		ax[1].plot(taus_gtempo, Gtau_gtempo, ls = (0, (5, 5)), color=colors[i], label=r"$\chi_2$=%s"%chi2)
	ax[1].set_ylim(top=0)

	axin = ax[1].inset_axes([0.15,0.1,0.4,0.4])
	for i in range(len(chis)):
		chi2 = chis[i]
		Gtau_gtempo = read_gtempo_orb3(chi2)
		taus_gtempo = numpy.linspace(0, beta, num=len(Gtau_gtempo)) * 2
		axin.plot(taus_gtempo, numpy.abs(Gtau_qmc-Gtau_gtempo), color=colors[i])
	axin.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
	axin.set_ylabel('Abs. Error', fontsize=fontsize2, labelpad=-0.1)
	axin.minorticks_on()
	axin.tick_params("both", which='both', direction='in', labelsize=fontsize2)
	axin.set_xlim((-0.1,20+0.1))
	axin.set_ylim(bottom=0)



	for i in range(2):
		ax[i].set_ylabel(r'$G(\tau)$', fontsize=fontsize)
		ax[i].set_xlabel(r'$D\tau$', fontsize=fontsize)
		ax[i].set_xlim([0, 20])

		ax[i].minorticks_on()
		ax[i].tick_params("both", which='both', direction='in', labelsize=fontsize)

		ax[i].legend(loc=(0.56, 0.1), fontsize=12)#, handlelength=1.0, handleheight=1,)

	trans = mtransforms.ScaledTranslation(-20/72, 7/72, fig.dpi_scale_trans)
	xpoint, ypoint = -0.2, 0.85
	labels = ["(a)", "(b)"]
	for i in range(2):
		ax[i].text(xpoint, ypoint, labels[i], transform=ax[i].transAxes + trans, fontsize=fontsize3, va='bottom')

	plt.tight_layout(h_pad=0.1, w_pad=0, pad=0.2)
	plt.savefig("pics/orb23_conv.pdf", bbox_inches='tight', dpi=300)
	plt.show()



if True:
	fontsize=20
	fontsize2=14
	fontsize3=22
	ns = [0, 1, 9]
	fig, ax = plt.subplots(1,2,figsize=(8,3))
	color_gtempo = ("r", "g", "b")
	color_qmc = ["tomato", "springgreen", "cyan"]
	taus = numpy.arange(0,10.05, 0.05) * 2
	for i in range(len(ns)):
		idx = ns[i]
		dat_gtempo = read_gtempo_orb2(200, idx+1)
		dat_qmc = read_bethe_triqs_orb2(idx)

		ax[0].plot(taus, dat_qmc, color=color_qmc[i])
		ax[0].plot(taus, dat_gtempo, color=color_gtempo[i], ls = (0, (5, 5)), label="%s-th"%(idx+1))
	
	ax[0].set_xlabel(r'$D\tau$', fontsize=fontsize)
	ax[0].set_ylabel(r'$G(\tau)$', fontsize=fontsize)

	axin = ax[0].inset_axes([0.19,0.1,0.4,0.4])
	for i in range(len(ns)):
		idx = ns[i]
		dat_gtempo = read_gtempo_orb2(200, idx+1)
		dat_qmc = read_bethe_triqs_orb2(idx)
		axin.plot(taus, numpy.abs(dat_gtempo-dat_qmc), color=color_gtempo[i], lw=1.2, label="%s-th"%(idx+1))
	axin.set_ylabel('Abs. Error', fontsize=fontsize2)
	axin.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
	axin.minorticks_on()
	axin.tick_params("both", which='both', direction='in', labelsize=fontsize2)
	axin.set_xlim((-0.1,20+0.1))
	axin.set_ylim(bottom=0)


	ns = [0, 1, 4]
	for i in range(len(ns)):
		idx = ns[i]
		res = get_bethe_orb3(idx)
		ax[1].plot(res[0]*2, res[1], color=color_qmc[i])
		ax[1].plot(res[0]*2, res[2], color=color_gtempo[i], ls = (0, (5, 5)), label="%s-th"%(idx+1))
	ax[1].set_xlabel(r'$D\tau$', fontsize=fontsize)
	ax[1].set_ylabel(r'$G(\tau)$', fontsize=fontsize)

	axin = ax[1].inset_axes([0.19,0.1,0.4,0.4])

	for i in range(len(ns)):
		idx = ns[i]
		res = get_bethe_orb3(idx)
		axin.plot(res[0]*2, res[1]-res[2], color=color_gtempo[i], lw=0.4, label="%s-th"%(idx+1))
	axin.set_ylabel('Abs. Error', fontsize=fontsize2)
	axin.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
	axin.minorticks_on()
	axin.tick_params("both", which='both', direction='in', labelsize=fontsize2)
	axin.set_xlim((-0.1,20+0.1))
	axin.set_ylim(bottom=0)


	
	line = [[(0, 0)]]
	lci = [mcol.LineCollection(2 * line, linestyles=[(0, (2, 2)), (0, ())], colors=[color_gtempo[i], color_qmc[i]]) for i in range(3)]
	ax[0].legend(lci, ['1st', "2nd", "10th"], handler_map={type(lci[0]): HandlerDashedLines()},
          handlelength=1.8, handleheight=1, loc=(0.64,0.1), fontsize=fontsize2)
	ax[1].legend(lci, ['1st', "2nd", "5th"], handler_map={type(lci[0]): HandlerDashedLines()},
          handlelength=1.8, handleheight=1, loc=(0.64,0.1), fontsize=fontsize2)

	for i in range(2):
		# ax[i].legend(loc=(0.65,0.2))
		ax[i].minorticks_on()
		ax[i].tick_params("both", which='both', direction='in', labelsize=fontsize)
		ax[i].set_xlim((-0.1,20.1))

	trans = mtransforms.ScaledTranslation(-20/72, 7/72, fig.dpi_scale_trans)
	xpoint, ypoint = -0.2, 0.85
	labels = ["(a)", "(b)"]
	for i in range(2):
		ax[i].text(xpoint, ypoint, labels[i], transform=ax[i].transAxes + trans, fontsize=fontsize3, va='bottom')

	plt.tight_layout(h_pad=0.1, w_pad=0, pad=0.2)
	plt.savefig("pics/orb23_DMFT.pdf")
	plt.show()




