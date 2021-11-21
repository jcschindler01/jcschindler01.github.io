
import numpy as np
import matplotlib.pyplot as plt


"""
Consider a chain of coupled oscillators each of mass m,
each with self-coupling freq ws and neighbor-coupling freq wc.

Let

w0 = sqrt( ws^2 + 2*wc^2 ),

E0 = (1/2)*m*w0^2,

x = 2*wc^2 / w0^2,

s = w0*t,

ydot = dy/ds.


The no coupling limit is x=0, the no self limit is x=1.
Note that ws^2 = (1-x)*w0^2 and wc^2 = (1/2)*x*w0^2


Lagrangian is 

L = E0 * SUM_{ij=1...N} ( ydoti * Tij * ydotj - yi * Vij * yj )

where

Tij = delta_{i,j}

and Vij depends on boundary conditions.

For closed bcs:

Vij = delta_{i,j} - (x/2) * (delta_{i,j-1} + delta_{i-1,j}).

Open and periodic follow from this with slight variations.
"""

######################## matrix builders #################

def diag_ones(N):
	"""
	N x N matrix with ones on the diagonal.
	"""
	return np.diag(np.ones(N))

def offd_ones(N):
	"""
	N x N matrix with ones above and below the diagonal.
	"""
	## initialize
	a = np.zeros((N,N))
	## fill upper
	for i in range(N):
		for j in range(N):
			if i==j-1:
				a[i,j] = 1.
	## fill lower
	for i in range(N):
		for j in range(N):
			if i-1==j:
				a[i,j] = 1.
	## return
	return a

def corner_ones(N):
	"""
	N x N matrix with ones in top right and bottom left.
	"""
	## initialize
	a = np.zeros((N,N))
	## fill values
	a[0,N-1] = 1.
	a[N-1,0] = 1.
	## return
	return a

def firstlast_ones(N):
	"""
	N x N matrix with ones in first and last diagonal elements.
	"""
	## initialize
	a = np.zeros((N,N))
	## fill values
	a[0,0] = 1.
	a[N-1,N-1] = 1.
	## return
	return a

##############################################################

############## potential builders ###############################

def Vij_closed(N, x=0.2):
	"""
	Potential matrix with closed boundary conditions (i.e. y_{0} = y_{N+1} = 0).
	"""
	## init
	Vij = diag_ones(N) - 0.5 * x * offd_ones(N)
	## return 
	return Vij


def Vij_periodic(N, x=0.2):
	"""
	Potential matrix with periodic boundary conditions (i.e. spring connects y_1 to y_N).
	"""
	## init
	Vij = Vij_closed(N,x=x)  - 0.5 * x * corner_ones(N)
	## return 
	return Vij


def Vij_open(N, x=0.2):
	"""
	Potential matrix with open boundary conditions (i.e. y_1 and y_N only connected on one side).
	"""
	## init
	Vij = Vij_closed(N,x=x)  - 0.5 * x * firstlast_ones(N)
	## return 
	return Vij


def Vij_lambda(N, x=0.2, lam_open=0., lam_per=0.):
	"""
	Potential matrix with parameters (lamda_open,lambda_periodic) such that
	(0,0) = closed
	(1,0) = open
	(0,1) = periodic
	"""
	## print
	print "(lambda_open, lambda_per) = (%s, %s)"%(lam_open, lam_per)
	## terms
	V_closed = diag_ones(N) - 0.5 * x * offd_ones(N)
	deltaV_open =  - 0.5 * x * firstlast_ones(N)
	deltaV_periodic = - 0.5 * x * corner_ones(N)
	## total
	Vij = V_closed + lam_open * deltaV_open + lam_per * deltaV_periodic
	## return 
	return Vij 




#####################################################################3


##################### other ######################################

def eigmap(N=8, Vfunc=Vij_closed, npoints=1001, eps=0., sty=dict()):
	"""
	"""
	## params
	x = np.linspace(eps,1.-eps,npoints)
	## initialize
	eig = np.zeros((N,len(x)))
	## fill values by getting eigenvalue
	for i in range(len(x)):
		Vij = Vfunc(N,x=x[i])
		eig[:,i] = np.linalg.eigvalsh(Vij)
	## figure
	plt.title('eigmap')
	plt.grid(1)
	plt.xlim(0,1)
	plt.ylim(0,2)
	plt.xlabel('coupling param x ( no coupling <------> coupling only)')
	plt.ylabel('squared eigenfreqs wi^2 \n(units w0^2 = ws^2 + 2 wc^2)')
	print 'eigmap(N=%d, Vfunc=%s)'%(N,Vfunc.__name__)
	## plot
	for i in range(N):
		style = dict(ls='-', lw=3,c='k',alpha=0.5,zorder=50)
		style.update(sty)
		plt.plot(x,eig[i],**style)
	## plot w0 and its mirrorfor reference
	style2 = dict(c='k',alpha=0.1,ls='--',lw=10,zorder=10)
	plt.plot(x,(1.-x), **style2)
	plt.plot(x,(1.+x), **style2)
	## return
	return x, eig


def eigmap_compare(N=4):
	"""
	Plot eigenvalues vs coupling parameter for closed, open, and periodic case.
	"""
	## plot eigenvalues for closed case with N=N and print a sample value
	x, eig = eigmap(Vfunc=Vij_closed, N=N, sty=dict(c='r',ls='-'))
	print eig[:,500]
	## plot eigenvalues for open case with N=N+1 and print a sample value
	NN = N+1
	x, eig = eigmap(Vfunc=Vij_open  , N=NN, sty=dict(c='b',ls='-.'))
	print eig[:,500]
	## plot eigenvalues for periodic case with N=2(N-1) (even) or N=2(N-1) + 1 (odd) or  and print a sample value
	NN = 2*(N-1) + N%2
	x, eig = eigmap(Vfunc=Vij_periodic  , N=NN, sty=dict(c='g'))
	print eig[:,500]
	## the big gray dashed lines are w=ws
	## show
	plt.show()


def eigvecs(Vfunc=Vij_closed, N=24, mu=.5):
	## matrix and eigensysten
	Vij = Vfunc(N, x=mu)
	vals, vecs = np.linalg.eigh(Vij)
	## x values along unit line
	x0 = np.linspace(0.,1.,10)
	x  = np.linspace(0.,1.,N+2)[1:-1]
	## transpose eigenvec matrix and normalize vals
	vecs = np.transpose(vecs)
	normvals = vals / np.max(vals)
	colvals = (vals - np.min(vals)) / (np.max(vals) - np.min(vals))
	## size and shape
	yvals = 1. * vals
	dy = 0.8 * (np.max(yvals) - np.min(yvals)) / float(N)
	## figure
	plt.title('eigvecs')
	plt.grid(1)
	plt.xlim(0,1)
	plt.ylim(0,2)
	plt.xticks(np.arange(0,1.000001,.25))
	plt.yticks(np.arange(0,2.000001,.1))
	plt.xlabel('Oscillator Position $x$')
	plt.ylabel(r'Zero Reference Lines: Mode Frequency $\omega^2/\omega_0^2$'+'\n' +
				'Relative Height: Eigenvector Value $y_i$ (all with same arbitrary normalization)')
	## odd and even colormaps
	cms = [plt.cm.winter, plt.cm.Wistia]
	## plot
	for i in range(N):
		cm = cms[i%2]
		cv = colvals[i]
		style1 = dict(c=cm(cv),ls='-',lw=1,alpha=.8)
		style2 = dict(c=cm(cv),marker='o',markersize=4,ls='-',lw=1,alpha=.8)
		plt.plot(x0, yvals[i]+0.*x0,      **style1)
		plt.plot(x,  yvals[i]+dy*vecs[i], **style2)

def eigvecs_iterator(N=24,mu=0.5):
	for Vfunc in [Vij_closed,Vij_open,Vij_periodic]:
		print Vfunc.__name__
		plt.figure(figsize=(4,100))
		eigvecs(Vfunc=Vfunc,N=N,mu=mu)
		plt.tight_layout()
		fname = ('eigvecs_N(%d)_mu(%s)_V(%s)'%(N,mu,Vfunc.__name__)).replace('.','d')
		plt.savefig('fig/%s.pdf'%(fname))

def eigvals_vs_lambda(N=6, x=0.5):
	## lambda params
	ss = np.linspace(0,1,101)
	lam_open = 0.*ss
	lam_per  = 1.*ss
	## iterate and plot
	for i in range(len(ss)):
		Vij = Vij_lambda(N, x=x, lam_open=lam_open[i], lam_per=lam_per[i])
		vals = np.linalg.eigvalsh(Vij)
		for k in range(N):
			cv = float(k)/float(N)
			cm = plt.cm.gist_rainbow
			style = dict(marker='.', ls='none', c=cm(cv), markersize=10)
			plt.plot([ss[i]], [vals[k]], **style)
	## figure
	plt.grid(1)
	plt.xlim(0,1)
	plt.ylim(0,2)
	plt.xticks(np.arange(0,1.000001,.25))
	plt.yticks(np.arange(0,2.000001,.5))
	plt.title('eigenvalues vs lambda')
	plt.xlabel('Parameter $s$  such that\n' + r'$(\lambda_{open}(s), \,\lambda_{per}(s)): \; (%s,\,%s) \to (%s,\,%s)$'%(lam_open[0], lam_per[0], lam_open[-1], lam_per[-1]))
	plt.ylabel(r'$\omega^2/\omega_0^2$')
	plt.tight_layout()
	## save
	fname = 'eigvals_vs_lambda_N(%s)_x(%s)_lopen(%s-to-%s)_lper(%s-to-%s)'%(N,x,lam_open[0],lam_open[-1],lam_per[0],lam_per[-1])
	plt.savefig('fig/%s.pdf'%(fname))


def eigvecs_vs_lambda(N=50, mu=0.5, levels=[0]):
	## lambda params
	ss = np.linspace(0,1,25)
	lam_open = 1.-1.*ss
	lam_per  = 1.*ss
	## iterate and plot
	for i in range(len(ss)):
		## matrix and eigensysten
		Vij = Vij_lambda(N, x=mu, lam_open=lam_open[i], lam_per=lam_per[i])
		vals, vecs = np.linalg.eigh(Vij)
		## x values along unit line
		x0 = np.linspace(0.,1.,10)
		x  = np.linspace(0.,1.,N+2)[1:-1]
		## transpose eigenvec matrix and normalize vals
		vecs = np.transpose(vecs)
		normvals = vals / np.max(vals)
		colvals = (vals - np.min(vals)) / (np.max(vals) - np.min(vals))
		## size and shape
		if i==0:
			## use energy levels from intial lambda values
			yvals = 1. * vals
			dy = 0.8 * (np.max(yvals) - np.min(yvals)) / float(N)
		## figure
		plt.title('eigvecs')
		plt.grid(1)
		plt.xlim(0,1)
		#plt.ylim(0,2)
		plt.xticks(np.arange(0,1.000001,.25))
		plt.yticks(np.arange(0,2.000001,.1))
		plt.xlabel('Oscillator Position $x$')
		plt.ylabel(r'Zero Reference Lines: Mode Frequency $\omega^2/\omega_0^2$'+'\n' +
					'Relative Height: Eigenvector Value $y_i$ (all with same arbitrary normalization)')
		#plt.tight_layout()
		## colormaps
		cms = [plt.cm.gnuplot2_r]
		## plot
		for k in levels:
			cm = cms[0]
			cv = float(i)/float(len(ss))
			av = np.linspace(.5,1.,len(ss))[i]
			style1 = dict(c=cm(cv),ls='-',lw=1,alpha=av)
			style2 = dict(c=cm(cv),marker='o',markersize=4,ls='-',lw=1,alpha=av)
			plt.plot(x0, yvals[k]+0.*x0,      **style1)
			plt.plot(x,  yvals[k]+dy*vecs[k], **style2)
		



##########################################################


############### main ##################################


eigvecs_iterator(N=51,mu=1.)


#############################################################




