
import numpy as np
import matplotlib.pyplot as plt
plt.style.use("classic")

def f(x,r):
	return r*x*(1.-x)

def generate_data(x0=.5,r=2.5,npoints=1e3):
	x = np.array([x0])
	for i in range(int(npoints-1)):
		x = np.append(x,f(x[-1],r))
	return x

def bifurcation():
	x0vals = np.random.rand(4)
	rvals  = np.linspace(3.55,3.6,801)
	for r in rvals:
		for x0 in x0vals:
			y = generate_data(x0=x0,r=r,npoints=3000)[2400:]
			rr = r*np.ones_like(y)
			plt.plot(rr,y,'k.', alpha=.01, markersize=2)
	#plt.savefig('bifurcation3.png',dpi=1200)


def cobweb(x,r,sty=dict(c='r',ls='-')):
	## outline
	xx = np.linspace(0,1,1000)
	plt.plot(xx,xx,'b-',lw=2, zorder=50)
	plt.plot(xx,f(xx,r),'k-',lw=3, zorder=100)
	#plt.gca().set_aspect('equal')
	#plt.xticks(np.arange(0,1.1,.25))
	#plt.yticks(np.arange(0,1.1,.25))
	#plt.grid()
	## first point
	a = x[0]
	sty0=dict(ls=':', c='b', lw=.5, zorder=1)
	plt.plot([a,a],[0,a],**sty0)
	## other points
	for i in range(len(x)-1):
		a, b = x[i], x[i+1]
		plt.plot([a,a,b],[a,b,b],**sty)
	## plot
	##plt.show()

def diverge_web(x1,x2,r,npoints=1000):
	y1 = generate_data(x0=x1,r=r,npoints=npoints)
	y2 = generate_data(x0=x2,r=r,npoints=npoints)
	cobweb(y1,r=r,sty=dict(c='r',ls='-',lw=1,alpha=1))
	cobweb(y2,r=r,sty=dict(c='b',ls='-.',lw=1,alpha=1))
	#plt.savefig('diverge-web.png',dpi=1000)
	plt.figure()
	plt.plot(y1,'r-')
	plt.plot(y2,'b-')
	plt.xlim(0,50)
	#plt.savefig('diverge.png',dpi=1000)

def volume(x1,r,int=.01,step=1000):
	x0 = np.linspace(x1,x1+int,100)
	x1 = np.zeros_like(x0)
	for i in range(len(x0)):
		x1[i] = generate_data(x0=x0[i],r=r,npoints=step)[-1]
	plt.plot(x0,x1,'kx')
	plt.xlim(0,1)
	plt.ylim(0,1)
	plt.show()

def lyap(x1,r,intt=.01):
	x0 = np.linspace(x1,x1+intt,100)
	x1 = np.zeros_like(x0)
	Q = np.zeros(100)
	for n in range(len(Q)):
		for i in range(len(x0)):
			x1[i] = generate_data(x0=x0[i],r=r,npoints=n)[-1]
		z = x1[1:-1] - np.roll(x1,1)[1:-1]
		Q[n] = np.sum(z**2)
	plt.plot(Q,'k-')
	#plt.savefig('lyap.png',dpi=1000)


x0 = np.random.rand()
r = 3.56 #np.random.rand()*4.
npoints = 500

x = generate_data(x0=x0,r=r,npoints=npoints)

plt.figure(1, figsize=(4,4))
s = .1
#plt.axes([1.*s,1.*s,1.-2.*s,1.-2.*s])
#plt.gca().axis('off')
plt.xlim(0,1)
plt.ylim(0,1)
plt.xticks([0,.5,1])
plt.yticks([0,.5,1])

plt.gca().set_aspect(1)

cobweb(x,r,sty=dict(c='r',ls='-',lw=.5, zorder=1000))

plt.annotate(s="\n".join(["r=%f"%r,"x0=%f"%x0]), xycoords="axes fraction", xy=(.01,.99), ha="left", va="top", fontsize=12)

plt.savefig("r%fx%fn%d.png"%(r,x0,npoints), dpi=300)

