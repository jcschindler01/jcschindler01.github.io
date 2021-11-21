
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as an
import shutil
import time
import psi as ps
plt.style.use("classic")


## params
mu = 20.
nmax = 100

## animation params
nperiods = 5. ## number of timescales T to animate
persteps = 2.  ## number of steps per T period
Ttime = .2 ## time in seconds of single T in animation

## initial data xvals
x0 = np.linspace(0,1,5001)

## initial wavefunction
psi0 = np.exp(-100.*(x0-.5)**2) * np.exp(1j*50.*x0)

## get psi
psi = ps.psi_from_wavefunction(x0, psi0, mu, nmax=nmax)

## vals
x = np.linspace(0,1,5001)

## energy
Etot1 = np.sum((1./float(len(x)-1))*psi.energy(x,0))
Etot2 = np.sum(0.5 * psi.w**2 * np.abs(psi.ck)**2)

"""
Set up plots.
"""

## figure
fig = plt.figure(1, figsize=(10,8))

## spectrum plot
ax1 = plt.subplot(223)
plt.grid(1)
plt.title("Spectrum")
plt.xlabel("$k/\pi$")
plt.ylabel("$|c_n|^2$")
plt.annotate(s=r"$m/\pi\approx%.1f$"%(psi.mm/np.pi), xy=(.02,.96), xycoords='axes fraction', 
	bbox=dict(fc='w', ec='none', alpha=.7), ha="left", va="top", fontsize=11)

## classical field plot
ax2 = plt.subplot(222)
plt.xlim(0,1)
plt.xticks(np.arange(0,1.1,.1))
plt.ylim(-10,10)
plt.grid(1)
plt.title("Classical Field")
plt.xlabel("$x$")
plt.ylabel("$y(x,t)$")


## wavefunction plot
ax3 = plt.subplot(224)
plt.xlim(0,1)
plt.ylim(-10,10)
plt.xticks(np.arange(0,1.1,.1))
plt.grid(1)
plt.title("Wavefunction")
plt.xlabel("$x$")
plt.ylabel("$\psi(x,t)$")

## spectrum plot 2
ax4 = plt.subplot(221)
plt.grid(1)
plt.title("Spectrum")
plt.xlabel("$k/\pi$")
plt.ylabel("$\omega/\pi$")
plt.annotate(s=r"$m/\pi\approx%.1f$"%(psi.mm/np.pi), xy=(.02,.96), xycoords='axes fraction', 
	bbox=dict(fc='w', ec='none', alpha=.7), ha="left", va="top", fontsize=11)


## layout
plt.tight_layout()

"""
Static plots.
"""


## spectrum plot 1
plt.sca(ax1)
plt.plot([psi.mm/np.pi,psi.mm/np.pi],[0.,1.5*np.max(np.abs(psi.ck)**2)], 'r--', label='$k=m$')
plt.plot(psi.k/np.pi, np.abs(psi.ck)**2, 'ko-', label='$|c_k|^2$')
plt.ylim(0.,1.5*np.max(np.abs(psi.ck)**2))
plt.legend(loc='upper right', fontsize=8)


## spectrum plot 2
plt.sca(ax4)
kk = np.linspace(0,np.max(psi.k),1001)
ww = np.sqrt(psi.mm**2 + kk**2)
plt.plot(kk/np.pi, kk/np.pi, 'b--', label='$\omega=k$')
plt.plot(kk/np.pi, (psi.mm + 0.*kk)/np.pi, 'r--', label='$\omega=m$')
plt.plot(kk/np.pi, ww/np.pi, 'k-', lw=2, label='$\omega(k)=\sqrt{m^2+k^2}$')
plt.legend(loc='upper right', fontsize=8)


"""
Time dependent plots.
"""

t = 0.


## classical field plot
plt.sca(ax2)
l0, = plt.plot([],[], 'g-', lw=1, alpha=.6, label="$y^{\prime}(x)/m$")
l1, = plt.plot([],[], 'r-', lw=1, alpha=.6, label='$\dot{y}(x)/m$')
l2, = plt.plot([],[], 'k-', lw=4, alpha=.5, label='$E(x)/m^2$')
l3, = plt.plot([],[], 'b-', lw=3, alpha=1., label='$y(x)$')
plt.legend(loc='upper right', fontsize=8)

ss = "\n".join([
	r"$2E_{tot}/m^2 \approx %.3f$"%(2.*Etot2/psi.mm**2), 
	r"$T = 2\pi \, (2E_{tot})^{-1/2}$",
	])
plt.annotate(s=ss, xy=(.02,.04), xycoords='axes fraction', 
	bbox=dict(fc='w', ec='none', alpha=.9), ha="left", va="bottom", fontsize=11)

a1 = plt.annotate(s="", xy=(.02,.96), xycoords='axes fraction', 
	bbox=dict(fc='w', ec='none', alpha=.9), ha="left", va="top", fontsize=11)

## wavefunction plot
plt.sca(ax3)
l4, = plt.plot([],[], 'r-', lw=1, alpha=.6, label='$Im(\psi(x))$')
l5, = plt.plot([],[], 'k-', lw=4, alpha=.5, label='$|\psi(x)|^2$')
l6, = plt.plot([],[], 'b-', lw=3, alpha=1., label='$Re(\psi(x))$')
plt.legend(loc='upper right', fontsize=8)

ss = "\n".join([
	r"$2E_{tot}/m^2 \approx %.3f$"%(2.*Etot2/psi.mm**2), 
	r"$T = 2\pi \, (2E_{tot})^{-1/2}$",
	])
plt.annotate(s=ss, xy=(.02,.04), xycoords='axes fraction', 
	bbox=dict(fc='w', ec='none', alpha=.9), ha="left", va="bottom", fontsize=11)

a2 = plt.annotate(s="", xy=(.02,.96), xycoords='axes fraction', 
	bbox=dict(fc='w', ec='none', alpha=.9), ha="left", va="top", fontsize=11)


"""
Animation.
"""

## timescale
T = 2.*np.pi/np.sqrt(2.*Etot2)

## tvals
tvals = T * np.linspace(0., nperiods, int(nperiods*persteps+1))
interval = int( 1000. * Ttime / persteps )

## init
def init():
	## return
	return [l0,l1,l2,l3,l4,l5,l6,a1,a2]



## update
def update(t):
	## print
	print("t/T = %7.3f"%(t/T))
	## classical plot
	l0.set_data(x, psi.yprime(x,t)/psi.mm)
	l1.set_data(x, psi.ydot(x,t)/psi.mm)
	l2.set_data(x, psi.energy(x,t)/psi.mm**2)
	l3.set_data(x, psi.y(x,t))
	a1.set_text(r"$t/T \approx %7.3f$"%(t/T))
	## wavefunction plot
	l4.set_data(x, psi.psi(x,t).imag)
	l5.set_data(x, np.abs(psi.psi(x,t))**2)
	l6.set_data(x, psi.psi(x,t).real )
	a2.set_text(r"$t/T \approx %7.3f$"%(t/T))
	## return
	return [l0,l1,l2,l3,l4,l5,l6,a1,a2]
	

## animate
print( "animate...")
anim = an.FuncAnimation(fig, update, frames=tvals, interval=interval, init_func=init, blit=True)
print( "animate done")

## save
if True:
	## get path with timestamp
	ts = str(time.time()).replace(".","")
	path = "anim_out/anim_%s"%(ts)
	#anim.save("anim_out/anim_%s.mp4"%(ts), writer='ffmpeg', fps=None, codec=None, dpi=None, bitrate=6000)
	## save
	print( "save...")
	anim.save("%s.gif"%(path), writer='imagemagick', dpi=None)
	print( "save done")
	## copy to temp
	print( "copy...")
	shutil.copy("%s.gif"%(path), "anim_out/000_temp.gif")
	print( "copy done")

## show
plt.show()

