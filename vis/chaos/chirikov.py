## Iterate and plot orbits of the chirikov standard map given by
##      P(n+1) = P(n) - K sin(x(n))
##      x(n+1) = x(n) + P(n+1)
## where -pi < x <= pi and -inf < P < inf.


## preamble
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm






## subs

def modpi(x):   # x + 2pi = x so keep x in the range (-pi,pi]
    if x <= -np.pi:
        x += 2.*np.pi
    if x > np.pi:
        x += -2.*np.pi
    if x<=-np.pi or x>np.pi:
        x = modpi(x)
    return x
    
def arraymodpi(xx):
    yy = np.zeros_like(xx)
    for i in range(len(xx)):
        yy[i] = modpi(xx[i])
    return yy
          
def nextvals(vals,K):
    [P,x] = vals
    P = P - K * np.sin(x)
    x = x + P
    x = modpi(x)
    return np.array([P,x])

def extend1(PP,xx,K):
    if len(xx) != len(PP):
        print 'ERROR: len(xx) != len(PP)'
    else:
        [P,x] = nextvals([PP[-1:],xx[-1]],K)
        PP = np.append(PP,P)
        xx = np.append(xx,x)
    return [PP,xx]

def extendN(PP,xx,K,N):
    for n in range(N):
        [PP,xx]=extend1(PP,xx,K)
    return [PP,xx]

def rand_orbit(K,N):
    [PP,xx] = [arraymodpi(np.random.normal(0.0,np.pi/2.,1)),arraymodpi(np.random.normal(0.0,np.pi/2.,1))]
    [PP,xx] = extendN(PP,xx,K,N)
    return [PP,xx]

def orbit(K,N,init_point):
    [PP,xx] = init_point
    [PP,xx] = extendN(PP,xx,K,N)
    return [PP,xx]


def orbitplot(K,Nsteps,init_points,newfig=True):
    Norbits = len(init_points)
    if newfig:
        plt.figure('Kicked Rotor: K = %s, Nsteps=%s, Norbits=%s' %(K,Nsteps,Norbits))
    plt.axis([-3.15,3.15,-3.5,3.5])
    colors = iter(cm.cool(np.linspace(0,0.6,Norbits)))
    for i in range(Norbits):
        [PP,xx] = orbit(K,Nsteps,init_points[i])
        print PP,xx
        col=next(colors)
#        col='b'
        plt.scatter(xx,PP,marker='.',s=1.5,edgecolor='none',c=col)
#    plt.title('K=%s'%K)
#    plt.savefig('../figures/orbitplot000.pdf')

def random_initpoints(Norbits):
    pts = []
    for i in range(Norbits):
        pts += [[arraymodpi(np.random.normal(0.0,np.pi/2.,1)),arraymodpi(np.random.normal(0.0,np.pi/2.,1))]]
    return pts

def Plinear_initpoints(Norbits):
    Pmin = -4.
    Pmax = 4
    Ppts = np.linspace(Pmin,Pmax,Norbits)
    initpoints = []
    for i in range(len(Ppts)):
        initpoints += [[np.array([Ppts[i]]),np.array([0.0])]]
    return initpoints
    

def random_orbitplot(K,Nsteps,Norbits):
    init = random_initpoints(Norbits)
    orbitplot(K,Nsteps,init)
   
def sixplots():
    Norbits = 50
    Nsteps = 2000
    plt.figure('Kicked Rotor: Norbits=%s, Nsteps=%s' %(Norbits,Nsteps),figsize=(12,6))
    init = Plinear_initpoints(Norbits)
    K = np.linspace(0.1,.999,6)
    for i in range(6):
        plt.subplot(2,3,i+1)
        plt.annotate(s='K=%s'%(K[i]),xy=(0.0,-3.0))
        plt.gca().xaxis.set_visible(False)
        plt.gca().yaxis.set_visible(False)
        orbitplot(K[i],Nsteps,init,newfig=False)
    plt.savefig('../figures/sixplots000.jpg')

def rotate(xy,theta):
    [x,y]=xy
    th = np.pi*theta/180.
    for i in range(len(x)):
        [x[i],y[i]] = [np.cos(th)*x[i]+np.sin(th)*y[i],np.cos(th)*y[i]-np.sin(th)*x[i]]
    xy = [x,y]
    return xy
        
def randab(a,b):
    return (b-a)*np.random.rand(1) + a




############################################# art ####################
    
def artA(newplot=True):
    kvals = np.linspace(0.001,.99,8)
    if newplot:
        plt.ioff()
        plt.figure('art',figsize=(10,10))
        plt.axis([-3.15,3.15,-3.5,3.5])
        plt.gca().xaxis.set_visible(False)
        plt.gca().yaxis.set_visible(False)
        plt.gca().set_axis_bgcolor('k')
        plt.gca().set_aspect(.75)
    for k in kvals:
        Norbits = 25
        Nsteps = 500
        init = Plinear_initpoints(Norbits)
        col='white'
        for i in range(Norbits):
            [PP,xx] = orbit(k,Nsteps,init[i])
            print PP,xx
            plt.scatter(xx,PP,marker='.',s=1.5,edgecolor='none',c=col)
    if newplot:
        plt.savefig('../figures/art/artA.jpg',dpi=600,facecolor='black')
        print 'art = complete'
###############################  used^^^^ ####################

def artB(newplot=True):
    kvals = np.linspace(0.0001,1.2,8)
    bgcol1 = cm.seismic(0.05)
    bgcol2 = cm.seismic(0)
    if newplot:
        plt.ioff()
        plt.figure('art',figsize=(10,10))
        plt.axis([-3.15,3.15,-3.5,3.5])
        plt.gca().xaxis.set_visible(False)
        plt.gca().yaxis.set_visible(False)
        plt.gca().set_axis_bgcolor(bgcol1)
        for spine in ['top','bottom','left','right']:
            plt.gca().spines[spine].set_color(bgcol2)
#        plt.gca().axis('off')
        plt.gca().set_aspect(.75)
    colors = iter(cm.cool(np.linspace(0,.25,len(kvals))))
    for k in kvals:
        Norbits = 25
        Nsteps = 500
        init = Plinear_initpoints(Norbits)
        col=next(colors)
        for i in range(Norbits):
            [PP,xx] = rotate(orbit(k,Nsteps,init[i]),150.*k)
            print PP,xx
            plt.scatter(xx,PP,marker='.',s=1.8,edgecolor='none',c=col)
#            plt.draw()
    if newplot:
        plt.savefig('../figures/art/artB.jpg',facecolor=bgcol2,dpi=600)
        print 'art = complete'

################################ used ^^^^^^^^  #################

def artD(newplot=True):
    Ncurves = 900
    colset = cm.autumn
    bgcol = 'midnightblue'
    if newplot:
        plt.ioff()
        plt.figure('art',figsize=(12,10))
        plt.axis([-3.15,3.15,-3.5,3.5])
        plt.gca().xaxis.set_visible(False)
        plt.gca().yaxis.set_visible(False)
        plt.gca().set_axis_bgcolor(bgcol)
#        for spine in ['top','bottom','left','right']:
#            plt.gca().spines[spine].set_color(bgcol)
        plt.gca().axis('off')
        plt.gca().set_aspect(.75)
    for i in range(Ncurves):
        k = randab(0.00001,2.0)
        Nsteps = randab(200,2000)
        [p0,x0] = [randab(-4.,4.),randab(-np.pi,np.pi)]
#        col = (np.random.rand(),np.random.rand(),np.random.rand(),0.5*(np.random.rand()+1.))
        col = colset(randab(.6,1))
        theta = randab(0,280)
        [pp,xx] = rotate(orbit(k,Nsteps,[p0,x0]),theta)
        print k,pp,xx
        plt.scatter(xx,pp,marker='.',s=1.5,edgecolor='none',c=col)
#        plt.draw()
    if newplot:
        plt.savefig('../figures/art/artD.jpg',facecolor=bgcol,dpi=800)
        print 'art = complete'
    
########################################


def fullmapA(k=.97,newplot=True):
    # params
    Norbits = 180
    Nsteps = 1000
    # set up axes
    cmap = cm.hot
    bgcol1 = 'w'
    bgcol2 = 'snow'
    if newplot:
        plt.ioff()
        plt.figure('art',figsize=(10,10))
        plt.axis([-3.15,3.15,-3.5,3.5])
        plt.gca().xaxis.set_visible(False)
        plt.gca().yaxis.set_visible(False)
        plt.gca().set_axis_bgcolor(bgcol1)
        for spine in ['top','bottom','left','right']:
            plt.gca().spines[spine].set_color('k')
        plt.gca().set_aspect(.75)
    colors = iter(cm.cool(np.linspace(0,.72,Norbits)))
#    colors = iter(cmap(np.linspace(.5,.15,Norbits)))
    # create plot
    init = Plinear_initpoints(Norbits)
    for i in range(Norbits):
        [PP,xx] = orbit(k,Nsteps,init[i])
        print PP,xx
        col=next(colors)
        plt.scatter(xx,PP,marker='.',s=2,edgecolor='none',c=col)
    if newplot:
        plt.savefig('../figures/art/fullmapA.jpg',facecolor=bgcol2,dpi=600)
        print 'art = complete'






## main

def main():
    fullmapA(.97)


    




main()
















