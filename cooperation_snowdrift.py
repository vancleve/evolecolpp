import numpy as np
import matplotlib.pyplot as plt
import cooperation_snowdrift as snow

def evolve(rng, ngens, N, z0, mu, scale, b1, b2, c1, c2, grain=1):
    if type(z0) == float:
        P = snow.pop_init(N,z0)
    elif type(z0) == list and len(z0) == 2*N:
        P = snow.pop_init_dist(N,z0)
    else:
        print("Initial trait value z0 must be single value of list of length 2N")
        return False
    
    muts = list((snow.mutations(P),))
    for gen in range(ngens):
        X = snow.evolve_step(rng, P, mu, scale, b1, b2, c1, c2)
        if gen % grain == 0:
            muts.append(snow.mutations(P))
    
    return muts

def plotEvolve(vals, bins=100, interpolation='bicubic', vmin=None, vmax=None, tmult=1):
    parr = np.array(vals)
    nruns = parr.shape[0]

    phist = np.zeros((nruns,bins))
    for i in xrange(nruns):
        ph, xe = np.histogram(parr[i].ravel(), bins=bins, range=[0,1])
        phist[i] = ph

    extent=(min(xe),max(xe),0,nruns)
    plt.imshow(phist, aspect='auto', interpolation=interpolation,
               origin='lower', cmap=plt.cm.binary, extent=(min(xe),max(xe),0,nruns),
               vmin=vmin, vmax=vmax)
    yticks(around(linspace(0,nruns,11),-1), around(linspace(0,(nruns-1)*tmult+1,11),-1).astype('int'))
    plt.show()

# fix the seed and initialize the random number generator
SEED=314
gsl_rng = snow.GSLrng(SEED)


try: # run nothing by default if loaded in ipython
    __IPYTHON__
except: # run a sim and plot (if not loaded into ipython)
    # Generate Figure 1A from Doebeli, Hauert, and Killingback (2004, Science)
    # (noisier than Doebeli et al. due to smaller population size and higher mutation)
    w=10 # 'strength' of selection
    muts = evolve(gsl_rng, 40000, 1000, 0.01, 0.1, 0.005, 6.0*w, -1.4*w, 4.56*w, -1.6*w, grain=10)

    plotEvolve(muts)
