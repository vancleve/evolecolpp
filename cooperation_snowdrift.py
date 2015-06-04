import numpy as np
import matplotlib.pyplot as plt
import cooperation_snowdrift

def evolve(rng, ngens, N, z0, mu, scale, b1, b2, c1, c2, grain=1):
    if type(z0) == float:
        P = cooperation_snowdrift.pop_init(N,z0)
    elif type(z0) == list and len(z0) == 2*N:
        P = cooperation_snowdrift.pop_init_dist(N,z0)
    else:
        print("Initial trait value z0 must be single value of list of length 2N")
        return False
    
    muts = list((cooperation_snowdrift.mutations(P),))
    for gen in range(ngens):
        X = cooperation_snowdrift.evolve_step(rng, P, mu, scale, b1, b2, c1, c2)
        if gen % grain == 0:
            muts.append(cooperation_snowdrift.mutations(P))
    
    return muts

def plotEvolve(vals, bins=100, interpolation='bicubic'):
    parr = np.array(vals)
    times = np.repeat(np.arange(parr.shape[0]), parr.shape[1])
    plist = parr.ravel()
    
    phist, xe, ye = np.histogram2d(plist, times, bins=bins, range=[[0,1], [0, len(vals)]])

    plt.imshow(phist.T, aspect='auto', interpolation=interpolation,
           origin='lower', cmap=cm.binary, extent=(min(xe),max(xe),min(ye),max(ye)))



# fix the seed and initialize the random number generator
SEED=314
gsl_rng = cooperation_snowdrift.GSLrng(SEED)

# run a sim and plot (if not loaded into ipython)
if '__IPYTHON__' not in vars(__builtin__):
    # Generate Figure 1A from Doebeli, Hauert, and Killingback (2004, Science)
    # (noisier than Doebeli et al. due to smaller population size and higher mutation)
    w=10 # 'strength' of selection
    muts = evolve(gsl_rng, 40000, 1000, 0.01, 0.1, 0.005, 6.0*w, -1.4*w, 4.56*w, -1.6*w, grain=10)    

    plotEvolve(muts)
