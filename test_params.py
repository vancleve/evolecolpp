import social_evol_inf_alleles

SEED=10
gsl_rng=social_evol_inf_alleles.GSLrng(SEED)


def evolve(rng, ngens, N, z0, mu, scale, b1, b2, c1, c2, grain=1):
    if type(z0) == float:
        P = social_evol_inf_alleles.pop_init(N,z0)
    elif type(z0) == list and len(z0) == 2*N:
        P = social_evol_inf_alleles.pop_init_dist(N,z0)
    else:
        print("Initial trait value z0 must be single value of list of length 2N")
        return False
    
    phenos = list((social_evol_inf_alleles.phenotypes(P),))
    for gen in range(ngens):
        X = social_evol_inf_alleles.evolve_step(rng, P, mu, scale, b1, b2, c1, c2)
        if gen % grain == 0:
            phenos.append(social_evol_inf_alleles.phenotypes(P))
    
    return phenos

def plotEvolve(phenos,bins=100, interpolation='bilinear'):

    parr = array(phenos)
    times = repeat(arange(parr.shape[0]), parr.shape[1])
    plist = parr.ravel()
    
    phist, xe, ye = histogram2d(plist, times, bins=bins, range=[[0,1], [0, len(phenos)]])

    imshow(phist.T, aspect='auto', interpolation=interpolation,
           origin='lower', cmap=cm.binary, extent=(min(xe),max(xe),min(ye),max(ye)))
    colorbar()
