import social_evol_inf_alleles

SEED=10
gsl_rng=social_evol_inf_alleles.GSLrng(SEED)


def evolve(rng, ngens, N, z0, mu, scale, b1, b2, c1, c2):
    P = social_evol_inf_alleles.pop_init(N,z0)
    phenos = list((social_evol_inf_alleles.phenotypes(P),))
    for gen in range(ngens):
        X = social_evol_inf_alleles.evolve_step(rng, P, mu, scale, b1, b2, c1, c2)
        phenos.append(social_evol_inf_alleles.phenotypes(P))

    return phenos

def plotEvolve(phenos,bins=10):
    phist = zeros((len(phenos),bins))
    for i in range(len(phenos)):
        h,b = histogram(phenos[i], range=(0,1.0), bins=bins)
        phist[i] = h

    imshow(array(phist),aspect='auto', origin='lower', cmap='binary'); colorbar(); xticks(linspace(0,bins-1,10), map(lambda x: '{:.2f}'.format(x), linspace(0.5/bins,1-0.5/bins,10)))
