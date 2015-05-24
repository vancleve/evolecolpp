import social_evol_inf_alleles

SEED=10
gsl_rng=social_evol_inf_alleles.GSLrng(SEED)
NSAM=20
N=500
GENS=10000 # failed assertion for 737 vs successful completion for 738
THETA=50
RHO=THETA
mu=THETA/(4.*N)
mu_del = 0.1*mu
r=RHO/(4.*N)


def evolve(rng, ngens, N, z0, mu, scale, b1, b2, c1, c2):
    P = social_evol_inf_alleles.pop_init(N,0.05)
    phenos = list((social_evol_inf_alleles.phenotypes(P),))
    for gen in range(ngens):
        X = social_evol_inf_alleles.evolve_step(gsl_rng,P,mu,0.05,6.0, -1.4, 4.56, -1.6)
        phenos.append(social_evol_inf_alleles.phenotypes(P))

    return phenos

def plotEvolve(phenos):
    phist = list()
    for i in range(len(phenos)):
        h,b = histogram(phenos[i], range=(0,1.0), bins=10)
        phist.append(h)

    imshow(array(phist),aspect='auto', origin='lower', cmap='binary'); colorbar(); xticks(range(10), arange(0.05,1,0.1))
