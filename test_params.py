import social_evol

SEED=10
gsl_rng=social_evol.GSLrng(SEED)
NSAM=20
N=500
GENS=737 # failed assertion for 737 vs successful completion for 738
THETA=50
RHO=THETA
mu=THETA/(4.*N)
mu_del = 0.1*mu
r=RHO/(4.*N)

