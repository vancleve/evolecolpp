#Simulate a population of N=500 diploids for 10N generations
#with theta = rho = 50.  Take a sample of size n = 20 chromosomes,
#and get the mean SFS
import social_evol
import sys

SEED=10
gsl_rng=social_evol.GSLrng(SEED)
NSAM=20
N=500
GENS=int(sys.argv[1]) # failed assertion for 737 vs successful completion for 738
THETA=50
RHO=THETA
mu=THETA/(4.*N)
mu_del = 0.1*mu
r=RHO/(4.*N)

XX = social_evol.poptype(N)
social_evol.evolve(gsl_rng,XX,N,GENS,mu,mu_del,r,6.0, -1.4, 4.56, -1.6)
XX.sfs=social_evol.sfs(gsl_rng,XX,NSAM)

print(XX.sfs)
