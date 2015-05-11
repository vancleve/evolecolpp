#Simulate a population of N=500 diploids for 10N generations
#with theta = rho = 50.  Take a sample of size n = 20 chromosomes,
#and get the mean SFS
import fwdpp_boost_python
import random
import sys

SEED=int(sys.argv[1])
gsl_rng=fwdpp_boost_python.GSLrng(SEED)
NREPS=100
NSAM=20
N=500
GENS=10*N
THETA=50
RHO=THETA
mu=THETA/(4.*N)
r=RHO/(4.*N)
SFS=[0]*(NSAM-1)
for i in range(0,NREPS):
    XX=fwdpp_boost_python.evolve(gsl_rng,N,GENS,mu,r)
    XX.sfs=fwdpp_boost_python.sfs(gsl_rng,XX,NSAM)
    for i in range(0,(NSAM-1)):
        SFS[i]=SFS[i]+XX.sfs[i]


#get mean of SFS
for i in range(0,(NSAM-1)):
    SFS[i]=float(SFS[i])/float(NREPS)
print(SFS)
