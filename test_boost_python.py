import fwdpp_boost_python
import sys

SEED=10
gsl_rng=fwdpp_boost_python.GSLrng(SEED)
NSAM=20
N=500
# Assertion failed: (this->m_holder.m_size > n), function operator[], file /usr/local/include/boost/container/vector.hpp, line 1488
GENS=int(sys.argv[1])
# Success
#GENS=971
THETA=50
RHO=THETA
mu=THETA/(4.*N)
r=RHO/(4.*N)

XX=fwdpp_boost_python.evolve(gsl_rng,N,GENS,mu,r)
XX.sfs=fwdpp_boost_python.sfs(gsl_rng,XX,NSAM)

print(XX.sfs)
