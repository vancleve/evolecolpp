/*
In addition to the usual fwdpp depdencies, we need boost.python.
*/

#include <limits>
#include <algorithm>
#include <boost/python.hpp>
#include <boost/container/list.hpp>
#include <boost/container/vector.hpp>
#include <boost/pool/pool_alloc.hpp>
#include <boost/unordered_set.hpp>
#include <boost/functional/hash.hpp>

// Include custom mutation header
#include<mutation_inf_alleles.hpp>
// Main fwdpp library header
#include <fwdpp/diploid.hh>
// Include the necessary "sugar" components
// We need to get the 'deep' version of singlepop, as we need to make a custom singlepop_serialized_t for our sim
#include <fwdpp/sugar/singlepop/singlepop.hpp>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/sugar/serialization.hpp>

/*
  We will use a gsl_rng_mt19937 as our RNG.
  This type is implicitly convertible to gsl_rng *,
  and auto-handles the gsl_rng_free steps, etc.
 */
using GSLrng = KTfwd::GSLrng_t<KTfwd::GSL_RNG_MT19937>;

//FWDPP-related types
using mtype = KTfwd::mutation;
using mlist_t = boost::container::list<mtype,boost::pool_allocator<mtype> >;
using gamete_t = KTfwd::gamete_base<mtype,mlist_t>;
using glist_t = boost::container::list<gamete_t, boost::pool_allocator<gamete_t>>;

//We need to define a custom diploid genotype for our model
struct diploid_t : public KTfwd::tags::custom_diploid_t
{
  using first_type = glist_t::iterator;
  using second_type = glist_t::iterator;
  first_type first;
  second_type second;
  unsigned i;
  //constructors, etc.
  diploid_t() : first(first_type()),second(second_type()),i(std::numeric_limits<unsigned>::max()) {}
  //"perfect forwarding" constructor does not work with iterator from boost containers...
  //diploid_t(first_type && g1, first_type && g2) : first(std::forward(g1)),second(std::forward(g2)),i(numeric_limits<unsigned>::max()) {}
  diploid_t(first_type g1, first_type g2) : first(g1),second(g2),i(std::numeric_limits<unsigned>::max()) {}
  //The following constructors SHOULD be generated automagically by your compiler, so you don't have to:
  //(no idea what, if any, performance effect this may have.  Worst case is prob. the move constructor doesn't get auto-generated...
  //diploid_t( const diploid_t & ) = default;
  //diploid_t( diploid_t && ) = default;
  //diploid_t & operator=(const diploid_t &) = default;
};

//Define our our population type via KTfwd::sugar 
using poptype = KTfwd::sugar::singlepop_serialized<mtype,
						   KTfwd::mutation_writer,
						   KTfwd::mutation_reader<mtype>,
						   mlist_t,
						   glist_t,
						   boost::container::vector< diploid_t >,
						   boost::container::vector<mtype>,
						   boost::container::vector<unsigned>,
						   boost::unordered_set<double,boost::hash<double>,
									KTfwd::equal_eps>
						   >;

/*! \brief Continuous Snowdrift Game from Doebeli, Hauert, and Killingback (2004, Science, 306:859--862)
    as implemented by Wakano and Lehmann (2014, J Theor Biol, 351:83--95)
    \param phenotypes phenotypes in population
    \param fitnesses fitnesses in population
    \param b1 linear benefit term
    \param b2 quadratic benefit term
    \param c1 linear cost term
    \param c2 quadratic cost term
    \return Mean payoff (fitness) from pairwise interactions in snowdrift game
    \ingroup fitness
  */
struct snowdrift_diploid
{
  using result_type = double;
  inline result_type operator()( const poptype::dipvector_t::const_iterator dip,
				 const std::vector<double> & phenotypes, 
				 const double & b1, const double & b2,
				 const double & c1, const double & c2) const
  {
    unsigned N = phenotypes.size();
    double zself = phenotypes[dip->i];
    result_type fitness = 0;
    for (unsigned j = 0; j < N; ++j)	  
      {
	if (dip->i != j)
	  {
	    double zpair = zself+phenotypes[j];
	    // Payoff function from Fig 1
	    double a = b1*zpair + b2*zpair*zpair - c1*zself - c2*zself*zself;
	    fitness += 1 + std::max(a, 0.0);
	  }
      }
    return fitness/double(N-1);
  }
};

// Truncated Laplace distribution using rejection method (inefficient)
double ran_trunc_laplace( gsl_rng * r, double mean, double scale, double lo, double hi )
{
  double rv = mean + gsl_ran_laplace(r, scale);
  while (rv > hi && rv < lo)
    {
      rv = mean + gsl_ran_laplace(r, scale);
    }
  return mean + gsl_ran_laplace(r, scale);
}

// Calculate phnenotype for sinlge "locus"
double phenotypef_1locus( const diploid_t & dip )
{
  gamete_t::mcont_iterator m1 = dip.first->smutations.begin();
  gamete_t::mcont_iterator m2 = dip.second->smutations.begin();

  double pval = ((*m1)->h * (*m1)->s + (*m2)->h * (*m2)->s) / ((*m1)->h + (*m2)->h);
  return pval;
}

// Initialize population to be homogeneous with mutation value mutval0
poptype pop_init( const unsigned & N,
		  const double & mutval0)
{
  poptype pop(N); // population initially is single gamete in 2*N copies
  mtype mut0(mutval0, mutval0, 2*N); // initial mutation also in 2*N copies

  // insert mutation into population mutation list
  poptype::mlist_t::iterator mitr = pop.mutations.insert(pop.mutations.end(), std::move(mut0));

  // insert mutation in mutation list of gamete
  poptype::glist_t::iterator g0 = pop.gametes.begin();
  g0->smutations.insert( g0->smutations.begin(), mitr );

  return pop;
}

// Evolve population one generation using infinite alleles mutation model
double evolve_step( GSLrng & rng,
		    poptype & pop,
		    const unsigned & generation,
		    const double & mu, const double scale,
		    const double & b1, const double & b2, const double & c1, const double & c2)
{
  std::vector<double> phenotypes(pop.N);
  //Fill phenotypes
  unsigned i = 0;

  for( auto & dip : pop.diploids ) 
    { 
      dip.i = i; 
      phenotypes[i++] = phenotypef_1locus(dip); 
    }

  double wbar =
    KTfwd::sample_diploid(rng,
			  &pop.gametes,
			  &pop.diploids,
			  &pop.mutations,
			  pop.N,
			  mu,

			  // infinite alleles mutation model where mutations have
			  // truncated laplace distribution on (0,1)
			  std::bind(inf_alleles(),
				    rng,
				    &pop.mut_lookup,
				    std::placeholders::_1,
				    [&scale](gsl_rng * r, double c){
				      return ran_trunc_laplace(r, c, scale, 0.0, 1.0);
				    }),

			  // Recombination policy: do nothing, return zero breaks
			  [](poptype::glist_t::iterator & g1,
			     poptype::glist_t::iterator & g2){return 0;},

			  // Mutation insertion policy
			  std::bind(KTfwd::insert_at_end<poptype::mutation_t,
				    poptype::mlist_t>,
				    std::placeholders::_1,std::placeholders::_2),

			  // Gamete insertion policy
			  std::bind(KTfwd::insert_at_end<poptype::gamete_t,poptype::glist_t>,
				    std::placeholders::_1,std::placeholders::_2),
			  
			  // Fitness function
			  std::bind(snowdrift_diploid(),
				    std::placeholders::_1,
				    std::cref(phenotypes),
				    b1, b2, c1, c2),
			  
			  // Mutation removal from gamete policy: never remove
			  []( mlist_t::iterator & mit ) { return false; } 
			  );
  
  KTfwd::remove_lost(&pop.mutations, &pop.mut_lookup);

  return wbar;
}

using namespace boost::python;

// Calculate phenotypes
boost::python::list phenotypes(const poptype & pop)
{
  boost::python::list phenos;
  
  for( auto & dip : pop.diploids ) 
    { 
      phenos.append(phenotypef_1locus(dip)); 
    }

  return phenos;
}

//Now, we can expose the stuff to python
BOOST_PYTHON_MODULE(social_evol)
{
  //Expose the type based on fwdpp's "sugar" layer
  class_<poptype>("poptype",init<unsigned>())
    .def("clear",&poptype::clear)
    ;
  //Expose the GSL wrapper
  class_<GSLrng>("GSLrng",init<unsigned>())
    ;
  //Expose the function to run the model
  def("pop_init",pop_init);
  def("evolve_step",evolve_step);
  def("phenotypes", phenotypes);
}
