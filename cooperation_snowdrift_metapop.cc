#include <limits>
#include <algorithm>
#include <cmath>
#include <boost/python.hpp>
#include <boost/container/list.hpp>
#include <boost/container/vector.hpp>
#include <boost/pool/pool_alloc.hpp>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#include <boost/functional/hash.hpp>

// Include custom mutation header
#include<mutation_inf_alleles.hpp>

// Main fwdpp library header
#include <fwdpp/diploid.hh>
// Include the necessary "sugar" components
#include <fwdpp/sugar/metapop/metapop.hpp>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/sugar/serialization.hpp>
// GSL discrete distributions
#include <fwdpp/internal/gsl_discrete.hpp>

/*
  We will use a gsl_rng_mt19937 as our RNG.
  This type is implicitly convertible to gsl_rng *,
  and auto-handles the gsl_rng_free steps, etc.
 */
using GSLrng = KTfwd::GSLrng_t<KTfwd::GSL_RNG_MT19937>;
using lookup_t = KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr;

//FWDPP-related types
using mtype = KTfwd::mutation;
using mlist_t = boost::container::list<mtype,boost::pool_allocator<mtype> >;
using gamete_t = KTfwd::gamete_base<mtype,mlist_t>;
using glist_t = boost::container::list<gamete_t, boost::pool_allocator<gamete_t>>;

// Custom diploid genotype
struct diploid_t : public KTfwd::tags::custom_diploid_t
{
  using first_type = glist_t::iterator;
  using second_type = glist_t::iterator;
  first_type first;
  second_type second;
  unsigned i;
  //constructors, etc.
  diploid_t() : first(first_type()),second(second_type()),i(std::numeric_limits<unsigned>::max()) {}
  diploid_t(first_type g1, first_type g2) : first(g1),second(g2),i(std::numeric_limits<unsigned>::max()) {}
};

//Define our our population type via KTfwd::sugar 
using poptype = KTfwd::sugar::metapop_serialized<mtype,
						 KTfwd::mutation_writer,
						 KTfwd::mutation_reader<mtype>,
						 mlist_t,
						 glist_t,
						 boost::container::vector<diploid_t>,
						 boost::container::vector<glist_t>,
						 boost::container::vector<boost::container::vector<diploid_t>>,
						 boost::container::vector<mtype>,
						 boost::container::vector<unsigned>,
						 boost::unordered_set<double,boost::hash<double>,
								      KTfwd::equal_eps>>;

/*! \brief Continuous Snowdrift Game from Doebeli, Hauert, and Killingback (2004, Science, 306:859--862)
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
  inline result_type operator()( gsl_rng * r,
				 const poptype::dipvector_t::const_iterator dip,
				 const std::vector<double> & phenotypes, 
				 const double & b1, const double & b2,
				 const double & c1, const double & c2) const
  {
    unsigned N = phenotypes.size();
    double zself = phenotypes[dip->i];
    result_type fitness = 0;

    // pick random partner
    unsigned partner = gsl_rng_uniform_int(r, N);
    while (partner == dip->i)
      {
	partner = gsl_rng_uniform_int(r, N);
      }
	   
    double zpair = zself+phenotypes[partner];
    // Payoff function from Fig 1
    fitness += 1 + b1*zpair + b2*zpair*zpair - c1*zself - c2*zself*zself;
    
    return std::max(fitness, 0.0);
  }
};

// Truncated Laplace distribution using rejection method (inefficient)
double ran_trunc_laplace( gsl_rng * r, const double & mean, const double & scale, const double & lo, const double & hi )
{
  double rv = mean + gsl_ran_laplace(r, scale);
  while (rv > hi || rv < lo)
    {
      rv = mean + gsl_ran_laplace(r, scale);
    }
  return rv;
}

struct inf_alleles
{
  using result_type = mtype;
  mtype operator()( poptype::gamete_t & g,
		    gsl_rng * r,
		    poptype::lookup_table_t * lookup,
		    const double & scale, const double & lo, const double & hi) const
  {
    auto mitr = g.smutations.begin();

    // Make sure mutation list is not emppty. This requires proper initialization of the population.
    assert(mitr!=g.smutations.end());

    double mutval = ran_trunc_laplace(r,(*mitr)->s,scale,lo,hi);
    while(lookup->find(mutval) != lookup->end())
      {
	mutval = ran_trunc_laplace(r,(*mitr)->s,scale,lo,hi);
      }

    lookup->insert(mutval);
    return mtype(mutval,mutval,1);
  }
};

// Calculate phnenotype for sinlge "locus"
double phenotypef_1locus( const diploid_t & dip )
{
  gamete_t::mcont_iterator m1 = dip.first->smutations.begin();
  gamete_t::mcont_iterator m2 = dip.second->smutations.begin();

  double pval = ((*m1)->h * (*m1)->s + (*m2)->h * (*m2)->s) / ((*m1)->h + (*m2)->h);
  return pval;
}

// Effective "haploid" phenotype: just use the first gamete
double phenotypef_1locus_haploid( const diploid_t & dip )
{
  gamete_t::mcont_iterator m1 = dip.first->smutations.begin();

  double pval = (*m1)->s;
  return pval;
}

// migration for island model with hard selection.
// returns vector of lookup tables, one for each deme
std::vector<lookup_t> mig_island_hard(const double & m, std::vector<double> & wtot)
{
  unsigned n = wtot.size();
  std::vector<double> migrants(n, 0.0);
  std::vector<lookup_t> lookups;

  // Migrants from each deme
  for ( unsigned d = 0 ; d < n ; ++d )
    {
      for ( unsigned i = 0 ; i < n ; ++i )
	{
	  migrants[i] = ((1 - m) * (i == d) + m * (i != d)) * wtot[i];
	}
      lookups.emplace_back(lookup_t(gsl_ran_discrete_preproc(n, &migrants[0])));
    }

  return lookups;
}

// Evolve population one generation using infinite alleles mutation model
std::vector<double> evolve_step( GSLrng & rng,
				 poptype & pop,
				 const double & mu, const double & scale,
				 const double & m,
				 const double & b1, const double & b2, const double & c1, const double & c2)
{
  unsigned n = pop.Ns.size();
  std::vector<unsigned>::iterator maxNit = std::max_element(pop.Ns.begin(), pop.Ns.end());
  std::vector<std::vector<double>> phenotypes(n,
					      std::vector<double>(*maxNit));
  //Fill phenotypes
  for ( unsigned i = 0; i != n; i++ ) 
    {
      for ( unsigned j = 0; j != pop.Ns[i]; j++ )
	{
	  pop.diploids[i][j].i = j; 
	  phenotypes[i][j] = phenotypef_1locus(pop.diploids[i][j]); 
	}
    }

  // vector of fitness functions, one for each deme
  std::vector<std::function<double (poptype::dipvector_t::const_iterator)> > vff;
  for ( unsigned i = 0; i != n; i++ ) 
    {
      vff.push_back(std::bind(snowdrift_diploid(),
			      rng,
			      std::placeholders::_1,
			      std::cref(phenotypes[i]),
			      b1, b2, c1, c2));
    }

  // Calculate total fitness (dispersing offspring) for hard selection
  std::vector<double> wtot(n, 0.0);
  for ( unsigned i = 0; i != n; i++ ) 
    {
      for ( auto dptr = pop.diploids[i].begin() ; dptr != pop.diploids[i].end() ; ++dptr )
	{
	  wtot[i] += vff[i](dptr);
	}
    }

  // "migration policy" for island model and hard selection
  std::vector<lookup_t> lookups = mig_island_hard(m, wtot);
  
  std::vector<double> wbars =
    KTfwd::sample_diploid(rng,
			  &pop.gametes,
			  &pop.diploids,
			  &pop.mutations,
			  &pop.Ns[0],
			  mu,
			  // infinite alleles mutation model where mutations have
			  // truncated laplace distribution on (0,1)
			  std::bind(inf_alleles(),
				    std::placeholders::_1,
				    rng,
				    &pop.mut_lookup,
				    scale,
				    0.,
				    1.),
			  
			  // "null" recombination policy
			  [](poptype::glist_t::iterator & ,   
			     poptype::glist_t::iterator & ,
			     poptype::glist_t *){return 0;},

			  // Mutation insertion policy
			  std::bind(KTfwd::insert_at_end<poptype::mutation_t,
				    poptype::mlist_t>,
				    std::placeholders::_1,std::placeholders::_2),

			  // Gamete insertion policy
			  std::bind(KTfwd::insert_at_end<poptype::gamete_t,poptype::glist_t>,
				    std::placeholders::_1,std::placeholders::_2),
			  
			  // Vector of fitness functions
			  vff,
			  
			  // Mutation removal from gamete policy: never remove
			  []( mlist_t::iterator &  ) { return false; },

			  // Migration policy
			  [&rng, &lookups](const size_t & dest_deme)
			  { return gsl_ran_discrete(rng,lookups[dest_deme].get()); }
			  );
  
  KTfwd::remove_lost(&pop.mutations, &pop.mut_lookup);

  return wbars;
}

using namespace boost::python;

// Calculate phenotypes
boost::python::list phenotypes(const poptype & pop)
{
  boost::python::list phenos;

  for ( auto & deme : pop.diploids )
    {
      boost::python::list plist;
      for ( auto & dip : deme ) 
	{ 
	  plist.append(phenotypef_1locus(dip));
	}
      phenos.append(plist);
    }

  return phenos;
}

// Return list of mutations
boost::python::list mutations(const poptype & pop)
{
  boost::python::list mutations;
  gamete_t::mcont_iterator m1, m2;

  for ( auto & deme : pop.diploids )
    {
      boost::python::list mlist;
      for( auto & dip : deme ) 
	{
	  m1 = dip.first->smutations.begin();
	  m2 = dip.second->smutations.begin();

	  mlist.append((*m1)->s);
	  mlist.append((*m2)->s);
	}
      mutations.append(mlist);
    }

  return mutations;
}

// Initialize population to be homogeneous with mutation value mutval0
poptype pop_init( const boost::python::list & pyNs,
		  const double & mutval0)
{
  std::vector<unsigned> Ns = boost::python::extract<std::vector<unsigned>>(pyNs);
  unsigned n = Ns.size();

  // total population size
  unsigned nN = 1;
  for (auto & Nval : Ns)
    {
      nN *= Nval;
    }
  
  poptype pop(&Ns[0], n); // each of n demes of size N in vector Ns is initially single gamete in 2*N copies
  mtype mut0(mutval0, mutval0, 2*nN); // initial mutation in 2*nN copies
  pop.mut_lookup.insert(mutval0); // insert mutation value into lookup table
  
  // insert mutation into population mutation list
  poptype::mlist_t::iterator mitr = pop.mutations.insert(pop.mutations.end(), std::move(mut0));

  // insert mutation into gamete mutation list in each deme
  for (auto & deme_glist : pop.gametes)
    {
      // insert mutation in mutation list of gamete
      poptype::glist_t::iterator g0 = deme_glist.begin();
      g0->smutations.insert( g0->smutations.begin(), mitr );
    }

  return pop;
}

// // Initialize population to given distribution of mutation values
// poptype pop_init_dist( const unsigned & N,
// 		       const boost::python::list & mutvals)
// {
//   double mutval;
//   poptype pop(N); // population initially is single gamete in 2*N copies
  
//   // set initial gamete to zero count since it will be replaced
//   poptype::glist_t::iterator g0 = pop.gametes.begin();
//   g0->n = 0;

//   unsigned i = 0;
//   boost::unordered_map<double, poptype::glist_t::iterator, boost::hash<double>, KTfwd::equal_eps> new_gametes;
//   // loop over all diploids creating a new gametes and new mutations
//   // first gamete in diploid
//   for ( auto & dip : pop.diploids ) 
//     {
//       mutval = boost::python::extract<double>(mutvals[i]); i++; // convert mutation value
//       auto search = new_gametes.find(mutval);
//       if (search == new_gametes.end()) // gamete with mutval hasn't been seen yet
// 	{
// 	  poptype::gamete_t g(1); // new gamete
// 	  mtype mut(mutval, mutval, 1); // create new mutation object
// 	  pop.mut_lookup.insert(mutval); // save mutation values in lookup table

// 	  // insert new gametes into population gamete list
// 	  poptype::glist_t::iterator gitr = pop.gametes.insert(pop.gametes.end(), std::move(g));
      
// 	  // insert mutation into population mutation list
// 	  poptype::mlist_t::iterator mitr = pop.mutations.insert(pop.mutations.end(), std::move(mut));

// 	  // insert mutation in mutation list of gamete
// 	  gitr->smutations.insert( gitr->smutations.begin(), mitr );

// 	  // insert gamete into map and assign to diploid
// 	  new_gametes[mutval] = gitr;
// 	  dip.first  = gitr;
// 	}
//       else // gamete has been seen
// 	{
// 	  auto gam = search->second;
// 	  auto mut = gam->smutations.begin();

// 	  gam->n++; // increment gamete
// 	  (*mut)->n++; // increment mutation

// 	  dip.first = gam;
// 	}
//     }
//   // second gamete
//   for ( auto & dip : pop.diploids ) 
//     {
//       mutval = boost::python::extract<double>(mutvals[i]); i++; // convert mutation value
//       auto search = new_gametes.find(mutval);
//       if (search == new_gametes.end()) // gamete with mutval hasn't been seen yet
// 	{
// 	  poptype::gamete_t g(1); // new gamete
// 	  mtype mut(mutval, mutval, 1); // create new mutation object
// 	  pop.mut_lookup.insert(mutval); // save mutation values in lookup table

// 	  // insert new gametes into population gamete list
// 	  poptype::glist_t::iterator gitr = pop.gametes.insert(pop.gametes.end(), std::move(g));
      
// 	  // insert mutation into population mutation list
// 	  poptype::mlist_t::iterator mitr = pop.mutations.insert(pop.mutations.end(), std::move(mut));

// 	  // insert mutation in mutation list of gamete
// 	  gitr->smutations.insert( gitr->smutations.begin(), mitr );

// 	  // insert gamete into map and assign to diploid
// 	  new_gametes[mutval] = gitr;
// 	  dip.second  = gitr;
// 	}
//       else // gamete has been seen
// 	{
// 	  auto gam = search->second;
// 	  auto mut = gam->smutations.begin();

// 	  gam->n++; // increment gamete
// 	  (*mut)->n++; // increment mutation

// 	  dip.second = gam;
// 	}
//     }
  
//   return pop;
// }

//Now, we can expose the stuff to python
BOOST_PYTHON_MODULE(cooperation_snowdrift)
{
  //Expose the type based on fwdpp's "sugar" layer
  class_<poptype>("poptype",init<unsigned *, size_t>())
    .def("clear",&poptype::clear)
    ;
  //Expose the GSL wrapper
  class_<GSLrng>("GSLrng",init<unsigned>())
    ;
  //Expose the function to run the model
  def("pop_init",pop_init);
  // def("pop_init_dist",pop_init_dist);
  def("evolve_step",evolve_step);
  def("phenotypes", phenotypes);
  def("mutations", mutations);
}

