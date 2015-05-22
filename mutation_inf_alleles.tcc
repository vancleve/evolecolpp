//  -*- C++ -*- 
#ifndef _MUTATION_TCC_
#define _MUTATION_TCC_

#include <type_traits>
#include <algorithm>
#include <numeric>
#include <vector>

#include <gsl/gsl_randist.h>

#if defined(HAVE_BOOST_VECTOR) && !defined(USE_STANDARD_CONTAINERS)
#include <boost/container/vector.hpp>
#endif

namespace KTfwd
{
  // Overload mutate_gametes for infinite alleles mutations
  template< typename iterator_type,
	    typename mutation_model,
	    typename mutation_insertion_policy,
	    typename gamete_insertion_policy,
	    typename list_type_allocator,
	    typename list_type_allocator2,
	    template<typename,typename> class list_type,
	    template<typename,typename> class list_type2>
  iterator_type mutate_gamete( gsl_rng * r,
			       const double & mu,
			       list_type< typename iterator_type::value_type,list_type_allocator > * gametes,
			       list_type2<typename iterator_type::value_type::mutation_type,list_type_allocator2 > * mutations, 
			       iterator_type & g,
			       const mutation_model & mmodel,
			       const mutation_insertion_policy & mpolicy,
			       const gamete_insertion_policy & gpolicy)
  {
    assert( g != gametes->end() );
    unsigned newmut = gsl_ran_flat(r,0,1) < mu;
    if ( newmut )
      {
	assert( g->n > 0 );
	g->n--;
	typename iterator_type::value_type ng(1); // create new gamete in a single copy with no mutations

	// generate new mutation
	auto m = mmodel(g); 

	// add new mutation to mutation list
	auto mitr = mpolicy(std::move(m),mutations);

	// add new mutation to gamete mutation list
	if (m.neutral)
	  {
	    ng.mutations.insert( ng.mutations.end(), mitr );
	  }
	else
	  {
	    ng.smutations.insert( ng.mutations.end(), mitr );
	  }

	// add new gamete to gamete list
	return gpolicy(std::move(ng),gametes);
      }
    return g;
  }
}

struct inf_alleles
{
  template< typename iterator_type,
	    typename mlist_t,
	    typename lookup_table_t,
	    typename mutation_model>  
  typename mlist_t::value_type operator()( gsl_rng * r,
					   lookup_table_t * lookup,
					   iterator_type & g,
					   const mutation_model & mmodel) const
  {
    double mutval;
    typename iterator_type::value_type::mcont_iterator mitr = g->smutations.begin();
    
    mutval = mmodel(r, (*mitr)->s);
    while(lookup->find(mutval) != lookup->end())
      {
	mutval = mmodel(r, (*mitr)->s);
      }
    
    lookup->insert(mutval);
    return typename mlist_t::value_type(mutval,mutval,1);
  }
};
#endif /* _MUTATION_TCC_ */
