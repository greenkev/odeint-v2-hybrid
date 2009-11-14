/*
 boost header: numeric/odeint/integrator_constant_step.hpp

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 Copyright 2009 Andre Bergner

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_NUMERIC_ODEINT_INTEGRATOR_CONSTANT_STEP_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_INTEGRATOR_CONSTANT_STEP_HPP_INCLUDED

#include <boost/numeric/odeint/observer.hpp>

namespace boost {
namespace numeric {
namespace odeint {


    template< 
        class Stepper ,
        class DynamicalSystem ,
        class Observer
        >
    size_t integrate(
	Stepper &stepper ,
	DynamicalSystem &system ,
	typename Stepper::time_type start_time ,
	typename Stepper::time_type dt ,
	typename Stepper::container_type &state ,
	typename Stepper::time_type end_time ,
	Observer &observer
	)
    {
	size_t iteration = 0;
        while( start_time < end_time )
	{
	    observer( start_time , state , system );
            stepper.next_step( system , state , start_time , dt );
            start_time += dt;
	    ++iteration;
        }
	observer( start_time , state , system );

	return iteration;
    }



    template<
        class Stepper ,
        class DynamicalSystem
        >
    size_t integrate(
	Stepper &stepper ,
	DynamicalSystem &system ,
	typename Stepper::time_type start_time ,
	typename Stepper::time_type dt ,
	typename Stepper::container_type &state ,
	typename Stepper::time_type end_time
	)
    {
	return integrate(
	    stepper , system , start_time , dt , state , end_time ,
	    do_nothing_observer<
	    typename Stepper::time_type ,
	    typename Stepper::container_type ,
	    DynamicalSystem >
	    );
    }



    template<
	class Stepper , 
        class DynamicalSystem ,
        class Observer
        >
    typename Stepper::time_type integrate_steps(
	Stepper &stepper ,
	DynamicalSystem &system ,
	typename Stepper::time_type start_time ,
	typename Stepper::time_type dt ,
	typename Stepper::container_type &state ,
	size_t num_of_steps ,
	Observer &observer
	)
    {
	size_t iteration = 0;
        while( iteration < num_of_steps )
	{
	    observer( start_time , state , system );
            stepper.next_step( system , state , start_time , dt );
            start_time += dt;
	    ++iteration;
        }
	observer( start_time , state , system );

	return start_time;
    }


    template<
	class Stepper , 
        class DynamicalSystem
	>
    typename Stepper::time_type integrate_steps(
	Stepper &stepper ,
	DynamicalSystem &system ,
	typename Stepper::time_type start_time ,
	typename Stepper::time_type dt ,
	typename Stepper::container_type &state ,
	size_t num_of_steps 
	)
    {
	return integrate_steps(
	    stepper , system , start_time , dt , state , num_of_steps ,
	    do_nothing_observer<
	    typename Stepper::time_type ,
	    typename Stepper::container_type ,
	    DynamicalSystem >
	    );
    }

    
    
   

} // odeint
} // numeric
} // boost

#endif //BOOST_NUMERIC_ODEINT_INTEGRATOR_CONSTANT_STEP_HPP_INCLUDED