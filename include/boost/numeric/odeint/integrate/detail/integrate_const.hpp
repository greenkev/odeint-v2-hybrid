/*
 [auto_generated]
 boost/numeric/odeint/integrate/detail/integrate_const.hpp

 [begin_description]
 integrate const implementation
 [end_description]

 Copyright 2009-2012 Karsten Ahnert
 Copyright 2009-2012 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef BOOST_NUMERIC_ODEINT_INTEGRATE_DETAIL_INTEGRATE_CONST_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_INTEGRATE_DETAIL_INTEGRATE_CONST_HPP_INCLUDED

#include <boost/numeric/odeint/util/unwrap_reference.hpp>
#include <boost/numeric/odeint/stepper/stepper_categories.hpp>
#include <boost/numeric/odeint/util/unit_helper.hpp>
// #include <boost/math/tools/roots.hpp>
#include <boost/numeric/odeint/integrate/detail/integrate_adaptive.hpp>

#include <boost/numeric/odeint/util/detail/less_with_sign.hpp>

namespace boost {
namespace numeric {
namespace odeint {
namespace detail {

// For performing bisection search for event locations
template<class Stepper, class State, class Time, class Event>
struct EventFunctor {
	Stepper &stepper;
	State &state;
	typename boost::numeric::odeint::unwrap_reference< Event >::type &event;

	EventFunctor(Stepper &stepper, State &state, Event event) : stepper(stepper), state(state), event(event) {
	}

	inline Time operator()(Time t) {
		stepper.calc_state( t, state );
		return event(state, t);
	}
};


template<class T>
struct tolerance {
	T tol;
	tolerance(T tol) : tol(tol) {}

	bool operator()(T &a, T &b) {
		return std::abs<T>(a - b) < tol;
	}
};

template<class T>
tolerance<T> make_tolerance(T tol) {
	return tolerance<T>(tol);
}

// forward declaration
template< class Stepper , class System , class State , class Time , class Observer >
size_t integrate_adaptive(
        Stepper stepper , System system , State &start_state ,
        Time &start_time , Time end_time , Time &dt ,
        Observer observer , controlled_stepper_tag
);

// Basic Stepper, No Event Function
template< class Stepper , class System , class State , class Time , class Observer >
size_t integrate_const(
        Stepper stepper , System system , State &start_state ,
        Time start_time , Time end_time , Time dt ,
        Observer observer , stepper_tag 
)
{
    typename odeint::unwrap_reference< Observer >::type &obs = observer;
    
    Time time = start_time;
    int step = 0;
    
    while( less_eq_with_sign( time+dt , end_time , dt ) )
    {
        obs( start_state , time );
        stepper.do_step( system , start_state , time , dt );
        // direct computation of the time avoids error propagation happening when using time += dt
        // we need clumsy type analysis to get boost units working here
        ++step;
        time = start_time + static_cast< typename unit_value_type<Time>::type >(step) * dt;
    }
    obs( start_state , time );

    return step;
}

// Controlled Stepper, No Event Fucntion
template< class Stepper , class System , class State , class Time , class Observer >
size_t integrate_const(
        Stepper stepper , System system , State &start_state ,
        Time start_time , Time end_time , Time dt ,
        Observer observer , controlled_stepper_tag 
)
{
    typename odeint::unwrap_reference< Observer >::type &obs = observer;
    
    Time time = start_time;
    const Time time_step = dt;
    int step = 0;
    
    while( less_eq_with_sign( time+time_step , end_time , dt ) )
    {
        obs( start_state , time );
        detail::integrate_adaptive( stepper , system , start_state , time , time+time_step , dt ,
                                    null_observer() , controlled_stepper_tag() );
        // direct computation of the time avoids error propagation happening when using time += dt
        // we need clumsy type analysis to get boost units working here
        ++step;
        time = start_time + static_cast< typename unit_value_type<Time>::type >(step) * time_step;
    }
    obs( start_state , time );
    
    return step;
}

// Dense Stepper, No Event Fucntion
template< class Stepper , class System , class State , class Time , class Observer >
size_t integrate_const(
        Stepper stepper , System system , State &start_state ,
        Time start_time , Time end_time , Time dt ,
        Observer observer , dense_output_stepper_tag 
)
{
    typename odeint::unwrap_reference< Observer >::type &obs = observer;
    
    Time time = start_time;
    
    stepper.initialize( start_state , time , dt );
    obs( start_state , time );
    time += dt;

    int obs_step( 1 );
    int real_step( 0 );
    
    while( less_with_sign( time+dt , end_time , dt ) )
    {
        while( less_eq_with_sign( time , stepper.current_time() , dt ) )
        {
            stepper.calc_state( time , start_state );
            obs( start_state , time );
            ++obs_step;
            // direct computation of the time avoids error propagation happening when using time += dt
            // we need clumsy type analysis to get boost units working here
            time = start_time + static_cast< typename unit_value_type<Time>::type >(obs_step) * dt;
        }
        // we have not reached the end, do another real step
        if( less_with_sign( stepper.current_time()+stepper.current_time_step() ,
                            end_time ,
                            stepper.current_time_step() ) )
        {
            while( less_eq_with_sign( stepper.current_time() , time , dt ) )
            {
                stepper.do_step( system );
                ++real_step;
            }
        }
        else if( less_with_sign( stepper.current_time() , end_time , stepper.current_time_step() ) )
        { // do the last step ending exactly on the end point
            stepper.initialize( stepper.current_state() , stepper.current_time() , end_time - stepper.current_time() );
            stepper.do_step( system );
            ++real_step;
        }
        
    }
    // last observation, if we are still in observation interval
    if( less_eq_with_sign( time , end_time , dt ) )
    {
        stepper.calc_state( time , start_state );
        obs( start_state , time );
    }
    
    return real_step;
}

// Basic Stepper, With Event Function; Not Implemented throws informative error
template< class Stepper , class System , class State , class Time , class Event, class Observer >
size_t integrate_const(
        Stepper stepper , System system , State &start_state ,
        Time start_time , Time end_time , Time dt , Event event,
        Observer observer , stepper_tag 
)
{
    throw std::runtime_error("Attempted to call an unimplemented integrate_const function.\nBasic and Error stepping with event function not supported.");
    // typename odeint::unwrap_reference< Observer >::type &obs = observer;
	// typename boost::numeric::odeint::unwrap_reference< Event >::type &evt = event;
    
    // Time time = start_time;
    // int step = 0;
    
    // while( less_eq_with_sign( time+dt , end_time , dt ) )
    // {
    //     obs( start_state , time );
    //     stepper.do_step( system , start_state , time , dt );
    //     // direct computation of the time avoids error propagation happening when using time += dt
    //     // we need clumsy type analysis to get boost units working here

    //     // Check if the event occured
    //     if( evt( stepper.current_state(), time ) <= 0 ) {
    //         EventFunctor<Stepper, State, Time, Event> f(stepper, start_state, event);
    //         std::pair<Time, Time> root;
    //         root = boost::math::tools::bisect(f, start_time + static_cast< typename unit_value_type<Time>::type >(step-1) * dt,
    //                                           time, make_tolerance( 1e-10 ) );

    //         // Step back
    //         stepper.initialize( stepper.current_state(), stepper.current_time(), (root.first + root.second)/2 - stepper.current_time() );
    //         stepper.do_step(system);

    //         // Now return
    //         obs(stepper.current_state(), stepper.current_time() );
    //         // overwrite start_state with the final point
    //         boost::numeric::odeint::copy( stepper.current_state() , start_state );
    //         return stepper.current_time();
    //     }

    //     ++step;
    //     time = start_time + static_cast< typename unit_value_type<Time>::type >(step) * dt;
    // }
    // obs( start_state , time );

    // return step;
}

// Controlled Stepper, With Event Fucntion:  Not Implemented
template< class Stepper , class System , class State , class Time , class Event, class Observer >
size_t integrate_const(
        Stepper stepper , System system , State &start_state ,
        Time start_time , Time end_time , Time dt , Event event,
        Observer observer , controlled_stepper_tag 
)
{
    throw std::runtime_error("Attempted to call an unimplemented integrate_const function.\n Controlled stepping with event function not supported.");
    // typename odeint::unwrap_reference< Observer >::type &obs = observer;
    
    // Time time = start_time;
    // const Time time_step = dt;
    // int step = 0;
    
    // while( less_eq_with_sign( time+time_step , end_time , dt ) )
    // {
    //     obs( start_state , time );
    //     detail::integrate_adaptive( stepper , system , start_state , time , time+time_step , dt ,
    //                                 null_observer() , controlled_stepper_tag() );
    //     // direct computation of the time avoids error propagation happening when using time += dt
    //     // we need clumsy type analysis to get boost units working here
    //     ++step;
    //     time = start_time + static_cast< typename unit_value_type<Time>::type >(step) * time_step;
    // }
    // obs( start_state , time );
    
    // return step;
}

// Dense Stepper, With Event Fucntion:  Not Implemented
template< class Stepper , class System , class State , class Time , class Event, class Observer >
size_t integrate_const(
        Stepper stepper , System system , State &start_state ,
        Time start_time , Time end_time , Time dt , Event event,
        Observer observer , dense_output_stepper_tag 
)
{
    // throw std::runtime_error("Attempted to call an unimplemented integrate_const function.\n Dense stepping with event function not supported.");
    typename odeint::unwrap_reference< Observer >::type &obs = observer;
	typename odeint::unwrap_reference< Event >::type &evt = event;
    
    Time time = start_time; //This always represent the 'next' observation point

	Time t_old;
    
    stepper.initialize( start_state , time , dt );
    obs( start_state , time );
    time += dt;

    int obs_step( 1 );
    int real_step( 0 );
    
    //Loop as long as the next observation will not be past the ending time
    while( less_with_sign( time+dt , end_time , dt ) )
    {
        //Run observations until the next observation is ahead of the integrator
        while( less_eq_with_sign( time , stepper.current_time() , dt ) )
        {
            stepper.calc_state( time , start_state );
            obs( start_state , time );
            ++obs_step;
            // direct computation of the time avoids error propagation happening when using time += dt
            // we need clumsy type analysis to get boost units working here
            time = start_time + static_cast< typename unit_value_type<Time>::type >(obs_step) * dt;
        }
        // we have not reached the end, do another real step
        if( less_with_sign( stepper.current_time()+stepper.current_time_step() ,
                            end_time ,
                            stepper.current_time_step() ) )
        {
            // real step until we need a new observation
            while( less_eq_with_sign( stepper.current_time() , time , dt ) )
            {
			    t_old = stepper.current_time();
                stepper.do_step( system );
                // std::cout << "Normal Step" << std::endl;
                ++real_step;

                // Did the event occur?
                if( evt( stepper.current_state(), stepper.current_time() ) <= 0 ) {
                    EventFunctor<Stepper, State, Time, Event> f(stepper, start_state, event);
                    std::pair<Time, Time> root;
                    root = boost::math::tools::bisect(f, t_old, stepper.current_time(), make_tolerance( 1e-10 ) );

                    //Run observations until the next observation is ahead of the event
                    while( less_eq_with_sign( time ,  (root.first + root.second)/2 , dt ) )
                    {
                        stepper.calc_state( time , start_state );
                        obs( start_state , time );
                        ++obs_step;
                        // direct computation of the time avoids error propagation happening when using time += dt
                        // we need clumsy type analysis to get boost units working here
                        time = start_time + static_cast< typename unit_value_type<Time>::type >(obs_step) * dt;
                    }
                    // std::cout << "stepper.current_time() " << stepper.current_time() 
                                // << "\troot.first " << root.first
                                // << "\troot.second " << root.second
                                // << "\ttime of obs" << time << std::endl;
                    // Step backward to the event bisection
                    stepper.initialize( stepper.current_state(), stepper.current_time(), (root.first + root.second)/2 - stepper.current_time() );
                    stepper.do_step(system);
                    // std::cout << "Step to Event" << std::endl;

                    // Now return
                    obs(stepper.current_state(), stepper.current_time() );
                    // overwrite start_state with the final point
                    boost::numeric::odeint::copy( stepper.current_state() , start_state );
                    return real_step;
                }
            }
        }
        else if( less_with_sign( stepper.current_time() , end_time , stepper.current_time_step() ) )
        { // do the last step ending exactly on the end point
            stepper.initialize( stepper.current_state() , stepper.current_time() , end_time - stepper.current_time() );
			t_old = stepper.current_time();
            stepper.do_step( system );

            // std::cout << "Final Step to Time Endind" << std::endl;
            ++real_step;
            // Did the event occur?
            if( evt( stepper.current_state(), stepper.current_time() ) <= 0 ) {
                EventFunctor<Stepper, State, Time, Event> f(stepper, start_state, event);
                std::pair<Time, Time> root;
                root = boost::math::tools::bisect(f, t_old, stepper.current_time(), make_tolerance( 1e-10 ) );

                //Run observations until the next observation is ahead of the event
                while( less_eq_with_sign( time ,  (root.first + root.second)/2 - stepper.current_time() , dt ) )
                {
                    stepper.calc_state( time , start_state );
                    obs( start_state , time );
                    ++obs_step;
                    // direct computation of the time avoids error propagation happening when using time += dt
                    // we need clumsy type analysis to get boost units working here
                    time = start_time + static_cast< typename unit_value_type<Time>::type >(obs_step) * dt;
                }

                // Step backward to the event bisection
                stepper.initialize( stepper.current_state(), stepper.current_time(), (root.first + root.second)/2 - stepper.current_time() );
                stepper.do_step(system);

                // Now return
                obs(stepper.current_state(), stepper.current_time() );
                // overwrite start_state with the final point
                boost::numeric::odeint::copy( stepper.current_state() , start_state );
                return real_step;
            }
        }
        
    }
    // last observation, if we are still in observation interval
    if( less_eq_with_sign( time , end_time , dt ) )
    {
        stepper.calc_state( time , start_state );
        obs( start_state , time );
    }
    
    return real_step;
}

} } } }

#endif
