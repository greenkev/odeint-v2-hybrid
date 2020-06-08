/*
 [auto_generated]
 boost/numeric/odeint/integrate/integrate_until.hpp

 [begin_description]
 Integration of ODEs with a stopping condition
 [end_description]

 Copyright 2009-2011 Karsten Ahnert
 Copyright 2009-2011 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */


#ifndef BOOST_NUMERIC_ODEINT_INTEGRATE_INTEGRATE_UNTIL_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_INTEGRATE_INTEGRATE_UNTIL_HPP_INCLUDED

#include <boost/numeric/odeint/stepper/runge_kutta_dopri5.hpp>
#include <boost/numeric/odeint/stepper/controlled_runge_kutta.hpp>
#include <boost/numeric/odeint/integrate/null_observer.hpp>
#include <boost/numeric/odeint/integrate/detail/integrate_until.hpp>

namespace boost {
namespace numeric {
namespace odeint {

template< class System , class State , class Time , class Event, class Observer >
Time integrate_until( System system , State &start_state , Time start_time , Time end_time , Time dt , Event event, Observer observer )
{
	typedef boost::numeric::odeint::runge_kutta_dopri5<State> base_stepper;
	typedef typename boost::numeric::odeint::result_of::make_dense_output< base_stepper >::type dense_stepper_type;
	dense_stepper_type stepper = make_dense_output( 1.0e-6 , 1.0e-6 , base_stepper() );

	return integrate_until( stepper , system , start_state , start_time , end_time , dt , event, observer );
	return 0;
}

template< class System , class State , class Time , class Event >
Time integrate_until( System system , State &start_state , Time start_time , Time end_time , Time dt, Event event )
{
	return integrate_until(system, start_state, start_time, end_time, dt, event, null_observer() );
}

// template<class System, class State, class Time, class Stepper, class Event, class Observer>
// Time integrate_until( Stepper stepper , System system , State &start_state ,
//         Time start_time , Time end_time , Time dt , Event event , Observer observer ) {

// }

/**
 * \fn integrate_until( System system , State &start_state , Time start_time , Time end_time , Time dt , Event event, Observer observer )
 * \brief Integrates the ODE.
 *
 * Integrates the ODE given by system from start_time either to end_time or until event returns a value below zero
 * with start_state as initial condition and dt as initial time step.
 * This function uses a dense output dopri5 stepper and performs an adaptive
 * integration with step size control, thus dt changes during the integration.
 * This method uses standard error bounds of 1E-6.
 * After each step, the observer and event is called.
 * The precise location of the event is determined using the dense output and bisection.
 *
 * \param system The system function to solve, hence the r.h.s. of the
 * ordinary differential equation.
 * \param start_state The initial state.
 * \param start_time Start time of the integration.
 * \param end_time End time of the integration.
 * \param dt Initial step size, will be adjusted during the integration.
 * \param event The event function which will halt the integration.
 * \param observer Observer that will be called after each time step.
 * \return The time of the event
 */

/**
 * \fn integrate_until( System system , State &start_state , Time start_time , Time end_time , Time dt , Event event )
 * \brief Integrates the ODE.
 *
 * Integrates the ODE given by system from start_time either to end_time or until event returns a value below zero
 * with start_state as initial condition and dt as initial time step.
 * This function uses a dense output dopri5 stepper and performs an adaptive
 * integration with step size control, thus dt changes during the integration.
 * This method uses standard error bounds of 1E-6.
 * After each step, the observer and event is called.
 * The precise location of the event is determined using the dense output and bisection.
 *
 * \param system The system function to solve, hence the r.h.s. of the
 * ordinary differential equation.
 * \param start_state The initial state.
 * \param start_time Start time of the integration.
 * \param end_time End time of the integration.
 * \param dt Initial step size, will be adjusted during the integration.
 * \param event The event function which will halt the integration.
 * \return The time of the event
 */

} // namespace odeint
} // namespace numeric
} // namespace boost

#endif // BOOST_NUMERIC_ODEINT_INTEGRATE_INTEGRATE_UNTIL_HPP_INCLUDED
