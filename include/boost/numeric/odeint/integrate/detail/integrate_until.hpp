/*
 [auto_generated]
 boost/numeric/odeint/integrate/detail/integrate_until.hpp

 [begin_description]
 Default integrate until implementation.
 [end_description]

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/



#ifndef INTEGRATE_UNTIL_HPP_
#define INTEGRATE_UNTIL_HPP_

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/util/unwrap_reference.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/numeric/odeint/util/detail/less_with_sign.hpp>

namespace boost {
namespace numeric {
namespace odeint {

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

template<class System, class State, class Time, class Stepper, class Event, class Observer>
Time integrate_until( Stepper stepper , System system , State &start_state ,
        Time start_time , Time end_time , Time dt , Event event , Observer observer ) {

	using namespace boost::numeric::odeint::detail;
	typename boost::numeric::odeint::unwrap_reference< Observer >::type &obs = observer;
	typename boost::numeric::odeint::unwrap_reference< Event >::type &evt = event;

	size_t count = 0;
	Time t_old;
	stepper.initialize( start_state , start_time , dt );

	while( less_with_sign( stepper.current_time() , end_time , stepper.current_time_step() ) )
	{
		while( less_eq_with_sign( stepper.current_time() + stepper.current_time_step() ,
			   end_time ,
			   stepper.current_time_step() ) )
		{   //make sure we don't go beyond the end_time
			obs( stepper.current_state() , stepper.current_time() );
			t_old = stepper.current_time();
			stepper.do_step( system );

			// Did the event occur?
			if( evt( stepper.current_state(), stepper.current_time() ) <= 0 ) {
				EventFunctor<Stepper, State, Time, Event> f(stepper, start_state, event);
				std::pair<Time, Time> root;
				root = boost::math::tools::bisect(f, t_old, stepper.current_time(), make_tolerance( 1e-10 ) );

				// Step back
				stepper.initialize( stepper.current_state(), stepper.current_time(), (root.first + root.second)/2 - stepper.current_time() );
				stepper.do_step(system);

				// Now return
				obs(stepper.current_state(), stepper.current_time() );
				// overwrite start_state with the final point
				boost::numeric::odeint::copy( stepper.current_state() , start_state );
				return stepper.current_time();
			}

			++count;
		}
		stepper.initialize( stepper.current_state() , stepper.current_time() , end_time - stepper.current_time() );
	}
	obs( stepper.current_state() , stepper.current_time() );
	// overwrite start_state with the final point
	boost::numeric::odeint::copy( stepper.current_state() , start_state );
	return stepper.current_time();
}

} // namespace odeint
} // namespace numeric
} // namespace boost

#endif /* INTEGRATE_UNTIL_HPP_ */
