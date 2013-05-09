/*
  [auto_generated]
  boost/numeric/odeint/stepper/extrapolation_stepper.hpp

  [begin_description]
  extrapolation stepper
  [end_description]

  Copyright 2009-2013 Karsten Ahnert
  Copyright 2009-2013 Mario Mulansky

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE_1_0.txt or
  copy at http://www.boost.org/LICENSE_1_0.txt)
*/


#ifndef BOOST_NUMERIC_ODEINT_STEPPER_EXTRAPOLATION_STEPPER_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_STEPPER_EXTRAPOLATION_STEPPER_HPP_INCLUDED


#include <iostream>

#include <algorithm>

#include <boost/config.hpp> // for min/max guidelines
#include <boost/static_assert.hpp>

#include <boost/numeric/odeint/util/bind.hpp>
#include <boost/numeric/odeint/util/unwrap_reference.hpp>

#include <boost/numeric/odeint/stepper/controlled_runge_kutta.hpp>
#include <boost/numeric/odeint/stepper/modified_midpoint.hpp>
#include <boost/numeric/odeint/stepper/controlled_step_result.hpp>
#include <boost/numeric/odeint/algebra/range_algebra.hpp>
#include <boost/numeric/odeint/algebra/default_operations.hpp>
#include <boost/numeric/odeint/algebra/algebra_dispatcher.hpp>
#include <boost/numeric/odeint/algebra/operations_dispatcher.hpp>

#include <boost/numeric/odeint/util/state_wrapper.hpp>
#include <boost/numeric/odeint/util/is_resizeable.hpp>
#include <boost/numeric/odeint/util/resizer.hpp>
#include <boost/numeric/odeint/util/unit_helper.hpp>
#include <boost/numeric/odeint/util/detail/less_with_sign.hpp>

namespace boost {
namespace numeric {
namespace odeint {

template<
    unsigned short Order ,
    class State ,
    class Value = double ,
    class Deriv = State ,
    class Time = Value ,
    class Algebra = typename algebra_dispatcher< State >::algebra_type ,
    class Operations = typename operations_dispatcher< State >::operations_type ,
    class Resizer = initially_resizer
    >
class extrapolation_stepper : public explicit_error_stepper_base<
    extrapolation_stepper< Order , State , Value , Deriv , Time , Algebra , Operations , Resizer > ,
    Order , Order , Order-2 , State , Value , Deriv , Time , Algebra , Operations , Resizer >
{

 private:
    // check for Order being odd
    BOOST_STATIC_ASSERT_MSG( (Order % 2) == 1 , "extrapolation_stepper requires odd Order" );

 public:

    typedef State state_type;
    typedef Value value_type;
    typedef Deriv deriv_type;
    typedef Time time_type;
    typedef Algebra algebra_type;
    typedef Operations operations_type;
    typedef Resizer resizer_type;
#ifndef DOXYGEN_SKIP
    typedef state_wrapper< state_type > wrapped_state_type;
    typedef state_wrapper< deriv_type > wrapped_deriv_type;
    typedef error_stepper_tag stepper_category;

    typedef extrapolation_stepper< Order , State , Value , Deriv , Time , Algebra , Operations , Resizer > stepper_type;
    typedef explicit_error_stepper_base< stepper_type , Order , Order , Order-2 , State , Value , Deriv , Time , Algebra , Operations , Resizer > base_stepper_type;
    
    typedef std::vector< value_type > value_vector;
    typedef std::vector< value_vector > value_matrix;
    typedef std::vector< size_t > int_vector;
    typedef std::vector< wrapped_state_type > state_table_type;
    
#endif //DOXYGEN_SKIP

    typedef unsigned short order_type;
    static const order_type order_value = base_stepper_type::order_value;
    static const order_type stepper_order_value = base_stepper_type::stepper_order_value;
    static const order_type error_order_value = base_stepper_type::error_order_value;
    
    const static size_t m_k_max = (order_value-1)/2;
    
    extrapolation_stepper()
        : m_midpoint() , m_interval_sequence( m_k_max+1 ) , m_coeff( m_k_max+1 ) , m_table( m_k_max )
    {
        for( unsigned short i = 0; i < m_k_max+1; i++ )
        {
            m_interval_sequence[i] = 2 * (i+1);
            m_coeff[i].resize(i);
            for( size_t k = 0 ; k < i ; ++k  )
            {
                const value_type r = static_cast< value_type >( m_interval_sequence[i] ) / static_cast< value_type >( m_interval_sequence[k] );
                m_coeff[i][k] = static_cast< value_type >( 1 ) / ( r*r - static_cast< value_type >( 1 ) ); // coefficients for extrapolation
            }
        }
    }
    
    
    template< class System , class StateIn , class DerivIn , class StateOut , class Err >
    void do_step_impl( System system , const StateIn &in , const DerivIn &dxdt ,
                       time_type t , StateOut &out , time_type dt , Err &xerr )
    {
        // normal step
        do_step_impl( system , in , dxdt , t , out , dt );
        
        // additionally, perform the error calculation
        
    }


    template< class System , class StateIn , class DerivIn , class StateOut >
    void do_step_impl( System system , const StateIn &in , const DerivIn &dxdt ,
                       time_type t , StateOut &out , time_type dt )
    {
        m_resizer.adjust_size( in , detail::bind( &stepper_type::template resize_impl< StateIn > , detail::ref( *this ) , detail::_1 ) );

        // ToDo: better implementation for that?
        if( &in == &out )
        {
            // special care for  in == out
            m_xout_resizer.adjust_size( in , detail::bind( &stepper_type::template resize_m_xout< StateIn > , detail::ref( *this ) , detail::_1 ) );

            size_t k = 0;
            m_midpoint.set_steps( m_interval_sequence[k] );
            m_midpoint.do_step( system , in , dxdt , t , m_xout.m_v , dt );
            for( k = 1 ; k <= m_k_max ; ++k )
            {
                m_midpoint.set_steps( m_interval_sequence[k] );
                m_midpoint.do_step( system , in , dxdt , t , m_table[k-1].m_v , dt );
                extrapolate( k , m_table , m_coeff , m_xout.m_v );
            }
            boost::numeric::odeint::copy( m_xout.m_v , out );
        } else {
            size_t k = 0;
            m_midpoint.set_steps( m_interval_sequence[k] );
            m_midpoint.do_step( system , in , dxdt , t , out , dt );
            for( k = 1 ; k <= m_k_max ; ++k )
            {
                m_midpoint.set_steps( m_interval_sequence[k] );
                m_midpoint.do_step( system , in , dxdt , t , m_table[k-1].m_v , dt );
                extrapolate( k , m_table , m_coeff , out );
            }
        }
    }


    template< class StateIn >
    void adjust_size( const StateIn &x )
    {
        resize_impl( x );
        m_midpoint.adjust_size( x );
    }

private:

    template< class StateIn >
    bool resize_impl( const StateIn &x )
    {
        bool resized( false );
        for( size_t i = 0 ; i < m_k_max ; ++i )
            resized |= adjust_size_by_resizeability( m_table[i] , x , typename is_resizeable<state_type>::type() );
        return resized;
    }

    template< class StateIn >
    bool resize_m_xout( const StateIn &x )
    {
        return adjust_size_by_resizeability( m_xout , x , typename is_resizeable<state_type>::type() );
    }

    

    template< class StateInOut >
    void extrapolate( size_t k , state_table_type &table , const value_matrix &coeff , StateInOut &xest )
    /* polynomial extrapolation, see http://www.nr.com/webnotes/nr3web21.pdf
       uses the obtained intermediate results to extrapolate to dt->0 
    */
    {
        static const value_type val1 = static_cast< value_type >( 1.0 );
        for( int j=k-1 ; j>0 ; --j )
        {
            m_algebra.for_each3( table[j-1].m_v , table[j].m_v , table[j-1].m_v ,
                                 typename operations_type::template scale_sum2< value_type , value_type >( val1 + coeff[k][j] , -coeff[k][j] ) );
        }
        m_algebra.for_each3( xest , table[0].m_v , xest ,
                             typename operations_type::template scale_sum2< value_type , value_type >( val1 + coeff[k][0] , -coeff[k][0]) );
    }


 private:
    modified_midpoint< state_type , value_type , deriv_type , time_type , algebra_type , operations_type , resizer_type > m_midpoint;
    
    algebra_type m_algebra;
    
    resizer_type m_resizer;
    resizer_type m_xout_resizer;
    
    int_vector m_interval_sequence; // stores the successive interval counts
    value_matrix m_coeff;

    wrapped_state_type m_xout;
    state_table_type m_table; // sequence of states for extrapolation

};


} } }
#endif