/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "ConservativeVariables1DXIC.h"

template<>
InputParameters validParams<ConservativeVariables1DXIC>()
{
  InputParameters params = validParams<InitialCondition>();
  params.addRequiredParam<FunctionName>("area", "function to compute the cross section");
    // Initial conditions:
  params.addRequiredParam<Real>("pressure_init_left", "Initial pressure on the left");
  params.addRequiredParam<Real>("pressure_init_right", "Initial pressure on the right");
  params.addRequiredParam<Real>("vel_init_left", "Initial velocity on the left");
  params.addRequiredParam<Real>("vel_init_right", "Inital velocity on the right");
  params.addRequiredParam<Real>("temp_init_left", "Initial value of the temperature");
  params.addRequiredParam<Real>("temp_init_right", "Initial value of the temperature");
    // Membrane position:
  params.addParam<Real>("membrane", 0.5, "The value of the membrane");
  params.addParam<Real>("length", 0.01, "To smooth the IC over a given length");
    // Equation of state
  params.addRequiredParam<UserObjectName>("eos", "parameters for eos.");
  return params;
}

ConservativeVariables1DXIC::ConservativeVariables1DXIC(const std::string & name,
                     InputParameters parameters) :
    InitialCondition(name, parameters),
    // Function
    _area(getFunction("area")),
	// IC parameters
    _p_left(getParam<Real>("pressure_init_left")),
    _p_right(getParam<Real>("pressure_init_right")),
    _v_left(getParam<Real>("vel_init_left")),
    _v_right(getParam<Real>("vel_init_right")),
    _t_left(getParam<Real>("temp_init_left")),
    _t_right(getParam<Real>("temp_init_right")),
    // Position of the membrane:
    _membrane(getParam<Real>("membrane")),
    _length(getParam<Real>("length")),
  	// User Objects
    _eos(getUserObject<EquationOfState>("eos"))
{}

Real
ConservativeVariables1DXIC::value(const Point & p)
{
// Define and compute parameters used to smooth the initial condition if wished
Real _x1 = _membrane - 0.5 * _length;
Real _x2 = _x1 + _length;
Real _a_p, _b_p, _a_vel, _b_vel, _a_t, _b_t;
_a_p = ( _p_left - _p_right) / ( _x1 - _x2 );
_b_p = ( _x1*_p_right - _x2*_p_left ) / ( _x1 - _x2 );
_a_vel = ( _v_left - _v_right) / ( _x1 - _x2 );
_b_vel = ( _x1*_v_right - _x2*_v_left ) / ( _x1 - _x2 );
_a_t = ( _t_left - _t_right) / ( _x1 - _x2 );
_b_t = ( _x1*_t_right - _x2*_t_left ) / ( _x1 - _x2 );
// Get the name of the variable this object acts on
std::string _name_var = _var.name();
// Compute the pressure, velocity and temperature values
Real _pressure = 0.;
Real _temp = 0.;
Real _vel = 0.;
  if ( p(0) <= _x1 )
	{
	 _pressure = _p_left;
         _vel = _v_left;
         _temp = _t_left;
	}
  else if ( p(0) >= _x2 )
	{
	 _pressure = _p_right;
         _vel = _v_right;
         _temp = _t_right;
	}
  else 
	{
        _pressure = ( _a_p * p(0) + _b_p );
        _vel = ( _a_vel * p(0) + _b_vel );
        _temp = ( _a_t * p(0) + _b_t );
	}
// Compute the conservative variables
Real _density = (_pressure + _eos.Pinf()) / (_eos.Cv()*(_eos.gamma()-1)*_temp);
Real _int_energy = (_pressure+_eos.gamma()*_eos.Pinf())/(_density*(_eos.gamma()-1)) + _eos.qcoeff();
Real _tot_energy = _density*(_int_energy + 0.5*_vel*_vel);
// Value of the area:
    Real _A = _area.value(0., p);
// Return the value of the initial condition. Identify the name of the variable
// Density: rhoA
if ( _name_var == "rhoA" )
	return ( _density*_A );
// Momentum: rhouA and rhovA
else if ( _name_var == "rhouA" )
{
	return ( _density*_vel*_A);
}
else if ( _name_var == "rhovA" )
{
	return ( _density*_vel*_A);
}
// total energy: rhoEA
else if ( _name_var == "rhoEA" )
	return ( _tot_energy*_A );
else
	return 0;
}
