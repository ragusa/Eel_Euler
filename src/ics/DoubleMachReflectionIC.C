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

#include "DoubleMachReflectionIC.h"

template<>
InputParameters validParams<DoubleMachReflectionIC>()
{
    InputParameters params = validParams<InitialCondition>();
    params.addRequiredParam<FunctionName>("area", "function to compute the cross section");
    // Initial conditions:
    params.addRequiredParam<Real>("pressure_init_left", "Initial pressure on the left");
    params.addRequiredParam<Real>("pressure_init_right", "Initial pressure on the right");
    params.addRequiredParam<Real>("vel_x_init_left", "Initial x velocity on the left");
    params.addRequiredParam<Real>("vel_x_init_right", "Inital x velocity on the right");
    params.addRequiredParam<Real>("vel_y_init_left", "Initial y velocity on the left");
    params.addRequiredParam<Real>("vel_y_init_right", "Inital y velocity on the right");
    params.addRequiredParam<Real>("rho_init_left", "Initial value of the temperature");
    params.addRequiredParam<Real>("rho_init_right", "Initial value of the temperature");
    // Membrane position:
    params.addRequiredParam<Real>("x_point_source", "Position of the point source: x");
    params.addRequiredParam<Real>("y_point_source", "Position of the point source: y");
    // Equation of state
    params.addRequiredParam<UserObjectName>("eos", "parameters for eos.");
    return params;
}

DoubleMachReflectionIC::DoubleMachReflectionIC(const std::string & name,
                     InputParameters parameters) :
    InitialCondition(name, parameters),
    // Function
    _area(getFunction("area")),
	// IC parameters
    _p_left(getParam<Real>("pressure_init_left")),
    _p_right(getParam<Real>("pressure_init_right")),
    _v_x_left(getParam<Real>("vel_x_init_left")),
    _v_x_right(getParam<Real>("vel_x_init_right")),
    _v_y_left(getParam<Real>("vel_y_init_left")),
    _v_y_right(getParam<Real>("vel_y_init_right")),
    _rho_left(getParam<Real>("rho_init_left")),
    _rho_right(getParam<Real>("rho_init_right")),
    // Position of the membrane:
    _x_pt_source(getParam<Real>("x_point_source")),
    _y_pt_source(getParam<Real>("y_point_source")),
  	// User Objects
    _eos(getUserObject<EquationOfState>("eos"))
{}

Real
DoubleMachReflectionIC::value(const Point & p)
{
// Get the name of the variable this object acts on
std::string _name_var = _var.name();
// Compute the pressure, velocity and temperature values
Real _pressure = 0.;
Real _density = 0.;
Real _vel_x = 0.;
Real _vel_y = 0.;
Real _limit = _x_pt_source + (p(1)+0.5)*std::sqrt(_y_pt_source);
if ( p(0) <= _limit ) {
    _pressure = _p_left;
    _vel_x = _v_x_left;
    _vel_y = _v_y_left;
    _density = _rho_left;
	}
else {
    _pressure = _p_right;
    _vel_x = _v_x_right;
    _vel_y = _v_y_right;
    _density = _rho_right;
	}
// Compute the conservative variables
Real _int_energy = (_pressure+_eos.gamma()*_eos.Pinf())/(_density*(_eos.gamma()-1)) + _eos.qcoeff();
Real _tot_energy = _density*(_int_energy + 0.5*(_vel_x*_vel_x+_vel_y*_vel_y));
// Value of the area:
    Real _A = _area.value(0., p);
// Return the value of the initial condition. Identify the name of the variable
// Density: rhoA
if ( _name_var == "rhoA" )
	return ( _density*_A );
// Momentum: rhouA and rhovA
else if ( _name_var == "rhouA" )
{
	return ( _density*_vel_x*_A);
}
else if ( _name_var == "rhovA" )
{
	return ( _density*_vel_y*_A);
}
// total energy: rhoEA
else if ( _name_var == "rhoEA" )
	return ( _tot_energy*_A );
else
	return 0;
}
