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

#include "ConservativeVariables2DIC.h"

template<>
InputParameters validParams<ConservativeVariables2DIC>()
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
    params.addRequiredParam<Real>("x_point_source", "Position of the point source: x");
    params.addRequiredParam<Real>("y_point_source", "Position of the point source: y");
    params.addParam<Real>("length", 0.05, "distance from the point source.");
    params.addParam<Real>("smoothing", 0., "parameter for smoothing.");
    // Equation of state
    params.addRequiredParam<UserObjectName>("eos", "parameters for eos.");
    return params;
}

ConservativeVariables2DIC::ConservativeVariables2DIC(const std::string & name,
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
    _x_pt_source(getParam<Real>("x_point_source")),
    _y_pt_source(getParam<Real>("y_point_source")),
    _length(getParam<Real>("length")),
    _smoothing(getParam<Real>("smoothing")),
  	// User Objects
    _eos(getUserObject<EquationOfState>("eos"))
{}

Real
ConservativeVariables2DIC::value(const Point & p)
{
// Get the name of the variable this object acts on
std::string _name_var = _var.name();
    
// Compute r1 and r2 (radius) to smooth data:
Real _r1 = _length - 0.5 * _smoothing;
Real _r2 = _length + 0.5 * _smoothing;
//    std::cout<<"r1="<<_r1<<std::endl;
//    std::cout<<"r2="<<_r2<<std::endl;
    
// Compute parameters for smoothing of initial data:
Real _a_p, _b_p, _a_vel, _b_vel, _a_t, _b_t;
_a_p = ( _p_left - _p_right) / ( _r1 - _r2 );
_b_p = ( _r1*_p_right - _r2*_p_left ) / ( _r1 - _r2 );
_a_vel = ( _v_left - _v_right) / ( _r1 - _r2 );
_b_vel = ( _r1*_v_right - _r2*_v_left ) / ( _r1 - _r2 );
_a_t = ( _t_left - _t_right) / ( _r1 - _r2 );
_b_t = ( _r1*_t_right - _r2*_t_left ) / ( _r1 - _r2 );
//    std::cout<<"a_p="<<_a_p<<std::endl;
//    std::cout<<"b_p="<<_b_p<<std::endl;
//    std::cout<<"a_vel="<<_a_vel<<std::endl;
//    std::cout<<"b_vel="<<_b_vel<<std::endl;
//    std::cout<<"a_t="<<_a_t<<std::endl;
//    std::cout<<"b_t="<<_b_t<<std::endl;
    
// Compute the pressure, velocity and temperature values
Real _pressure = 0.;
Real _temp = 0.;
Real _vel = 0.;
Real _radius = (p(0)-_x_pt_source)*(p(0)-_x_pt_source) + (p(1)-_y_pt_source)*(p(1)-_y_pt_source);
//    std::cout<<"radius="<<_radius<<std::endl;
    
if ( std::sqrt(_radius) <= _r1 ) {
    _pressure = _p_left;
    _vel = _v_left;
    _temp = _t_left;
}
else if ( std::sqrt(_radius) > _r2 ) {
    _pressure = _p_right;
    _vel = _v_right;
    _temp = _t_right;
}
else {
    _pressure = ( _a_p * std::sqrt(_radius) + _b_p );
    _vel = ( _a_vel * std::sqrt(_radius) + _b_vel );
    _temp = ( _a_t * std::sqrt(_radius) + _b_t );
//    std::cout<<"pressure="<<_pressure<<std::endl;
}
//    std::cout<<"pressure="<<_pressure<<std::endl;
//    std::cout<<"vel="<<_vel<<std::endl;
//    std::cout<<"temperature="<<_temp<<std::endl;
    
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
