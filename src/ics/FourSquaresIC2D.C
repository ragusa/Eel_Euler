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

#include "FourSquaresIC2D.h"

template<>
InputParameters validParams<FourSquaresIC2D>()
{
    InputParameters params = validParams<InitialCondition>();
    // Input parameters:
    params.addRequiredParam<FunctionName>("area", "function to compute the cross section");
    params.addRequiredParam<Real>("press_upper_left_corner", "Initial pressure value of the upper left corner");
    params.addRequiredParam<Real>("press_upper_right_corner", "Initial pressure value of the upper right corner");
    params.addRequiredParam<Real>("press_bottom_left_corner", "Initial pressure value of the bottom left corner");
    params.addRequiredParam<Real>("press_bottom_right_corner", "Initial pressure value of the bottom right corner");
    params.addRequiredParam<Real>("dens_upper_left_corner", "Initial density value of the upper left corner");
    params.addRequiredParam<Real>("dens_upper_right_corner", "Initial density value of the upper right corner");
    params.addRequiredParam<Real>("dens_bottom_left_corner", "Initial density value of the bottom left corner");
    params.addRequiredParam<Real>("dens_bottom_right_corner", "Initial density value of the bottom right corner");
    params.addRequiredParam<Real>("x_vel_upper_left_corner", "Initial x-velocity value of the upper left corner");
    params.addRequiredParam<Real>("x_vel_upper_right_corner", "Initial x-velocity value of the upper right corner");
    params.addRequiredParam<Real>("x_vel_bottom_left_corner", "Initial x-velocity value of the bottom left corner");
    params.addRequiredParam<Real>("x_vel_bottom_right_corner", "Initial x-velocity value of the bottom right corner");
    params.addRequiredParam<Real>("y_vel_upper_left_corner", "Initial y-velocity value of the upper left corner");
    params.addRequiredParam<Real>("y_vel_upper_right_corner", "Initial y-velocity value of the upper right corner");
    params.addRequiredParam<Real>("y_vel_bottom_left_corner", "Initial y-velocity value of the bottom left corner");
    params.addRequiredParam<Real>("y_vel_bottom_right_corner", "Initial y-velocity value of the bottom right corner");
    // Coordinates of the common node to all squares:
    params.addParam<Real>("x_node", 0.5, "x coord of the common node to all squares.");
    params.addParam<Real>("y_node", 0.5, "y coord of the common node to all squares.");
    // Parameter to smooth the initial discontinuity:
    params.addParam<Real>("length", 0., "length to smooth the data over.");
    // Equation of state
    params.addRequiredParam<UserObjectName>("eos", "Equation of state");
    return params;
}

FourSquaresIC2D::FourSquaresIC2D(const std::string & name,
                     InputParameters parameters) :
    InitialCondition(name, parameters),
    // Function
    _area(getFunction("area")),
	// Input parameters:
    _press_upper_left_corner(getParam<Real>("press_upper_left_corner")),
    _press_upper_right_corner(getParam<Real>("press_upper_right_corner")),
    _press_bottom_left_corner(getParam<Real>("press_bottom_left_corner")),
    _press_bottom_right_corner(getParam<Real>("press_bottom_right_corner")),
    _dens_upper_left_corner(getParam<Real>("dens_upper_left_corner")),
    _dens_upper_right_corner(getParam<Real>("dens_upper_right_corner")),
    _dens_bottom_left_corner(getParam<Real>("dens_bottom_left_corner")),
    _dens_bottom_right_corner(getParam<Real>("dens_bottom_right_corner")),
    _x_vel_upper_left_corner(getParam<Real>("x_vel_upper_left_corner")),
    _x_vel_upper_right_corner(getParam<Real>("x_vel_upper_right_corner")),
    _x_vel_bottom_left_corner(getParam<Real>("x_vel_bottom_left_corner")),
    _x_vel_bottom_right_corner(getParam<Real>("x_vel_bottom_right_corner")),
    _y_vel_upper_left_corner(getParam<Real>("y_vel_upper_left_corner")),
    _y_vel_upper_right_corner(getParam<Real>("y_vel_upper_right_corner")),
    _y_vel_bottom_left_corner(getParam<Real>("y_vel_bottom_left_corner")),
    _y_vel_bottom_right_corner(getParam<Real>("y_vel_bottom_right_corner")),
    // Coordinate of the common node to all squares:
    _x_node(getParam<Real>("x_node")),
    _y_node(getParam<Real>("y_node")),
    // Length
    _length(getParam<Real>("length")),
    // Equation of state
    _eos(getUserObject<EquationOfState>("eos"))
{}

Real
FourSquaresIC2D::value(const Point & p)
{
    // If statement on the node coordinates:
    Real _p_ic = 0.; Real _rho_ic = 0.;
    Real _x_vel_ic = 0.; Real _y_vel_ic = 0.;
    
    if (p(0) <= _x_node && p(1) >= _y_node) {
        _p_ic = _press_upper_left_corner;
        _rho_ic = _dens_upper_left_corner;
        _x_vel_ic = _x_vel_upper_left_corner;
        _y_vel_ic = _y_vel_upper_left_corner;
    }
    else if (p(0) > _x_node && p(1) >= _y_node ) {
        _p_ic = _press_upper_right_corner;
        _rho_ic = _dens_upper_right_corner;
        _x_vel_ic = _x_vel_upper_right_corner;
        _y_vel_ic = _y_vel_upper_right_corner;
    }
    else if (p(0) <= _x_node && p(1) < _y_node) {
        _p_ic = _press_bottom_left_corner;
        _rho_ic = _dens_bottom_left_corner;
        _x_vel_ic = _x_vel_bottom_left_corner;
        _y_vel_ic = _y_vel_bottom_left_corner;
    }
    else if (p(0) > _x_node && p(1) < _y_node) {
        _p_ic = _press_bottom_right_corner;
        _rho_ic = _dens_bottom_right_corner;
        _x_vel_ic = _x_vel_bottom_right_corner;
        _y_vel_ic = _y_vel_bottom_right_corner;
    }
    
    // Compute the x,y-momentum and the total energy:
    Real _x_rhou_ic = _rho_ic * _x_vel_ic;
    Real _y_rhou_ic = _rho_ic * _y_vel_ic;
    Real _e_ic = _eos.e_from_p_rho(_p_ic, _rho_ic);
    Real _rhoE_ic = _rho_ic*(_e_ic + 0.5*(_x_vel_ic*_x_vel_ic+_y_vel_ic*_y_vel_ic) );
    
    // Get the name of the variable this object acts on
    std::string _name_var = _var.name();
    
    // Value of the area:
    Real _A = _area.value(0., p);
    
    // Return the value:
    // Density
    if ( _name_var == "rhoA" )
        return ( _rho_ic*_A );
    // Momentum: rhouA and rhovA
    else if ( _name_var == "rhouA" )
        return ( _x_rhou_ic*_A);
    else if ( _name_var == "rhovA" )
        return ( _y_rhou_ic*_A);
    // total energy: rhoEA
    else if ( _name_var == "rhoEA" )
        return ( _rhoE_ic*_A );
    else
        mooseError("The variable name is not initialyze by this function.");
        return 0;
}
