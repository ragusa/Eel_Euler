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
/**
This function computes the pressure. It is dimension agnostic.
**/
#include "PressureAux.h"

template<>
InputParameters validParams<PressureAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("rhoA", "rhoA");
  params.addRequiredCoupledVar("rhouA_x", "rhouA_x");
  params.addCoupledVar("rhouA_y", "rhouA_y");
  params.addCoupledVar("rhouA_z", "rhouA_z");
  params.addRequiredCoupledVar("rhoEA", "rhoEA");
  params.addRequiredCoupledVar("area", "area");
  params.addRequiredParam<UserObjectName>("eos", "Equation of state");
  return params;
}

PressureAux::PressureAux(const std::string & name, InputParameters parameters) :
    AuxKernel(name, parameters),
    // Coupled variables
    _rhoA(coupledValue("rhoA")),
    _rhouA_x(coupledValue("rhouA_x")),
    _rhouA_y(_mesh.dimension()>=2 ? coupledValue("rhouA_y"): _zero),
    _rhouA_z(_mesh.dimension()==3 ? coupledValue("rhouA_z"): _zero),
    _rhoEA(coupledValue("rhoEA")),
    _area(coupledValue("area")),
    // User Objects for eos
    _eos(getUserObject<EquationOfState>("eos"))

{}

Real
PressureAux::computeValue()
{
    // Computes the density, the norm of the velocity and the total energy:
    Real _rho = _rhoA[_qp] / _area[_qp];
    Real _rhoE = _rhoEA[_qp] / _area[_qp];
    Real _vel_x = _rhouA_x[_qp] / _rhoA[_qp];
    Real _vel_y = _rhouA_y[_qp] / _rhoA[_qp];
    Real _vel_z = _rhouA_z[_qp] / _rhoA[_qp];
    Real _norm_vel = std::sqrt(_vel_x*_vel_x + _vel_y*_vel_y + _vel_z*_vel_z);
    
    // Computes the pressure
    return _eos.pressure(_rho, _norm_vel, _rhoE);
}
