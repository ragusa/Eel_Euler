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
This function compute the Mach number. It is dimension agnostic.
**/
#include "MachNumberAux.h"

template<>
InputParameters validParams<MachNumberAux>()
{
  InputParameters params = validParams<AuxKernel>();
    // Conservative coupled variables:
    params.addRequiredCoupledVar("rhoA", "rhoA");
    params.addRequiredCoupledVar("rhouA_x", "x component of momentum");
    params.addCoupledVar("rhouA_y", "y component of momentum");
    params.addCoupledVar("rhouA_z", "z component of momentum");
    // Aux coupled variables:
    params.addRequiredCoupledVar("pressure", "pressure");
    params.addRequiredCoupledVar("area", "area");
    params.addRequiredParam<UserObjectName>("eos", "Equation of state");
  return params;
}

MachNumberAux::MachNumberAux(const std::string & name, InputParameters parameters) :
    AuxKernel(name, parameters),
    // Coupled variables
    _rhoA(coupledValue("rhoA")),
    _rhouA_x(coupledValue("rhouA_x")),
    _rhouA_y(_mesh.dimension()>=2 ? coupledValue("rhouA_y") : _zero),
    _rhouA_z(_mesh.dimension()==3 ? coupledValue("rhouA_z") : _zero),
    // Aux coupled variables:
    _pressure(coupledValue("pressure")),
    _area(coupledValue("area")),
    // User Objects for eos
    _eos(getUserObject<EquationOfState>("eos"))
{}

Real
MachNumberAux::computeValue()
{
    // Compute the norm of velocity and density:
    Real _rho = _rhoA[_qp] / _area[_qp];
    Real _vel_x = _rhouA_x[_qp] / _rhoA[_qp];
    Real _vel_y = _rhouA_y[_qp] / _rhoA[_qp];
    Real _vel_z = _rhouA_z[_qp] / _rhoA[_qp];
    Real _norm_vel2 = _vel_x*_vel_x + _vel_y*_vel_y + _vel_z*_vel_z;
    
    // Computes the speed of sounds:
    Real _c2 = _eos.c2_from_p_rho(_rho, _pressure[_qp]);
//    std::cout<<"c2="<<_c2<<std::endl;
//    std::cout<<"pressure="<<_pressure[_qp]<<std::endl;

    // Return the value of the Mach number:
//    std::cout<<"Mach="<<std::sqrt(_norm_vel2 / _c2)<<std::endl;
    return std::sqrt(_norm_vel2 / _c2);
}
