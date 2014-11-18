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
This function computes the fluid internal energy 'rhoe' from the conservative variables. It is dimension agnostic.
**/
#include "InternalEnergyAux.h"

template<>
InputParameters validParams<InternalEnergyAux>()
{
  InputParameters params = validParams<AuxKernel>();
    /// Coupled variables
  params.addRequiredCoupledVar("rhoA", "fluid density: rho*A");
  params.addRequiredCoupledVar("rhouA_x", "fluid x momentum component");
  params.addCoupledVar("rhouA_y", "fluid x momentum component");
  params.addCoupledVar("rhouA_z", "fluid x momentum component");
  params.addRequiredCoupledVar("rhoEA", "rhoEA");
    /// Coupled aux variables
  params.addRequiredCoupledVar("area", "area");  
  return params;
}

InternalEnergyAux::InternalEnergyAux(const std::string & name, InputParameters parameters) :
    AuxKernel(name, parameters),
    _rhoA(coupledValue("rhoA")),
    _rhouA_x(coupledValue("rhouA_x")),
    _rhouA_y(_mesh.dimension()>=2 ? coupledValue("rhouA_y") : _zero),
    _rhouA_z(_mesh.dimension()==3 ? coupledValue("rhouA_z") : _zero),
    _rhoEA(coupledValue("rhoEA")),
    _area(coupledValue("area"))
{}

Real
InternalEnergyAux::computeValue()
{
    
    /// Compute density, norm of velocity and total energy:
    Real _rhoE = _rhoEA[_qp] / _area[_qp];
    Real _vel_x = _rhouA_x[_qp] / _rhoA[_qp];
    Real _vel_y = _rhouA_y[_qp] / _rhoA[_qp];
    Real _vel_z = _rhouA_z[_qp] / _rhoA[_qp];
    Real _rho = _rhoA[_qp] / _area[_qp];
    Real _norm_vel2 = _vel_x*_vel_x + _vel_y*_vel_y + _vel_z*_vel_z;
    
    /// Return internal energy:
    return _rhoE - 0.5*_rho*_norm_vel2;
}
