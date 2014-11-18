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
This function computes the total energy
**/
#include "TotalEnergyAux.h"

template<>
InputParameters validParams<TotalEnergyAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("rhoEA", "rhoEA");
  params.addRequiredCoupledVar("area", "area");
  return params;
}

TotalEnergyAux::TotalEnergyAux(const std::string & name, InputParameters parameters) :
    AuxKernel(name, parameters),
  // Coupled variables
    _rhoEA( coupledValue("rhoEA")),
  // Coupled Aux Variables
    _area(coupledValue("area"))
{}

Real
TotalEnergyAux::computeValue()
{
  return _rhoEA[_qp] / _area[_qp] ;
}
