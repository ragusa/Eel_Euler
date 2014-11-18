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
This function aims at computing the Temperature at the nodes and its gradient. This auxiliary variable is coupled
to rho_bar, m_bar and E_bar defined as the product of the usual density, momentum and energy, and the cross section
A computed by the function AreaFunction.
**/
#include "TemperatureAux.h"

template<>
InputParameters validParams<TemperatureAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("pressure", "pressure");
  params.addRequiredCoupledVar("density", "density: rho");;
  params.addRequiredParam<UserObjectName>("eos", "The parameters for the eos.");
  return params;
}

TemperatureAux::TemperatureAux(const std::string & name, InputParameters parameters) :
    AuxKernel(name, parameters),
  // Coupled variables
   _pressure(coupledValue("pressure")),
   _rho(coupledValue("density")),
  // User Objects for eos
   _eos(getUserObject<EquationOfState>("eos"))
{}

/** Overwrite the function computevalue() */
Real
TemperatureAux::computeValue()
{
    return _eos.temperature_from_p_rho(_pressure[_qp], _rho[_qp]);
}
