/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                               */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "EelPressureBasedVisc.h"

/**
This Kernel computes the pressure based viscosity: Jameson-Schmidt-Turkel (JST), Hassan-Morgan-Peraire (HMP) and Swanson-Turkel (ST).
*/
template<>
InputParameters validParams<EelPressureBasedVisc>()
{
  InputParameters params = validParams<Kernel>();
    // Aux variable:
    params.addRequiredCoupledVar("pressure", "pressure variables");
  return params;
}

EelPressureBasedVisc::EelPressureBasedVisc(const std::string & name,
                       InputParameters parameters) :
  Kernel(name, parameters),
    // Coupled aux variables
    _grad_press(coupledGradient("pressure"))
{
}

Real EelPressureBasedVisc::computeQpResidual()
{
    return _grad_press[_qp]*_grad_test[_i][_qp];
}

Real EelPressureBasedVisc::computeQpJacobian()
{
  return ( 0 );
}

Real EelPressureBasedVisc::computeQpOffDiagJacobian( unsigned int _jvar)
{ 
    return ( 0 );
}
