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
This function compute the norm of a vector. It is used for LAPIDUS since the gradient of the norm of the velocity vector is required.
**/
#include "NormVectorAux.h"

template<>
InputParameters validParams<NormVectorAux>()
{
  InputParameters params = validParams<AuxKernel>();
    // Conservative coupled variables:
    params.addRequiredCoupledVar("x_component", "x component of vector");
    params.addCoupledVar("y_component", "y component of vector");
    params.addCoupledVar("z_component", "z component of vector");
  return params;
}

NormVectorAux::NormVectorAux(const std::string & name, InputParameters parameters) :
    AuxKernel(name, parameters),
    // Coupled variables
    _x_comp(coupledValue("x_component")),
    _y_comp(isCoupled("y_component") ? coupledValue("y_component") : _zero),
    _z_comp(isCoupled("z_component") ? coupledValue("z_component") : _zero)
{}

Real
NormVectorAux::computeValue()
{
    // Return the norm of the vector:
    return std::sqrt(_x_comp[_qp]*_x_comp[_qp]+_y_comp[_qp]*_y_comp[_qp]+_z_comp[_qp]*_z_comp[_qp]);
}
