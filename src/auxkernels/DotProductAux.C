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
#include "DotProductAux.h"

template<>
InputParameters validParams<DotProductAux>()
{
  InputParameters params = validParams<AuxKernel>();
    // Vector 1:
    params.addRequiredCoupledVar("x_component", "x component of vector 1");
    params.addCoupledVar("y_component", "y component of vector 1");
    params.addCoupledVar("z_component", "z component of vector 1");
    // Vector 1:
    params.addRequiredCoupledVar("vector", "second vector for dot product: will take the gradient of it.");
  return params;
}

DotProductAux::DotProductAux(const std::string & name, InputParameters parameters) :
    AuxKernel(name, parameters),
    // Vector 1:
    _x_comp(coupledValue("x_component")),
    _y_comp(_mesh.dimension()>=2 ? coupledValue("y_component") : _zero),
    _z_comp(_mesh.dimension()==3 ? coupledValue("z_component") : _zero),
    // Vector 2:
    _grad_vector(coupledGradient("vector"))
{}

Real
DotProductAux::computeValue()
{
    // Initialyze and normalize vector:
    RealVectorValue _x_vec(_x_comp[_qp], _y_comp[_qp], _z_comp[_qp]);
    // Return the dot product between vectors 1 and 2:
    return _x_vec*_grad_vector[_qp]/_grad_vector[_qp].size();
}
