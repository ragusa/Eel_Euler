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

#include "AreaFunction2D.h"

template<>
InputParameters validParams<AreaFunction2D>()
{
  InputParameters params = validParams<Function>();
  params.addParam<Real>("Ao", "Coefficient for the function");
  params.addParam<Real>("Bo", "Coefficient for the function");
  return params;
}

AreaFunction2D::AreaFunction2D(const std::string & name, InputParameters parameters) :
    Function(name, parameters),
    _Ao(getParam<Real>("Ao")),
    _Bo(getParam<Real>("Bo"))
{}

Real
AreaFunction2D::value(Real /*t*/, const Point & p)
{
    return _Ao * ( 1 + 0.5*cos(2*libMesh::pi*p(0)) ) + _Bo;
}

RealVectorValue 
AreaFunction2D::gradient(Real /*t*/, const Point & p)
{
	return 0;
}
