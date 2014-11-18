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

#include "AreaFunction.h"

template<>
InputParameters validParams<AreaFunction>()
{
  InputParameters params = validParams<Function>();
  params.addParam<Real>("left", 0., "Value of the left bound");
  params.addParam<Real>("length", 1.,  "Length of the variable area");
  params.addParam<Real>("Ao", "Coefficient for the function");
  params.addParam<Real>("Bo", "Coefficient for the function");
  return params;
}

AreaFunction::AreaFunction(const std::string & name, InputParameters parameters) :
    Function(name, parameters),
    _left(getParam<Real>("left")),
    _length(getParam<Real>("length")),
    _Ao(getParam<Real>("Ao")),
    _Bo(getParam<Real>("Bo"))
{}

Real
AreaFunction::value(Real /*t*/, const Point & p)
{
  //std::cout << " p(0): " << p(0) << std::endl;
  Real _right; _right = _length + _left;
  //std::cout << "right: " << _right << std::endl;
  //std::cout << "left: " << _left << std::endl;
  if ( _Ao != 0) {
      if ( p(0) < _left )
      {
        //std::cout << 1 << std::endl;
          return _Ao;
      }
      else if ( p(0) >= _left && p(0) <= _right )
      {
          //std::cout << 2 << std::endl;
          //return _Ao * (1 - 0.5*sin(libMesh::pi*(p(0)-_left)/_length) ) + _Bo;
          return _Ao * ( 1 + 0.5*cos(2*libMesh::pi*p(0))) + _Bo;
      }
      else
      {
          //std::cout << 3 << std::endl;
          return _Ao;
      }
  }
  else
	return _Bo;
}

RealVectorValue 
AreaFunction::gradient(Real /*t*/, const Point & p)
{
	return 0;
}
