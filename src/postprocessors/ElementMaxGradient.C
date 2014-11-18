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

#include "ElementMaxGradient.h"

template<>
InputParameters validParams<ElementMaxGradient>()
{
  InputParameters params = validParams<ElementPostprocessor>();
    // Coupled aux variable:
    params.addRequiredCoupledVar("pressure", "pressure of the fluid.");
  return params;
}

ElementMaxGradient::ElementMaxGradient(const std::string & name, InputParameters parameters) :
    ElementPostprocessor(name, parameters),
    _grad_press(coupledGradient("pressure")),
    _value(-std::numeric_limits<Real>::max())
{}

void
ElementMaxGradient::initialize()
{
  _value = -std::numeric_limits<Real>::max();
}

void
ElementMaxGradient::execute()
{
    // Loop over quadrature points:
    Real _local_max = 0.;
    //std::cout << "nb quad=" << _qrule->n_points() << std::endl;
    for (int _qp=0; _qp < _qrule->n_points(); _qp++) {
        //std::cout << "qp=" << _qp << std::endl;
        //Real _norm_grad_press2 = _grad_press[_qp](0)*_grad_press[_qp](0) + _grad_press[_qp](1)*_grad_press[_qp](1) + _grad_press[_qp](2)*_grad_press[_qp](2);
        _local_max = std::max(std::fabs(_grad_press[_qp](0)), _local_max);
    }
    //std::cout << "local=" << _local_max << std::endl;
    _value = std::max(_value, _local_max);
    //std::cout<<"execute="<<_value<<std::endl;
}

void
ElementMaxGradient::finalize()
{
    //std::cout<<"finalize="<<_value<<std::endl;
    gatherMax(_value);
    //std::cout<<"finalize="<<_value<<std::endl;
}

Real
ElementMaxGradient::getValue()
{
    gatherMax(_value);
    //std::cout<<"getValue="<<_value<<std::endl;
    return _value;
}

void
ElementMaxGradient::threadJoin(const UserObject & y)
{
  ElementPostprocessor::threadJoin(y);
  const ElementMaxGradient & pps = dynamic_cast<const ElementMaxGradient &>(y);
  _value = std::max(_value, pps._value);
}
