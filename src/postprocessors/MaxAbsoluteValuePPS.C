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

#include "MaxAbsoluteValuePPS.h"
/* This pps computes the maximum absolute value at the quadrature points (infinite norm). */
template<>
InputParameters validParams<MaxAbsoluteValuePPS>()
{
  InputParameters params = validParams<ElementPostprocessor>();
    params.addRequiredCoupledVar("variable", "variable this pps is acting on.");
  return params;
}

MaxAbsoluteValuePPS::MaxAbsoluteValuePPS(const std::string & name, InputParameters parameters) :
    ElementPostprocessor(name, parameters),
    _u(coupledValue("variable")),
    _value(-std::numeric_limits<Real>::max())
{}

void
MaxAbsoluteValuePPS::initialize()
{
  _value = -std::numeric_limits<Real>::max();
}

void
MaxAbsoluteValuePPS::execute()
{
    // Loop over quadrature points:
    Real _local_max = 0.;
    for (int _qp=0; _qp < _qrule->n_points(); _qp++) {
        _local_max = std::max(std::fabs(_u[_qp]), _local_max);
    }
    _value = std::max(_value, _local_max);
}

void
MaxAbsoluteValuePPS::finalize()
{
    gatherMax(_value);
}

Real
MaxAbsoluteValuePPS::getValue()
{
    gatherMax(_value);
    return _value;
}

void
MaxAbsoluteValuePPS::threadJoin(const UserObject & y)
{
  ElementPostprocessor::threadJoin(y);
  const MaxAbsoluteValuePPS & pps = dynamic_cast<const MaxAbsoluteValuePPS &>(y);
  _value = std::max(_value, pps._value);
}
