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

#include "ElementMaxDuDtValue.h"

template<>
InputParameters validParams<ElementMaxDuDtValue>()
{
  InputParameters params = validParams<ElementPostprocessor>();
    // Coupled aux variable:
    params.addRequiredCoupledVar("variable", "variable this object is acting on.");
    params.addCoupledVar("variable2", "second variable that multiplies the time derivative.");
  return params;
}

ElementMaxDuDtValue::ElementMaxDuDtValue(const std::string & name, InputParameters parameters) :
    ElementPostprocessor(name, parameters),
    _var(coupledValue("variable")),
    _var_old(coupledValueOld("variable")),
    _var2(isCoupled("variable2") ? coupledValue("variable2") : _zero),
    _value(-std::numeric_limits<Real>::max())
{}

void
ElementMaxDuDtValue::initialize()
{
  _value = -std::numeric_limits<Real>::max();
}

void
ElementMaxDuDtValue::execute()
{
    // Loop over quadrature points:
    Real _local_max = 0.;
    for (int _qp=0; _qp < _qrule->n_points(); _qp++)
    {
        Real var2 = isCoupled("variable2") ? _var2[_qp] : 1.;
        _local_max = std::max( var2*std::fabs(_var[_qp] - _var_old[_qp]) / std::fabs(_var[_qp]), _local_max);
//        _local_max = std::max( var2*std::fabs(_var[_qp] - _var_old[_qp]) / (_dt * std::fabs(_var[_qp])), _local_max);
    }
    
    _value = std::max(_value, _local_max);
}

void
ElementMaxDuDtValue::finalize()
{
    gatherMax(_value);
}

Real
ElementMaxDuDtValue::getValue()
{
    gatherMax(_value);
    return _value;
}

void
ElementMaxDuDtValue::threadJoin(const UserObject & y)
{
  ElementPostprocessor::threadJoin(y);
  const ElementMaxDuDtValue & pps = dynamic_cast<const ElementMaxDuDtValue &>(y);
  _value = std::max(_value, pps._value);
}
