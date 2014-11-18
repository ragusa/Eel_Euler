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

#include "NodalMinMultipleValues.h"

#include <algorithm>
#include <limits>

template<>
InputParameters validParams<NodalMinMultipleValues>()
{
  InputParameters params = validParams<NodalVariablePostprocessor>();
    // Variable:
    params.addRequiredParam<VariableName>("variable", "The name of the variable that this object operates on");
    // Output type
    params.addRequiredParam<std::string>("output_type", "Output type: rho*c*c, rho*c*vel or rho*vel*vel.");
    // Conservative variables:
    params.addRequiredCoupledVar("rhoA", "rho*A");
    params.addRequiredCoupledVar("rhouA_x", "alpha*rho*u*A");
    params.addCoupledVar("rhouA_y", "rho*v*A");
    params.addRequiredCoupledVar("rhoEA", "rho*E*A");
    // Auxkernel variable:
    params.addRequiredCoupledVar("area", "area");
    // Equation of state:
    params.addRequiredParam<UserObjectName>("eos", "Equation of state");
  return params;
}

NodalMinMultipleValues::NodalMinMultipleValues(const std::string & name, InputParameters parameters) :
  NodalVariablePostprocessor(name, parameters),
    _value(-std::numeric_limits<Real>::max()),
    // Function Mach number:
    _output_name(getParam<std::string>("output_type")),
    _output_type("RHOVEL2, RHOCVEL, RHOC2, INVALID", _output_name),
    // Variable:
    _var(_subproblem.getVariable(_tid, parameters.get<VariableName>("variable"))),
    // Conservative variables
    _rhoA(coupledValue("rhoA")),
    _rhouA_x(coupledValue("rhouA_x")),
    _rhouA_y(_mesh.dimension()>=2 ? coupledValue("rhouA_y"): _zero),
    _rhoEA(coupledValue("rhoEA")),
    // Auxkernel variable:
    _area(coupledValue("area")),
    // User Objects for eos
    _eos(getUserObject<EquationOfState>("eos"))
{}

void
NodalMinMultipleValues::initialize()
{
  _value = -std::numeric_limits<Real>::max();
}

void
NodalMinMultipleValues::execute()
{
    // Compute the pressure from the equation of state and the conservative variables:
    Real rho = _rhoA[_qp] / _area[_qp];
    RealVectorValue vel(_rhouA_x[_qp]/_rhoA[_qp], _rhouA_y[_qp]/_rhoA[_qp]);
    Real rhoE = _rhoEA[_qp] / _area[_qp];
    Real pressure = _eos.pressure(rho, vel.size(), rhoE);

    // Compute the speed of sound:
    Real c2 = _eos.c2_from_p_rho(rho, pressure);
    
    // Return value:
    Real value_temp = 0.;
    switch (_output_type)
    {
        case RHOVEL2:
            value_temp = rho * vel.size() * vel.size();
            break;
        case RHOCVEL:
            value_temp = rho * std::sqrt(c2) * vel.size();
            break;
        case RHOC2:
            value_temp = rho * c2;
            break;
        default:
            mooseError("The output type is not supporter by this function.");
            break;
    }
    _value = std::max(_value, -value_temp);
}

Real
NodalMinMultipleValues::getValue()
{
  gatherMax(_value);
  return -_value;
}

void
NodalMinMultipleValues::threadJoin(const UserObject & y)
{
  const NodalMinMultipleValues & pps = static_cast<const NodalMinMultipleValues &>(y);
  _value = std::max(_value, pps._value);
}
