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

#include "RayleighFannoFlow.h"

/**
 * This function defines the valid parameters for
 * this Kernel and their default values
 */
template<>
InputParameters validParams<RayleighFannoFlow>()
{
  InputParameters params = validParams<ODEKernel>();
    params.addParam<Real>("friction", 0., "constant friction parameter.");
    params.addParam<Real>("Dh", 1., "constant hydraulic diameter parameter.");
    params.addRequiredParam<UserObjectName>("eos", "The name of equation of state object to use.");
  return params;
}

RayleighFannoFlow::RayleighFannoFlow(const std::string & name, InputParameters parameters) :
    // You must call the constructor of the base class first
    ODEKernel(name, parameters),
    _f(getParam<Real>("friction")),
    _Dh(getParam<Real>("Dh")),
    _eos(getUserObject<EquationOfState>("eos"))
{
}

Real
RayleighFannoFlow::computeQpResidual()
{
    
    // Compute M^2 and M^3:
    Real _M2 = _u[_i]*_u[_i];
    Real _M3 = _u[_i]*_u[_i]*_u[_i];
    
    // Return the value:
    return _eos.gamma()*_M3*(1+0.5*(_eos.gamma()-1)*_M2)*_f/(_Dh*(1-_M2));
}

Real
RayleighFannoFlow::computeQpJacobian()
{
    return 0.;
}

Real
RayleighFannoFlow::computeQpOffDiagJacobian(unsigned int jvar)
{
    return 0.;
}
