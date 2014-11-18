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

#include "EelFannoFlow.h"

/**
This Kernel computes the convection flux of the continuity equation :
rho*u*A where A is the area of the geometry.
*/
template<>
InputParameters validParams<EelFannoFlow>()
{
  InputParameters params = validParams<Kernel>();
    params.addParam<Real>("friction", 0., "constant friction parameter.");
    params.addParam<Real>("Dh", 1., "constant hydraulic diameter parameter.");
    params.addRequiredParam<UserObjectName>("eos", "The name of equation of state object to use.");
  return params;
}

EelFannoFlow::EelFannoFlow(const std::string & name,
                       InputParameters parameters) :
  Kernel(name, parameters),
    _f(getParam<Real>("friction")),
    _Dh(getParam<Real>("Dh")),
    _eos(getUserObject<EquationOfState>("eos"))
{}

Real EelFannoFlow::computeQpResidual()
{
    // Compute M^2 and M^3:
    Real _M2 = _u[_qp]*_u[_qp];
    Real _M3 = _u[_qp]*_u[_qp]*_u[_qp];
    //std::cout<<"M2="<<_M2<<std::endl;
    //std::cout<<"M3="<<_M3<<std::endl;
    //std::cout<<"friction="<<_f<<std::endl;
    //std::cout<<"Dh="<<_Dh<<std::endl;
    // Return the value:
    return -_grad_u[_qp](0) -_eos.gamma()*_M3*(1+0.5*(_eos.gamma()-1)*_M2)*_f/(_Dh*(1-_M2));
}

Real EelFannoFlow::computeQpJacobian()
{
  return ( 0 );
}

Real EelFannoFlow::computeQpOffDiagJacobian( unsigned int _jvar)
{ 
    return ( 0 );
}
