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

#include "EelMass.h"

/**
This Kernel computes the convection flux of the continuity equation :
rho*u*A where A is the area of the geometry.
*/
template<>
InputParameters validParams<EelMass>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredCoupledVar("rhouA_x", "x component of rhouA");
  params.addCoupledVar("rhouA_y", "y component of rhouA");
  params.addCoupledVar("rhouA_z", "z component of rhouA");
  return params;
}

EelMass::EelMass(const std::string & name,
                       InputParameters parameters) :
  Kernel(name, parameters),
    // Coupled aux variables
    _rhouA_x(coupledValue("rhouA_x")),
    _rhouA_y(_mesh.dimension()>=2 ? coupledValue("rhouA_y") : _zero ),
    _rhouA_z(_mesh.dimension()==3 ? coupledValue("rhouA_z") : _zero ),
    // Parameters for jacobian
    _rhouA_x_nb(coupled("rhouA_x")),
    _rhouA_y_nb(isCoupled("rhouA_y") ? coupled("rhouA_y") : -1),
    _rhouA_z_nb(isCoupled("rhouA_z") ? coupled("rhouA_z") : -1)
{}

Real EelMass::computeQpResidual()
{
    // Compute convective part of the continuity equation:
    RealVectorValue _conv(_rhouA_x[_qp], _rhouA_y[_qp], _rhouA_z[_qp]);
    
    // Return the total expression for the continuity equation:
    return -_conv * _grad_test[_i][_qp];
}

Real EelMass::computeQpJacobian()
{
  return ( 0 );
}

Real EelMass::computeQpOffDiagJacobian( unsigned int _jvar)
{
    if (_jvar == _rhouA_x_nb)
        return -_phi[_j][_qp]*_grad_test[_i][_qp](0);
    else if (_jvar == _rhouA_y_nb && _mesh.dimension() >= 2)
        return -_phi[_j][_qp]*_grad_test[_i][_qp](1);
    else if (_jvar == _rhouA_z_nb && _mesh.dimension() == 3)
        return -_phi[_j][_qp]*_grad_test[_i][_qp](2);
    else
        return 0.;
}
