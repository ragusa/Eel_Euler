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

#include "EelMomentum.h"
/**
This function computes the x, y and z momentum equationS. It is dimension agnostic. 
 */
template<>
InputParameters validParams<EelMomentum>()
{
  InputParameters params = validParams<Kernel>();
    params.addRequiredCoupledVar("rhouA_x", "x component of momentum");
    params.addCoupledVar("rhouA_y", "y component of momentum");
    params.addCoupledVar("rhouA_z", "z component of momentum");
    params.addCoupledVar("rhoEA", "energy: used in the jacobian matrix");
    params.addRequiredCoupledVar("pressure", "pressure");
    params.addCoupledVar("rhoA", "rho*A: only used in the friction term.");
    params.addRequiredCoupledVar("area", "area");
    params.addRequiredParam<UserObjectName>("eos", "Equation of state");
    params.addParam<int>("component", 0, "component of the momentum equation to compute (0,1,2)->(x,y,z)");
    params.addParam<Real>("friction", 0., "friction coefficient for wall friction term.");
    params.addParam<Real>("Dh", 1., "Hydraulic diameter for the friction term.");
    params.addParam<RealVectorValue>("gravity", (0., 0., 0.), "Gravity vector.");
  return params;
}

EelMomentum::EelMomentum(const std::string & name,
                       InputParameters parameters) :
  Kernel(name, parameters),
    // Coupled auxilary variables
    _rhouA_x(coupledValue("rhouA_x")),
    _rhouA_y(_mesh.dimension()>=2 ? coupledValue("rhouA_y") : _zero),
    _rhouA_z(_mesh.dimension()==3 ? coupledValue("rhouA_z") : _zero),
    _pressure(coupledValue("pressure")),
    _rhoA(coupledValue("rhoA")),
    _rhoEA(coupledValue("rhoEA")),
    _area(coupledValue("area")),
    _grad_area(coupledGradient("area")),
    // Equation of state:
    _eos(getUserObject<EquationOfState>("eos")),
    // Parameters:
    _component(getParam<int>("component")),
    _friction(getParam<Real>("friction")),
    _Dh(getParam<Real>("Dh")),
    _gravity(getParam<RealVectorValue>("gravity")),
    // Parameters for jacobian:
    _rhoA_nb(coupled("rhoA")),
    _rhouA_x_nb(coupled("rhouA_x")),
    _rhouA_y_nb(isCoupled("rhouA_y") ? coupled("rhouA_y") : -1),
    _rhouA_z_nb(isCoupled("rhouA_z") ? coupled("rhouA_z") : -1),
    _rhoEA_nb(isCoupled("rhoEA") ? coupled("rhoEA") : -1)


{
    if ( _component > 2 )
        mooseError("ERROR: the integer variable 'component' can only take three values: 0, 1 and 2 that correspond to x, y and z momentum components, respectively.");
    if ( isCoupled("friction") != isCoupled("density") )
        std::cout<<"WARNING: the density variable is only used in the wall friction term. When running a simulation with wall friction, both the friction factor and the density variables have to be supplied in the input file."<<std::endl;
}

Real EelMomentum::computeQpResidual()
{
  // Convection term: _u = rho*vel*vel*A
    RealVectorValue _vector_vel( _rhouA_x[_qp]/_rhoA[_qp], _rhouA_y[_qp]/_rhoA[_qp], _rhouA_z[_qp]/_rhoA[_qp] );
    RealVectorValue _advection = _u[_qp] * _vector_vel;
    
  // Pressure term:
    Real _press = _pressure[_qp]*_area[_qp];
    
  // Source term: P*dA/dx
    Real _PdA = _pressure[_qp]*_grad_area[_qp](_component);
    
  // Wall friction term:
    Real _wall_friction = 0.5 * _friction * _rhoA[_qp] * _vector_vel.size() * _vector_vel(_component) / _Dh;
    
  // Gravity force:
    Real _gravity_force = _gravity(_component) * _rhoA[_qp];
    
  // Return the kernel value:
    return -( _advection*_grad_test[_i][_qp] + _press*_grad_test[_i][_qp](_component) + (_PdA - _wall_friction - _gravity_force)*_test[_i][_qp] );
}

Real EelMomentum::computeQpJacobian()
{
    // Compute the momentum vector rhouA:
    RealVectorValue _rhouA_vec(_rhouA_x[_qp], _rhouA_y[_qp], _rhouA_z[_qp]);
    
    // Compute the velocity vector:
    RealVectorValue _vector_vel( _rhouA_x[_qp]/_rhoA[_qp], _rhouA_y[_qp]/_rhoA[_qp], _rhouA_z[_qp]/_rhoA[_qp] );
    _vector_vel(_component) *= 2.;
    
    // Compute the derivative of \partial_(x,y,z) (AP) - P \partial_(x,y,z) A:
    Real _press_term = _eos.dAp_drhouA(_rhoA[_qp], _u[_qp], _rhoEA[_qp])*(_grad_area[_qp](_component)/_area[_qp]*_test[_i][_qp] + _grad_test[_i][_qp](_component));
    
    // Return the value of the jacobian:
    return -_phi[_j][_qp] * ( _vector_vel * _grad_test[_i][_qp] + _press_term );
}

Real EelMomentum::computeQpOffDiagJacobian( unsigned int _jvar)
{
    // Compute rhouA_vec:
    RealVectorValue _rhouA_vec(_rhouA_x[_qp], _rhouA_y[_qp], _rhouA_z[_qp]);
    
    // Compute the velocity vector:
    RealVectorValue _vector_vel( _rhouA_x[_qp]/_rhoA[_qp], _rhouA_y[_qp]/_rhoA[_qp], _rhouA_z[_qp]/_rhoA[_qp] );
    
    // density (rho*A):
    if (_jvar == _rhoA_nb) {
        Real _press_term = _eos.dAp_drhoA(_rhoA[_qp], _rhouA_vec.size(), _rhoEA[_qp])*(_grad_area[_qp](_component)/_area[_qp]*_test[_i][_qp]+_grad_test[_i][_qp](_component));
        return _phi[_j][_qp] * (_u[_qp] / _rhoA[_qp] * _vector_vel * _grad_test[_i][_qp] - _press_term );
    }
    
    // x-momentum component:
    else if (_jvar == _rhouA_x_nb ) {
        Real _press_term = _eos.dAp_drhouA(_rhoA[_qp], _rhouA_x[_qp], _rhoEA[_qp])*(_grad_area[_qp](_component)/_area[_qp]*_test[_i][_qp]+_grad_test[_i][_qp](_component));
        return -_phi[_j][_qp] * ( _grad_test[_i][_qp](0)*_u[_qp]/_rhoA[_qp] + _press_term );
    }
    // y-momentum component:
    else if (_jvar == _rhouA_y_nb && _mesh.dimension()>=2 ) {
        Real _press_term = _eos.dAp_drhouA(_rhoA[_qp], _rhouA_y[_qp], _rhoEA[_qp])*(_grad_area[_qp](_component)/_area[_qp]*_test[_i][_qp]+_grad_test[_i][_qp](_component));
        return -_phi[_j][_qp] * ( _grad_test[_i][_qp](1)*_u[_qp]/_rhoA[_qp] + _press_term );
    }
    // z-momentum component:
    else if (_jvar == _rhouA_z_nb && _mesh.dimension()==3 ) {
        Real _press_term = _eos.dAp_drhouA(_rhoA[_qp], _rhouA_z[_qp], _rhoEA[_qp])*(_grad_area[_qp](_component)/_area[_qp]*_test[_i][_qp]+_grad_test[_i][_qp](_component));
        return -_phi[_j][_qp] * ( _grad_test[_i][_qp](2)*_u[_qp]/_rhoA[_qp] + _press_term );
    }
    
    // energy (rho*E*A):
    else if (_jvar == _rhoEA_nb) {
        return -_phi[_j][_qp] * _eos.dAp_drhoEA(_rhoA[_qp], _rhouA_vec.size(), _rhoEA[_qp]) * (_grad_area[_qp](_component)/_area[_qp]*_test[_i][_qp]+_grad_test[_i][_qp](_component));
    }
    else
        return 0.;
}
