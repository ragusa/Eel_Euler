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

#include "EelArtificialVisc.h"
/**
This function computes the dissipative terms for all of the equations. It is dimension agnostic.
 */
template<>
InputParameters validParams<EelArtificialVisc>()
{
  InputParameters params = validParams<Kernel>();
    // Equation and diffusion names:
    params.addParam<std::string>("equation_name", "INVALID", "Name of the equation.");
    params.addParam<std::string>("diffusion_name", "PARABOLIC", "Name of the diffusion.");
    // Coupled aux variables
    params.addRequiredCoupledVar("density", "density of the fluid");
    params.addRequiredCoupledVar("velocity_x", "x component of the velocity");
    params.addCoupledVar("velocity_y", "y component of the velocity");
    params.addCoupledVar("velocity_z", "z component of the velocity");
    params.addRequiredCoupledVar("internal_energy", "internal energy of the fluid");
    params.addRequiredCoupledVar("area", "area of the geometry");
    params.addRequiredCoupledVar("norm_velocity", "norm of the velocity vector");
  return params;
}

EelArtificialVisc::EelArtificialVisc(const std::string & name,
                       InputParameters parameters) :
  Kernel(name, parameters),
    // Declare equation types
    _equ_name(getParam<std::string>("equation_name")),
    _diff_name(getParam<std::string>("diffusion_name")),
    _equ_type("CONTINUITY, XMOMENTUM, YMOMENTUM, ZMOMENTUM, ENERGY, INVALID", _equ_name),
    _diff_type("ENTROPY, PARABOLIC, NONE, INVALID",_diff_name),
    // Coupled auxilary variables
    _rho(coupledValue("density")),
    _grad_rho(coupledGradient("density")),
    _vel_x(coupledValue("velocity_x")),
    _vel_y(_mesh.dimension()>=2 ? coupledValue("velocity_y") : _zero),
    _vel_z(_mesh.dimension()==3 ? coupledValue("velocity_z") : _zero),
    _grad_vel_x(coupledGradient("velocity_x")),
    _grad_vel_y(_mesh.dimension()>=2 ? coupledGradient("velocity_y") : _grad_zero),
    _grad_vel_z(_mesh.dimension()==3 ? coupledGradient("velocity_z") : _grad_zero),
    _grad_rhoe(coupledGradient("internal_energy")),
    _area(coupledValue("area")),
    _norm_vel(coupledValue("norm_velocity")),
    _grad_norm_vel(coupledGradient("norm_velocity")),
    // Material property: viscosity coefficient.
    _mu(getMaterialProperty<Real>("mu")),
    _kappa(getMaterialProperty<Real>("kappa"))
{
//    _equ_type = _equ_name;
//    _diff_type = _diff_name;
}

Real EelArtificialVisc::computeQpResidual()
{
    // Determine if cell is on boundary or not and then compute a unit vector 'l=grad(norm(vel))/norm(grad(norm(vel)))':
    Real isonbnd = 1.;
    if (_mesh.isBoundaryNode(_current_elem->node(_i))==true) {
        isonbnd = 0.;
    }

    // If statement on diffusion type:
    if (_diff_type == 1) {
        switch (_equ_type) {
            case CONTINUITY:
                return _area[_qp]*_kappa[_qp]*_grad_rho[_qp]*_grad_test[_i][_qp];
                break;
            case XMOMENTUM:
                return _area[_qp]*_kappa[_qp]*(_rho[_qp]*_grad_vel_x[_qp]+_vel_x[_qp]*_grad_rho[_qp])*_grad_test[_i][_qp];
                break;
            case YMOMENTUM:
                return _area[_qp]*_kappa[_qp]*(_rho[_qp]*_grad_vel_y[_qp]+_vel_y[_qp]*_grad_rho[_qp])*_grad_test[_i][_qp];
                break;
            case ZMOMENTUM:
                return _area[_qp]*_kappa[_qp]*(_rho[_qp]*_grad_vel_z[_qp]+_vel_z[_qp]*_grad_rho[_qp])*_grad_test[_i][_qp];
                break;
            case ENERGY:
                return _area[_qp]*_kappa[_qp]*(_grad_rhoe[_qp]+_rho[_qp]*_norm_vel[_qp]*_grad_norm_vel[_qp]+0.5*_norm_vel[_qp]*_norm_vel[_qp]*_grad_rho[_qp])*_grad_test[_i][_qp];
                break;
            default:
                mooseError("INVALID equation name.");
        }
    }
    else if (_diff_type == 0) {
        // Compute 0.5*rho*mu*grad(vel)_symmetric: (get a symmetric tensor)
        TensorValue<Real> grad_vel_tensor(_grad_vel_x[_qp], _grad_vel_y[_qp], _grad_vel_z[_qp]);
        grad_vel_tensor = ( grad_vel_tensor + grad_vel_tensor.transpose() );
        grad_vel_tensor *= 0.5 * _rho[_qp] * _mu[_qp];
        
        // Compute f = kappa * grad(rho):
        RealVectorValue f(_kappa[_qp]*_grad_rho[_qp](0), _kappa[_qp]*_grad_rho[_qp](1), _kappa[_qp]*_grad_rho[_qp](2));
        
        // Compute velocity vector and its norm:
        RealVectorValue vel_vector(_vel_x[_qp], _vel_y[_qp], _vel_z[_qp]);
        Real norm_vel2 = vel_vector.size_sq();
        
        // Compute h = kappa * grad(rho*e):
        RealVectorValue h(_kappa[_qp]*_grad_rhoe[_qp](0), _kappa[_qp]*_grad_rhoe[_qp](1), _kappa[_qp]*_grad_rhoe[_qp](2));
        
        // return the dissipative terms:
            switch (_equ_type) {
                case CONTINUITY: // div(kappa grad(rho))
                    return isonbnd*_area[_qp] * f * _grad_test[_i][_qp];
                    break;
                case XMOMENTUM:
                    return isonbnd*_area[_qp]*( _vel_x[_qp]*f + grad_vel_tensor.row(0) ) * _grad_test[_i][_qp];
                    break;
                case YMOMENTUM:
                    return isonbnd*_area[_qp]*( _vel_y[_qp]*f + grad_vel_tensor.row(1) ) * _grad_test[_i][_qp];
                    break;
                case ZMOMENTUM:
                    return isonbnd*_area[_qp]*( _vel_z[_qp]*f + grad_vel_tensor.row(2) ) * _grad_test[_i][_qp];
                    break;
                case ENERGY:
                    return isonbnd*_area[_qp]*( h + 0.5*f*norm_vel2 + grad_vel_tensor*vel_vector )*_grad_test[_i][_qp];
                    break;
                default:
                    mooseError("INVALID equation name.");
            }
    }
    else {
        mooseError("INVALID dissipation terms.");
    }
}

Real EelArtificialVisc::computeQpJacobian()
{
    // We assumed that the all of the above regularization can be approximated by the parabolic regularization:
    Real _mu_jac = std::max(_mu[_qp], _kappa[_qp]);
    return _mu_jac*_grad_phi[_j][_qp]*_grad_test[_i][_qp];
}

Real EelArtificialVisc::computeQpOffDiagJacobian( unsigned int _jvar)
{
    // With above assumption, there is no contribution from the dissipative terms to the off diagonal terms of the jacobian matrix.
    return 0.*_jvar;
}
