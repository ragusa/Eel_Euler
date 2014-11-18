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

#include "LowMachPreconditioner.h"
/**
This function computes the x, y and z momentum equationS. It is dimension agnostic. 
 */
template<>
InputParameters validParams<LowMachPreconditioner>()
{
  InputParameters params = validParams<Kernel>();
    params.addRequiredCoupledVar("rhoA", "rho*A");
    params.addRequiredCoupledVar("rhouA_x", "x component of momentum");
    params.addCoupledVar("rhouA_y", "y component of momentum");
    params.addCoupledVar("rhouA_z", "z component of momentum");
    params.addRequiredCoupledVar("pressure", "pressure");
    params.addRequiredCoupledVar("area", "area");
    params.addRequiredParam<UserObjectName>("eos", "Equation of state");
    params.addParam<Real>("Mach_nb_ref", 1., "Reference Mach number.");
  return params;
}

LowMachPreconditioner::LowMachPreconditioner(const std::string & name,
                       InputParameters parameters) :
  Kernel(name, parameters),
    // Coupled auxilary variables
    _rhoA(coupledValue("rhoA")),
    _rhouA_x(coupledValue("rhouA_x")),
    _rhouA_y(_mesh.dimension()>=2 ? coupledValue("rhouA_y") : _zero),
    _rhouA_z(_mesh.dimension()==3 ? coupledValue("rhouA_z") : _zero),
    _pressure_old(coupledValueOld("pressure")),
    _pressure_older(coupledValueOlder("pressure")),
    _area(coupledValue("area")),
    // Equation of state:
    _eos(getUserObject<EquationOfState>("eos")),
    // Parameters:
    _Mach_ref(getParam<Real>("Mach_nb_ref")),
    // Parameters for jacobian:
    _rhoA_nb(coupled("rhoA")),
    _rhouA_x_nb(coupled("rhouA_x")),
    _rhouA_y_nb(isCoupled("rhouA_y") ? coupled("rhouA_y") : -1),
    _rhouA_z_nb(isCoupled("rhouA_z") ? coupled("rhouA_z") : -1)


{}

Real LowMachPreconditioner::computeQpResidual()
{
// Compute the weight for BDF2:
    Real weight0 = (2.*_dt+_dt_old)/(_dt*(_dt+_dt_old));
    Real weight1 = -(_dt+_dt_old)/(_dt*_dt_old);
    Real weight2 = _dt/(_dt_old*(_dt+_dt_old));
    
  // Compute density, velocity and total energy:
    Real rho = _rhoA[_qp] / _area[_qp];
    RealVectorValue vector_vel( _rhouA_x[_qp]/_rhoA[_qp], _rhouA_y[_qp]/_rhoA[_qp], _rhouA_z[_qp]/_rhoA[_qp] );
    Real rhoE = _u[_qp] / _area[_qp];
    
    // Compute the pressure:
    Real press = _eos.pressure(rho, vector_vel.size(), rhoE);
    
    // Reference Mach number term:
    Real Mach_term = (1. - _Mach_ref*_Mach_ref) / (_Mach_ref*_Mach_ref);
//    std::cout<<Mach_term<<std::endl;
    
  // Return the kernel value:
    if (_t_step <= 1) {
        return Mach_term * (press - _pressure_old[_qp]) / _dt * _test[_i][_qp];
    }
    else
        return Mach_term * (weight0*press + weight1*_pressure_old[_qp] + weight2*_pressure_older[_qp]) * _test[_i][_qp];
}

Real LowMachPreconditioner::computeQpJacobian()
{
    // Compute weight:
    Real weight0 = (2.*_dt+_dt_old)/(_dt*(_dt+_dt_old));
    
    // Reference Mach number term:
    Real Mach_term = (1. - _Mach_ref*_Mach_ref) / (_Mach_ref*_Mach_ref);
    
    // Compute the momentum vector rhouA:
    RealVectorValue rhouA_vec(_rhouA_x[_qp], _rhouA_y[_qp], _rhouA_z[_qp]);
    
    // Compute the derivative of pressure with respect to rhoEA
    Real press_term = _eos.dAp_drhoEA(_rhoA[_qp], rhouA_vec.size(), _u[_qp]);
    
    // Return the value of the jacobian:
    if (_t_step <= 1) {
        return _phi[_j][_qp] * Mach_term * press_term / _dt * _test[_i][_qp];
    }
    else
        return _phi[_j][_qp] * Mach_term * press_term * weight0 * _test[_i][_qp];
}

Real LowMachPreconditioner::computeQpOffDiagJacobian( unsigned int _jvar)
{
    // Compute weight:
    Real weight0 = (2.*_dt+_dt_old)/(_dt*(_dt+_dt_old));
    
    // Reference Mach number term:
    Real Mach_term = (1. - _Mach_ref*_Mach_ref) / (_Mach_ref*_Mach_ref);
    
    // Compute the momentum vector rhouA:
    RealVectorValue rhouA_vec(_rhouA_x[_qp], _rhouA_y[_qp], _rhouA_z[_qp]);
    
    Real press_term = 0.;
    
    // density (rho*A):
    if (_jvar == _rhoA_nb) {
        press_term = _eos.dAp_drhoA(_rhoA[_qp], rhouA_vec.size(), _u[_qp]);
    }
    // x-momentum component:
    else if (_jvar == _rhouA_x_nb ) {
        press_term = _eos.dAp_drhouA(_rhoA[_qp], _rhouA_x[_qp], _u[_qp]);
    }
    // y-momentum component:
    else if (_jvar == _rhouA_y_nb && _mesh.dimension()>=2 ) {
        press_term = _eos.dAp_drhouA(_rhoA[_qp], _rhouA_y[_qp], _u[_qp]);
    }
    // z-momentum component:
    else if (_jvar == _rhouA_z_nb && _mesh.dimension()==3 ) {
        press_term = _eos.dAp_drhouA(_rhoA[_qp], _rhouA_z[_qp], _u[_qp]);
    }
    else
        press_term = 0.;
    
    // Return the value of the jacobian:
    if (_t_step <= 1) {
        return _phi[_j][_qp] * Mach_term * press_term / _dt * _test[_i][_qp];
    }
    else
        return _phi[_j][_qp] * Mach_term * press_term * weight0 * _test[_i][_qp];
}
