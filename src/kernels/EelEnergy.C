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

#include "EelEnergy.h"
/**
This function computes the convective part of the total energy equation.
 */
template<>
InputParameters validParams<EelEnergy>()
{
  InputParameters params = validParams<Kernel>();
    params.addRequiredCoupledVar("rhoA", "density");
    params.addRequiredCoupledVar("rhouA_x", "x component of rhouA");
    params.addCoupledVar("rhouA_y", "y component of rhouA");
    params.addCoupledVar("rhouA_z", "z component of rhouA");
    params.addRequiredCoupledVar("pressure", "pressure");
    params.addRequiredCoupledVar("area", "area");
    params.addRequiredParam<UserObjectName>("eos", "Equation of state");
    params.addParam<std::string>("Hw_fn_name", "Function name for the wall heat transfer.");
    params.addParam<std::string>("Tw_fn_name", "Function name for the wall temperature.");
    params.addParam<Real>("Hw", 0., "Wall heat transfer.");
    params.addParam<Real>("Tw", 0., "Wall temperature.");
    params.addParam<Real>("aw", 0., "Wall heat surface.");
    params.addParam<RealVectorValue>("gravity", (0., 0., 0.), "Gravity vector.");
  return params;
}

EelEnergy::EelEnergy(const std::string & name,
                       InputParameters parameters) :
  Kernel(name, parameters),
    // Coupled variables
    _rhoA(coupledValue("rhoA")),
    _rhouA_x(coupledValue("rhouA_x")),
    _rhouA_y(_mesh.dimension()>=2 ? coupledValue("rhouA_y") : _zero),
    _rhouA_z(_mesh.dimension()==3 ? coupledValue("rhouA_z") : _zero),
    _pressure(coupledValue("pressure")),
    _area(coupledValue("area")),
    // Parameters:
    _Hw_fn_name(isParamValid("Hw_fn_name") ? getParam<std::string>("Hw_fn_name") : std::string(" ")),
    _Tw_fn_name(isParamValid("Tw_fn_name") ? getParam<std::string>("Tw_fn_name") : std::string(" ")),
    _Hw(getParam<Real>("Hw")),
    _Tw(getParam<Real>("Tw")),
    _aw(getParam<Real>("aw")),
    _gravity(getParam<RealVectorValue>("gravity")),
    // Equation of state:
    _eos(getUserObject<EquationOfState>("eos")),
    // Parameters for jacobian:
    _rhoA_nb(coupled("rhoA")),
    _rhouA_x_nb(coupled("rhouA_x")),
    _rhouA_y_nb(isCoupled("rhouA_y") ? coupled("rhouA_y") : -1),
    _rhouA_z_nb(isCoupled("rhouA_z") ? coupled("rhouA_z") : -1)
{
}

Real EelEnergy::computeQpResidual()
{
    // Compute convective part of the energy equation:
    RealVectorValue _conv;
    _conv(0) = _rhouA_x[_qp] * ( _u[_qp] + _pressure[_qp]*_area[_qp] ) / _rhoA[_qp];
    _conv(1) = _rhouA_y[_qp] * ( _u[_qp] + _pressure[_qp]*_area[_qp] ) / _rhoA[_qp];
    _conv(2) = _rhouA_z[_qp] * ( _u[_qp] + _pressure[_qp]*_area[_qp] ) / _rhoA[_qp];
    
    // Gravity work:
    RealVectorValue _vector_vel(_rhouA_x[_qp]/_rhoA[_qp], _rhouA_y[_qp]/_rhoA[_qp], _rhouA_z[_qp]/_rhoA[_qp]);
    Real _gravity_work = _rhoA[_qp]*_gravity*_vector_vel;
    
    // Wall heat tranfer (WHT):
    Real rho = _rhoA[_qp] / _area[_qp];
    
    Real Hw_val = isParamValid("Hw_fn_name") ? getFunctionByName(_Hw_fn_name).value(_t, _q_point[_qp]) : _Hw;
    Real Tw_val = isParamValid("Tw_fn_name") ? getFunctionByName(_Tw_fn_name).value(_t, _q_point[_qp]) : _Tw;
    Real WHT = Hw_val * _aw * ( _eos.temperature_from_p_rho(_pressure[_qp], rho) - Tw_val );

    // Returns the residual
    return -_conv * _grad_test[_i][_qp] + (WHT+_gravity_work)*_test[_i][_qp];
}

Real EelEnergy::computeQpJacobian()
{
    // Compute the velocity and momentum vectors:
    RealVectorValue _vector_vel(_rhouA_x[_qp]/_rhoA[_qp], _rhouA_y[_qp]/_rhoA[_qp], _rhouA_z[_qp]/_rhoA[_qp]);
    RealVectorValue _rhouA_vec(_rhouA_x[_qp], _rhouA_y[_qp], _rhouA_z[_qp]);
    
    // Return the value of the jacobian:
    return -_phi[_j][_qp]*_vector_vel*(1+_eos.dAp_drhoEA(_rhoA[_qp], _rhouA_vec.size(), _u[_qp]))*_grad_test[_i][_qp];
}

Real EelEnergy::computeQpOffDiagJacobian( unsigned int _jvar)
{
    // Compute the velocity and momentum vectors:
    RealVectorValue _vector_vel(_rhouA_x[_qp]/_rhoA[_qp], _rhouA_y[_qp]/_rhoA[_qp], _rhouA_z[_qp]/_rhoA[_qp]);
    RealVectorValue _rhouA_vec(_rhouA_x[_qp], _rhouA_y[_qp], _rhouA_z[_qp]);
    
    // jacobian term from the density (rho*A):
    if (_jvar == _rhoA_nb) {
        return _phi[_j][_qp]*( _vector_vel*(_u[_qp]+_area[_qp]*_pressure[_qp])/_rhoA[_qp] - _vector_vel*_eos.dAp_drhoA(_rhoA[_qp], _rhouA_vec.size(), _u[_qp]) )*_grad_test[_i][_qp];
    }
    // x-momentum components:
    else if (_jvar == _rhouA_x_nb) {
        Real _press_term = _vector_vel*_grad_test[_i][_qp]*_eos.dAp_drhouA(_rhoA[_qp], _rhouA_x[_qp], _u[_qp]);
        return -_phi[_j][_qp]*( (_u[_qp]+_area[_qp]*_pressure[_qp])/_rhoA[_qp]*_grad_test[_i][_qp](0) + _press_term );
    }
    // y-momentum components:
    else if (_jvar == _rhouA_y_nb && _mesh.dimension() >= 2) {
        Real _press_term = _vector_vel*_grad_test[_i][_qp]*_eos.dAp_drhouA(_rhoA[_qp], _rhouA_y[_qp], _u[_qp]);
        return -_phi[_j][_qp]*( (_u[_qp]+_area[_qp]*_pressure[_qp])/_rhoA[_qp]*_grad_test[_i][_qp](1) + _press_term );
    }
    // z-momentum components:
    else if (_jvar == _rhouA_z_nb && _mesh.dimension() == 3) {
        Real _press_term = _vector_vel*_grad_test[_i][_qp]*_eos.dAp_drhouA(_rhoA[_qp], _rhouA_z[_qp], _u[_qp]);
        return -_phi[_j][_qp]*( (_u[_qp]+_area[_qp]*_pressure[_qp])/_rhoA[_qp]*_grad_test[_i][_qp](2) + _press_term );
    }
    else
        return 0.;
}
