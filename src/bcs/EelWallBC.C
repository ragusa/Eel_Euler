#include "EelWallBC.h"

template<>
InputParameters validParams<EelWallBC>()
{
  InputParameters params = validParams<IntegratedBC>();

  params.addRequiredParam<std::string>("equation_name", "The name of the equation this BC is acting on");
    // Coupled variables:
    params.addRequiredCoupledVar("area", "area");
    // Equation of state
    params.addRequiredParam<UserObjectName>("eos", "Equation of state");
    // Variables used in the jacobian matrix:
    params.addRequiredCoupledVar("rhoA", "density: rho*A");
    params.addRequiredCoupledVar("rhouA_x", "x-momentum: rho*u*A");
    params.addCoupledVar("rhouA_y", "y-momentum: rho*v*A");
    params.addRequiredCoupledVar("rhoEA", "energy: rho*E*A");
  return params;
}

EelWallBC::EelWallBC(const std::string & name, InputParameters parameters) :
    IntegratedBC(name, parameters),
    // Name of the equation:
    _eqn_name(getParam<std::string>("equation_name")),
    _eqn_type("CONTINUITY, XMOMENTUM, YMOMENTUM, ENERGY, INVALID", "INVALID"),
    // Coupled variables:
    _rhoA(coupledValue("rhoA")),
    _rhouA_x(coupledValue("rhouA_x")),
    _rhouA_y(isCoupled("rhouA_y") ? coupledValue("rhouA_y") : _zero),
    _rhoEA(coupledValue("rhoEA")),
    _area(coupledValue("area")),
    // Equation of state:
    _eos(getUserObject<EquationOfState>("eos")),
    // Parameters for jacobian matrix:
    _rhoA_nb(isCoupled("rhoA") ? coupled("rhoA") : -1),
    _rhouA_x_nb(isCoupled("rhouA_x") ? coupled("rhouA_x") : -1),
    _rhouA_y_nb(isCoupled("rhouA_y") ? coupled("rhouA_y") : -1),
    _rhoEA_nb(isCoupled("rhoEA") ? coupled("rhoEA") : -1)
{
  _eqn_type = _eqn_name;
}

Real
EelWallBC::computeQpResidual()
{
    // Compute the pressure:
    RealVectorValue vel_vec(_rhouA_x[_qp]/_rhoA[_qp], _rhouA_y[_qp]/_rhoA[_qp], 0.);
    Real pressure = _eos.pressure(_rhoA[_qp]/_area[_qp], vel_vec.size(), _rhoEA[_qp]/_area[_qp]);
    
    // Switch statement on the
    switch (_eqn_type) {
        case CONTINUITY:
            return 0.;
//            break;
        case XMOMENTUM:
            return _area[_qp]*pressure*_normals[_qp](0)*_test[_i][_qp];
//            break;
        case YMOMENTUM:
            return _area[_qp]*pressure*_normals[_qp](1)*_test[_i][_qp];
//            break;
        case ENERGY:
            return 0.;
//            break;
        default:
            mooseError("The equation with name: \"" << _eqn_name << "\" is not supported in the \"EelWallBC\" type of boundary condition.");
//            return 0.;
//            break;
    }
}

Real
EelWallBC::computeQpJacobian()
{
    switch (_eqn_type) {
        case CONTINUITY:
            return 0.;
//            break;
        case XMOMENTUM:
            return _phi[_j][_qp]*_eos.dAp_drhouA(_rhoA[_qp],_u[_qp],_rhoEA[_qp])*_normals[_qp](0)*_test[_i][_qp];
//            break;
        case YMOMENTUM:
            return _phi[_j][_qp]*_eos.dAp_drhouA(_rhoA[_qp],_u[_qp],_rhoEA[_qp])*_normals[_qp](1)*_test[_i][_qp];
//            break;
        case ENERGY:
            return 0.;
//            break;
        default:
            mooseError("The equation with name: \"" << _eqn_name << "\" is not supported in the \"EelWallBC\" type of boundary condition.");
//            return 0.;
//            break;
    }
}

Real
EelWallBC::computeQpOffDiagJacobian(unsigned _jvar)
{
    RealVectorValue _rhouA_vec(_rhouA_x[_qp], _rhouA_y[_qp], 0.);
    switch (_eqn_type) {
        case CONTINUITY:
            return 0.;
//            break;
        case XMOMENTUM:
            if (_jvar == _rhoA_nb)
                return _phi[_j][_qp]*_eos.dAp_drhoA(_rhoA[_qp],_rhouA_vec.size(),_rhoEA[_qp])*_normals[_qp](0)*_test[_i][_qp];
            else if (_jvar == _rhouA_y_nb)
                return _phi[_j][_qp]*_eos.dAp_drhouA(_rhoA[_qp],_rhouA_y[_qp],_rhoEA[_qp])*_normals[_qp](0)*_test[_i][_qp];
            else if (_jvar == _rhoEA_nb)
                return _phi[_j][_qp]*_eos.dAp_drhoEA(_rhoA[_qp],_rhouA_vec.size(),_rhoEA[_qp])*_normals[_qp](0)*_test[_i][_qp];
            else
                return 0.;
//            break;
        case YMOMENTUM:
            if (_jvar == _rhoA_nb)
                return _phi[_j][_qp]*_eos.dAp_drhoA(_rhoA[_qp],_rhouA_vec.size(),_rhoEA[_qp])*_normals[_qp](0)*_test[_i][_qp];
            else if (_jvar == _rhouA_x_nb)
                return _phi[_j][_qp]*_eos.dAp_drhouA(_rhoA[_qp],_rhouA_x[_qp],_rhoEA[_qp])*_normals[_qp](0)*_test[_i][_qp];
            else if (_jvar == _rhoEA_nb)
                return _phi[_j][_qp]*_eos.dAp_drhoEA(_rhoA[_qp],_rhouA_vec.size(),_rhoEA[_qp])*_normals[_qp](0)*_test[_i][_qp];
            else
                return 0.;
//            break;
        case ENERGY:
            return 0.;
//            break;
        default:
            mooseError("The equation with name: \"" << _eqn_name << "\" is not supported in the \"EelWallBC\" type of boundary condition.");
//            return 0.;
//            break;
    }
}
