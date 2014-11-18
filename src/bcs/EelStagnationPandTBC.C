#include "EelStagnationPandTBC.h"

template<>
InputParameters validParams<EelStagnationPandTBC>()
{
  InputParameters params = validParams<IntegratedBC>();

  params.addRequiredParam<std::string>("equation_name", "The name of the equation this BC is acting on");

    // Coupled conservative variables:
    params.addCoupledVar("rhoA", "density: rhoA");
    params.addCoupledVar("rhouA_x", "x component of the momentum: rhouA_x");
    params.addCoupledVar("rhouA_y", "y component of the momentum: rhouA_y");
    params.addCoupledVar("rhoEA", "total energy: rho*E*A");
    // Coupled aux variables:
    params.addCoupledVar("area", "Coupled area variable");
    // Input parameters
    params.addRequiredParam<Real>("p0_bc", "Stagnation pressure at the boundary");
    params.addRequiredParam<Real>("T0_bc", "Stagnation temperature at the boundary");
    params.addParam<Real>("gamma0_bc", 0., "Stagnation angle");
    // Make the name of the EOS function a required parameter.
    params.addRequiredParam<UserObjectName>("eos", "The name of equation of state object to use.");

  return params;
}

EelStagnationPandTBC::EelStagnationPandTBC(const std::string & name, InputParameters parameters) :
    IntegratedBC(name, parameters),
    // Type of equation:
    _eqn_name(getParam<std::string>("equation_name")),
    _eqn_type("VOIDFRACTION, CONTINUITY, XMOMENTUM, YMOMENTUM, ZMOMENTUM, ENERGY, INVALID", _eqn_name),
    // Coupled aux variables:
    _rhoA(coupledValue("rhoA")),
    _rhouA_x(coupledValue("rhouA_x")),
    _rhouA_y(isCoupled("rhouA_y") ? coupledValue("rhouA_y") : _zero),
    _rhoEA(isCoupled("rhoEA") ? coupledValue("rhoEA") : _zero),
    // Coupled aux variables:
    _area(coupledValue("area")),
    // Stagnation variables:
    _p0_bc(getParam<Real>("p0_bc")),
    _T0_bc(getParam<Real>("T0_bc")),
    _gamma0_bc(getParam<Real>("gamma0_bc")),
    // Equation of state:
    _eos(getUserObject<EquationOfState>("eos")),
    // Parameters for jacobian matrix:
    _rhoA_nb(coupled("rhoA")),
    _rhouA_x_nb(coupled("rhouA_x")),
    _rhouA_y_nb(isCoupled("rhouA_y") ? coupled("rhouA_x") : -1),
    _rhoEA_nb(coupled("rhoEA"))
{
    //std::cout<<"bug1"<<std::endl;
    _rho0_bc = _eos.rho_from_p_T(_p0_bc, _T0_bc);
    _H0_bc = _eos.e_from_p_rho(_p0_bc, _rho0_bc) + _p0_bc / _rho0_bc;
    //std::cout<<"rho0="<<_rho0_bc<<std::endl;
    _K = (_p0_bc + _eos.Pinf()) / std::pow(_rho0_bc, _eos.gamma());
    _H_bar = _eos.gamma() * (_p0_bc + _eos.Pinf()) / _rho0_bc / (_eos.gamma() - 1);
    
//    std::cout<<"K="<<_K<<std::endl;
//    std::cout<<"H="<<_H_bar<<std::endl;
//    std::cout<<"H0="<<_H0_bc<<std::endl;
//    std::cout<<"rho0="<<_rho0_bc<<std::endl;
}

Real
EelStagnationPandTBC::computeQpResidual()
{
    // Compute u_star and v_star:
    Real u_star = _rhouA_x[_qp]/_rhoA[_qp];
    Real v_star = u_star * std::tan(_gamma0_bc);
    RealVectorValue _vel_star(u_star, v_star);
    Real _norm_vel_star2 = _vel_star.size_sq();
    // Compute rho_star and static pressure:
    Real rho_star = std::pow((_H_bar - 0.5*_norm_vel_star2)*(_eos.gamma()-1)/(_eos.gamma())/_K, 1./(_eos.gamma()-1));
    Real p_bc = _K * std::pow(rho_star, _eos.gamma()) - _eos.Pinf();
    
  switch (_eqn_type)
  {
  case CONTINUITY:
        return _area[_qp]*rho_star*( u_star*_normals[_qp](0) + v_star*_normals[_qp](1) ) * _test[_i][_qp];
//        break;
  case XMOMENTUM:
        return _area[_qp]*(u_star*rho_star*(_vel_star*_normals[_qp])+p_bc*_normals[_qp](0)) * _test[_i][_qp];
//        break;
  case YMOMENTUM:
        return _area[_qp]*(v_star*rho_star*(_vel_star*_normals[_qp])+p_bc*_normals[_qp](1)) * _test[_i][_qp];
//        break;
  case ENERGY:
        return _area[_qp]*rho_star*_H0_bc*(_vel_star*_normals[_qp]) * _test[_i][_qp];
//        break;
  default:
    mooseError("The equation with name: \"" << _eqn_name << "\" is not supported in the \"OneD7EqnStagnationPandTBC\" type of boundary condition.");
  }
}

Real
EelStagnationPandTBC::computeQpJacobian()
{
    // Compute u_star and v_star:
    Real u_star = _rhouA_x[_qp]/_rhoA[_qp];
    Real v_star = u_star * std::tan(_gamma0_bc);
    // Declare parameter:
    Real _press_term = 0.;
    
    switch (_eqn_type)
    {
        case CONTINUITY:
            return 0.;
//            break;
        case XMOMENTUM:
            _press_term = _eos.dAp_drhouA(_rhoA[_qp], _rhouA_x[_qp], _rhoEA[_qp]);
            return _phi[_j][_qp]*((2*u_star+_press_term)*_normals[_qp](0)+v_star*_normals[_qp](1))*_test[_i][_qp];
//            break;
        case YMOMENTUM:
            _press_term = _eos.dAp_drhouA(_rhoA[_qp], _rhouA_y[_qp], _rhoEA[_qp]);
            return _phi[_j][_qp]*(u_star*_normals[_qp](0)+(2*v_star+_press_term)*_normals[_qp](1))*_test[_i][_qp];
//            break;
        case ENERGY:
            return 0.;
//            break;
        default:
            mooseError("The equation with name: \"" << _eqn_name << "\" is not supported in the \"OneD7EqnStagnationPandTBC\" type of boundary condition.");
    }
}

Real
EelStagnationPandTBC::computeQpOffDiagJacobian(unsigned _jvar)
{
    // Compute u_star and v_star:
    Real u_star = _rhouA_x[_qp]/_rhoA[_qp];
    Real v_star = u_star * std::tan(_gamma0_bc);
    RealVectorValue _vel_star(u_star, v_star);
    // Declare parameter:
    Real _press_term = 0.;
    RealVectorValue _rhouA_vec(_rhouA_x[_qp], _rhouA_y[_qp], 0.);
    
    switch (_eqn_type)
    {
        case CONTINUITY:
            if (_jvar == _rhouA_x_nb)
                return _phi[_j][_qp]*_normals[_qp](0)*_test[_i][_qp];
            else if (_jvar == _rhouA_y_nb)
                return _phi[_j][_qp]*_normals[_qp](0)*_test[_i][_qp];
            else
                return 0.;
        case XMOMENTUM:
            if(_jvar == _rhoA_nb) {
                _press_term = _eos.dAp_drhoA(_rhoA[_qp], _rhouA_vec.size(), _rhoEA[_qp]);
                return -_phi[_j][_qp]*(_rhouA_x[_qp]*_vel_star*_normals[_qp]+_press_term*_normals[_qp](0))*_test[_i][_qp];
            }
            else if (_jvar == _rhouA_y_nb) {
                _press_term = _eos.dAp_drhouA(_rhoA[_qp], _rhouA_y[_qp], _rhoEA[_qp]);
                return _phi[_j][_qp]*(_rhouA_x[_qp]/_rhoA[_qp]*_normals[_qp](1)+_press_term*_normals[_qp](0))*_test[_i][_qp];
            }
            else if (_jvar == _rhoEA_nb) {
                _press_term = _eos.dAp_drhoEA(_rhoA[_qp], _rhouA_vec.size(), _rhoEA[_qp]);
                return _phi[_j][_qp]*_press_term*_test[_i][_qp]*_normals[_qp](0);
            }
            else
                return 0.;
//            break;
        case YMOMENTUM:
            if(_jvar == _rhoA_nb) {
                _press_term = _eos.dAp_drhoA(_rhoA[_qp], _rhouA_vec.size(), _rhoEA[_qp]);
                return -_phi[_j][_qp]*(_rhouA_y[_qp]*_vel_star*_normals[_qp]+_press_term*_normals[_qp](1))*_test[_i][_qp];
            }
            else if (_jvar == _rhouA_x_nb) {
                _press_term = _eos.dAp_drhouA(_rhoA[_qp], _rhouA_x[_qp], _rhoEA[_qp]);
                return _phi[_j][_qp]*(_rhouA_y[_qp]/_rhoA[_qp]*_normals[_qp](0)+_press_term*_normals[_qp](1))*_test[_i][_qp];
            }
            else if (_jvar == _rhoEA_nb) {
                _press_term = _eos.dAp_drhoEA(_rhoA[_qp], _rhouA_vec.size(), _rhoEA[_qp]);
                return _phi[_j][_qp]*_press_term*_test[_i][_qp]*_normals[_qp](1);
            }
            else
                return 0.;
//            break;
        case ENERGY:
            return 0.;
//            break;
        default:
            mooseError("The equation with name: \"" << _eqn_name << "\" is not supported in the \"OneD7EqnStagnationPandTBC\" type of boundary condition.");
    }
}
