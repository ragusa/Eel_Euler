#include "EelStaticPandTBC.h"

template<>
InputParameters validParams<EelStaticPandTBC>()
{
  InputParameters params = validParams<IntegratedBC>();

    params.addRequiredParam<std::string>("equation_name", "The name of the equation this BC is acting on");
    // Coupled variables:
    params.addRequiredCoupledVar("rhoA", "density: rho*A");
    params.addRequiredCoupledVar("rhouA_x", "x component of the momentum");
    params.addCoupledVar("rhouA_y", "y component of the momentum");
    params.addRequiredCoupledVar("rhoEA", "total energy: rho*E*A");
    params.addRequiredCoupledVar("area", "area aux variable");
    // Input parameters:
    params.addRequiredParam<Real>("p_bc", "Static pressure at the boundary");
    params.addParam<Real>("T_bc", 0.0, "Static temperature at the boundary");
    params.addParam<Real>("gamma_bc", 0.0, "inflow angle for inlet BC, [-], ignored for outlet condition");
    // Equation of state:
    params.addRequiredParam<UserObjectName>("eos", "The name of equation of state object to use.");
  return params;
}

EelStaticPandTBC::EelStaticPandTBC(const std::string & name, InputParameters parameters) :
    IntegratedBC(name, parameters),
    // Name of the equation:
    _eqn_name(getParam<std::string>("equation_name")),
    _eqn_type("CONTINUITY, XMOMENTUM, YMOMENTUM, ENERGY, INVALID", "INVALID"),
    // Coupled variables:
    _rhoA(coupledValue("rhoA")),
    _rhouA_x(coupledValue("rhouA_x")),
    _rhouA_y(_mesh.dimension()>=2 ? coupledValue("rhouA_y") : _zero),
    _rhoEA(coupledValue("rhoEA")),
    _rhoA_old(coupledValueOld("rhoA")),
    _rhouA_x_old(coupledValueOld("rhouA_x")),
    _rhouA_y_old(_mesh.dimension()>=2 ? coupledValueOld("rhouA_y") : _zero),
    _rhoEA_old(coupledValueOld("rhoEA")),
    _area(coupledValue("area")),
    // Boundary condition parameters:
    _p_bc(getParam<Real>("p_bc")),
    _T_bc(getParam<Real>("T_bc")),
    _gamma_bc(getParam<Real>("gamma_bc")),
    // Equation of state:
    _eos(getUserObject<EquationOfState>("eos")),
    // Parameter for jacobian matrix:
    _rhoA_nb(coupled("rhoA")),
    _rhouA_x_nb(coupled("rhouA_x")),
    _rhouA_y_nb(isCoupled("rhouA_y") ? coupled("rhouA_y") : -1),
    _rhoEA_nb(coupled("rhoEA"))
{
  _eqn_type = _eqn_name;
}

Real
EelStaticPandTBC::computeQpResidual()
{
    // Compute v dot n:
    RealVectorValue _vel_vec(_rhouA_x[_qp]/_rhoA[_qp], _rhouA_y[_qp]/_rhoA[_qp], 0.);
    RealVectorValue _vel_vec_old(_rhouA_x_old[_qp]/_rhoA_old[_qp], _rhouA_y_old[_qp]/_rhoA_old[_qp], 0.);
    Real _v_dot_n = _vel_vec * _normals[_qp];
    if ( _v_dot_n <0 ) // Inlet
    {
        Real _rho_bc = _eos.rho_from_p_T(_p_bc, _T_bc);
        Real _vel_x = _rhouA_x[_qp]/_rhoA[_qp];
        Real _e_bc = 0; Real _rhoE_bc = 0; Real _vel_y_bc = 0;
        Real _norm_vel2 = 0;
        if (_gamma_bc != 0)
            _vel_y_bc = _rhouA_x[_qp]/_rhoA[_qp]*std::tan(_gamma_bc);
        RealVectorValue _vel_bc(_vel_x, _vel_y_bc, 0.);
        
        switch (_eqn_type) {
            case CONTINUITY:
                return _area[_qp]*_rho_bc*(_vel_bc*_normals[_qp])*_test[_i][_qp];
//                break;
            case XMOMENTUM:
                return _area[_qp]*( _rho_bc*_vel_x*(_vel_bc*_normals[_qp]) +_p_bc*_normals[_qp](0))*_test[_i][_qp];
//                break;
            case YMOMENTUM:
                return _area[_qp]*( _rho_bc*_vel_y_bc*(_vel_bc*_normals[_qp]) +_p_bc*_normals[_qp](1) )*_test[_i][_qp];
//                break;
            case ENERGY:
                _e_bc = _eos.e_from_p_rho(_p_bc, _rho_bc);
                _norm_vel2 = _vel_bc.size_sq();
                _rhoE_bc = _rho_bc*(_e_bc + 0.5*_norm_vel2);
                return _area[_qp]*( (_rhoE_bc+_p_bc)*(_vel_bc*_normals[_qp]))*_test[_i][_qp];
//                break;
            default:
                mooseError("The equation with name: \"" << _eqn_name << "\" is not supported in the \"EelStaticPandTBC\" type of boundary condition.");
        }
    }
    else // outlet
    {
        //Real _Mach = _vel_vec.size() / std::sqrt(_eos.c2_from_p_rho(_rho_aux[_qp], _press[_qp]));
        Real press_old = _eos.pressure(_rhoA_old[_qp]/_area[_qp], _vel_vec_old.size(), _rhoEA_old[_qp]/_area[_qp]);
        Real _Mach = _vel_vec_old.size() / std::sqrt(_eos.c2_from_p_rho(_rhoA_old[_qp]/_area[_qp], press_old));
        //if (_t > 0.001) {
        //    std::cout<<"Mach="<<_Mach<<std::endl;
        //    std::cout<<"pressure="<<_press[_qp]<<std::endl;
        //}
        Real _p_bc2 = _p_bc;
        //std::cout<<"pbc2="<<_p_bc2<<std::endl;
        if (_Mach > 1.) {
            _p_bc2 = _eos.pressure(_rhoA[_qp]/_area[_qp], _vel_vec.size(), _rhoEA[_qp]/_area[_qp]);
        }
        switch (_eqn_type) {
            case CONTINUITY:
                return (_rhouA_x[_qp]*_normals[_qp](0)+_rhouA_y[_qp]*_normals[_qp](1))*_test[_i][_qp];
//                break;
            case XMOMENTUM:
                return (_u[_qp]*(_vel_vec*_normals[_qp])+_area[_qp]*_p_bc2*_normals[_qp](0))*_test[_i][_qp];
//                break;
            case YMOMENTUM:
                return (_u[_qp]*(_vel_vec*_normals[_qp])+_area[_qp]*_p_bc2*_normals[_qp](1))*_test[_i][_qp];
//                break;
            case ENERGY:
                return (_vel_vec*_normals[_qp])*(_u[_qp] + _area[_qp]*_p_bc2)*_test[_i][_qp];
//                break;
            default:
                mooseError("The equation with name: \"" << _eqn_name << "\" is not supported in the \"EelStaticPandTBC\" type of boundary condition.");
        }
    }
}

Real
EelStaticPandTBC::computeQpJacobian()
{
    // Declare some variables used in the switch statement:
    RealVectorValue _rhouA_vec(_rhouA_x[_qp], _rhouA_y[_qp], 0.);
    RealVectorValue _vel_vec = _rhouA_vec / _rhoA[_qp];
    RealVectorValue _vel_vec_old(_rhouA_x_old[_qp]/_rhoA_old[_qp], _rhouA_y_old[_qp]/_rhoA_old[_qp], 0.);
    Real _press_term = 0.;
    
    // Compute Mach number:
    Real _rho = _rhoA[_qp]/_area[_qp];
    Real _rhoE = _rhoEA[_qp]/_area[_qp];
    //Real _pressure = _eos.pressure(_rho, _vel_vec.size(), _rhoE);
    Real press_old = _eos.pressure(_rhoA_old[_qp]/_area[_qp], _vel_vec_old.size(), _rhoEA_old[_qp]/_area[_qp]);
    //Real _Mach = _vel_vec.size() / std::sqrt(_eos.c2_from_p_rho(_rho, _pressure));
    Real _Mach = _vel_vec_old.size() / std::sqrt(_eos.c2_from_p_rho(_rhoA_old[_qp]/_area[_qp], press_old));
    bool _mach_bool = _Mach >= 1 ? true : false;
    
    // Return the value of the jacobian matrix:
    switch (_eqn_type) {
        case CONTINUITY:
            return 0.;
//            break;
        case XMOMENTUM:
            if (_vel_vec*_normals[_qp] < 0) { // Inlet
                return _phi[_j][_qp]/_rhoA[_qp]*(2*_rhouA_x[_qp]*_normals[_qp](0)+_rhouA_y[_qp]*_normals[_qp](1))*_test[_i][_qp];
            }
            else { // Outlet
                _press_term = (double)_mach_bool*_eos.dAp_drhouA(_rhoA[_qp], _rhouA_x[_qp], _rhoEA[_qp]);
                return _phi[_j][_qp]/_rhoA[_qp]*((2*_rhouA_x[_qp]+_press_term)*_normals[_qp](0)+_rhouA_y[_qp]*_normals[_qp](1))*_test[_i][_qp];
            }
//            break;
        case YMOMENTUM:
            if (_vel_vec*_normals[_qp] < 0) { // Inlet
                return _phi[_j][_qp]/_rhoA[_qp]*(_rhouA_x[_qp]*_normals[_qp](0)+2*_rhouA_y[_qp]*_normals[_qp](1))*_test[_i][_qp];
            }
            else { // Outlet
                _press_term = (double)_mach_bool*_eos.dAp_drhouA(_rhoA[_qp], _rhouA_y[_qp], _rhoEA[_qp]);
                return _phi[_j][_qp]/_rhoA[_qp]*(_rhouA_x[_qp]*_normals[_qp](0)+(2*_rhouA_y[_qp]+_press_term)*_normals[_qp](1))*_test[_i][_qp];
            }
//            break;
        case ENERGY:
            if (_vel_vec*_normals[_qp]) { // inlet
                return _phi[_j][_qp]*_vel_vec*_normals[_qp]*_test[_i][_qp];
            }
            else { // outlet
                _press_term = (double)_mach_bool*_eos.dAp_drhoEA(_rhoA[_qp], _rhouA_vec.size(), _rhoEA[_qp]);
                return _phi[_j][_qp]*_vel_vec*_normals[_qp]*(1+_press_term)*_test[_i][_qp];
            }
//            break;
        default:
            mooseError("The equation with name: \"" << _eqn_name << "\" is not supported in the \"EelStaticPandTBC\" type of boundary condition.");
            //return 0.;
    }
}

Real
EelStaticPandTBC::computeQpOffDiagJacobian(unsigned _jvar)
{
    // Declare some variables used in the switch statement:
    RealVectorValue _rhouA_vec(_rhouA_x[_qp], _rhouA_y[_qp], 0.);
    RealVectorValue _vel_vec = _rhouA_vec / _rhoA[_qp];
    RealVectorValue _vel_vec_old(_rhouA_x_old[_qp]/_rhoA_old[_qp], _rhouA_y_old[_qp]/_rhoA_old[_qp], 0.);
    Real _press_term = 0.;
    
    // Compute Mach number:
    Real _rho = _rhoA[_qp]/_area[_qp];
    Real _rhoE = _rhoEA[_qp]/_area[_qp];
    Real _pressure = _eos.pressure(_rho, _vel_vec.size(), _rhoE);
    Real press_old = _eos.pressure(_rhoA_old[_qp]/_area[_qp], _vel_vec_old.size(), _rhoEA_old[_qp]/_area[_qp]);
    //Real _Mach = _vel_vec.size() / std::sqrt(_eos.c2_from_p_rho(_rho, _pressure));
    Real _Mach = _vel_vec_old.size() / std::sqrt(_eos.c2_from_p_rho(_rhoA_old[_qp]/_area[_qp], press_old));
    bool _mach_bool = _Mach >= 1 ? true : false;
 
    // Return the value of the jacobian matrix:
    switch (_eqn_type) {
        case CONTINUITY:
            if ((_vel_vec*_normals[_qp] < 0)) { // Inlet
                if (_jvar == _rhouA_x_nb)
                    return _phi[_j][_qp]*_normals[_qp](0)*_test[_i][_qp];
                else if (_jvar == _rhouA_y_nb)
                    return _phi[_j][_qp]*_normals[_qp](1)*_test[_i][_qp];
                else
                    return 0.;
            }
            else { // Outlet
                if (_jvar == _rhouA_x_nb)
                    return _phi[_j][_qp]*_normals[_qp](0)*_test[_i][_qp];
                else if (_jvar == _rhouA_y_nb)
                    return _phi[_j][_qp]*_normals[_qp](1)*_test[_i][_qp];
                else
                    return 0.;
            }
//            break;
        case XMOMENTUM:
            if (_vel_vec*_normals[_qp] < 0) { // Inlet
                if (_jvar == _rhoA_nb) {
                    return -_phi[_j][_qp]*_rhouA_x[_qp]*_vel_vec*_normals[_qp]/_rhoA[_qp]*_test[_i][_qp];
                }
                else if (_jvar == _rhouA_y_nb) {
                    return _phi[_j][_qp]*_rhouA_x[_qp]/_rhoA[_qp]*_normals[_qp](1)*_test[_i][_qp];
                }
                else if (_jvar == _rhoEA_nb) {
                    return 0.;
                }
                else
                    return 0.;
            }
            else { // Outlet
                if (_jvar == _rhoA_nb) {
                    _press_term = (double)_mach_bool*_eos.dAp_drhoA(_rhoA[_qp], _rhouA_vec.size(), _rhoEA[_qp]);
                    return _phi[_j][_qp]*(-_vel_vec*_normals[_qp]/_rhoA[_qp] + _press_term)*_test[_i][_qp];
                }
                else if (_jvar == _rhouA_y_nb) {
                    _press_term = (double)_mach_bool*_eos.dAp_drhouA(_rhoA[_qp], _rhouA_y[_qp], _rhoEA[_qp]);
                    return _phi[_j][_qp]*(_rhouA_x[_qp]/_rhoA[_qp]*_normals[_qp](1) + _press_term*_normals[_qp](0))*_test[_i][_qp];
                }
                else if (_jvar == _rhoEA_nb) {
                    _press_term = (double)_mach_bool*_eos.dAp_drhoEA(_rhoA[_qp], _rhouA_vec.size(), _rhoEA[_qp]);
                    return _phi[_j][_qp]*_press_term*_normals[_qp](0)*_test[_i][_qp];
                }
                else
                    return 0.;
            }
//            break;
        case YMOMENTUM:
            if (_vel_vec*_normals[_qp] < 0) { // Inlet
                if (_jvar == _rhoA_nb) {
                    return -_phi[_j][_qp]*_rhouA_y[_qp]*_vel_vec*_normals[_qp]/_rhoA[_qp]*_test[_i][_qp];
                }
                else if (_jvar == _rhouA_x_nb) {
                    return _phi[_j][_qp]*_rhouA_y[_qp]/_rhoA[_qp]*_normals[_qp](0)*_test[_i][_qp];
                }
                else if (_jvar == _rhoEA_nb) {
                    return 0.;
                }
                else
                    return 0.;
            }
            else { // Outlet
                if (_jvar == _rhoA_nb) {
                    _press_term = (double)_mach_bool*_eos.dAp_drhoA(_rhoA[_qp], _rhouA_vec.size(), _rhoEA[_qp]);
                    return _phi[_j][_qp]*(-_vel_vec*_normals[_qp]/_rhoA[_qp] + _press_term)*_test[_i][_qp];
                }
                else if (_jvar == _rhouA_x_nb) {
                    _press_term = (double)_mach_bool*_eos.dAp_drhouA(_rhoA[_qp], _rhouA_x[_qp], _rhoEA[_qp]);
                    return _phi[_j][_qp]*(_rhouA_y[_qp]/_rhoA[_qp]*_normals[_qp](0) + _press_term*_normals[_qp](1))*_test[_i][_qp];
                }
                else if (_jvar == _rhoEA_nb) {
                    _press_term = (double)_mach_bool*_eos.dAp_drhoEA(_rhoA[_qp], _rhouA_vec.size(), _rhoEA[_qp]);
                    return _phi[_j][_qp]*_press_term*_normals[_qp](1)*_test[_i][_qp];
                }
                else
                    return 0.;
            }
//            break;
        case ENERGY:
            if (_vel_vec*_normals[_qp]) { // inlet
                if (_jvar == _rhoA_nb) {
                    return -_phi[_j][_qp]*_vel_vec*_normals[_qp]/_rhoA[_qp]*(_rhoEA[_qp]+_pressure)*_test[_i][_qp];
                }
                else if (_jvar == _rhouA_x_nb) {
                    return _phi[_j][_qp]*_normals[_qp](0)*(_rhoEA[_qp]+_pressure)/_rhoA[_qp]*_test[_i][_qp];
                }
                else if (_jvar == _rhouA_y_nb) {
                    return _phi[_j][_qp]*_normals[_qp](1)*(_rhoEA[_qp]+_pressure)/_rhoA[_qp]*_test[_i][_qp];
                }
                else
                    return 0.;
            }
            else { // Outlet
                if (_jvar == _rhoA_nb) {
                    _press_term = (double)_mach_bool*_eos.dAp_drhoA(_rhoA[_qp], _rhouA_vec.size(), _rhoEA[_qp]);
                    return _phi[_j][_qp]*_test[_i][_qp];
                }
                else if (_jvar == _rhouA_x_nb) {
                    _press_term = (double)_mach_bool*_eos.dAp_drhouA(_rhoA[_qp], _rhouA_x[_qp], _rhoEA[_qp]);
                    return _phi[_j][_qp]*( _normals[_qp](0)*(_rhoEA[_qp]+_pressure)/_rhoA[_qp] + _vel_vec*_normals[_qp]*_press_term )*_test[_i][_qp];
                }
                else if (_jvar == _rhouA_y_nb) {
                    _press_term = (double)_mach_bool*_eos.dAp_drhoEA(_rhoA[_qp], _rhouA_y[_qp], _rhoEA[_qp]);
                    return _phi[_j][_qp]*( _normals[_qp](1)*(_rhoEA[_qp]+_pressure)/_rhoA[_qp] + _vel_vec*_normals[_qp]*_press_term )*_test[_i][_qp];
                }
                else
                    return 0.;
            }
//            break;
        default:
            mooseError("The equation with name: \"" << _eqn_name << "\" is not supported in the \"EelStaticPandTBC\" type of boundary condition.");
    }
}
