#include "EelFluxBC.h"

template<>
InputParameters validParams<EelFluxBC>()
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
    // Make the name of the EOS function a required parameter.
    params.addRequiredParam<UserObjectName>("eos", "The name of equation of state object to use.");

  return params;
}

EelFluxBC::EelFluxBC(const std::string & name, InputParameters parameters) :
    IntegratedBC(name, parameters),
    // Type of equation:
    _eqn_name(getParam<std::string>("equation_name")),
    _eqn_type("CONTINUITY, XMOMENTUM, YMOMENTUM, ZMOMENTUM, ENERGY, INVALID", _eqn_name),
    // Coupled aux variables:
    _rhoA(coupledValue("rhoA")),
    _rhouA_x(coupledValue("rhouA_x")),
    _rhouA_y(isCoupled("rhouA_y") ? coupledValue("rhouA_y") : _zero),
    _rhoEA(coupledValue("rhoEA")),
    _grad_rhoA(coupledGradient("rhoA")),
    _grad_rhouA(coupledGradient("rhouA_x")),
    _grad_rhovA(coupledGradient("rhouA_y")),
    _grad_rhoEA(coupledGradient("rhoEA")),
    // Coupled aux variables:
    _area(coupledValue("area")),
    // Material property:
    _kappa(getMaterialProperty<Real>("kappa")),
    _mu(getMaterialProperty<Real>("mu")),
    // Equation of state:
    _eos(getUserObject<EquationOfState>("eos")),
    // Parameters for jacobian matrix:
    _rhoA_nb(coupled("rhoA")),
    _rhouA_x_nb(coupled("rhouA_x")),
    _rhouA_y_nb(isCoupled("rhouA_y") ? coupled("rhouA_y") : -1),
    _rhoEA_nb(coupled("rhoEA"))
{
}

Real
EelFluxBC::computeQpResidual()
{
    // Compute the velocity vector:
    RealVectorValue _vel_vec(_rhouA_x[_qp]/_rhoA[_qp], _rhouA_y[_qp]/_rhoA[_qp], 0.);
    Real _rho = _rhoA[_qp]/_area[_qp];
    Real _rhoE = _rhoEA[_qp]/_area[_qp];
    Real _pressure = _eos.pressure(_rho, _vel_vec.size(), _rhoE);
    //std::cout<<"press="<<_pressure<<std::endl;
    //std::cout<<"vel="<<_vel_vec<<std::endl;
    
    // Switch statement on equation type:
    switch (_eqn_type)
    {
        case CONTINUITY:
//            return ( _u[_qp]*_vel_vec*_normals[_qp] - _kappa[_qp]*_grad_rhoA[_qp]*_normals[_qp] ) * _test[_i][_qp];
            return _u[_qp]*_vel_vec*_normals[_qp] * _test[_i][_qp];
//            break;
        case XMOMENTUM:
//            return ( _u[_qp]*( _vel_vec*_normals[_qp] ) + _area[_qp]*_pressure*_normals[_qp](0) - _kappa[_qp]*_grad_rhouA[_qp]*_normals[_qp]) * _test[_i][_qp];
            return ( _u[_qp]*( _vel_vec*_normals[_qp] ) + _area[_qp]*_pressure*_normals[_qp](0) ) * _test[_i][_qp];
//            break;
        case YMOMENTUM:
//            return ( _u[_qp]*( _vel_vec*_normals[_qp] ) + _area[_qp]*_pressure*_normals[_qp](1) - _kappa[_qp]*_grad_rhovA[_qp]*_normals[_qp]) * _test[_i][_qp];
            return ( _u[_qp]*( _vel_vec*_normals[_qp] ) + _area[_qp]*_pressure*_normals[_qp](1) ) * _test[_i][_qp];
//            break;
        case ENERGY:
//            return  (_vel_vec*_normals[_qp]*(_u[_qp]+_area[_qp]*_pressure) - _kappa[_qp]*_grad_rhoEA[_qp]*_normals[_qp]) * _test[_i][_qp];
            return  _vel_vec*_normals[_qp]*(_u[_qp]+_area[_qp]*_pressure) * _test[_i][_qp];
//            break;
        default:
            mooseError("The equation with name: \"" << _eqn_name << "\" is not supported in the \"OneD7EqnStagnationPandTBC\" type of boundary condition.");
    }
}

Real
EelFluxBC::computeQpJacobian()
{
    // Parameters used in the jacobian matrix:
    RealVectorValue _vel_vec(0., 0., 0.);
    RealVectorValue _rhouA_vec(0., 0., 0.);
    Real _press_term = 0.;
    
    // Comtribution from dissipative terms:
    Real _diss_jac = 0;//_kappa[_qp]*_grad_phi[_j][_qp]*_normals[_qp]*_test[_i][_qp];
    
    // Switch statement on equation type:
    switch (_eqn_type)
    {
        case CONTINUITY:
            return -_diss_jac;
//            break;
        case XMOMENTUM:
            // Local parameter:
            _vel_vec(0) = 2*_rhouA_x[_qp]/_rhoA[_qp];
            _vel_vec(1) = _rhouA_y[_qp]/_rhoA[_qp];
            _press_term = _eos.dAp_drhouA(_rhoA[_qp], _u[_qp], _rhoEA[_qp]);
            // Return values:
            return _phi[_j][_qp] * ( ( _vel_vec*_normals[_qp] ) + _press_term*_normals[_qp](0) ) * _test[_i][_qp] - _diss_jac;
//            break;
        case YMOMENTUM:
            // Local parameter:
            _vel_vec(0) = _rhouA_x[_qp]/_rhoA[_qp];
            _vel_vec(1) = 2*_rhouA_y[_qp]/_rhoA[_qp];
            _press_term = _eos.dAp_drhouA(_rhoA[_qp], _u[_qp], _rhoEA[_qp]);
            // Return values:
            return _phi[_j][_qp] * ( ( _vel_vec*_normals[_qp] ) + _press_term*_normals[_qp](1) ) * _test[_i][_qp] - _diss_jac;
//            break;
        case ENERGY:
            // Local parameter:
            _rhouA_vec(0) = _rhouA_x[_qp];
            _rhouA_vec(1) = _rhouA_y[_qp];
            _press_term = _eos.dAp_drhoEA(_rhoA[_qp], _rhouA_vec.size(), _rhoEA[_qp]);
            // Return value:
            return  _phi[_j][_qp] * (_vel_vec*_normals[_qp])*(1+_press_term) * _test[_i][_qp] - _diss_jac;
//            break;
        default:
            mooseError("The equation with name: \"" << _eqn_name << "\" is not supported in the \"OneD7EqnStagnationPandTBC\" type of boundary condition.");
    }
}

Real
EelFluxBC::computeQpOffDiagJacobian(unsigned _jvar)
{
    // Parameters used in the jacobian matrix:
    RealVectorValue _rhouA_vec(_rhouA_x[_qp], _rhouA_y[_qp], 0.);
    RealVectorValue _vel_vec = _rhouA_vec / _rhoA[_qp];
    Real _rho = _rhoA[_qp]/_area[_qp];
    Real _rhoE = _rhoEA[_qp]/_area[_qp];
    Real _pressure = _eos.pressure(_rho, _vel_vec.size(), _rhoE);
    Real _press_term = 0.;
    
    // Switch statement on equation type:
    switch (_eqn_type)
    {
        case CONTINUITY:
            if (_jvar == _rhouA_x_nb)
                return _phi[_j][_qp]*_normals[_qp](0)*_test[_i][_qp];
            if (_jvar == _rhouA_y_nb)
                return _phi[_j][_qp]*_normals[_qp](1)*_test[_i][_qp];
            else
                return 0.;
//            break;
        case XMOMENTUM:
            if (_jvar == _rhoA_nb) {
                // Local parameters:
                _press_term = _eos.dAp_drhoA(_rhoA[_qp], _rhouA_vec.size(), _rhoEA[_qp]);
                // Return value:
                return _phi[_j][_qp]*(-_rhouA_x[_qp]*_vel_vec*_normals[_qp]/_rhoA[_qp]+_press_term*_normals[_qp](0))*_test[_i][_qp];
            }
            else if (_jvar == _rhouA_y_nb) {
                // Local parameters:
                _press_term = _eos.dAp_drhouA(_rhoA[_qp], _rhouA_y[_qp], _rhoEA[_qp]);
                // Return value:
                return _phi[_j][_qp]*(_rhouA_x[_qp]*_normals[_qp](1)+_press_term*_normals[_qp](0))*_test[_i][_qp];
            }
            else if (_jvar == _rhoEA_nb) {
                return _phi[_j][_qp]*_eos.dAp_drhoEA(_rhoA[_qp], _rhouA_vec.size(), _rhoEA[_qp])*_normals[_qp](0)*_test[_i][_qp];
            }
            else
                return 0.;
//            break;
        case YMOMENTUM:
            if (_jvar == _rhoA_nb) {
                // Local parameters:
                _press_term = _eos.dAp_drhoA(_rhoA[_qp], _rhouA_vec.size(), _rhoEA[_qp]);
                // Return value:
                return _phi[_j][_qp]*(-_rhouA_y[_qp]*_vel_vec*_normals[_qp]/_rhoA[_qp]+_press_term*_normals[_qp](1))*_test[_i][_qp];
            }
            else if (_jvar == _rhouA_x_nb) {
                // Local parameters:
                _press_term = _eos.dAp_drhouA(_rhoA[_qp], _rhouA_x[_qp], _rhoEA[_qp]);
                // Return value:
                return _phi[_j][_qp]*(_rhouA_y[_qp]*_normals[_qp](0)+_press_term*_normals[_qp](1))*_test[_i][_qp];
            }
            else if (_jvar == _rhoEA_nb) {
                return _phi[_j][_qp]*_eos.dAp_drhoEA(_rhoA[_qp], _rhouA_vec.size(), _rhoEA[_qp])*_normals[_qp](1)*_test[_i][_qp];
            }
            else
                return 0.;
//            break;
        case ENERGY:
            // Local parameter:
            _rhouA_vec(0) = _rhouA_x[_qp];
            _rhouA_vec(1) = _rhouA_y[_qp];
            if (_jvar == _rhoA_nb) {
                // Local parameters:
                _press_term = _eos.dAp_drhoA(_rhoA[_qp], _rhouA_vec.size(), _rhoEA[_qp]);
                // Return value:
                return _phi[_j][_qp]*_vel_vec*_normals[_qp]*(-(_u[_qp]+_pressure)/_rhoA[_qp]+_press_term)*_test[_i][_qp];
            }
            else if (_jvar == _rhouA_x_nb) {
                // Local parameters:
                _press_term = _eos.dAp_drhouA(_rhoA[_qp], _rhouA_x[_qp], _rhoEA[_qp]);
                // Return value:
                return _phi[_j][_qp]*(_normals[_qp](0)/_rhoA[_qp]*(_rhoEA[_qp]+_pressure)+_vel_vec*_normals[_qp]*_press_term)*_test[_i][_qp];
            }
            else if (_jvar == _rhouA_y_nb) {
                // Local parameters:
                _press_term = _eos.dAp_drhouA(_rhoA[_qp], _rhouA_y[_qp], _rhoEA[_qp]);
                // Return value:
                return _phi[_j][_qp]*(_normals[_qp](1)/_rhoA[_qp]*(_rhoEA[_qp]+_pressure)+_vel_vec*_normals[_qp]*_press_term)*_test[_i][_qp];
            }
            else
                return 0.;
//            break;
        default:
            mooseError("The equation with name: \"" << _eqn_name << "\" is not supported in the \"OneD7EqnStagnationPandTBC\" type of boundary condition.");
    }
}
