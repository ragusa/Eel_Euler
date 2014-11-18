#include "EelMassInflowBC.h"

template<>
InputParameters validParams<EelMassInflowBC>()
{
  InputParameters params = validParams<IntegratedBC>();

  params.addRequiredParam<std::string>("equation_name", "The name of the equation this BC is acting on");
  // Coupled variables:
  params.addRequiredCoupledVar("pressure", "pressure of the fluid");
  params.addCoupledVar("area", "area aux variable");
  // Input parameters:
  params.addRequiredParam<Real>("rhou0_bc", "Specified inflow momentum.");
  params.addRequiredParam<Real>("T0_bc", "Specified inflow temperature.");
  params.addParam<Real>("alpha0_bc", 0., "Inflow angle.");
  // Equation of state:
  params.addRequiredParam<UserObjectName>("eos", "The name of equation of state object to use.");

  return params;
}

EelMassInflowBC::EelMassInflowBC(const std::string & name, InputParameters parameters) :
    IntegratedBC(name, parameters),
   // Name of the equation:
    _eqn_name(getParam<std::string>("equation_name")),
    _eqn_type("CONTINUITY, XMOMENTUM, YMOMENTUM, ENERGY, INVALID", "INVALID"),
   // Coupled variables:
    _pressure(coupledValue("pressure")),
    _area(coupledValue("area")),
   // Boundary condition parameters:
    _rhou0_bc(getParam<Real>("rhou0_bc")),
    _T0_bc(getParam<Real>("T0_bc")),
    _alpha0_bc(getParam<Real>("alpha0_bc")),
    // Equation of state:
    _eos(getUserObject<EquationOfState>("eos"))
{
  _eqn_type = _eqn_name;
}

Real
EelMassInflowBC::computeQpResidual()
{
    /*std::cout << "m_dot_bc=" << _m_dot_bc << std::endl;
    std::cout << "T_bc=" << _T_bc << std::endl;
    std::cout << "gamma_bc=" << _gamma_bc << std::endl;*/
    
    // Compute the density, internal energy, total energy, x-component and y-component of the momentum and velocity vectors, using the bcs parameters:
    Real _rho_bc = _eos.rho_from_p_T(_pressure[_qp], _T0_bc);
    Real _u_bc = cos(_alpha0_bc)*_rhou0_bc/_rho_bc;
    Real _v_bc = sin(_alpha0_bc)*_rhou0_bc/_rho_bc;
    Real _e_bc = _eos.e_from_p_rho(_pressure[_qp], _rho_bc);
    RealVectorValue _vel_bc(_u_bc, _v_bc, 0.);
    Real _E_bc = _e_bc + 0.5*_vel_bc.size_sq();
    
    /*std::cout <<"##########################" << _eqn_name << std::endl;
    std::cout<<"normal(0)="<<_normals[_qp](0)<<std::endl;
    std::cout<<"normal(1)="<<_normals[_qp](1)<<std::endl;
    std::cout<<"normal(2)="<<_normals[_qp](2)<<std::endl;
    std::cout<<"vel x="<<_vel_x_bc<<std::endl;
    std::cout<<"vel y="<<_vel_y_bc<<std::endl;
    std::cout<<"vel x="<<_vel_bc(0)<<std::endl;
    std::cout<<"vel y="<<_vel_bc(1)<<std::endl;
    std::cout<<"vel dot normal="<<_vel_bc*_normals[_qp]<<std::endl;*/
    
    switch (_eqn_type) {
        case CONTINUITY:
            return _area[_qp]*_rhou0_bc*_normals[_qp](0)*_test[_i][_qp];
            break;
        case XMOMENTUM:
            return _area[_qp]*(_rho_bc*_u_bc*_vel_bc*_normals[_qp] +_pressure[_qp]*_normals[_qp](0))*_test[_i][_qp];
            break;
        case YMOMENTUM:
            return _area[_qp]*(_rho_bc*_v_bc*_vel_bc*_normals[_qp] +_pressure[_qp]*_normals[_qp](1))*_test[_i][_qp];
            break;
        case ENERGY:
            return _area[_qp]*(_rho_bc*_E_bc+_pressure[_qp])*(_vel_bc*_normals[_qp])*_test[_i][_qp];
            break;
        default:
            mooseError("The equation with name: \"" << _eqn_name << "\" is not supported in the \"EelMassInflowBC\" type of boundary condition.");
            return 0.;
            break;
        }
}

Real
EelMassInflowBC::computeQpJacobian()
{
  // TODO
  return 0;
}

Real
EelMassInflowBC::computeQpOffDiagJacobian(unsigned jvar)
{
  // TODO
  return 0;
}
