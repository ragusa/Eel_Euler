#include "StiffenedGasEquationOfState.h"
#include "MooseError.h"

template<>
InputParameters validParams<StiffenedGasEquationOfState>()
{
  InputParameters params = validParams<EquationOfState>();
    params.addParam<Real>("gamma", 0, "  gamma");
    params.addParam<Real>("Pinf", 0, "  P infinity");
    params.addParam<Real>("q", 0, "  q coefficient");
    params.addParam<Real>("q_prime", 0, "  q' coefficient");
    params.addParam<Real>("Cv", 0, "  heat capacity at constant volume");
    return params;
}

StiffenedGasEquationOfState::StiffenedGasEquationOfState(const std::string & name, InputParameters parameters) :
  EquationOfState(name, parameters),
    _gamma(getParam<Real>("gamma")),
    _Pinf(getParam<Real>("Pinf")),
    _qcoeff(getParam<Real>("q")),
    _qcoeff_prime(getParam<Real>("q_prime")),
    _Cv(getParam<Real>("Cv"))
{}

StiffenedGasEquationOfState::~StiffenedGasEquationOfState()
{
  // Destructor, empty
}

void
StiffenedGasEquationOfState::destroy()
{
}

Real StiffenedGasEquationOfState::pressure(Real rho, Real vel_norm, Real rhoE) const
{
  Real _e = ( rhoE - 0.5 * rho*vel_norm*vel_norm ) / rho;
  return ( (_gamma-1) * ( _e - _qcoeff) * rho - _gamma * _Pinf );
}

Real StiffenedGasEquationOfState::rho_from_p_T(Real pressure, Real temperature) const
{
    return ( (pressure + _Pinf) / ((_gamma-1)*_Cv*temperature) );
}

Real StiffenedGasEquationOfState::e_from_p_rho(Real pressure, Real rho) const
{
    return ( (pressure + _gamma*_Pinf) / ((_gamma-1)*rho) + _qcoeff );
}

Real StiffenedGasEquationOfState::temperature_from_p_rho(Real pressure, Real rho) const
{
    return ( (pressure + _Pinf) / ((_gamma-1)*_Cv*rho) );
}

Real StiffenedGasEquationOfState::c2_from_p_rho(Real rho, Real pressure) const
{
    return ( _gamma * ( pressure + _Pinf ) / rho );
}

Real StiffenedGasEquationOfState::dAp_drhoA(Real rhoA, Real rhouA_norm, Real rhoEA) const
{
    // Compute the norm of the velocity vector:
    Real _vel_norm = rhouA_norm / rhoA;
    
    // Return the value:
    return 0.5*(_gamma-1)*_vel_norm*_vel_norm - _qcoeff;
}

Real StiffenedGasEquationOfState::dAp_drhouA(Real rhoA, Real rhouA_component, Real rhoEA) const
{
    // Compute the velocity component:
    Real _vel = rhouA_component / rhoA;
    
    // Return the value:
    return -(_gamma-1)*_vel;
}

Real StiffenedGasEquationOfState::dAp_drhoEA(Real rhoA, Real rhouA_norm, Real rhoEA) const
{
    return (_gamma-1);
}
