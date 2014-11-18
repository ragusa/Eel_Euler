#include "TaitEOS.h"
#include "MooseError.h"

/*
 Tait equation of state:
 P = P0 * [ (\rho / \rho0)^gamma - 1 ] + P1 and e = Cv * (T - T0) + e0
 */
template<>
InputParameters validParams<TaitEOS>()
{
  InputParameters params = validParams<EquationOfState>();
    params.addParam<Real>("gamma", 0, "  gamma");
    params.addParam<Real>("P0", 0, "  P0 ");
    params.addParam<Real>("P1", 0, "  P1 ");
    params.addParam<Real>("rho0", 0, "  rho0 ");
    params.addParam<Real>("e0", 0, "  e0 ");
    params.addParam<Real>("T0", 0, "  T0 ");
    params.addParam<Real>("Cv", 0, "  heat capacity at constant volume");
    return params;
}

TaitEOS::TaitEOS(const std::string & name, InputParameters parameters) :
  EquationOfState(name, parameters),
    _gamma(getParam<Real>("gamma")),
    _P0(getParam<Real>("P0")),
    _P1(getParam<Real>("P1")),
    _rho0(getParam<Real>("rho0")),
    _e0(getParam<Real>("e0")),
    _T0(getParam<Real>("T0")),
    _Cv(getParam<Real>("Cv"))
{}

TaitEOS::~TaitEOS()
{
}

void
TaitEOS::destroy()
{
}

Real TaitEOS::pressure(Real rho, Real vel_norm, Real rhoE) const
{
    return ( _P0*(std::pow(rho / _rho0, _gamma) - 1.) + _P1);
}

Real TaitEOS::rho_from_p_T(Real pressure, Real temperature) const
{
    return ( _rho0 * std::pow( (pressure-_P1) / _P0 + 1., 1./_gamma) );
}

Real TaitEOS::e_from_p_rho(Real pressure, Real rho) const
{
    mooseError("Not implemented yet.");
    return 0.;
}

Real TaitEOS::temperature_from_p_rho(Real pressure, Real rho) const
{
    mooseError("Not implemented yet.");
    return 0.;
}

Real TaitEOS::c2_from_p_rho(Real rho, Real pressure) const
{
    return ( _gamma * pressure / rho );
}

Real TaitEOS::dAp_drhoA(Real rhoA, Real rhouA_norm, Real rhoEA) const
{
    // Return the value:
    return _gamma * _P0 * std::pow(rhoA / _rho0, _gamma-1.);
}

Real TaitEOS::dAp_drhouA(Real rhoA, Real rhouA_component, Real rhoEA) const
{
    return 0.;
}

Real TaitEOS::dAp_drhoEA(Real rhoA, Real rhouA_norm, Real rhoEA) const
{
    return 0.;
}
