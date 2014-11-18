#include "EquationOfState.h"
#include "MooseError.h"

template<>
InputParameters validParams<EquationOfState>()
{
  InputParameters params = validParams<UserObject>();
    params.addParam<Real>("gamma", 0, "  gamma");
    params.addParam<Real>("Pinf", 0, "  P infinity");
    params.addParam<Real>("q", 0, "  q coefficient");
    params.addParam<Real>("q_prime", 0, "  q' coefficient");
    params.addParam<Real>("Cv", 0, "  heat capacity at constant volume");
    return params;
}

EquationOfState::EquationOfState(const std::string & name, InputParameters parameters) :
  GeneralUserObject(name, parameters),
    _gamma(getParam<Real>("gamma")),
    _Pinf(getParam<Real>("Pinf")),
    _qcoeff(getParam<Real>("q")),
    _qcoeff_prime(getParam<Real>("q_prime")),
    _Cv(getParam<Real>("Cv"))
{}

EquationOfState::~EquationOfState()
{
  // Destructor, empty
}

void
EquationOfState::destroy()
{
}

Real EquationOfState::pressure(Real rho, Real vel_norm, Real rhoE) const
{
    this->error_not_implemented("pressure");
    return 0.;
}

Real EquationOfState::rho_from_p_T(Real pressure, Real temperature) const
{
    this->error_not_implemented("density");
    return 0.;
}

Real EquationOfState::e_from_p_rho(Real pressure, Real rho) const
{
    this->error_not_implemented("internal energy");
    return 0.;
}

Real EquationOfState::temperature_from_p_rho(Real pressure, Real rho) const
{
    this->error_not_implemented("temperature");
    return 0.;
}

Real EquationOfState::c2_from_p_rho(Real rho, Real pressure) const
{
    this->error_not_implemented("speed of sound");
    return 0.;
}

Real EquationOfState::dAp_drhoA(Real rhoA, Real rhouA_norm, Real rhoEA) const
{
    this->error_not_implemented("derivative of pressure with respect to density");
    return 0.;
}

Real EquationOfState::dAp_drhouA(Real rhoA, Real rhouA_component, Real rhoEA) const
{
    this->error_not_implemented("derivative of pressure with respect to momentum");
    return 0.;
}

Real EquationOfState::dAp_drhoEA(Real rhoA, Real rhouA_norm, Real rhoEA) const
{
    this->error_not_implemented("derivative of pressure with respect to total energy");
    return 0.;
}

Real EquationOfState::gamma() const
{
    return _gamma;
}

Real EquationOfState::Pinf() const
{
    return _Pinf;
}

Real EquationOfState::qcoeff() const
{
    return _qcoeff;
}

Real EquationOfState::qcoeff_prime() const
{
    return _qcoeff_prime;
}

Real EquationOfState::Cv() const
{
    return _Cv;
}

void EquationOfState::error_not_implemented(std::string method_name) const
{}
