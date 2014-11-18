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

#include "ExactSolAreaVariable.h"

template<>
InputParameters validParams<ExactSolAreaVariable>()
{
  InputParameters params = validParams<Function>();
    params.addRequiredParam<std::string>("variable_name", "The name of the variable");
    params.addRequiredParam<Real>("p0_bc", "Inlet stagnation pressure");
    params.addRequiredParam<Real>("T0_bc", "Inlet stagnation temperature");
    params.addRequiredParam<Real>("p_bc", "Outlet static pressure");
    params.addRequiredParam<Real>("length", "Length of the computational domain.");
//    params.addRequiredParam<UserObjectName>("eos", "The name of equation of state object to use.");
  return params;
}

ExactSolAreaVariable::ExactSolAreaVariable(const std::string & name, InputParameters parameters) :
    Function(name, parameters),
    _var_name(getParam<std::string>("variable_name")),
    _var_type("DENSITY, VELOCITY, PRESSURE, INVALID", _var_name),
    _Po(getParam<Real>("p0_bc")),
    _To(getParam<Real>("T0_bc")),
    _Pout(getParam<Real>("p_bc")),
    _length(getParam<Real>("length"))
//    _eos(getUserObject<EquationOfState>("eos"))
{}

Real
ExactSolAreaVariable::value(Real /*t*/, const Point & p)
{
    // Equation of state parameters:
    Real Pinf = 1.e9;
    Real gamma = 2.35;
    Real qcoeff = -1167.e3;
    Real Cv = 1816.;
    
    // Compute the value of the cross section at point p:
    Real area = 1. + 0.5*std::cos(2*libMesh::pi*p(0));
    Real area_out = 1. + 0.5*std::cos(2*libMesh::pi*_length);
 
    // Compute some stagnation variables:
//    Real rho_o = _eos.rho_from_p_T(_Po, _To);
//    Real eo = _eos.e_from_p_rho(_Po, rho_o);
    Real rho_o = (_Po + Pinf) / ((gamma-1)*Cv*_To);
    Real eo = (_Po + gamma*Pinf) / ((gamma-1)*rho_o);// + qcoeff;
    Real Ho = eo + _Po / rho_o;
    
    // Compute the outflow density and velocity:
//    Real rho_out = rho_o*std::pow( (_Pout+_eos.Pinf() )/( _Po+_eos.Pinf() ), 1./_eos.gamma());
//    Real vel_out = std::sqrt( 2.*( Ho-_eos.gamma()*(_Pout+_eos.Pinf()) )/( (_eos.gamma()-1.)*rho_out ) );
    Real rho_out = rho_o*std::pow( (_Pout+Pinf )/( _Po+Pinf ), 1./gamma);
    Real vel_out = std::sqrt( 2.*( Ho-gamma*(_Pout+Pinf)/( (gamma-1.)*rho_out ) ) );
    Real m_out = rho_out*vel_out*area_out;
    
//    std::cout<<"rho out="<<rho_out<<std::endl;
//    std::cout<<"vel out="<<vel_out<<std::endl;
//    std::cout<<"m out="<<m_out<<std::endl;
    
    // Compute the pressure (non-linear solve)
    Real pressure = _Pout;
    Real rho = rho_out;
    Real vel = vel_out;
    Real diff = 1.;
    while (std::fabs(diff) > 1.e-6) {
        diff = pressure;
//        rho = rho_o*std::pow( (pressure+_eos.Pinf() )/( _Po+_eos.Pinf() ), 1./_eos.gamma());
        rho = rho_o*std::pow( (pressure+Pinf )/( _Po+Pinf ), 1./gamma);
        Real rhs = -0.5 * (m_out/(rho*area)) * (m_out/(rho*area)) + Ho;
//        pressure = (_eos.gamma()-1.)*rho*rhs/_eos.gamma() - _eos.Pinf();
        pressure = (gamma-1.)*rho*rhs/gamma - Pinf;
        diff = diff - pressure;
    }
//    std::cout<<"pressure="<<pressure<<std::endl;
    
    // Compute the density and velocity values:
//    rho = rho_o*std::pow( (pressure+_eos.Pinf() )/( _Po+_eos.Pinf() ), 1./_eos.gamma());
    rho = rho_o*std::pow( (pressure+Pinf )/( _Po+Pinf ), 1./gamma);
    vel = m_out / (rho * area);
    
    // Return the value:
    switch (_var_type) {
        case DENSITY:
            return rho;
        case VELOCITY:
            return vel;
        case PRESSURE:
            return pressure;
        default:
            mooseError("The variable with name: \"" << _var_name << "\" is not supported in the \"EaxctSolAreaVariable\" type of function.");
    }
}

RealVectorValue 
ExactSolAreaVariable::gradient(Real /*t*/, const Point & p)
{
	return 0;
}
