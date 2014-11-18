#include "EelHRhoUBC.h"

template<>
InputParameters validParams<EelHRhoUBC>()
{
  InputParameters params = validParams<IntegratedBC>();

  params.addRequiredParam<std::string>("equation_name", "The name of the equation this BC is acting on");

    // Coupled conservative variables:
    params.addCoupledVar("rhoA", "density: rhoA");
    params.addCoupledVar("rhouA_x", "x component of the momentum: rhouA_x");
    params.addCoupledVar("rhoEA", "total energy: rho*E*A");
    // Coupled aux variables:
    params.addCoupledVar("area", "Coupled area variable");
    // Input parameters
    params.addRequiredParam<Real>("rhou", "Specified momentum: rhou");
    params.addRequiredParam<Real>("H", "Specified enthalpy");
    // Make the name of the EOS function a required parameter.
    params.addRequiredParam<UserObjectName>("eos", "The name of equation of state object to use.");

  return params;
}

EelHRhoUBC::EelHRhoUBC(const std::string & name, InputParameters parameters) :
    IntegratedBC(name, parameters),
    // Type of equation:
    _eqn_name(getParam<std::string>("equation_name")),
    _eqn_type("VOIDFRACTION, CONTINUITY, XMOMENTUM, ENERGY, INVALID", _eqn_name),
    // Coupled aux variables:
    _rhoA(coupledValue("rhoA")),
    _rhouA_x(coupledValue("rhouA_x")),
    _rhoEA(coupledValue("rhoEA")),
    // Coupled aux variables:
    _area(coupledValue("area")),
    // Stagnation variables:
    _rhou(getParam<Real>("rhou")),
    _H(getParam<Real>("H")),
    // Equation of state:
    _eos(getUserObject<EquationOfState>("eos"))
    // Parameters for jacobian matrix:
//    _rhoA_nb(coupled("rhoA")),
//    _rhouA_x_nb(coupled("rhouA_x")),
//    _rhouA_y_nb(isCoupled("rhouA_y") ? coupled("rhouA_x") : -1),
//    _rhoEA_nb(coupled("rhoEA"))
{}

Real
EelHRhoUBC::computeQpResidual()
{
    Real rho, P;
  switch (_eqn_type)
  {
  case CONTINUITY:
        return _area[_qp] * _rhou * _normals[_qp](0) * _test[_i][_qp];
  case XMOMENTUM:
        rho = _rhoA[_qp] / _area[_qp];
        P = _eos.pressure(rho, _rhouA_x[_qp]/_rhoA[_qp], _rhoEA[_qp]/_area[_qp]);
        return _area[_qp] * ( _rhou * _rhou / rho + P ) * _normals[_qp](0) * _test[_i][_qp];
  case ENERGY:
        return _area[_qp] * _rhou * _H * _normals[_qp](0) * _test[_i][_qp];
  default:
    mooseError("The equation with name: \"" << _eqn_name << "\" is not supported in the \"OneD7EqnStagnationPandTBC\" type of boundary condition.");
  }
}

Real
EelHRhoUBC::computeQpJacobian()
{
    return 0.;
}

Real
EelHRhoUBC::computeQpOffDiagJacobian(unsigned _jvar)
{
    return 0.;
}
