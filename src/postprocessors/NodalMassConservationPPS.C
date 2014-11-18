#include "NodalMassConservationPPS.h"

template<>
InputParameters validParams<NodalMassConservationPPS>()
{
  InputParameters params = validParams<NodalPostprocessor>();
  params.addCoupledVar("nx", "x-component of the normal");
  params.addCoupledVar("ny", "y-component of the normal");
  params.addCoupledVar("nz", "z-component of the normal");

  params.set<std::vector<VariableName> >("nx") = std::vector<VariableName>(1, "nodal_normal_x");
  params.set<std::vector<VariableName> >("ny") = std::vector<VariableName>(1, "nodal_normal_y");
  params.set<std::vector<VariableName> >("nz") = std::vector<VariableName>(1, "nodal_normal_z");

  params.addRequiredCoupledVar("rho_u", "x-component of the momentum");
  params.addCoupledVar("rho_v", "y-component of the momentum");
  params.addCoupledVar("rho_w", "z-component of the momentum");

  return params;
}

NodalMassConservationPPS::NodalMassConservationPPS(const std::string & name, InputParameters parameters) :
    NodalPostprocessor(name, parameters),
    _value(0.),
    _nx(coupledValue("nx")),
    _ny(coupledValue("ny")),
    _nz(coupledValue("nz")),
    _rho_u(coupledValue("rho_u")),
    _rho_v(isCoupled("rho_v") ? coupledValue("rho_v") : _zero),
    _rho_w(isCoupled("rho_w") ? coupledValue("rho_w") : _zero)
{
}

NodalMassConservationPPS::~NodalMassConservationPPS()
{
}

void
NodalMassConservationPPS::initialize()
{
  _value = 0.;
}

void
NodalMassConservationPPS::execute()
{
  RealVectorValue mom_vec(_rho_u[_qp], _rho_v[_qp], _rho_w[_qp]);
  RealVectorValue n(_nx[_qp], _ny[_qp], _nz[_qp]);
  _value += mom_vec * n;
}

PostprocessorValue
NodalMassConservationPPS::getValue()
{
  gatherSum(_value);
  return _value;
}

void
NodalMassConservationPPS::threadJoin(const UserObject & uo)
{
  const NodalMassConservationPPS & m = dynamic_cast<const NodalMassConservationPPS &>(uo);
  _value += m._value;
}
