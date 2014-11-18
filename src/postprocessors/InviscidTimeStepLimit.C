#include "InviscidTimeStepLimit.h"
#include "EquationOfState.h"

template<>
InputParameters validParams<InviscidTimeStepLimit>()
{
  InputParameters params = validParams<ElementPostprocessor>();
  // Coupled variables
  params.addRequiredCoupledVar("vel_mag", "Velocity magnitude");
  params.addRequiredCoupledVar("c", "Sound speed");
  params.addParam<Real>("beta", 0.8, "User supplied constant");

  return params;
}

InviscidTimeStepLimit::InviscidTimeStepLimit(const std::string & name, InputParameters parameters) :
    ElementPostprocessor(name, parameters),
    _dim(_mesh.dimension()),
    _vel_mag(coupledValue("vel_mag")),
    _c(coupledValue("c")),
    _beta(getParam<Real>("beta"))
{
}

InviscidTimeStepLimit::~InviscidTimeStepLimit()
{
}

void
InviscidTimeStepLimit::initialize()
{
  _value = std::numeric_limits<Real>::max();
}

void
InviscidTimeStepLimit::execute()
{
  Real h_min = _current_elem->hmin();
  for (unsigned qp = 0; qp < _qrule->n_points(); ++qp)
    _value = std::min(_value, _beta * h_min / (_vel_mag[qp] + _c[qp]));
}

Real
InviscidTimeStepLimit::getValue()
{
  _communicator.min(_value);
  return _value;
}

void
InviscidTimeStepLimit::threadJoin(const UserObject & uo)
{
  const InviscidTimeStepLimit & pps = dynamic_cast<const InviscidTimeStepLimit &>(uo);
  _value = std::min(_value, pps._value);
}
