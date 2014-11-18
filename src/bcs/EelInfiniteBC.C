#include "EelInfiniteBC.h"

template<>
InputParameters validParams<EelInfiniteBC>()
{
  InputParameters params = validParams<IntegratedBC>();

  params.addRequiredParam<std::string>("equation_name", "The name of the equation this BC is acting on");
  // Coupled variables:
    params.addRequiredCoupledVar("pressure", "fluid pressure");
    params.addRequiredCoupledVar("velocity_x", "x component of the velocity");
    params.addCoupledVar("velocity_y", "y component of the velocity");
    params.addCoupledVar("velocity_z", "z component of the velocity");
    params.addRequiredCoupledVar("area", "area");

  return params;
}

EelInfiniteBC::EelInfiniteBC(const std::string & name, InputParameters parameters) :
    IntegratedBC(name, parameters),
   // Name of the equation:
    _eqn_name(getParam<std::string>("equation_name")),
    _eqn_type("CONTINUITY, XMOMENTUM, YMOMENTUM, ZMOMENTUM, ENERGY, INVALID", "INVALID"),
   // Coupled variables:
    _pressure(coupledValue("pressure")),
    _vel_x(coupledValue("velocity_x")),
    _vel_y(_mesh.dimension()>=2 ? coupledValue("velocity_y") : _zero),
    _vel_z(_mesh.dimension()==3 ? coupledValue("velocity_z") : _zero),
    _area(coupledValue("area"))
{
  _eqn_type = _eqn_name;
}

Real
EelInfiniteBC::computeQpResidual()
{
    RealVectorValue _vector(_vel_x[_qp], _vel_y[_qp], _vel_z[_qp]);
    RealVectorValue _ones(1., 1., 1.);
    switch (_eqn_type) {
        case CONTINUITY:
            return _vector*_normals[_qp]*_u[_qp]*_test[_i][_qp];
            break;
        case XMOMENTUM:
            return ( _u[_qp]*_vector*_ones + _area[_qp]*_pressure[_qp] )*_normals[_qp](0)*_test[_i][_qp];
            break;
        case YMOMENTUM:
            return ( _u[_qp]*_vector*_ones + _area[_qp]*_pressure[_qp] )*_normals[_qp](1)*_test[_i][_qp];
            break;
        case ZMOMENTUM:
            return ( _u[_qp]*_vector*_ones + _area[_qp]*_pressure[_qp] )*_normals[_qp](2)*_test[_i][_qp];
            break;
        case ENERGY:
            return _vector*_normals[_qp]*( _u[_qp] + _area[_qp]*_pressure[_qp] )*_test[_i][_qp];
            break;
        default:
            mooseError("The equation with name: \"" << _eqn_name << "\" is not supported in the \"EelInfiniteBC\" type of boundary condition.");
            return 0.;
            break;
    }
}

Real
EelInfiniteBC::computeQpJacobian()
{
  // TODO
  return 0;
}

Real
EelInfiniteBC::computeQpOffDiagJacobian(unsigned jvar)
{
  // TODO
  return 0;
}
