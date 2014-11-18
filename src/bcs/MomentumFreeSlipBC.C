#include "MomentumFreeSlipBC.h"

template<>
InputParameters validParams<MomentumFreeSlipBC>()
{
  InputParameters params = validParams<NodalNormalBC>();
  params.addRequiredCoupledVar("rho_u", "x-component of velocity");
  params.addCoupledVar("rho_v", "y-component of velocity");
  params.addCoupledVar("rho_w", "z-component of velocity");

  return params;
}

MomentumFreeSlipBC::MomentumFreeSlipBC(const std::string & name, InputParameters parameters) :
    NodalNormalBC(name, parameters),
    _mesh_dimension(_mesh.dimension()),
    _rho_u(coupledValue("rho_u")),
    _rho_v(_mesh_dimension >= 2 ? coupledValue("rho_v") : _zero),
    _rho_w(_mesh_dimension >= 3 ? coupledValue("rho_w") : _zero)
{
}

MomentumFreeSlipBC::~MomentumFreeSlipBC()
{
}

bool
MomentumFreeSlipBC::shouldApply()
{
  // this prevents zeroing out the row
  return !_sys.currentlyComputingJacobian();
}

void
MomentumFreeSlipBC::computeResidual(NumericVector<Number> & residual)
{
  NonlinearSystem & nl = _fe_problem.getNonlinearSystem();
  NumericVector<Number> & re = nl.residualVector(Moose::KT_NONTIME);
  NumericVector<Number> & ret = nl.residualVector(Moose::KT_TIME);

  _qp = 0;
  if (_var.isNodalDefined())
  {
    if (_mesh_dimension == 1)
    {
      mooseError("Not implemented yet");
    }
    else if (_mesh_dimension == 2)
    {
      MooseVariable & rho_u_var = *getVar("rho_u", 0);
      unsigned int & rho_u_dof_idx = rho_u_var.nodalDofIndex();

      MooseVariable & rho_v_var = *getVar("rho_v", 0);
      unsigned int & rho_v_dof_idx = rho_v_var.nodalDofIndex();

      Real Re_u = residual(rho_u_dof_idx);
      Real Re_v = residual(rho_v_dof_idx);

      Real rho_u_val = ( Re_u * _ny[_qp] * _ny[_qp] - Re_v * _nx[_qp] * _ny[_qp]);
      Real rho_v_val = (-Re_u * _nx[_qp] * _ny[_qp] + Re_v * _nx[_qp] * _nx[_qp]);
//        std::cout<<"nx="<<_nx[_qp]<<std::endl;
//        std::cout<<"ny="<<_ny[_qp]<<std::endl;

      // NOTE: we have to handle all components at the same time, otherwise we'd work with the modified residual and we do not want that
      residual.set(rho_u_dof_idx, rho_u_val);
      residual.set(rho_v_dof_idx, rho_v_val);
    }
    else if (_mesh_dimension == 3)
    {
      mooseError("Not implemented yet");
    }
  }
}

Real
MomentumFreeSlipBC::computeQpResidual()
{
  return 0.;
}
