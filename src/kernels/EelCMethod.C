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

#include "EelCMethod.h"
/**
This function is based on the C-method theory. It computes the 
 */
template<>
InputParameters validParams<EelCMethod>()
{
  InputParameters params = validParams<Kernel>();
    // Coupled aux variables:
    params.addRequiredCoupledVar("pressure", "pressure");
    // Parameter for diffusion term:
    params.addParam<double>("kappa", 1., "Parameter for diffusion term: kappa.");
    // name of the pps computing max of eigenvalues:
    params.addRequiredParam<std::string>("max_eig_pps", "pps computing the max of eigenvalues.");
    params.addRequiredParam<std::string>("max_grad_pps", "pps computing the max of ||grad(P)||.");
  return params;
}

EelCMethod::EelCMethod(const std::string & name,
                       InputParameters parameters) :
  Kernel(name, parameters),
    // Coupled aux variables:
    _grad_press(coupledGradient("pressure")),
    // Parameters for diffution term:
    _kappa(getParam<double>("kappa")),
    // name of the pps:
    _max_eig_pps_name(getParam<std::string>("max_eig_pps")),
    _max_grad_press_pps_name(getParam<std::string>("max_grad_pps"))
{
}

Real EelCMethod::computeQpResidual()
{
    // Initialize some variables:
    Real _hmax = _current_elem->hmax();
    // pps:
    //std::cout<<_max_grad_press_pps_name<<std::endl;
    Real _Smax = getPostprocessorValue(_max_eig_pps_name);
    Real _gradP_max = getPostprocessorValue(_max_grad_press_pps_name);
    //std::cout << "kernel gradP_max=" << _gradP_max << std::endl;
    if (std::fabs(_gradP_max) < 1e-10)
        _gradP_max = 1.;
    //std::cout << "kernel gradP_max=" << _gradP_max << std::endl;
    // Compute the norm of grad(pressure):
    Real _norm_grad_press2 = _grad_press[_qp](0)*_grad_press[_qp](0); // + _grad_press[_qp](1)*_grad_press[_qp](1) + _grad_press[_qp](2)*_grad_press[_qp](2);
    //std::cout<<"kernel gradP="<<std::sqrt(_norm_grad_press2)<<std::endl;
    // Compute source term:
    Real _source = _Smax * std::fabs(_grad_press[_qp](0)) / _gradP_max;
    //std::cout << _Smax << std::endl;
    //std::cout << std::fabs(_grad_press[_qp](0)) / _gradP_max << std::endl;
    /// Returns the residual
    return (_Smax*_u[_qp]/_hmax - _source)*_test[_i][_qp] + _Smax*_hmax*_kappa*_grad_u[_qp](0)*_grad_test[_i][_qp](0);
}

Real EelCMethod::computeQpJacobian()
{
    return 0;
}

Real EelCMethod::computeQpOffDiagJacobian( unsigned int _jvar)
{
    return 0;
}
