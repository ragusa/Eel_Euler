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

#ifndef EELCMETHOD_H
#define EELCMETHOD_H

#include "Kernel.h"

// Forward Declarations
class EelCMethod;

template<>
InputParameters validParams<EelCMethod>();

class EelCMethod : public Kernel
{
public:

  EelCMethod(const std::string & name,
             InputParameters parameters);

protected:

  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian( unsigned int _jvar);

private:
    /// Coupled aux variables:
    VariableGradient & _grad_press;
    // Parameter for diffusion term:
    double _kappa;
    // Name of pps computing max(eigenvalues):
    std::string _max_eig_pps_name;
    // Name of pps computing max of ||grad(pressure)||
    std::string _max_grad_press_pps_name;
};

#endif // EELCMETHOD_H
