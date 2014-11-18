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

#ifndef EELFANNOFLOW_H
#define EELFANNOFLOW_H

#include "Kernel.h"
#include "EquationOfState.h"

class EelFannoFlow;

template<>
InputParameters validParams<EelFannoFlow>();
class EelFannoFlow : public Kernel
{
public:

  EelFannoFlow(const std::string & name,
             InputParameters parameters);

protected:
 
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian( unsigned int jvar );

private:
    // Parameters:
    const Real & _f;
    const Real & _Dh;
    // Equation of state:
    const EquationOfState & _eos;
};

#endif // EelFannoFlow_H
