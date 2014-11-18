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

#ifndef RAYLEIGHFANNOFLOW_H
#define RAYLEIGHFANNOFLOW_H


#include "ODEKernel.h"
#include "EquationOfState.h"

// Forward Declarations
class RayleighFannoFlow;

template<>
InputParameters validParams<RayleighFannoFlow>();

class RayleighFannoFlow : public ODEKernel
{
public: 
  RayleighFannoFlow(const std::string & name, InputParameters parameters);

protected:
    
    virtual Real computeQpResidual();

    virtual Real computeQpJacobian();

    virtual Real computeQpOffDiagJacobian(unsigned int jvar);

    // Friction parameter:
    const Real & _f;
    // Hydraulic parameter:
    const Real & _Dh;
    // Equation of state:
    const EquationOfState & _eos;
};


#endif /* RayleighFannoFlow_H */
