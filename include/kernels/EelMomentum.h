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

#ifndef EELMOMENTUM_H
#define EELMOMENTUM_H

#include "Kernel.h"
#include "EquationOfState.h"

// Forward Declarations
class EelMomentum;

template<>
InputParameters validParams<EelMomentum>();

class EelMomentum : public Kernel
{
public:

  EelMomentum(const std::string & name,
             InputParameters parameters);

protected:

  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian( unsigned int jvar );

private:
    // Aux variables:
    VariableValue & _rhouA_x;
    VariableValue & _rhouA_y;
    VariableValue & _rhouA_z;
    VariableValue & _pressure;
    VariableValue & _rhoA;
    VariableValue & _rhoEA;
    VariableValue & _area;
    VariableGradient & _grad_area;
    
    // Equation of state:
    const EquationOfState & _eos;
    
    // Parameters:
    int _component;
    Real _friction;
    Real _Dh;
    RealVectorValue _gravity;
    
    // Parameters for jacobian:
    unsigned int _rhoA_nb;
    unsigned int _rhouA_x_nb;
    unsigned int _rhouA_y_nb;
    unsigned int _rhouA_z_nb;
    unsigned int _rhoEA_nb;
};

#endif // EELMOMENTUM_H
