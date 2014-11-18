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

#ifndef EELENERGY_H
#define EElENERGY_H

#include "Kernel.h"
#include "EquationOfState.h"
#include "Function.h"

// Forward Declarations
class EelEnergy;

template<>
InputParameters validParams<EelEnergy>();

class EelEnergy : public Kernel
{
public:

  EelEnergy(const std::string & name,
             InputParameters parameters);

protected:

  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian( unsigned int _jvar);

private:
    // Coupled variables
    VariableValue & _rhoA;
    VariableValue & _rhouA_x;
    VariableValue & _rhouA_y;
    VariableValue & _rhouA_z;
    VariableValue & _pressure;
    VariableValue & _area;
    
    // Component:
    std::string _Hw_fn_name;
    std::string _Tw_fn_name;
    const Real & _Hw;
    const Real & _Tw;
    const Real & _aw;
    const RealVectorValue & _gravity;
    
    // Equation of state:
    const EquationOfState & _eos;
    
    // Parameters for jacobian:
    unsigned int _rhoA_nb;
    unsigned int _rhouA_x_nb;
    unsigned int _rhouA_y_nb;
    unsigned int _rhouA_z_nb;
};

#endif // EELENERGY_H
