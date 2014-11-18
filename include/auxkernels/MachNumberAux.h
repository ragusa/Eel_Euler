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

#ifndef MACHNUMBERAUX_H
#define MACHNUMBERAUX_H

#include "AuxKernel.h"
#include "EquationOfState.h"


//Forward Declarations
class MachNumberAux;

template<>
InputParameters validParams<MachNumberAux>();

/**
 * Coupled auxiliary value
 */
class MachNumberAux : public AuxKernel
{
public:

  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  MachNumberAux(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeValue();

    VariableValue & _rhoA;
    VariableValue & _rhouA_x;
    VariableValue & _rhouA_y;
    VariableValue & _rhouA_z;
    VariableValue & _pressure;
    VariableValue & _area;
    const EquationOfState & _eos;
};

#endif //MachNumberAux_H
