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

#ifndef TOTALENERGYAUX_H
#define TOTALENERGYAUX_H

#include "AuxKernel.h"


//Forward Declarations
class TotalEnergyAux;

template<>
InputParameters validParams<TotalEnergyAux>();

/**
 * Coupled auxiliary value
 */
class TotalEnergyAux : public AuxKernel
{
public:

  TotalEnergyAux(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeValue();

  VariableValue & _rhoEA;
  VariableValue & _area;
};

#endif //TOTALENERGYAUX_H
