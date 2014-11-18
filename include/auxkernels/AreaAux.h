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

#ifndef AREAAUX_H
#define AREAAUX_H

#include "AuxKernel.h"
#include "AreaFunction.h"

//Forward Declarations
class AreaAux;

template<>
InputParameters validParams<AreaAux>();

/**
 * Coupled auxiliary value
 */
class AreaAux : public AuxKernel
{
public:

  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  AreaAux(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeValue();
  Function & _area;
};

#endif //AREAAUX_H
