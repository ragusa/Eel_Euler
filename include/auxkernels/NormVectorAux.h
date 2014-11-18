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

#ifndef NORMVECTORAUX_H
#define NORMVECTORAUX_H

#include "AuxKernel.h"

//Forward Declarations
class NormVectorAux;

template<>
InputParameters validParams<NormVectorAux>();

class NormVectorAux : public AuxKernel
{
public:

  NormVectorAux(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeValue();

    VariableValue & _x_comp;
    VariableValue & _y_comp;
    VariableValue & _z_comp;
};

#endif //NORMVECTORAUX_H
