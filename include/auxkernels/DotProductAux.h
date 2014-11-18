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

#ifndef DOTPRODUCTAUX_H
#define DOTPRODUCTAUX_H

#include "AuxKernel.h"

//Forward Declarations
class DotProductAux;

template<>
InputParameters validParams<DotProductAux>();

class DotProductAux : public AuxKernel
{
public:

  DotProductAux(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeValue();
    // Vector 1:
    VariableValue & _x_comp;
    VariableValue & _y_comp;
    VariableValue & _z_comp;
    // Vector 2:
    VariableGradient & _grad_vector;
};

#endif //DOTPRODUCTAUX_H
