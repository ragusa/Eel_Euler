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

#ifndef AREAFUNCTION_H
#define AREAFUNCTION_H

#include "Function.h"

class AreaFunction;

template<>
InputParameters validParams<AreaFunction>();

class AreaFunction : public Function
{
public:
  AreaFunction(const std::string & name, InputParameters parameters);

  virtual Real value(Real t, const Point & p);

  virtual RealVectorValue gradient(Real t, const Point & p);

protected:
  Real _left;
  Real _length;
  Real _Ao;
  Real _Bo;
};

#endif //AREAFUNCTION_H
