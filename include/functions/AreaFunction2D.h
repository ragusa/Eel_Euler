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

#ifndef AreaFunction2D2D_H
#define AreaFunction2D2D_H

#include "Function.h"

class AreaFunction2D;

template<>
InputParameters validParams<AreaFunction2D>();

class AreaFunction2D : public Function
{
public:
  AreaFunction2D(const std::string & name, InputParameters parameters);

  virtual Real value(Real t, const Point & p);

  virtual RealVectorValue gradient(Real t, const Point & p);

protected:
  Real _Ao;
  Real _Bo;
};

#endif //AreaFunction2D2D_H
