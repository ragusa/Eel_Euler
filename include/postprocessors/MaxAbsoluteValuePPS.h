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

#ifndef MAXABSOLUTEVALUEPPS_H
#define MAXABSOLUTEVALUEPPS_H

#include "ElementPostprocessor.h"

//Forward Declarations
class MaxAbsoluteValuePPS;

template<>
InputParameters validParams<MaxAbsoluteValuePPS>();

class MaxAbsoluteValuePPS : public ElementPostprocessor
{
public:
  MaxAbsoluteValuePPS(const std::string & name, InputParameters parameters);

  virtual void initialize();
  virtual void execute();
  virtual void finalize();
  virtual Real getValue();
  virtual void threadJoin(const UserObject & y);

protected:
    // Variable this pps is acting on:
    VariableValue & _u;
    
    // Value storing the max:
    Real _value;
};

#endif
