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

#ifndef CONSERVATIVEVARIABLES2DIC_H
#define CONSERVATIVEVARIABLES2DIC_H

// MOOSE Includes
#include "InitialCondition.h"
#include "EquationOfState.h"
#include "AreaFunction.h"

// Forward Declarations
class ConservativeVariables2DIC;

template<>
InputParameters validParams<ConservativeVariables2DIC>();

/**
 * ConservativeVariables2DIC just returns a constant value.
 */
class ConservativeVariables2DIC : public InitialCondition
{
public:

  /**
   * Constructor: Same as the rest of the MOOSE Objects
   */
  ConservativeVariables2DIC(const std::string & name,
            InputParameters parameters);

  /**
   * The value of the variable at a point.
   *
   * This must be overriden by derived classes.
   */
  virtual Real value(const Point & p);

private:
    // Function are
    Function & _area;
    
    // Initial conditions for left and right values:
    Real _p_left;
    Real _p_right;
    Real _v_left;
    Real _v_right;
    Real _t_left;
    Real _t_right;
    
    // Position of the point source:
    Real _x_pt_source;
    Real _y_pt_source;
    Real _length;
    Real _smoothing;
    
    // Name of the variable:
    std::string _name_var;
    
    // Equation of state:
    const EquationOfState & _eos;

};

#endif //CONSERVATIVEVARIABLES2DIC_H
