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

#ifndef EELARTIFICIALVISC_H
#define EELARTIFICIALVISC_H

#include "Kernel.h"

// Forward Declarations
class EelArtificialVisc;

template<>
InputParameters validParams<EelArtificialVisc>();

class EelArtificialVisc : public Kernel
{
public:

  EelArtificialVisc(const std::string & name,
             InputParameters parameters);

protected:

  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian(unsigned int _jvar);
    
private:
    // Equations types
    enum EquationType
    {
        CONTINUITY = 0,
        XMOMENTUM = 1,
        YMOMENTUM = 2,
        ZMOMENTUM = 3,
        ENERGY = 4
    };
    // Diffusion types
    enum DiffusionType
    {
        ENTROPY = 0,
        PARABLOIC = 1,
        NONE = 2
    };
    // Diffusion name
    std::string _equ_name;
    std::string _diff_name;
    // Diffusion type
    MooseEnum _equ_type;
    MooseEnum _diff_type;
    // Coupled aux variables:
    VariableValue & _rho;
    VariableGradient & _grad_rho;
    VariableValue & _vel_x;
    VariableValue & _vel_y;
    VariableValue & _vel_z;
    VariableGradient & _grad_vel_x;
    VariableGradient & _grad_vel_y;
    VariableGradient & _grad_vel_z;
    VariableGradient & _grad_rhoe;
    VariableValue & _area;
    VariableValue & _norm_vel;
    VariableGradient & _grad_norm_vel;
    // Material property: viscosity coefficient.
    MaterialProperty<Real> & _mu;
    MaterialProperty<Real> & _kappa;
};

#endif // EELARTIFICIALVISC_H
