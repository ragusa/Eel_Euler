#ifndef EELFLUXBC_H
#define EELFLUXBC_H

#include "IntegratedBC.h"
#include "EquationOfState.h"

// Forward Declarations
class EelFluxBC;
class EquationOfState;

template<>
InputParameters validParams<EelFluxBC>();


/**
 * The boundary condition with specified stagnation pressure and temperature
 * A void fraction boundary has also to be included to close the boundary condition for 7eqn system
 */
class EelFluxBC : public IntegratedBC
{

public:
  EelFluxBC(const std::string & name, InputParameters parameters);

  virtual ~EelFluxBC(){}

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned jvar);

    enum EFlowEquationType
    {
    CONTINUITY = 0,
    XMOMENTUM = 1,
    YMOMENTUM = 2,
    ZMOMENTUM = 3,
    ENERGY = 4
    };

    // Eqn. name to be read from input file
    std::string _eqn_name;
    // which equation (mass/momentum/energy) this BC is acting on
    MooseEnum _eqn_type;

    // Coupled variables
    VariableValue & _rhoA;
    VariableValue & _rhouA_x;
    VariableValue & _rhouA_y;
    VariableValue & _rhoEA;
    
    // Gradient of coupled varibles:
    VariableGradient & _grad_rhoA;
    VariableGradient & _grad_rhouA;
    VariableGradient & _grad_rhovA;
    VariableGradient & _grad_rhoEA;
    
    // Coupled aux variables:
    VariableValue & _area;
    
    // Material property
    MaterialProperty<Real> & _kappa;
    MaterialProperty<Real> & _mu;

    // Equation of state:
    const EquationOfState & _eos;
    
    // Parameters for jacobian matrix:
    unsigned int _rhoA_nb;
    unsigned int _rhouA_x_nb;
    unsigned int _rhouA_y_nb;
    unsigned int _rhoEA_nb;
};

#endif // EelFluxBC_H

