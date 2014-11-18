#ifndef EELWALLBC_H
#define EELWALLBC_H

#include "IntegratedBC.h"
#include "EquationOfState.h"

class EelWallBC;

template<>
InputParameters validParams<EelWallBC>();

class EelWallBC : public IntegratedBC
{

public:
  EelWallBC(const std::string & name, InputParameters parameters);

  virtual ~EelWallBC(){}

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned jvar);

  enum EFlowEquationType
  {
    CONTINUITY = 0,
    XMOMENTUM = 1,
    YMOMENTUM = 2,
    ENERGY = 3
  };

    // Eqn. name to be read from input file
    std::string _eqn_name;
    
    // which equation (mass/momentum/energy) this BC is acting on
    MooseEnum _eqn_type;
    
    // Coupled aux variables
    VariableValue & _rhoA;
    VariableValue & _rhouA_x;
    VariableValue & _rhouA_y;
    VariableValue & _rhoEA;
    VariableValue & _area;
    
    // Equation of state for jacobian matrix:
    const EquationOfState & _eos;
    
    // Parameters for jacobian matrix:
    unsigned int _rhoA_nb;
    unsigned int _rhouA_x_nb;
    unsigned int _rhouA_y_nb;
    unsigned int _rhoEA_nb;
};

#endif // EELWALLBC_H

