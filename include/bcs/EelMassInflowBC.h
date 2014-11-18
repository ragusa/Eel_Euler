#ifndef EELMASSINFLOWBC_H
#define EELMASSINFLOWBC_H

#include "IntegratedBC.h"
#include "EquationOfState.h"

// Forward Declarations
class EelMassInflowBC;
class EquationOfState;

template<>
InputParameters validParams<EelMassInflowBC>();

class EelMassInflowBC : public IntegratedBC
{

public:
  EelMassInflowBC(const std::string & name, InputParameters parameters);

  virtual ~EelMassInflowBC(){}

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
    VariableValue & _pressure;
    VariableValue & _area;
    
    // Specified velocity values:
    Real _rhou0_bc;
    Real _T0_bc;
    Real _alpha0_bc;
    
    // Equation of state
    const EquationOfState & _eos;
};

#endif // EELMASSINFLOWBC_H

