#ifndef EELSTAGNATIONPANDTBC_H
#define EELSTAGNATIONPANDTBC_H

#include "IntegratedBC.h"
#include "EquationOfState.h"

// Forward Declarations
class EelStagnationPandTBC;
class EquationOfState;

template<>
InputParameters validParams<EelStagnationPandTBC>();


/**
 * The boundary condition with specified stagnation pressure and temperature
 * A void fraction boundary has also to be included to close the boundary condition for 7eqn system
 */
class EelStagnationPandTBC : public IntegratedBC
{

public:
  EelStagnationPandTBC(const std::string & name, InputParameters parameters);

  virtual ~EelStagnationPandTBC(){}

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned jvar);

    enum EFlowEquationType
    {
    CONTINUITY = 1,
    XMOMENTUM = 2,
    YMOMENTUM = 3,
    ZMOMENTUM = 4,
    ENERGY = 5
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
    
    // Coupled aux variables:
    VariableValue & _area;

    // Specified stagnation variables:
    Real _p0_bc;
    Real _T0_bc;
    Real _gamma0_bc;

    // Calculated rho_0, K, etc. on the boundary:
    Real _rho0_bc;
    Real _H0_bc;
    Real _K;
    Real _H_bar;

    // Equation of state:
    const EquationOfState & _eos;
    
    // Parameters for jacobian matrix:
    unsigned int _rhoA_nb;
    unsigned int _rhouA_x_nb;
    unsigned int _rhouA_y_nb;
    unsigned int _rhoEA_nb;
};

#endif // EELSTAGNATIONPANDTBC_H

