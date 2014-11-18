#ifndef EELHRHOUBC_H
#define EELHRHOUBC_H

#include "IntegratedBC.h"
#include "EquationOfState.h"

// Forward Declarations
class EelHRhoUBC;
class EquationOfState;

template<>
InputParameters validParams<EelHRhoUBC>();


/**
 * The boundary condition with specified stagnation pressure and temperature
 * A void fraction boundary has also to be included to close the boundary condition for 7eqn system
 */
class EelHRhoUBC : public IntegratedBC
{

public:
  EelHRhoUBC(const std::string & name, InputParameters parameters);

  virtual ~EelHRhoUBC(){}

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned jvar);

    enum EFlowEquationType
    {
    CONTINUITY = 1,
    XMOMENTUM = 2,
    ENERGY = 3
    };

    // Eqn. name to be read from input file
    std::string _eqn_name;
    // which equation (mass/momentum/energy) this BC is acting on
    MooseEnum _eqn_type;

    // Coupled variables
    VariableValue & _rhoA;
    VariableValue & _rhouA_x;
    VariableValue & _rhoEA;
    
    // Coupled aux variables:
    VariableValue & _area;

    // Specified stagnation variables:
    Real _rhou;
    Real _H;

    // Equation of state:
    const EquationOfState & _eos;
    
    // Parameters for jacobian matrix:
//    unsigned int _rhoA_nb;
//    unsigned int _rhouA_x_nb;
//    unsigned int _rhouA_y_nb;
//    unsigned int _rhoEA_nb;
};

#endif // EELHRHOUBC_H

