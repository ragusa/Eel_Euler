#ifndef EELINFINITEBC_H
#define EELINFINITEBC_H

#include "IntegratedBC.h"

class EelInfiniteBC;

template<>
InputParameters validParams<EelInfiniteBC>();

class EelInfiniteBC : public IntegratedBC
{

public:
  EelInfiniteBC(const std::string & name, InputParameters parameters);

  virtual ~EelInfiniteBC(){}

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

    /// Eqn. name to be read from input file
    std::string _eqn_name;
    /// which equation (mass/momentum/energy) this BC is acting on
    MooseEnum _eqn_type;
    /// Coupled aux variables
    VariableValue & _pressure;
    VariableValue & _vel_x;
    VariableValue & _vel_y;
    VariableValue & _vel_z;
    VariableValue & _area;
};

#endif // EELINFINITEBC_H

