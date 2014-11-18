#ifndef EELDBC_H
#define EELDBC_H

#include "IntegratedBC.h"

class EelDBC;

template<>
InputParameters validParams<EelDBC>();

class EelDBC : public IntegratedBC
{

public:
  EelDBC(const std::string & name, InputParameters parameters);

  virtual ~EelDBC(){}

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned jvar);

};

#endif // EELDBC_H

