#ifndef MOMENTUMFREESLIPBC_H
#define MOMENTUMFREESLIPBC_H

#include "NodalNormalBC.h"

class MomentumFreeSlipBC;

template<>
InputParameters validParams<MomentumFreeSlipBC>();

/**
 *
 */
class MomentumFreeSlipBC : public NodalNormalBC
{
public:
  MomentumFreeSlipBC(const std::string & name, InputParameters parameters);
  virtual ~MomentumFreeSlipBC();

  virtual void computeResidual(NumericVector<Number> & residual);
  virtual bool shouldApply();

protected:
  virtual Real computeQpResidual();

  const unsigned int _mesh_dimension;
  
  VariableValue & _rho_u;
  VariableValue & _rho_v;
  VariableValue & _rho_w;
};


#endif /* MOMENTUMFREESLIPBC_H_ */
