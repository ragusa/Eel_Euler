#ifndef NODALMASSCONSERVATIONPPS_H
#define NODALMASSCONSERVATIONPPS_H

#include "NodalPostprocessor.h"

class NodalMassConservationPPS;

template<>
InputParameters validParams<NodalMassConservationPPS>();

/**
 * Computes \int{\Gamma} \rho \vec u \hat n \d\Gamma at nodes
 *
 * NOTES:
 * - nodal normals has to be turned on
 * - will work only with nodal shape functions (i.e. like LAGRANGE)
 */
class NodalMassConservationPPS : public NodalPostprocessor
{
public:
  NodalMassConservationPPS(const std::string & name, InputParameters parameters);
  virtual ~NodalMassConservationPPS();

  virtual void initialize();
  virtual void execute();
  virtual PostprocessorValue getValue();
  virtual void threadJoin(const UserObject & uo);

protected:
  Real _value;
  /// Components of nodal normals
  VariableValue & _nx;
  VariableValue & _ny;
  VariableValue & _nz;
  /// Components of momentum vector
  VariableValue & _rho_u;
  VariableValue & _rho_v;
  VariableValue & _rho_w;
};

#endif /* NODALMASSCONSERVATIONPPS_H */
