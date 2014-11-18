#ifndef INVISCIDTIMESTEPLIMIT_H
#define INVISCIDTIMESTEPLIMIT_H

#include "ElementPostprocessor.h"

class InviscidTimeStepLimit;

template<>
InputParameters validParams<InviscidTimeStepLimit>();

/**
 * The inviscid time step stability limit:
 *
 * h_e \over {|\vec u| + c}
 */
class InviscidTimeStepLimit : public ElementPostprocessor
{
public:
  InviscidTimeStepLimit(const std::string & name, InputParameters parameters);
  virtual ~InviscidTimeStepLimit();

  virtual void initialize();
  virtual void execute();
  virtual Real getValue();
  virtual void threadJoin(const UserObject & uo);

protected:
  unsigned int _dim;
  /// The value of dt (NOTE: _dt member variable is already defined)
  Real _value;
  /// Velocity magnitude.  Hint: Use VectorMagnitudeAux in Moose for this
  VariableValue & _vel_mag;
  /// Sound Speed
  VariableValue & _c;
  Real _beta;
};


#endif /* INVISCIDTIMESTEPLIMIT_H */
