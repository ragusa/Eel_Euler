#ifndef EEL_EULERAPP_H
#define EEL_EULERAPP_H

#include "MooseApp.h"

class EelEulerApp;

template<>
InputParameters validParams<EelEulerApp>();

class EelEulerApp : public MooseApp
{
public:
  EelEulerApp(const std::string & name, InputParameters parameters);
  virtual ~EelEulerApp();

  static void registerApps();
  static void registerObjects(Factory & factory);
  static void associateSyntax(Syntax & syntax, ActionFactory & action_factory);
};

#endif /* EEL_EULERAPP_H */
