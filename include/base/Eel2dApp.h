#ifndef EEL_2DAPP_H
#define EEL_2DAPP_H

#include "MooseApp.h"

class Eel2dApp;

template<>
InputParameters validParams<Eel2dApp>();

class Eel2dApp : public MooseApp
{
public:
  Eel2dApp(const std::string & name, InputParameters parameters);
  virtual ~Eel2dApp();

  static void registerApps();
  static void registerObjects(Factory & factory);
  static void associateSyntax(Syntax & syntax, ActionFactory & action_factory);
};

#endif /* EEL_2DAPP_H */
