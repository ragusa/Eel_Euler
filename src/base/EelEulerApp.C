#include "EelEulerApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"

template<>
InputParameters validParams<EelEulerApp>()
{
  InputParameters params = validParams<MooseApp>();
  return params;
}

EelEulerApp::EelEulerApp(const std::string & name, InputParameters parameters) :
    MooseApp(name, parameters)
{
  srand(processor_id());

  Moose::registerObjects(_factory);
  ModulesApp::registerObjects(_factory);
  EelEulerApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  ModulesApp::associateSyntax(_syntax, _action_factory);
  EelEulerApp::associateSyntax(_syntax, _action_factory);
}

EelEulerApp::~EelEulerApp()
{
}

void
EelEulerApp::registerApps()
{
  registerApp(EelEulerApp);
}

void
EelEulerApp::registerObjects(Factory & factory)
{
}

void
EelEulerApp::associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
}
