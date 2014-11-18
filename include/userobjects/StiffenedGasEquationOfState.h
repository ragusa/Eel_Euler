#ifndef STIFENEDGASEQUATIONOFSTATE_H
#define STIFENEDGASEQUATIONOFSTATE_H

#include "EquationOfState.h"

// Forward Declarations
class StiffenedGasEquationOfState;

template<>
InputParameters validParams<StiffenedGasEquationOfState>();

class StiffenedGasEquationOfState : public EquationOfState
{
public:
  // Constructor
  StiffenedGasEquationOfState(const std::string & name, InputParameters parameters);

  // Destructor  
  virtual ~StiffenedGasEquationOfState(); 

  /**
   * Called when this object needs to compute something.
   */
  virtual void execute() {}

  /**
   * Called before execute() is ever called so that data can be cleared.
   */
  virtual void initialize(){}

  
  virtual void destroy();

  virtual void finalize() {};

    // The interface for derived StiffenedGasEquationOfState objects to implement...
    virtual Real pressure(Real rho=0., Real mom_norm=0., Real rhoE=0.) const;
    
    // The interface for derived StiffenedGasEquationOfState objects to implement...
    virtual Real rho_from_p_T(Real pressure=0., Real temperature=0.) const;
    
    // The interface for derived StiffenedGasEquationOfState objects to implement...
    virtual Real e_from_p_rho(Real pressure=0., Real rho=0.) const;
    
    // The interface for derived StiffenedGasEquationOfState objects to implement...
    virtual Real temperature_from_p_rho(Real pressure=0., Real rho=0.) const;
    
    // The interface for derived StiffenedGasEquationOfState objects to implement...
    virtual Real c2_from_p_rho(Real rho=0., Real pressure=0.) const;
    
    // Derivatives of the pressure:
    virtual Real dAp_drhoA(Real rhoA=0., Real rhouA_norm=0., Real rhoEA=0.) const;
    
    virtual Real dAp_drhouA(Real rhoA=0., Real rhouA_component=0., Real rhoEA=0.) const;
    
    virtual Real dAp_drhoEA(Real rhoA=0., Real rhouA_norm=0., Real rhoEA=0.) const;
    
//  Real gamma() const { return _gamma; }
//    
//  Real Pinf() const { return _Pinf; }
//    
//  Real qcoeff() const { return _qcoeff; }
//    
//  Real qcoeff_prime() const { return _qcoeff_prime; }
//    
//  Real Cv() const { return _Cv; }

protected:
    Real _gamma;
    Real _Pinf;
    Real _qcoeff;
    Real _qcoeff_prime;
    Real _Cv;
};


#endif // STIFENEDGASEQUATIONOFSTATE_H
