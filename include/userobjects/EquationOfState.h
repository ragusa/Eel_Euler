#ifndef EQUATIONOFSTATE_H
#define EQUATIONOFSTATE_H

#include "GeneralUserObject.h"

// Forward Declarations
class EquationOfState;

template<>
InputParameters validParams<EquationOfState>();

class EquationOfState : public GeneralUserObject
{
public:
  // Constructor
  EquationOfState(const std::string & name, InputParameters parameters);

  // Destructor  
  virtual ~EquationOfState(); 

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

    // The interface for derived EquationOfState objects to implement...
    virtual Real pressure(Real rho=0., Real mom_norm=0., Real rhoE=0.) const;
    
    // The interface for derived EquationOfState objects to implement...
    virtual Real rho_from_p_T(Real pressure=0., Real temperature=0.) const;
    
    // The interface for derived EquationOfState objects to implement...
    virtual Real e_from_p_rho(Real pressure=0., Real rho=0.) const;
    
    // The interface for derived EquationOfState objects to implement...
    virtual Real temperature_from_p_rho(Real pressure=0., Real rho=0.) const;
    
    // The interface for derived EquationOfState objects to implement...
    virtual Real c2_from_p_rho(Real rho=0., Real pressure=0.) const;
    
    // Derivatives of the pressure:
    virtual Real dAp_drhoA(Real rhoA=0., Real rhouA_norm=0., Real rhoEA=0.) const;
    
    virtual Real dAp_drhouA(Real rhoA=0., Real rhouA_component=0., Real rhoEA=0.) const;
    
    virtual Real dAp_drhoEA(Real rhoA=0., Real rhouA_norm=0., Real rhoEA=0.) const;

    Real gamma() const;
    
    Real Pinf() const;
    
    Real qcoeff() const;
    
    Real qcoeff_prime() const;
    
    Real Cv() const;
    
protected:
    Real _gamma;
    Real _Pinf;
    Real _qcoeff;
    Real _qcoeff_prime;
    Real _Cv;
    // Prints an error message for non-implemented functions
    void error_not_implemented(std::string method_name) const;
};


#endif // EQUATIONOFSTATE_H
