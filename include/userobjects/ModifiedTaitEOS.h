#ifndef MODIFIEDTAITEOS_H
#define MODIFIEDTAITEOS_H

#include "GeneralUserObject.h"

// Forward Declarations
class ModifiedTaitEOS;

template<>
InputParameters validParams<ModifiedTaitEOS>();

class ModifiedTaitEOS : public GeneralUserObject
{
public:
  // Constructor
  ModifiedTaitEOS(const std::string & name, InputParameters parameters);

  // Destructor  
  virtual ~ModifiedTaitEOS(); 

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

    // The interface for derived ModifiedTaitEOS objects to implement...
    virtual Real pressure(Real rho=0., Real mom_norm=0., Real rhoE=0.) const;
    
    // The interface for derived ModifiedTaitEOS objects to implement...
    virtual Real rho_from_p_T(Real pressure=0., Real temperature=0.) const;
    
    // The interface for derived ModifiedTaitEOS objects to implement...
    virtual Real e_from_p_rho(Real pressure=0., Real rho=0.) const;
    
    // The interface for derived ModifiedTaitEOS objects to implement...
    virtual Real temperature_from_p_rho(Real pressure=0., Real rho=0.) const;
    
    // The interface for derived ModifiedTaitEOS objects to implement...
    virtual Real c2_from_p_rho(Real rho=0., Real pressure=0.) const;
    
    // Derivatives of the pressure:
    virtual Real dAp_drhoA(Real rhoA=0., Real rhouA_norm=0., Real rhoEA=0.) const;
    
    virtual Real dAp_drhouA(Real rhoA=0., Real rhouA_component=0., Real rhoEA=0.) const;
    
    virtual Real dAp_drhoEA(Real rhoA=0., Real rhouA_norm=0., Real rhoEA=0.) const;
    
  Real gamma() const { return _gamma; }
    
  Real P0() const { return _P0; }
    
  Real P1() const { return _P1; }
    
  Real rho0() const { return _rho0; }
    
  Real e0() const { return _e0; }
    
  Real T0() const { return _T0; }
    
  Real Cv() const { return _Cv; }

protected:
    Real _gamma;
    Real _P0;
    Real _P1;
    Real _rho0;
    Real _e0;
    Real _T0;
    Real _Cv;
};


#endif // MODIFIEDTAITEOS_H
