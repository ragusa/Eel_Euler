#ifndef COMPUTEVISCCOEFF_H
#define COMPUTEVISCCOEFF_H

#include "Material.h"
#include "MaterialProperty.h"
#include "EquationOfState.h"

//Forward Declarations
class ComputeViscCoeff;

template<>
InputParameters validParams<ComputeViscCoeff>();

class ComputeViscCoeff : public Material
{
public:
  ComputeViscCoeff(const std::string & name, InputParameters parameters);

protected:
  virtual void computeQpProperties();

private:
    // Viscosity types
    enum ViscosityType
    {
        LAPIDUS = 0,
        FIRST_ORDER = 1,
        FIRST_ORDER_MACH = 2,
        ENTROPY = 3,
        PRESSURE_BASED = 4
    };
    std::string _visc_name;
    MooseEnum _visc_type;
    
    // Pressure-based viscosity:
    enum PBType
    {
        JST = 0,
        HMP = 1,
        ST = 2
    };
    VariableValue & _PBVisc;
    std::string _norm_pbs_name;
    MooseEnum _norm_pbs_type;
    
    // Boolean for jump
    bool _isJumpOn;
    bool _isShock;
    
    // Coupled aux variables: velocity
    VariableValue & _vel_x;
    VariableValue & _vel_y;
    VariableValue & _vel_z;
    VariableGradient & _grad_vel_x;
    VariableGradient & _grad_vel_y;
    VariableGradient & _grad_vel_z;
    
    // Coupled aux variables: pressure
    VariableValue & _pressure;
    VariableValue & _pressure_old;
    VariableValue & _pressure_older;
    VariableGradient & _grad_press;
    
    // Coupled aux variable: density
    VariableValue & _rho;
    VariableValue & _rho_old;
    VariableValue & _rho_older;
    VariableGradient & _grad_rho;
    
    // Coupled aux variable: norm of velocity
    VariableValue & _norm_vel;
    VariableGradient & _grad_norm_vel;
    
    // Jump of pressure gradient:
    VariableValue & _jump_grad_press;
    VariableValue & _jump_grad_dens;
    
    // Jump cross section:
    VariableValue & _area;
    VariableGradient & _grad_area;
    
    // Material properties
    MaterialProperty<Real> & _mu;
    MaterialProperty<Real> & _mu_max;
    MaterialProperty<Real> & _kappa;
    MaterialProperty<Real> & _kappa_max;
    MaterialProperty<RealVectorValue> & _l;
    
//    MaterialProperty<Real> & _residual;
    // Wall heat transfer
    std::string _Hw_fn_name;
    std::string _Tw_fn_name;
    const Real & _Hw;
    const Real & _Tw;
    const Real & _aw;
    
    // Multiplicative coefficient for viscosity:
    double _Ce;
    double _Cjump;
    double _Cmax;
    
    // UserObject: equation of state
    const EquationOfState & _eos;
    
    // Name of the posprocessors for pressure, velocity and void fraction:
    std::string _rhov2_pps_name;
    std::string _rhoc2_pps_name;
    std::string _press_pps_name;
};

#endif //ComputeViscCoeff_H
