/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "SmoothFunction.h"

/* This function is called to compute the jump of the gradient of a given quantity when using CONTINUOUS finite element. This function acts on the sides of the cell.*/
template<>
InputParameters validParams<SmoothFunction>()
{
  InputParameters params = validParams<InternalSideUserObject>();
    params.addRequiredCoupledVar("variable", "the variable name this userobject is acting on.");
    params.addRequiredParam<std::string>("var_name", "the name of the variable that will store the smoothed variable.");
  return params;
}

SmoothFunction::SmoothFunction(const std::string & name, InputParameters parameters) :
    InternalSideUserObject(name, parameters),
    _aux(_fe_problem.getAuxiliarySystem()),
    _u(coupledValue("variable")),
    _u_neighbor(coupledNeighborValue("variable")),
    _var_name(getParam<std::string>("var_name")),
    _value(0.)
{
}

SmoothFunction::~SmoothFunction()
{
}

void
SmoothFunction::initialize()
{
    //_value = 0.;
    NumericVector<Number> & sln = _aux.solution();
    _aux.system().zero_variable(sln, _aux.getVariable(_tid, _var_name).number());
}

void
SmoothFunction::execute()
{
    // Get the dimension of the geometry:
    int _dim = _fe_problem.mesh().dimension();
    
    // Compute the total perimeter/area of the elements:
    Real _perim_elem = 0.;
    Real _perim_nghb_elem = 0.;
    Real _size_side = 1.;
    if (_mesh.dimension() != 1) {
        for (unsigned int _jvar=0; _jvar<_current_elem->n_sides(); _jvar++) {
            _perim_elem += _current_elem->side(_jvar)->volume();
        }
        for (unsigned int _jvar=0; _jvar<_neighbor_elem->n_sides(); _jvar++) {
            _perim_nghb_elem += _neighbor_elem->side(_jvar)->volume();
        }
        _size_side = _current_side_volume;
    }
    else {
        _perim_elem = 1.;
        _perim_nghb_elem = 1.;
    }
    
    // Compute the weights function:
    Real _weight_elem = _size_side / _perim_elem;
    Real _weight_nghb_elem = _size_side / _perim_nghb_elem;
    
    // Initialyze some parameters:
    dof_id_type dof_nb_aux = _current_elem->n_dofs(_aux.number(), _fe_problem.getVariable(_tid, _var_name).number());
    dof_id_type dof_nb = 0.;
    dof_id_type dof_nb_neighbor = 0.;
    
    // Do the job only if the coupled variable 'variable' is defined on the nodes:
    if (dof_nb_aux != 0) {
        _value = 0.;
        NumericVector<Number> & sln = _aux.solution();
        
        // Determine the maximum value for smoothing:
        for (unsigned int qp = 0; qp < _q_point.size(); ++qp) {
            Real _value_temp = 0.5*(_u[qp]+_u_neighbor[qp]);
            _value = std::max(_value_temp, _value);
        }
        
        // Determine the degree of freedom:
        dof_nb = _current_elem->dof_number(_aux.number(), _fe_problem.getVariable(_tid, _var_name).number(), 0);
        dof_nb_neighbor = _neighbor_elem->dof_number(_aux.number(), _fe_problem.getVariable(_tid, _var_name).number(), 0);
        
        // Set the value:
        sln.add(dof_nb, _value*_weight_elem);
        sln.add(dof_nb_neighbor, _value*_weight_nghb_elem);
    }
}

void
SmoothFunction::destroy()
{
}

void
SmoothFunction::finalize()
{
    _aux.solution().close();
    
}

void
SmoothFunction::threadJoin(const UserObject & uo)
{
}
