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

#include "JumpGradientInterface.h"

/* This function is called to compute the jump of the gradient of a given quantity when using CONTINUOUS finite element. This function acts on the sides of the cell.*/
template<>
InputParameters validParams<JumpGradientInterface>()
{
  InputParameters params = validParams<InternalSideUserObject>();
    params.addRequiredCoupledVar("variable", "the variable name this userobject is acting on.");
    params.addRequiredParam<std::string>("jump_name", "the name of the variable that will store the jump");
  return params;
}

JumpGradientInterface::JumpGradientInterface(const std::string & name, InputParameters parameters) :
    InternalSideUserObject(name, parameters),
    _aux(_fe_problem.getAuxiliarySystem()),
    _grad_u(coupledGradient("variable")),
    _grad_u_neighbor(coupledNeighborGradient("variable")),
    _jump_name(getParam<std::string>("jump_name")),
    _value(0.)
{
}

JumpGradientInterface::~JumpGradientInterface()
{
}

void
JumpGradientInterface::initialize()
{
    //_value = 0.;
    NumericVector<Number> & sln = _aux.solution();
    _aux.system().zero_variable(sln, _aux.getVariable(_tid, _jump_name).number());
    
}

void
JumpGradientInterface::execute()
{
    // Get the dimension of the geometry:
    int _dim = _fe_problem.mesh().dimension();
    
    // Compute the total perimeter/area of the elements:
    Real _perim_elem = 0.;
    Real _perim_nghb_elem = 0.;
    Real _size_side = 0.5;
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
    dof_id_type dof_nb_aux = _current_elem->n_dofs(_aux.number(), _fe_problem.getVariable(_tid, _jump_name).number());
    dof_id_type dof_nb = 0.;
    dof_id_type dof_nb_neighbor = 0.;
    
    // Do the job only if the coupled variable 'variable' is defined on the nodes:
    if (dof_nb_aux != 0)
    {
        _value = 0.;
        NumericVector<Number> & sln = _aux.solution();
        
        // Compute the jump of the given variable:(grad(f_i) - grad(f_ip1))*_normals
        for (unsigned int qp = 0; qp < _q_point.size(); ++qp)
            _value = std::max(std::fabs(_grad_u[qp]*_normals[qp] - _grad_u_neighbor[qp]*_normals[qp]), _value);
        
        dof_nb = _current_elem->dof_number(_aux.number(), _fe_problem.getVariable(_tid, _jump_name).number(), 0);
        dof_nb_neighbor = _neighbor_elem->dof_number(_aux.number(), _fe_problem.getVariable(_tid, _jump_name).number(), 0);
        
        // Set the value:
        sln.add(dof_nb, _value*_weight_elem);
        sln.add(dof_nb_neighbor, _value*_weight_nghb_elem);
    }
}

void
JumpGradientInterface::destroy()
{
}

void
JumpGradientInterface::finalize()
{
    _aux.solution().close();
    
}

void
JumpGradientInterface::threadJoin(const UserObject & uo)
{
}
