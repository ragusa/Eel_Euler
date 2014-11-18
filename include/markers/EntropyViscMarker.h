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

#ifndef ENTROPYVISCMARKER_H
#define ENTROPYVISCMARKER_H

#include "IndicatorMarker.h"

class EntropyViscMarker;

template<>
InputParameters validParams<EntropyViscMarker>();

class EntropyViscMarker : public IndicatorMarker
{
public:
  EntropyViscMarker(const std::string & name, InputParameters parameters);
  virtual ~EntropyViscMarker(){};

  virtual void markerSetup();

protected:
  virtual MarkerValue computeElementMarker();

  Real _coarsen;
  Real _refine;

  Real _max;
  Real _min;
  Real _delta;
  Real _refine_cutoff;
  Real _coarsen_cutoff;
};

#endif /* ENTROPYVISCMARKER_H */
