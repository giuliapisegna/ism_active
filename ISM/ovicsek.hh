/*
 * ovicsek.hh -- Simulation class for Vicsek model (as originally formulated
 *
 */

#ifndef OVICSEK_HH
#define OVICSEK_HH

#include "glsim/random.hh"
#include "glsim/stochastic.hh"
#include "glsim/simulation.hh"
#include "glsim/olconfiguration.hh"
#include "glsim/observable.hh"
#include "glsim/md.hh"
#include "glsim/avevar.hh"

#include "ovicsek_env.hh"

/*****************************************************************************
 *
 * OVicsekSimulation
 * 
 */

class OVicsekSimulation : public glsim::Simulation {
public:
  OVicsekSimulation(OVicsekEnvironment& e,glsim::OLconfiguration &c);
  ~OVicsekSimulation();
  const char* name() const {return "Vicsek's model (original)";}

  void step();
  void log(),log_start_sim();

protected:
  void update_observables();

  OVicsekEnvironment&         env;
  glsim::OLconfiguration& conf;

private:
  void vnoise(double*);
  void update_velocities();

  glsim::NeighbourList_subcells *NN;
  double  (*confb)[3];
  double  rcsq,v0sq;

  glsim::Uniform_real *ranz,*ranphi;
} ;


/******************************************************************************
 *
 * Observable
 *
 */

class OVicsekObservable_parameters : public glsim::Parameters {
public:
  OVicsekObservable_parameters(const char* scope);
} ;

class OVicsekObservable :  public glsim::SBObservable {
public:
  OVicsekObservable(OVicsekEnvironment&,glsim::OLconfiguration&);

  void interval_and_file();
  void write_header();
  void observe();

private:
  OVicsekEnvironment  &env;
  glsim::OLconfiguration &conf;
  OVicsekObservable_parameters par;

  void update();
} ;

inline OVicsekObservable::OVicsekObservable(OVicsekEnvironment& e,glsim::OLconfiguration &c) :
  SBObservable(e),
  env(e),
  conf(c),
  par(e.scope())
{}


#endif /* OVICSEK_HH */
