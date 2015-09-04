/*
 * vici.hh -- Vicsek model integrator
 *
 * Integrates the equation of motion for the Vicsek model, understood
 * as the overdamped inertial spin model (with Brownian noise),
 * without specifiying the velocity-velocity interaction.
 *
 * Tomás S. Grigera <tgrigera@iflysib.unlp.edu.ar>
 *
 * La Plata, September 2015
 *
 */

#ifndef VICI_HH
#define VICI_HH

#include "vici.hh"

#include "glsim/random.hh"
#include "glsim/stochastic.hh"
#include "glsim/simulation.hh"
#include "glsim/olconfiguration.hh"
#include "glsim/observable.hh"
#include "glsim/md.hh"
#include "glsim/avevar.hh"

#include "social.hh"

/*****************************************************************************
 *
 * Environment
 * 
 */

class VicsekIntegratorParameters : public glsim::Parameters {
public:
  VicsekIntegratorParameters(const char* scope);
} ;

class VicsekEnvironment : public glsim::SimEnvironment {
public:
  VicsekEnvironment(const char* scope=glsim::Parameters::default_scope);

  long    VSsteps;
  double  time_step;
  bool    fixed_graph;
  double  temperature;
  // double  eta;
  bool    rescale_v0;

  // System info
  double  v0;
  double  social_mass[2];

  // Observables
  int     total_number;
  double  total_social_mass;
  double  social_total_energy,social_potential_energy,social_kinetic_energy;
  double  polarization,v0sqave;
  double  Vcm[3],total_spin[3],total_spinsq;
  glsim::AveVar<true> spinsqavar;

protected:
  void    init_local() {SimEnvironment::init_local(); common_init();}
  void    warm_init_local() {SimEnvironment::warm_init_local(); common_init();}
  void    update_observables();

private:
  glsim::StochasticEnvironment SE;
  VicsekIntegratorParameters par;

  void common_init();
  void vserial(oarchive_t &ar) {ar << *this;}
  void vserial(iarchive_t &ar) {ar >> *this;}
  template <typename Archive>
  void serialize(Archive &ar,const unsigned int version);
  friend class boost::serialization::access;

public:
  static const int class_version=1;
} ;

BOOST_CLASS_VERSION(VicsekEnvironment,VicsekEnvironment::class_version);


/*****************************************************************************
 *
 * VicsekSimulation
 * 
 */

class VicsekSimulation : public glsim::Simulation {
public:
  VicsekSimulation(VicsekEnvironment& e,glsim::OLconfiguration &c,VicsekInteraction *i);
  ~VicsekSimulation();
  const char* name() const {return "Vicsek's model (as an overdamped ISM)";}

  void step();
  void log(),log_start_sim();

protected:
  void update_observables();

  VicsekEnvironment&         env;
  glsim::OLconfiguration& conf;
  VicsekInteraction       *inter;

private:
  double  (*confb)[3];
  double  mass,xDt,Dt;
  double  v0sq,sv,sa,rho,c0,c1dt,c1mc2,c2dt;
  glsim::Gaussian_distribution* noise;
} ;


/******************************************************************************
 *
 * Observable
 *
 */

class VicsekObservable_parameters : public glsim::Parameters {
public:
  VicsekObservable_parameters(const char* scope);
} ;

class VicsekObservable :  public glsim::SBObservable {
public:
  VicsekObservable(VicsekEnvironment&,glsim::OLconfiguration&);

  void interval_and_file();
  void write_header();
  void observe();

private:
  VicsekEnvironment  &env;
  glsim::OLconfiguration &conf;
  VicsekObservable_parameters par;

  void update();
} ;

inline VicsekObservable::VicsekObservable(VicsekEnvironment& e,glsim::OLconfiguration &c) :
  SBObservable(e),
  env(e),
  conf(c),
  par(e.scope())
{}

#endif /* VICI_HH */
