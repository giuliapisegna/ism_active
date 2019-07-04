/*
 * isi_potential.hh -- The inertial spin integrator with external positional potential
 *
 * Integrates the equation of motion for an inertial spin model (with
 * Brownian noise), without specifiying the velocity-velocity
 * interaction.
 *
 * Tomás S. Grigera <tgrigera@iflysib.unlp.edu.ar>
 *
 * La Plata, July 2015
 *
 */

#ifndef ISI_HH
#define ISI_HH

#include "glsim/random.hh"
#include "glsim/stochastic.hh"
#include "glsim/simulation.hh"
#include "glsim/olconfiguration.hh"
#include "glsim/observable.hh"
#include "glsim/md.hh"
#include "glsim/avevar.hh"

#include "social_potential.hh"

/*****************************************************************************
 *
 * Environment
 * 
 */

class ISMParameters : public glsim::Parameters {
public:
  ISMParameters(const char* scope);
} ;

class ISMEnvironment : public glsim::SimEnvironment {
public:
  ISMEnvironment(const char* scope=glsim::Parameters::default_scope);

  long    ISsteps;
  double  time_step;
  bool    fixed_graph;
  double  temperature;
  double  eta;
  bool    rescale_v0;
  
  //Potential properties
  double  A;
  double  omega;

  // System info
  double  v0;
  double  social_mass[2];

  // Observables
  int     total_number;
  double  total_social_mass;
  double  social_total_energy,social_potential_energy,social_kinetic_energy;
  double  polarization,v0sqave;
  double  Vcm[3], Rcm[3], total_spin[3],total_spinsq;
  glsim::AveVar<true> spinsqavar;

//Landmark position
  double rl[3];
//derivatime in time of this position
  double drl_x;
//modifica!!
 //double R0;

protected:
  void    init_local() {SimEnvironment::init_local(); common_init();}
  void    warm_init_local() {SimEnvironment::warm_init_local(); common_init();}
  void    update_observables();

private:
  glsim::StochasticEnvironment SE;
  ISMParameters par;

  void common_init();
  void vserial(oarchive_t &ar) {ar << *this;}
  void vserial(iarchive_t &ar) {ar >> *this;}
  template <typename Archive>
  void serialize(Archive &ar,const unsigned int version);
  friend class boost::serialization::access;

public:
  static const int class_version=1;
} ;

BOOST_CLASS_VERSION(ISMEnvironment,ISMEnvironment::class_version);


/*****************************************************************************
 *
 * ISMSimulation
 * 
 */

class ISMSimulation : public glsim::Simulation {
public:
  ISMSimulation(ISMEnvironment& e,glsim::OLconfiguration &c,VicsekInteraction *i);
  ~ISMSimulation();
  const char* name() const {return "Inertial spin model";}

  void step();
  void log(),log_start_sim();

protected:
  void update_observables();

  ISMEnvironment&         env;
  glsim::OLconfiguration& conf;
  VicsekInteraction       *inter;

private:
  double  (*confb)[3];
  double  mass,xDt,Dt;
  double  v0sq,sv,sa,rho,c0,c1dt,c1mc2,c2dt;
  glsim::BivariateGaussian_distribution* noise;
} ;


/******************************************************************************
 *
 * Observable
 *
 */

class ISMObservable_parameters : public glsim::Parameters {
public:
  ISMObservable_parameters(const char* scope);
} ;

class ISMObservable :  public glsim::SBObservable {
public:
  ISMObservable(ISMEnvironment&,glsim::OLconfiguration&);

  void interval_and_file();
  void write_header();
  void observe();

protected:
  ISMEnvironment  &env;
  glsim::OLconfiguration &conf;
  ISMObservable_parameters par;

  void update();
} ;

inline ISMObservable::ISMObservable(ISMEnvironment& e,glsim::OLconfiguration &c) :
  SBObservable(e),
  env(e),
  conf(c),
  par(e.scope())
{}

#endif /* ISI_HH */
