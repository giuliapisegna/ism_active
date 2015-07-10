/*
 * isi.hh -- The inertial spin integrator
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
#include "glsim/md.hh"

#include "social.hh"

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
  double  v0;
  double  temperature;
  double  eta;

  // Observables
  int     total_number;
  double  total_social_mass;
  double  social_total_energy,social_potential_energy,social_kinetic_energy;
  double  polarization;

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
  ISMSimulation(ISMEnvironment& e,glsim::OLconfiguration &c,SocialInteractions *i);
  ~ISMSimulation();
  const char* name() const {return "Inertial spin model";}

  void step();
  void log(),log_start_sim();

protected:
  void update_observables();

  ISMEnvironment&         env;
  glsim::OLconfiguration& conf;
  SocialInteractions      *inter;

private:
  double  (*confb)[3];
  double  mass,xDt,Dt;
  double  v0sq,sv,sa,rho,c0,c1dt,c1mc2,c2dt;
  glsim::BivariateGaussian_distribution* noise;
} ;


#endif /* ISI_HH */

