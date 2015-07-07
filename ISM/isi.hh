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

class ISMSimulation : public glsim::MDSimulation {
public:
  ISMSimulation(ISMEnvironment& e,OLconfiguration &c,glsim::Interactions *i);
  ~ISMSimulation();
  const char* name() const {return "Inertial spin model";}

  void step();

private:
  glsim::ISMEnvironment&   env;
  OLconfiguration&         conf;
  glsim::Interactions     *inter;

private:
  double Dt;
  double c0,c1dt,c1mc2,c2dt;
  BivariateGaussian_distribution* noise;
} ;
