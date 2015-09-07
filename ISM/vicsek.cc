/*
 * vicsek.cc -- The Vicsek model
 *
 * Tomás S. Grigera <tgrigera@iflysib.unlp.edu.ar>
 *
 * La Plata, September 2015
 *
 */

#include "glsim/histogram.hh"
#include "vici.hh"


/*****************************************************************************/

class Energy_histogram :  public glsim::SBObservable {
public:
  Energy_histogram(VicsekEnvironment&,int N,double J);
  ~Energy_histogram();
  void interval_and_file();
  void write_header() {}
  void observe();

private:
  VicsekEnvironment  &env;
  glsim::Histogram   histo;
} ;

Energy_histogram::Energy_histogram(VicsekEnvironment& e,int N,double J) :
  SBObservable(e),
  env(e),
  histo(200,-J*N,J*N)
{
  init_local();
}

void Energy_histogram::interval_and_file()
{
  obs_interval=1;
  obs_file_prefix="ehisto";
}

void Energy_histogram::observe()
{
  histo.push(env.social_potential_energy);
}

Energy_histogram::~Energy_histogram()
{
  if (!of) return;
  std::ostringstream o;
  o << histo;
  fputs(o.str().c_str(),of);
}

/*****************************************************************************/


void wmain(int argc, char *argv[])
{
  VicsekEnvironment  env;
  VicsekParameters VP;
  glsim::OLconfiguration conf;
  VicsekObservable obs(env,conf);
  // glsim::Trajectory traj(env,conf,
  // 			 glsim::OLconfig_file::options().r_frame());

  glsim::SimulationCL CL("vicsek (Vicsek's model)","(C) 2015 Tomas S. Grigera",env.scope());
  CL.parse_command_line(argc,argv);
  glsim::prepare(CL,env,conf);

  VicsekInteraction inter(VP,conf);
  VicsekSimulation sim(env,conf,&inter);
  Energy_histogram ehis(env,conf.N,inter.coupling());

// traj.observe_first();
  ehis.observe_first();
  obs.observe_first();
  sim.run();
  env.save();
  conf.save(env.configuration_file_fin);
}

int main(int argc, char *argv[])
{
  return glsim::StandardEC(argc,argv,wmain);
}
