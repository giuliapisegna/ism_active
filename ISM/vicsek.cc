/*
 * vicsek.cc -- The Vicsek model
 *
 * Tomás S. Grigera <tgrigera@iflysib.unlp.edu.ar>
 *
 * La Plata, September 2015
 *
 */

#include "vici.hh"

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
  // // traj.observe_first();
  obs.observe_first();
  sim.run();
  env.save();
  conf.save(env.configuration_file_fin);
}

int main(int argc, char *argv[])
{
  return glsim::StandardEC(argc,argv,wmain);
}
