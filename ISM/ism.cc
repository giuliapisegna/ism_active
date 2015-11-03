/*
 * ism.cc -- The Inertial Spin Model
 *
 * Tomás S. Grigera <tgrigera@iflysib.unlp.edu.ar>
 *
 * La Plata, July 2015
 *
 */

#include <limits.h>
#include <math.h>

#include "isi.hh"

std::ostream& operator<<(std::ostream&o,double* d)
{
  o << '(' << d[0] << ", " << d[1] << ", " << d[2] << ')';
  return o;
}

/*****************************************************************************/


void wmain(int argc, char *argv[])
{
  ISMEnvironment  env;
  VicsekParameters VP;
  glsim::OLconfiguration conf;
  ISMObservable obs(env,conf);
  // glsim::Trajectory traj(env,conf,
  // 			 glsim::OLconfig_file::options().r_frame());
  glsim::SimulationCL CL("ism (inertial spin model)","(C) 2015 Tomas S. Grigera",env.scope());
  CL.parse_command_line(argc,argv);
  glsim::prepare(CL,env,conf);

  MetricVicsekInteraction inter(VP,conf);
  ISMSimulation sim(env,conf,&inter);
  // traj.observe_first();
  obs.observe_first();
  sim.run();
  env.save();
  conf.save(env.configuration_file_fin);
}

int main(int argc, char *argv[])
{
  return glsim::StandardEC(argc,argv,wmain);
}
