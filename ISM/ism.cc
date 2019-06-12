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

#include "glsim/offlattice.hh"
#include "isi.hh"

/*****************************************************************************/


void wmain(int argc, char *argv[])
{
  ISMEnvironment  env;
  VicsekParameters VP;

  glsim::OLconfiguration conf;
  ISMObservable obs(env,conf);
  glsim::Trajectory traj(env,conf);
  glsim::SimulationCL CL("ism (inertial spin model)","(C) 2015 Tomas S. Grigera",env.scope());
  CL.parse_command_line(argc,argv);
  glsim::prepare(CL,env,conf);
  VicsekInteraction *inter;

//Choice of type of integration: metric with magnetic field, metric without magnetic field, topological.

if(VP.value("Vicsek.metric").as<bool>() == 1 && VP.value("Vicsek.withfield").as<bool>() ==1){
        inter = new  MetricVicsekInteractionwField<>(VP,conf);
} else if(VP.value("Vicsek.metric").as<bool>() == 1) {
       inter = new   MetricVicsekInteraction<>(VP,conf);
} else{ inter = new  TopologicalVicsekInteraction(VP,conf);}

     //I_AM_HERE;
   /* //old algorithm
    VicsekInteraction *inter = VP.value("Vicsek.metric").as<bool>() ?
    (VicsekInteraction*) new MetricVicsekInteraction<>(VP,conf):
    (VicsekInteraction*) new TopologicalVicsekInteraction(VP,conf); */

  ISMSimulation sim(env,conf,inter);
I_AM_HERE;
  traj.observe_first();
  obs.observe_first();
I_AM_HERE;
  sim.run();
  env.save();
  conf.save(env.configuration_file_fin);
  delete inter;

}

int main(int argc, char *argv[])
{
  return glsim::StandardEC(argc,argv,wmain);
}
