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
#include "isi_potential.hh"

/*****************************************************************************/


void wmain(int argc, char *argv[])
{
  ISMEnvironment  env;
  VicsekParameters VP;

  I_AM_HERE;

  ISMEnvironment  env;
  VicsekParameters VP;
  glsim::OLconfiguration conf;
  ISMObservable obs(env,conf);
  glsim::Trajectory traj(env,conf);
  glsim::SimulationCL CL("ism (inertial spin model) with a harmonic potential in the positions","(C) 2019 Tomas S. Grigera e Giulia Pisegna",env.scope());
  CL.parse_command_line(argc,argv);
  I_AM_HERE;
  glsim::prepare(CL,env,conf);
  I_AM_HERE;
  VicsekInteraction *inter;

  //create one inter which contains a positional potential

  std::cout << "conf N " << conf.N << '\n';

  inter = new MetricVicsekInteractionwPotential<>(VP,conf);
  I_AM_HERE;

/*

if(VP.value("Vicsek.metric").as<bool>() == 1 && VP.value("Vicsek.withfield").as<bool>() ==1){
        inter = new  MetricVicsekInteractionwField<>(VP,conf);
} else if(VP.value("Vicsek.metric").as<bool>() == 1) {
       inter = new   MetricVicsekInteraction<>(VP,conf);
} else{ inter = new  TopologicalVicsekInteraction(VP,conf);} */

     //I_AM_HERE;
   /* //old algorithm
    VicsekInteraction *inter = VP.value("Vicsek.metric").as<bool>() ?
    (VicsekInteraction*) new MetricVicsekInteraction<>(VP,conf):
    (VicsekInteraction*) new TopologicalVicsekInteraction(VP,conf); */

  ISMSimulation sim(env,conf,inter);
  I_AM_HERE;
  traj.observe_first();
  I_AM_HERE;
  obs.observe_first();
 I_AM_HERE;
  sim.run();
I_AM_HERE;
  env.save();
I_AM_HERE;
  conf.save(env.configuration_file_fin);
I_AM_HERE;
  delete inter;
I_AM_HERE;


}

int main(int argc, char *argv[])
{
  return glsim::StandardEC(argc,argv,wmain);
}
