/*
 * ism-pl.cc -- Estimate persistence length in ISM
 *
 * Makes some assumptions about inital conf and .ini (see docs)
 *
 *  - N=2, perfectly polarized, coincident at r=0
 *  - will stop when mutual distance reaches cutoff
 *
 * Tomás S. Grigera <tgrigera@iflysib.unlp.edu.ar>
 *
 * La Plata, March 2018
 *
 */

#include <limits.h>
#include <math.h>

#include "glsim/offlattice.hh"
#include "isi.hh"

/*****************************************************************************/

class ISMplsimulation : public ISMSimulation {
public:
  ISMplsimulation(ISMEnvironment& e,glsim::OLconfiguration &c,VicsekInteraction *i) :
    ISMSimulation(e,c,i) {}
  void step();
} ;
  
void ISMplsimulation::step() {
  double rr[3];
  
  ISMSimulation::step();
  double dsq=conf.distancesq(0,1);
  env.run_completed = sqrt(dsq)>=inter->cutoff();
}
    
class ISMplObservable: public ISMObservable {
public:
  ISMplObservable(ISMEnvironment&e,glsim::OLconfiguration&c) :
    ISMObservable(e,c) {}

  void write_header();
  void observe();

} ;

void ISMplObservable::write_header()
{
  fprintf(of,"#   (1)| |     (2)| |     (3)| |     (4)| |     (5)| |     (6)| |     (7)| |     (8)| |     (9)| |     (10)| |    (11)| |    (12)| |    (13)| |    (14)| |    (15)| |    (16)| |    (17)|      (18)|\n");
  fprintf(of,"#- Step and time -| |------- Social energy --------| | Av v^2 | |--- Center of mass velocity --| | Polariz.| |---------- Total spin --------|  |---------- Spin (single conf) ----------|----------|\n");
  fprintf(of,"#   Step       Time  Potential    Kinetic      Total  <|v_i|^2>       VCMx       VCMy       VXMz         Phi         Sx         Sy         Sz    Average   Variance        Min        Max  MutalDsq\n");

}

void ISMplObservable::observe()
{
  update();
  double dsq=conf.distancesq(0,1);
  fprintf(of,"%8ld %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %11.4e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
	  env.steps_completed,env.time_completed,
	  env.social_potential_energy/conf.N,env.social_kinetic_energy/conf.N,env.social_total_energy/conf.N,
	  env.v0sqave,env.Vcm[0],env.Vcm[1],env.Vcm[2],
	  env.polarization,
	  env.total_spin[0],env.total_spin[1],env.total_spin[2],
	  env.spinsqavar.ave(),env.spinsqavar.var(),env.spinsqavar.min(),env.spinsqavar.max(),
	  dsq);
}

void wmain(int argc, char *argv[])
{
  ISMEnvironment  env;
  VicsekParameters VP;
  glsim::OLconfiguration conf;
  ISMplObservable obs(env,conf);
  glsim::Trajectory traj(env,conf);
  glsim::SimulationCL CL("Persistence length for ism","(C) 2018 Tomas S. Grigera",env.scope());
  CL.parse_command_line(argc,argv);
  glsim::prepare(CL,env,conf);

  conf.box_length[0]=conf.box_length[1]=conf.box_length[2]=4*VP.value("Vicsek.cutoff").as<double>();
  conf.a[0][0]=conf.a[0][1]=conf.a[0][2]=0.;
  conf.a[1][0]=conf.a[1][1]=conf.a[1][2]=0.;

  VicsekInteraction *inter = VP.value("Vicsek.metric").as<bool>() ?
    (VicsekInteraction*) new MetricVicsekInteraction<>(VP,conf) :
    (VicsekInteraction*) new TopologicalVicsekInteraction(VP,conf);
  ISMplsimulation sim(env,conf,inter);
  traj.observe_first();
  obs.observe_first();
  sim.run();
  env.save();
  conf.save(env.configuration_file_fin);
  delete inter;
}

int main(int argc, char *argv[])
{
  return glsim::StandardEC(argc,argv,wmain);
}
