#include <limits.h>
#include <math.h>
#include "glsim/offlattice.hh"
#include "vici_g.hh"


void wmain(int argc, char *argv[]){
//creo environment, creo una conf, un obs e una traj: la obs la devo scrivere io per dirgli chi deve stampare. Invece la traj e'gia' tutta dentro glsim e quindi sa come fare il tutto
  VicEnvironment env;
  glsim::OLconfiguration conf;
  VicObservable obs(env,conf);
  glsim::Trajectory traj(env,conf);
 
  glsim::SimulationCL CL("Vicsek 2d","(C) 2015 Tomas S. Grigera",env.scope());
//parse_command_line acquisisce tutto dalla riga di comando, anche il file.ini: quindi attenzione! devi creare tutti gli oggetti che servono prima di questa chiamata
  CL.parse_command_line(argc,argv);

//il prepare cerca di capire se ci sono gia' degli env od obs gia' salvati da cui cominciare o se deve crearne di nuovi. NB: sempre prima della chiamata a sim.
  glsim::prepare(CL,env,conf);


  //glsim::NeighbourList_subcells NN(env.rc);

//creo oggetto di simulazione
  VicSimulation sim(env,conf);

//queste funzioni _first sono gia' nella base class. Quando attacco una simulazione a un'altra evitano che la prima riga del file che continua sia uguale all'ultima del file precedente. 
  traj.observe_first();
  obs.observe_first();
  sim.run();

//questa e' la parte che stampa direttamente sul file
  env.save();
  conf.box_length[2]=0;
  conf.save(env.configuration_file_fin);
}

int main(int argc, char *argv[]){
 return glsim::StandardEC(argc,argv,wmain);
}
