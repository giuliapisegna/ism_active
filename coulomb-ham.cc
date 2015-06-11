/*
 * coulomb-ham.cc -- Two charged particles with spherical constraint
 *
 * This integrates the equations of motion directly in the hamiltonian
 * formulation.
 *
 * Tomás S. Grigera <tgrigera@iflysib.unlp.edu.ar>
 *
 * La Plata, June 2015
 *
 */

#include <limits.h>
#include <math.h>

#include "glsim/olconfiguration.hh"
#include "glsim/interactions.hh"
#include "glsim/mdenvironment.hh"
#include "glsim/mdobservable.hh"
#include "glsim/trajectory.hh"
#include "glsim/md.hh"

std::ostream& operator<<(std::ostream&o,double* d)
{
  o << '(' << d[0] << ", " << d[1] << ", " << d[2] << ')';
  return o;
}

class Coulomb {
public:
  Coulomb(const char*);
  void init(glsim::OLconfiguration &c);
  const char* name() const {return "Point charges";}
  bool within_cutoff(double dsq,int ta,int tb)
  {return true;}
  double mass(short t) {return qmass[t];}
  bool   has_efield() const {return true;}
  double cutoff() {return DBL_MAX;}

  double external_field(double *d,short t) {return 0;}
  double external_field(double *d,short t,double *f) { f[0]=f[1]=f[2]=0.; return 0.;}
  double pair_potential(double dsq,short t1,short t2);
  double pair_potential_r(double dsq,short t1,short t2) {return pair_potential(dsq,t1,t2);}
  double pair_potential(double dsq,short t1,short t2,double &vpsr);
  double pair_potential_r(double dsq,short t1,short t2,double &vpsr)
  {return pair_potential(dsq,t1,t2,vpsr);}

private:
  double Q[2],qmass[2];
} ;

Coulomb::Coulomb(const char *scope) :
  qmass{1.,1.},
  Q{5.,-5.}
  // Q{0.,0.}
{}

void Coulomb::init(glsim::OLconfiguration &c)
{
  if (c.type==0) {
    c.type=new short[c.N];
    memset(c.type,0,c.N*sizeof(short));
  }
}

double Coulomb::pair_potential(double dsq,short t1,short t2)
{
  return Q[t1]*Q[t2]/sqrt(dsq);
}

double Coulomb::pair_potential(double dsq,short t1,short t2,double &vpsr)
{
  double r=sqrt(dsq);
  vpsr = -Q[t1]*Q[t2]/(dsq*r);
  return Q[t1]*Q[t2]/r;
}

/*****************************************************************************
 *
 * CoulConf
 *
 */
class CoulConf : public glsim::OLconfiguration {
public:
  CoulConf();
  void init(const char* s);

  double RAD,RADSQ;
  double center[2][3];
} ;

CoulConf::CoulConf() :
  RAD(1.),
  center{ {2,2,2}, {7,2,2} }
  // center{ {0,0,0}, {0,0,0} }
{
  RADSQ=RAD*RAD;
}

void CoulConf::init(const char* s) {
  if (s) {
    load(s);
    return;
  }
  box_length[0]=100.;
  box_length[1]=100.;
  box_length[2]=100.;
  box_angles[0]=90.;
  box_angles[1]=90.;
  box_angles[2]=90.;
  N=2;
  if (r) delete[] r;
  r=new double[2][3];
  r[0][0]=center[0][0]+0;
  r[0][1]=center[0][1]+RAD;
  r[0][2]=center[0][2]+0;
  r[1][0]=center[1][0]+0;
  r[1][1]=center[1][1]+RAD;
  r[1][2]=center[1][2]+0;
  if (v) delete[] v;
  v=new double[2][3];
  v[0][0]=-2.;
  v[0][1]=0;
  v[0][2]=0;
  v[1][0]=1;
  v[1][1]=0;
  v[1][2]=0;
  if (id) delete[] id;
  id=new short[2];
  id[0]=0;
  id[1]=1;
  if (type) delete[] type;
  type=new short[2];
  type[0]=0;
  type[1]=1;
}


/*****************************************************************************
 *
 * Constrained MD
 *
 */

void vprod(double *prod,double *a,double *b)
{
  prod[0]=a[1]*b[2] - a[2]*b[1];
  prod[1]=a[2]*b[0] - a[0]*b[2];
  prod[2]=a[0]*b[1] - a[1]*b[0];
}

class ConstrainedMD : public glsim::MDSimulation {
public:
  ConstrainedMD(glsim::MDEnvironment& e,CoulConf &c,glsim::Interactions *i);
  const char* name() const {return "Constrainted";}

  void step();

private:
  glsim::MDEnvironment&   env;
  CoulConf&               conf;
  glsim::Interactions     *inter;

  void log();

  double            (*L)[3];
  double           Dt,Dt2,Dtsq;
} ;

ConstrainedMD::ConstrainedMD(glsim::MDEnvironment& e,CoulConf &c,
			     glsim::Interactions *i) :
  MDSimulation(e,c,i),
  env(e),
  conf(c),
  inter(i)
{
  Dt=env.time_step;
  Dt2=Dt/2;
  Dtsq=Dt*Dt;

  L=new double[conf.N][3];
  for (int i=0; i<conf.N; ++i) {
    double xmC[3];
    xmC[0]=conf.r[i][0]-conf.center[i][0];
    xmC[1]=conf.r[i][1]-conf.center[i][1];
    xmC[2]=conf.r[i][2]-conf.center[i][2];
    vprod(L[i],xmC,conf.v[i]);
    L[i][0]*=inter->mass(conf.type[i]);
  }
}


inline double dotp(double x[],double y[])
{
  return x[0]*y[0]+x[1]*y[1]+x[2]*y[2];
}

template <typename T> T sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

void rotate(double v[],double phi)
{
  double x,y; 

  x = cos(phi) * v[0] - sin(phi) * v[1];
  y = sin(phi) * v[0] + cos(phi) * v[1];
  v[0]=x;
  v[1]=y;
}

void ConstrainedMD::step()
{
  double l[3],dl[3];

  // Just 2 dimension for now
  for (int i=0; i<conf.N; ++i) {

    double xmC[3];
    xmC[0]=conf.r[i][0]-conf.center[i][0];
    xmC[1]=conf.r[i][1]-conf.center[i][1];
    xmC[2]=conf.r[i][2]-conf.center[i][2];

    // double check=dotp(L[i],xmC);
    // std::cerr << "L dot r " << check << "\n";

    double I = inter->mass(conf.type[i]) * conf.RADSQ;
    double phi = sqrt(dotp(L[i],L[i])) * Dt / I;
    rotate(xmC,phi);

    conf.r[i][0] = conf.center[i][0] + xmC[0];
    conf.r[i][1] = conf.center[i][1] + xmC[1];
    conf.r[i][2] = conf.center[i][2] + xmC[2];
  }


  // Forces at next step
  env.Epot=inter->force_and_potential_energy(conf);
  // Update L
  for (int i=0; i<conf.N; ++i) {
    double I = inter->mass(conf.type[i]) * conf.RADSQ;
    double xmC[3];
    xmC[0]=conf.r[i][0]-conf.center[i][0];
    xmC[1]=conf.r[i][1]-conf.center[i][1];
    xmC[2]=conf.r[i][2]-conf.center[i][2];

    vprod(dl,xmC,conf.a[i]);
    L[i][0] += Dt * dl[0];
    L[i][1] += Dt * dl[1];
    L[i][2] += Dt * dl[2];

    vprod(conf.v[i],L[i],xmC);
    conf.v[i][0]/=I;
    conf.v[i][1]/=I;
    conf.v[i][2]/=I;

    // double check=dotp(conf.v[i],xmC);
    // std::cerr << "v dot r " << check << "\n";

  }
  inter->fold_coordinates(conf);

  env.time_completed+=Dt;
  env.time_in_run+=Dt;
  env.run_completed = env.steps_in_run>=env.MDsteps;
  conf.time=env.time_completed;
  conf.step=env.steps_completed;
}

void ConstrainedMD::log()
{
  char buff[300];

  double xmC[3];
  xmC[0]=conf.r[0][0]-conf.center[0][0];
  xmC[1]=conf.r[0][1]-conf.center[0][1];
  xmC[2]=conf.r[0][2]-conf.center[0][2];

  double ldotr = dotp(L[0],xmC);
  double rsq = dotp(xmC,xmC);
  sprintf(buff,"L=(%8.4e, %8.4e, %8.4e) rsq %8.4e L dot r = %8.4e\n",
	  L[0][0],L[0][1],L[0][2],rsq,ldotr);
  glsim::logs(glsim::info) << buff;
}


/*****************************************************************************/


void wmain(int argc, char *argv[])
{
  glsim::MDEnvironment  env;
  CoulConf conf;
  Coulomb CP(env.scope());
  glsim::MDObservable obs(env,conf);
  glsim::Trajectory traj(env,conf,
			 glsim::OLconfig_file::options().r_frame());
  glsim::SimulationCL CL("GS_ljmd","(C) 2015 Tomas S. Grigera",env.scope());
  CL.parse_command_line(argc,argv);
  glsim::prepare(CL,env,conf);

  glsim::Interactions_isotropic_pairwise_naive<Coulomb> inter(CP,conf);
  traj.observe_first();
  ConstrainedMD sim(env,conf,&inter);
  obs.observe_first();
  sim.run();
  env.save();
  conf.save(env.configuration_file_fin);
}

int main(int argc, char *argv[])
{
  return glsim::StandardEC(argc,argv,wmain);
}
