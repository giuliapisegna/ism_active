
#include <limits.h>
#include <math.h>

#include "glsim/olconfiguration.hh"
#include "glsim/interactions.hh"
#include "glsim/mdenvironment.hh"
#include "glsim/vverletmd.hh"

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
  Q{-.5,.5}
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

  double RAD;
  double center[2][3];
} ;

CoulConf::CoulConf() :
  RAD(1.),
  center{ {0,0,0}, {5,0,0} }
{
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
  r[0][0]=center[0][0]+RAD;
  r[0][1]=center[0][1]+0;
  r[0][2]=center[0][2]+0;
  r[1][0]=center[1][0]+RAD;
  r[1][1]=center[1][1]+0;
  r[1][2]=center[1][2]+0;
  if (v) delete[] v;
  v=new double[2][3];
  v[0][0]=0;
  v[0][1]=RAD;
  v[0][2]=0;
  v[1][0]=0;
  v[1][1]=RAD;
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

class ConstrainedMD : public glsim::VVerletMD {
public:
  ConstrainedMD(glsim::MDEnvironment& e,CoulConf &c,glsim::Interactions *i);

  void step();

private:
  glsim::MDEnvironment&   env;
  CoulConf&               conf;
  glsim::Interactions     *inter;

  double           Dt,Dt2,Dtsq2;
} ;

ConstrainedMD::ConstrainedMD(glsim::MDEnvironment& e,CoulConf &c,
			     glsim::Interactions *i) :
  VVerletMD(e,c,i),
  env(e),
  conf(c),
  inter(i)
{
  Dt=env.time_step;
  Dt2=Dt/2;
  Dtsq2=Dt*Dt2;
}


void ConstrainedMD::step()
{
  for (int i=0; i<conf.N; ++i) {
    conf.r[i][0] += Dt*conf.v[i][0] + Dtsq2*conf.a[i][0];
    conf.r[i][1] += Dt*conf.v[i][1] + Dtsq2*conf.a[i][1];
    conf.r[i][2] += Dt*conf.v[i][2] + Dtsq2*conf.a[i][2];
    conf.v[i][0] += Dt2*conf.a[i][0];
    conf.v[i][1] += Dt2*conf.a[i][1];
    conf.v[i][2] += Dt2*conf.a[i][2];
  }
  // Shake-like corections for r
  for (int i=0; i<conf.N; ++i) {
    double drc[3];
    double *r = conf.r[i];
    drc[0] = r[0]-conf.center[i][0];
    drc[1] = r[1]-conf.center[i][1];
    drc[2] = r[2]-conf.center[i][2];
    double drcsq= drc[0]*drc[0] + drc[1]*drc[1] + drc[2]*drc[2];
    double gamma=conf.RAD/sqrt(drcsq);
    r[0]=conf.center[i][0]+gamma*drc[0];
    r[1]=conf.center[i][1]+gamma*drc[1];
    r[2]=conf.center[i][2]+gamma*drc[2];
  }
  Epot=inter->acceleration_and_potential_energy(conf);
  for (int i=0; i<conf.N; ++i) {
    conf.v[i][0] += Dt2*conf.a[i][0];
    conf.v[i][1] += Dt2*conf.a[i][1];
    conf.v[i][2] += Dt2*conf.a[i][2];
  }
  // Shake-like corections for v
  for (int i=0; i<conf.N; ++i) {
    double drc[3];
    double *r = conf.r[i];
    drc[0] = r[0]-conf.center[i][0];
    drc[1] = r[1]-conf.center[i][1];
    drc[2] = r[2]-conf.center[i][2];
    double gamma=(conf.v[i][0]*drc[0] + conf.v[i][1]*drc[1] + conf.v[i][2]*drc[2])
      /(conf.RAD*conf.RAD);
    conf.v[i][0]-=gamma*drc[0];
    conf.v[i][1]-=gamma*drc[1];
    conf.v[i][2]-=gamma*drc[2];
  }
  inter->fold_coordinates(conf);

  env.time_completed+=Dt;
  env.time_in_run+=Dt;
  env.run_completed = env.steps_in_run>=env.MDsteps;
  conf.time=env.time_completed;
  conf.step=env.steps_completed;
}




/*****************************************************************************/


void wmain(int argc, char *argv[])
{
  glsim::MDEnvironment  env;
  CoulConf conf;
  Coulomb CP(env.scope());
  glsim::SimulationCL CL("GS_ljmd","(C) 2015 Tomas S. Grigera",env.scope());
  CL.parse_command_line(argc,argv);
  glsim::prepare(CL,env,conf);

  glsim::Interactions_isotropic_pairwise_naive<Coulomb> inter(CP,conf);
  ConstrainedMD sim(env,conf,&inter);
  sim.run();
  env.save();
  conf.save(env.configuration_file_fin);
}

int main(int argc, char *argv[])
{
  return glsim::StandardEC(argc,argv,wmain);
}
