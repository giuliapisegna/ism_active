/*
 * ism_create_conf.cc -- create configuration for the ISM
 *
 */

#include <assert.h>

#include "glsim/parameters.hh"
#include "glsim/random.hh"
#include "glsim/cerrors.h"
#include "glsim/olconfiguration.hh"
#include "3dvecs.hh"

struct scomp {
  int    N;
  double v0;
  double s0;
  double boxl[3];
  bool   planar;
  bool   polarize;

  scomp() : N(0) {}
  ~scomp() {}
} ;

void confined_positions(glsim::OLconfiguration &conf,scomp &SC,double dconf)
{
  conf.N=SC.N;
  conf.step=0;
  conf.time=0;
  conf.box_angles[0]=conf.box_angles[1]=conf.box_angles[2]=90.;
  memcpy(conf.box_length,SC.boxl,3*sizeof(double));

  int i;
  conf.id=new short[conf.N];
  for (i=0; i<conf.N; conf.id[i]=i++) ;
  conf.type=new short[conf.N];
  for (int j=0; j<SC.N; ++j)
    conf.type[j]=0;
  conf.r=new double[conf.N][3];
  glsim::Uniform_real ranx(0,dconf);
  glsim::Uniform_real rany(0,dconf);
  glsim::Uniform_real ranz(0,dconf);
  for (i=0; i<conf.N; ++i) {
    conf.r[i][0]=ranx();
    conf.r[i][1]=rany();
    conf.r[i][2]=0;
  }
  if (!SC.planar)
    for (i=0; i<conf.N; ++i) conf.r[i][2]=ranz();
}

void random_positions(glsim::OLconfiguration &conf,scomp &SC)
{
  conf.N=SC.N;
  conf.step=0;
  conf.time=0;
  conf.box_angles[0]=conf.box_angles[1]=conf.box_angles[2]=90.;
  memcpy(conf.box_length,SC.boxl,3*sizeof(double));

  int i;
  conf.id=new short[conf.N];
  for (i=0; i<conf.N; conf.id[i]=i++) ;
  conf.type=new short[conf.N];
  for (int j=0; j<SC.N; ++j)
    conf.type[j]=0;
  conf.r=new double[conf.N][3];
  glsim::Uniform_real ranx(0,conf.box_length[0]);
  glsim::Uniform_real rany(0,conf.box_length[1]);
  glsim::Uniform_real ranz(0,conf.box_length[2]);
  for (i=0; i<conf.N; ++i) {
    conf.r[i][0]=ranx();
    conf.r[i][1]=rany();
    conf.r[i][2]=0;
  }
  if (!SC.planar)
    for (i=0; i<conf.N; ++i) conf.r[i][2]=ranz();
}

void fcc_positions(glsim::OLconfiguration &conf,scomp &SC)
{
  int   ix,iy,iz,i;
  int   nx,ny,nz;
  double d,x,y,z,vol;

  conf.N=SC.N;
  conf.step=0;
  conf.time=0;
  conf.box_angles[0]=conf.box_angles[1]=conf.box_angles[2]=90.;
  memcpy(conf.box_length,SC.boxl,3*sizeof(double));
  conf.id=new short[conf.N];
  for (i=0; i<conf.N; conf.id[i]=i++) ;
  conf.type=new short[conf.N];
  for (int j=0; j<SC.N; ++j)
    conf.type[j]=0;
  conf.r=new double[conf.N][3];

  vol=conf.box_length[0]*conf.box_length[1]*conf.box_length[2];
  d=powf(2*conf.N/vol,-1./3.);
  nx=ceil(conf.box_length[0]/d-1e-6);
  ny=ceil(conf.box_length[1]/d-1e-6);
  nz=ceil(conf.box_length[2]/d-1e-6);

  i=0;
  for (ix=0; ix<nx; ix++) {
    x=ix*d;
    for (iy=0; iy<ny; iy++) {
      y=iy*d;
      for (iz=0; iz<nz; iz++) {
	if ( (ix+iy+iz)%2 !=0 ) continue;
	if (i>=conf.N) goto exhausted;
	z=iz*d;
	conf.r[i][0]=x;
	conf.r[i][1]=y;
	conf.r[i][2]=z;
	i++;
      }
    }
  }
  goto go;
 exhausted:
  std::cout << "WARNING: box was not uniformly filled (use N=(m^3)/2)\n";

 go:
  std::cout << "First neighbour distance: " << sqrt(2.)*d << '\n';
}

void sc_positions(glsim::OLconfiguration &conf,scomp &SC)
{
  int   ix,iy,iz,i;
  int   nx,ny,nz;
  double d,x,y,z,vol;

  conf.N=SC.N;
  conf.step=0;
  conf.time=0;
  conf.box_angles[0]=conf.box_angles[1]=conf.box_angles[2]=90.;
  memcpy(conf.box_length,SC.boxl,3*sizeof(double));
  conf.id=new short[conf.N];
  for (i=0; i<conf.N; conf.id[i]=i++) ;
  conf.type=new short[conf.N];
  for (int j=0; j<SC.N; ++j)
    conf.type[j]=0;
  conf.r=new double[conf.N][3];

  vol=conf.box_length[0]*conf.box_length[1]*conf.box_length[2];
  d=powf(conf.N/vol,-1./3.);
  nx=ceil(conf.box_length[0]/d-1e-6);
  ny=ceil(conf.box_length[1]/d-1e-6);
  nz=ceil(conf.box_length[2]/d-1e-6);

  i=0;
  for (ix=0; ix<nx; ix++) {
    x=ix*d;
    for (iy=0; iy<ny; iy++) {
      y=iy*d;
      for (iz=0; iz<nz; iz++) {
	if (i>=conf.N) goto exhausted;
	z=iz*d;
	conf.r[i][0]=x;
	conf.r[i][1]=y;
	conf.r[i][2]=z;
	i++;
      }
    }
  }
  goto go;
 exhausted:
  std::cout << "WARNING: box was not uniformly filled (use N=m^3)\n";

 go:
  std::cout << "First neighbour distance: " << d << '\n';
}

void random_velocities(glsim::OLconfiguration &conf,double v0,std::vector<double> spin,bool polarize)
{
  delete[] conf.v;
  delete[] conf.a;
  conf.v = new double[conf.N][3];
  conf.a = new double[conf.N][3];
  glsim::Spherical3d_distribution sr;
  glsim::Uniform_real ran(0,2*M_PI);

  double u[3],v[3];
  double s0,s[3];

  if (spin.size()==1) s0=spin[0];
  else s0=sqrt(modsq(spin.data()));

  if (s0==0) {                      // just put random v, a=0
    for (int i=0; i<conf.N; ++i) {
      sr(u);
      conf.v[i][0]=v0*u[0];
      conf.v[i][1]=v0*u[1];
      conf.v[i][2]=v0*u[2];
    }
    return;
  }
    
  for (int i=0; i<conf.N; ++i) {
    if (spin.size()==1) {           // got modulus only, use random spin direction
      sr(u);
      s[0]=spin[0]*u[0];
      s[1]=spin[0]*u[1];
      s[2]=spin[0]*u[2];
    } else {
      s[0]=spin[0];
      s[1]=spin[1];
      s[2]=spin[2];
    }

    // Get something perpendicular to s
    u[0]=s[2];
    u[1]=s[0];
    u[2]=s[1];
    vprod(v,s,u);
    if (modsq(v)==0) {
      u[0]=0;
      u[1]=s[1];
      u[2]=s[2];
      vprod(v,s,u);
    }
    vprod(u,s,v);   // now u \perp v \perp s
    normalize(u);
    normalize(v);
    double theta= polarize ? 0 : ran();
    conf.v[i][0]=v0*(u[0]*cos(theta) + v[0]*sin(theta));
    conf.v[i][1]=v0*(u[1]*cos(theta) + v[1]*sin(theta));
    conf.v[i][2]=v0*(u[2]*cos(theta) + v[2]*sin(theta));

    vprod(conf.a[i],s,conf.v[i]);
  }
}


/*****************************************************************************
 *
 * options and main
 *
 */

static struct options_ {
  std::string         ofile;
  unsigned long       seed;
  int                 N;
  double              density;
  bool                fcc_lattice;
  bool                sc_lattice;
  double              confine;
  double              v0;
  std::vector<double> spin;
  bool                planar;
  bool                polarize;

  options_() :
    spin(3,0.)
  {}
  
} options;

class CLoptions : public glsim::UtilityCL {
public:
  CLoptions();
  void show_usage() const;
} ;

CLoptions::CLoptions() : UtilityCL("ism_greate_conf")
{
  hidden_command_line_options().add_options()
    ("out_file",po::value<std::string>(&options.ofile)->required(),"output file")
    ("Nparts",po::value<int>(&options.N)->required(),"total number of particles")
    ;
  
  command_line_options().add_options()
    ("density,d",po::value<double>(&options.density)->default_value(1.),"average density")
    ("seed,S",po::value<unsigned long>(&options.seed)->default_value(0),"random number seed")
    ("FCC,F",po::bool_switch(&options.fcc_lattice),
     "use FCC lattice for positions (default random positions)")
    ("SC,C",po::bool_switch(&options.sc_lattice),
     "use SC (simple cubic lattice) for positions")
    ("confined",po::value<double>(&options.confine)->default_value(-1),
     "confine particles in a cube of size arg with vertex at the origin (arg can be 0)")
    ("vmodulus,v",po::value<double>(&options.v0)->default_value(1.),
     "speed (equal for all particles)")
    ("spin,s",po::value<std::vector<double>>(&options.spin),"spin/chi (equal for all particles). Repeat 3 times for components, or give once for moduluse (random orientation")
    ("planar,P",po::bool_switch(&options.planar)->default_value(false),
     "keep positions and velocities in the XY plane, density is then surface density")
    ("polarize,p",po::bool_switch(&options.polarize)->default_value(false),
     "polarize in a direction perpendicular to the spin")
    ;

  positional_options().add("Nparts",1).add("out_file",1);
}

void CLoptions::show_usage() const
{
  std::cerr
    << "usage: " << progname << " [options] Nparts outfile\n\n"
    << "Create an off-lattice configuration with Nparts total particles and save to outfile.  Options:\n\n";
  show_command_line_options(std::cerr);

  std::cerr
    << "\nIf you give spin different from zero, the velocities will be given a random\n"
    << "orientation in the plane perpendicular to the spin, unless you also specify\n"
    << "polarize (-p),in which case all will be parallel.\n"
    << "\n";
}

void wmain(int argc,char *argv[])
{
  CLoptions opt;
  opt.parse_command_line(argc,argv);

  glsim::OLconfiguration conf;

  scomp SC;
  SC.N=options.N;
  SC.planar=options.planar;
  double volume=options.N/options.density;
  if (opt.count("spin"))
    options.spin=opt.value("spin").as<std::vector<double>>();
  SC.v0=options.v0;

  if (options.planar) {
    if (options.spin.size()==1) {
      options.spin.resize(3);
    } else {
      if (options.spin.size()!=3) options.spin.resize(3);
      if (options.spin[2]==0) options.spin[2]=1;
    }
    options.spin[0]=options.spin[1]=0;
    SC.boxl[2]=0.01;
    SC.boxl[0]=pow(volume,1./2.);
    SC.boxl[1]=SC.boxl[0];
  } else {
    SC.boxl[0]=pow(volume,1./3.);
    SC.boxl[2]=SC.boxl[1]=SC.boxl[0];
  }

  glsim::Random_number_generator RNG(glsim::gsl_rng_mt19937,options.seed);

  if (options.fcc_lattice)
    fcc_positions(conf,SC);
  else if (options.sc_lattice)
    sc_positions(conf,SC);
  else if (options.confine>=0)
    confined_positions(conf,SC,options.confine);
  else
    random_positions(conf,SC);
  conf.name="Created by ism_create_conf";
  random_velocities(conf,options.v0,options.spin,options.polarize);

  conf.save(options.ofile);
}

int main(int argc, char *argv[])
{
  return glsim::StandardEC(argc,argv,wmain);
}
