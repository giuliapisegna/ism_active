/*
 * cktiso.cc -- Compute the isotropic (averaged over directions of k) C(k,t) at fixed k
 *
 */


#include <gsl/gsl_sf_trig.h>

#include "glsim/blib.hh"
#include "glsim/offlattice.hh"

#include "3dvecs.hh"

double corr(glsim::OLconfiguration& conf0,glsim::OLconfiguration& conf,double k)
{
  k/=M_PI; // Divide by PI because of GSL's definition 
            // of the sinc function
  double c=0;
  for (int i=0; i<conf.N; ++i) {
    c+=dotp(conf0.v[i],conf.v[i])/conf.N;
  }
  for (int i=0; i<conf.N-1; ++i) {
    for (int j=i+1; j<conf.N; ++j) {
      double r=sqrt(conf.distancesq(i,j));
      c+=2*dotp(conf0.v[i],conf.v[j])*gsl_sf_sinc(k*r)/conf.N;
    }
  }
  return c;
}

class comp_corr {
public:
  double        c;
  
  comp_corr(glsim::OLconfiguration &conf0_,glsim::OLconfiguration &conf_,
	    double k_,double c_) :
    conf0(conf0_), conf(conf_),
    k(k_), c(c_) {}
  comp_corr(const comp_corr &cc) :
    conf0(cc.conf0), conf(cc.conf), k(cc.k), c(0) {}

  void operator()(int i,int j,double dsq);
  void reduce(comp_corr &cc) {c+=cc.c;}

private:
  double                  k;
  glsim::OLconfiguration& conf0,conf;
} ;

void comp_corr::operator()(int i,int j,double dsq)
{
  double r=sqrt(dsq);
  c+=2*dotp(conf0.v[i],conf.v[j])*gsl_sf_sinc(k*r)/conf.N;
}

double corr_mt(glsim::OLconfiguration& conf0,glsim::OLconfiguration& conf,double k)
{
  static glsim::MetricNearestNeighbours NN(1.1*sqrt(modsq(conf.box_length)));

  k/=M_PI; // Divide by PI because of GSL's definition of the sinc function
  double c=0;
  #pragma omp parallel for schedule(static) reduction(+:c)
  for (int i=0; i<conf.N; ++i) {
    c+=dotp(conf0.v[i],conf.v[i])/conf.N;
  }
  comp_corr CC(conf0,conf,k,c);
  for_each_pair_mt(NN,CC);
  return CC.c;
}

double find_k(std::string& find_k_file)
{
  FILE *f=fopen(find_k_file.c_str(),"r");
  if (f==0) throw glsim::Clib_error(HERE);
  char buf[5000];
  double kmax,Ckmax=-1;
  while (ungetc(fgetc(f),f)!=EOF) {
    double k,C;
    fgets(buf,5000,f);
    if (*buf=='#') continue;
    if (sscanf(buf,"%lg %lg",&k,&C)!=2)
      throw glsim::Clib_error(HERE);
    if (C>Ckmax) {Ckmax=C; kmax=k;}
  }
  return kmax;
}


/******************************************************************************
 *
 * Options
 *
 */

struct optlst {
public:
  double      k;
  bool        find_k;
  std::string find_k_file;
  bool        normalize;
  bool        multi_thread;

  std::vector<std::string> ifiles;

  optlst() : k(-1) {}
} options;

class CLoptions : public glsim::UtilityCL {
public:
  CLoptions();
  void show_usage() const;
} ;

CLoptions::CLoptions() : glsim::UtilityCL("cktiso")
{
  hidden_command_line_options().add_options()
    ("ifiles",po::value<std::vector<std::string> >(&options.ifiles)->required(),"input files")
    ;
  command_line_options().add_options()
    ("normalize,N",po::bool_switch(&options.normalize)->default_value(false),
     "Normalize correlation to 1 at t=0")
    ("wavevector,k",po::value<double>(&options.k),"wavevector modulus")
    ("find-k-in-file,f",po::value<std::string>(&options.find_k_file),
     "find k as that which maximizes C(k), to be read from the given file")
    ("multi-thread,m",po::bool_switch(&options.multi_thread)->default_value(false),
     "use multiple threads to loop over pairs of particles")
     ;

  positional_options().add("ifiles",-1);
}

void CLoptions::show_usage() const
{
  std::cerr
    << "usage: " << progname << "[options] kn ifile [ifile ....]\n\n"
    << "Computes isotropic C(k,t) (def2) from given trajectory files.\n"
    << "This computes\n\n"
    << "         C_R(k,t) = (1/N) \\sum_{ij} \\hat \\delta v_i \\hat\\delta v_k(t) exp(i k r_{ij}(t)),\n\n"
    << "where \\delta v_i(t) = (v_i(t)/v_0 - (1/N) \\sum_i v_i(t)/v_0) and\n"
    << "\\hat\\delta v_i(t)=\\delta v_i(t)/\\sqrt(\\sum\\delta v_k(t)/N).\n"
    << "\n"
    << " Options (MUST give -k or -f):\n";
  show_command_line_options(std::cerr);
}

/*****************************************************************************
 *
 * main and deltat
 *
 */

void normalize_vel(glsim::OLconfiguration& conf)
{
  /* Compute polarization */
  double V[3];
  memset(V,0,3*sizeof(double));
  for (int i=0; i<conf.N; ++i) {
    double v0=sqrt(modsq(conf.v[i]));
    conf.v[i][0]/=v0;
    conf.v[i][1]/=v0;
    conf.v[i][2]/=v0;
    V[0]+=conf.v[i][0];
    V[1]+=conf.v[i][1];
    V[2]+=conf.v[i][2];
  }
  V[0]/=conf.N;
  V[1]/=conf.N;
  V[2]/=conf.N;

  /* Substract polarization */
  for (int i=0; i<conf.N; ++i) {
    conf.v[i][0]-=V[0];
    conf.v[i][1]-=V[1];
    conf.v[i][2]-=V[2];
  }

  /* Normalize */
  double fn=0;
  for (int i=0; i<conf.N; ++i)
    fn+=modsq(conf.v[i]);
  fn=sqrt(fn/conf.N);
  for (int i=0; i<conf.N; ++i) {
    conf.v[i][0]/=fn;
    conf.v[i][1]/=fn;
    conf.v[i][2]/=fn;
  }
}

double get_deltat(glsim::H5_multi_file &ifs,glsim::OLconfiguration &conf)
{
  ifs.read();
  double t0=conf.time;
  ifs.read();
  ifs.rewind();
  return conf.time-t0;
}

void substract_cm(glsim::OLconfiguration &conf)
{
  std::vector<double> CM=conf.center_of_mass();
  for (int i=0; i<conf.N; ++i) {
    conf.r[i][0]-=CM[0];
    conf.r[i][1]-=CM[1];
    conf.r[i][2]-=CM[2];
  }
}

void wmain(int argc,char *argv[])
{
  CLoptions o;
  o.parse_command_line(argc,argv);
  
  glsim::OLconfiguration conf0,conf;
  glsim::OLconfig_file   cfile(&conf);
  glsim::H5_multi_file   ifs(options.ifiles,cfile);

  // Check options and find k
  if (options.k>=0 && options.find_k_file.size()>0)
    throw glsim::Runtime_error("Must give either -k or -f");
  if (options.find_k_file.size()>0)
    options.k=find_k(options.find_k_file);
  if (options.k<0)
    throw glsim::Runtime_error("Must give one of -k or -f");

  double (*corrf)(glsim::OLconfiguration& conf0,glsim::OLconfiguration& conf,double k);
  if (options.multi_thread) corrf=corr_mt;
  else corrf=corr;

  // All ready, go

  double deltat=get_deltat(ifs,conf);
  std::vector<double> ckt;
  hsize_t ntimes=ifs.size()/2;
  ckt.resize(ntimes);
  for (hsize_t rect0=0; rect0<ifs.size(); ++rect0) {
    ifs.seek(rect0);
    ifs.read();
    conf.unfold_coordinates_reset();
    conf0=conf;
    substract_cm(conf0);
    normalize_vel(conf0);
    do {
      conf.unfold_coordinates();
      normalize_vel(conf);
      substract_cm(conf);
      ckt[ifs.pos()-1-rect0]=corrf(conf0,conf,options.k);

      std::cout << "done  " << rect0 << " " << ifs.pos()-1 << '\n';

    } while (ifs.read() && ifs.pos()-1-rect0<ntimes);
  }
  for (int i=0; i<ntimes; ++i)
    ckt[i]/=ifs.size()-i;

  // Write
  std::cout << "# C(k,t) for k = " << options.k << '\n';
  std::cout << "# t  C(k,t)\n";
  double norm = options.normalize ? 1./ckt[0] : 1.;
  for (int i=0; i<ntimes; ++i)
    std::cout << i*deltat << "   " << norm*ckt[i] << '\n';
}

int main(int argc, char *argv[])
{
  return glsim::StandardEC(argc,argv,wmain);
}
