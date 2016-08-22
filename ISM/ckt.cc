/*
 * ckt.cc -- Compute C(k,t) (analogous to the intermediate scattering function)
 *           aka  Def2 in k space
 *
 */

#define HAVE_LIBFFTW3

#include "glsim/blib.hh"
#include "glsim/offlattice.hh"
#include "glsim/timecorr.hh"
#include "glsim/fft.hh"

#include "3dvecs.hh"

using glsim::vcomplex;
using glsim::dcomplex;

#ifdef HAVE_LIBFFTW3
typedef glsim::RealFFTW     rFFT;
typedef glsim::ComplexFFTW  cFFT;
#else
typedef glsim::RealFFT_gsl_2n     rFFT;
typedef glsim::ComplexFFT_gsl_2n  cFFT;
#endif /* HAVE_LIBFFTW3 */


/*****************************************************************************
 *
 *  Class for C(k,t).  I'll choose one direction of k, parallel to the
 *  cooordinate axes (set on construction), this way PBCs are
 *  correctly enforced.
 *
 */

class Ckt {
public:
  Ckt(double box_length[],double deltat_,int kn_,int kdir_,bool connect);

  Ckt& push_config(glsim::OLconfiguration &conf);
  Ckt& compute_Ckt();
  const vcomplex& Ckt_data() const {return Ckt_;} 

private:
  bool                              connect;
  int                               Npart,kn,kdir;
  double                            deltak_[3];
  std::vector<std::vector<double>>  k;
  double                            deltat;
  std::vector<vcomplex>             jkx_,jky_,jkz_;
  vcomplex                          Ckt_;
  
  void polarization(glsim::OLconfiguration& conf,double V[]);
  vcomplex j_k(glsim::OLconfiguration &,double k[],double V[]);

  friend std::ostream& operator<<(std::ostream&,const Ckt&); 
} ;

Ckt::Ckt(double box_length[],double deltat_,int kn_,int kdir_,bool connect_) :
  connect(connect_),
  deltat(deltat_),
  kn(kn_), kdir(kdir_)
{
  // choose deltak
  deltak_[0]=2*M_PI/box_length[0];
  deltak_[1]=2*M_PI/box_length[1];
  deltak_[2]=2*M_PI/box_length[2];

  if (kdir==0 || kdir==1 || kdir==2) {
    k.resize(1);
  } else {
    k.resize(6);
    kdir=2;
  }
  for (int i=0; i<k.size(); ++i) {
    k[i].resize(3,0);
    if (kdir<0)
      k[i][-kdir-1]=-deltak_[-kdir-1]*kn;
    else
      k[i][kdir]=deltak_[kdir]*kn;
    kdir--;
  }
  jkx_.resize(k.size());
  jky_.resize(k.size());
  jkz_.resize(k.size());
}

void Ckt::polarization(glsim::OLconfiguration& conf,double V[])
{
  memset(V,0,3*sizeof(double));
  for (int i=0; i<conf.N; ++i) {
    double v0=sqrt(modsq(conf.v[i]));
    V[0]+=conf.v[i][0]/v0;
    V[1]+=conf.v[i][1]/v0;
    V[2]+=conf.v[i][2]/v0;
  }
  V[0]/=conf.N;
  V[1]/=conf.N;
  V[2]/=conf.N;
}

vcomplex Ckt::j_k(glsim::OLconfiguration &conf,double *k,double V[])
{
  std::vector<dcomplex> rk(3,dcomplex(0));

  for (int i=0; i<conf.N; i++) {
    double v0=sqrt(modsq(conf.v[i]));
    double kr= k[0]*conf.r[i][0] + k[1]*conf.r[i][1] + k[2]*conf.r[i][2];
    rk[0]+=(conf.v[i][0]/v0-V[0])*exp(dcomplex(0,-1)*kr);
    rk[1]+=(conf.v[i][1]/v0-V[1])*exp(dcomplex(0,-1)*kr);
    rk[2]+=(conf.v[i][2]/v0-V[2])*exp(dcomplex(0,-1)*kr);
  }
  return rk;
}

Ckt& Ckt::push_config(glsim::OLconfiguration &conf)
{
  Npart=conf.N;
  double V[3]={0.,0.,0.};
  if (connect) polarization(conf,V);
  for (int i=0; i<k.size(); ++i) {
    std::vector<dcomplex> jk_;
    jk_=j_k(conf,k[i].data(),V);
    jkx_[i].push_back(jk_[0]);
    jky_[i].push_back(jk_[1]);
    jkz_[i].push_back(jk_[2]);
  }

  return *this;
}

Ckt& Ckt::compute_Ckt()
{
  dcomplex fac=1./dcomplex(k.size()*Npart,0);
  int Cktlen=jkx_[0].size()/2;
  Ckt_.clear();
  Ckt_.resize(Cktlen,dcomplex(0,0));
  cFFT FF(glsim::FFT::in_place);

  for (int i=0; i<k.size(); ++i) {
    correlation_1d_tti_fft(jkx_[i],FF,Cktlen);
    for (int j=0; j<FF.tdata_rw().size(); j++) Ckt_[j]+=fac*FF.tdatum(j);
    correlation_1d_tti_fft(jky_[i],FF,Cktlen);
    for (int j=0; j<FF.tdata_rw().size(); j++) Ckt_[j]+=fac*FF.tdatum(j);
    correlation_1d_tti_fft(jkz_[i],FF,Cktlen);
    for (int j=0; j<FF.tdata_rw().size(); j++) Ckt_[j]+=fac*FF.tdatum(j);
  }

  return *this;
}

struct optlst {
public:
  int  kn,kdir;
  bool normalize,connect;
  std::vector<std::string> ifiles;
} options;


std::ostream& operator<<(std::ostream& o,const Ckt& Ckt_)
{
  o << "# time   Ck'   Ck''\n";
  if (options.normalize) {
    dcomplex fac=1./Ckt_.Ckt_[0];
    for (int i=0; i<Ckt_.Ckt_.size(); i++) {
      dcomplex cc=Ckt_.Ckt_[i]*fac;
      o << i*Ckt_.deltat << "  " << cc.real()
	<< "  " << cc.imag() << '\n';
    }
  } else {
    for (int i=0; i<Ckt_.Ckt_.size(); i++)
      o << i*Ckt_.deltat << "  " << Ckt_.Ckt_[i].real()
	<< "  " << Ckt_.Ckt_[i].imag() << '\n';
  }
}

/******************************************************************************
 *
 * Options and main
 *
 */

class CLoptions : public glsim::UtilityCL {
public:
  CLoptions();
  void show_usage() const;
} ;

CLoptions::CLoptions() : glsim::UtilityCL("ckt")
{
  hidden_command_line_options().add_options()
    ("kn",po::value<int>(&options.kn)->required(),"wavevector")
    ("ifiles",po::value<std::vector<std::string> >(&options.ifiles)->required(),"input files")
    ;
  command_line_options().add_options()
    ("kdirection,d",po::value<int>(&options.kdir)->required(),
     "direction of k: 0=x, 1=y, 2=z, 3=average of 0,1,2 (not good for noncubic boxes)")
    ("connect,C",po::bool_switch(&options.connect)->default_value(false),
     "Computed connected C(k,t) (Roman style)")
    ("normalize,N",po::bool_switch(&options.normalize)->default_value(false),
     "Normalize correlation to 1 at t=0")
     ;

  positional_options().add("kn",1).add("ifiles",-1);
}

void CLoptions::show_usage() const
{
  std::cerr
    << "usage: " << progname << "[options] kn ifile [ifile ....]\n\n"
    << "Computes C(k,t) (def2) at wave vector kn (the desired multiple of Delta k).\n"
    << "Directions must be specified with option -d.\n"
    << "This computes\n\n"
    << "         C(k,t) = (1/N) \\sum_{ij} v_i(0) v_i(t) exp(i k r_{ij}(t)) / v0^2\n\n"
    << "or, if -C (connect) is given, the connected C(k), Roman style:\n\n"
    << "         C_R(k,t) = (1/N) \\sum_{ij} \\delta v_i \\delta v_k(t) exp(i k r_{ij}(t)),\n\n"
    << "where \\delta v_i(t) = v_i(t)/v_0 - (1/N) \\sum_i v_i(t)/v_0.  The La Plata style \n"
    << "connection can be obtained by computing C(k) and setting C(k=0) to 0.\n"
    << "\n"
    << " Options:\n";
  show_command_line_options(std::cerr);
}


/*****************************************************************************
 *
 * main and deltat
 *
 */

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
  
  glsim::OLconfiguration conf;
  glsim::OLconfig_file   cfile(&conf);
  glsim::H5_multi_file   ifs(options.ifiles,cfile);

  if (options.kdir<0 || options.kdir>3)
    throw glsim::Runtime_error("kdir must be 0,1,3");

  double deltat=get_deltat(ifs,conf);
  Ckt C(conf.box_length,deltat,options.kn,options.kdir,options.connect);
  while (ifs.read()) {
    conf.unfold_coordinates();
    substract_cm(conf);
    C.push_config(conf);
  }
  C.compute_Ckt();
  std::cout << C;
}

int main(int argc, char *argv[])
{
  return glsim::StandardEC(argc,argv,wmain);
}
