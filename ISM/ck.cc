/*
 * ck.cc -- Compute C(k) (analogous to the static structure factor)
 *
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
 *  Class for C(k).  I'll choose one direction of k, parallel to the
 *  cooordinate axes (set on construction), this way PBCs are
 *  correctly enforced.
 *
 */

class Ck {
public:
  Ck(double box_length[],int nk_,int kdir_,bool connected=false);

  Ck& push_config(glsim::OLconfiguration &conf);
  Ck& compute_Ck();
  const vcomplex& Ck_data() const {return Ck_;} 

private:
  bool                              connected;
  int                               nk,kdir,nsamp;
  double                            deltak_[3];
  std::vector<std::vector<double>>  k;
  double                            deltat;
  vcomplex                          Ck_;
  
  void polarization(glsim::OLconfiguration& conf,double V[]);

  friend std::ostream& operator<<(std::ostream&,const Ck&); 
} ;

Ck::Ck(double box_length[],int nk_,int kdir_,bool connected_) :
  connected(connected_), nk(nk_), kdir(kdir_), nsamp(0)
{
  // choose deltak
  deltak_[0]=2*M_PI/box_length[0];
  deltak_[1]=2*M_PI/box_length[1];
  deltak_[2]=2*M_PI/box_length[2];

  k.resize(3);
  for (int i=0; i<k.size(); ++i) {
    k[i].resize(3,0);
    k[i][i]=deltak_[i];
  }
  Ck_.resize(nk,0);
}

Ck& Ck::push_config(glsim::OLconfiguration &conf)
{
  dcomplex S;
  dcomplex vkx,vky,vkz;
  nsamp++;
  dcomplex unsamp=1./dcomplex(nsamp,0);
  for (int jk=0; jk<Ck_.size(); ++jk) {
    S=0;
    vkx=vky=vkz=0;
    for (int n=0; n<conf.N; ++n) {
      double kr= jk*k[kdir][0]*conf.r[n][0] + jk*k[kdir][1]*conf.r[n][1] + jk*k[kdir][2]*conf.r[n][2];
      vkx+=conf.v[n][0]*exp(dcomplex(0,-1)*kr);
      vky+=conf.v[n][1]*exp(dcomplex(0,-1)*kr);
      vkz+=conf.v[n][2]*exp(dcomplex(0,-1)*kr);
    }
    S=(vkx*conj(vkx)+vky*conj(vky)+vkz*conj(vkz))/dcomplex(conf.N,0.);

    // Add to the running average of C(k)
    dcomplex Q=S-Ck_[jk];
    Ck_[jk]+=Q*unsamp;
  }
  return *this;
}

std::ostream& operator<<(std::ostream& o,const Ck& Ck_)
{
  o << "# k   Ck'   Ck''\n";
  double dk = Ck_.kdir<3 && Ck_.kdir>=0 ? Ck_.deltak_[Ck_.kdir] : Ck_.deltak_[0];
  for (int i=0; i<Ck_.Ck_.size(); i++)
    o << i*dk << "  " << Ck_.Ck_[i].real() << "  " << Ck_.Ck_[i].imag() << '\n';
}

/******************************************************************************
 *
 * Options and main
 *
 */

struct optlst {
public:
  int  nk,kdir;
  bool connect;
  std::vector<std::string> ifiles;
} options;


class CLoptions : public glsim::UtilityCL {
public:
  CLoptions();
  void show_usage() const;
} ;

CLoptions::CLoptions() : glsim::UtilityCL("ck")
{
  hidden_command_line_options().add_options()
    ("nk",po::value<int>(&options.nk)->required(),"number of wavevectors to compute")
    ("ifiles",po::value<std::vector<std::string> >(&options.ifiles)->required(),"input files")
    ;
  command_line_options().add_options()
    ("kdirection,d",po::value<int>(&options.kdir)->required(),
     "direction of k: 0=x, 1=y, 2=z, 3=average over all cubic symmetry operations (not good for noncubic boxes)")
    ("connect,C",po::bool_switch(&options.connect)->default_value(false),
     "Computed connected C(k) (Roman style)")
     ;

  positional_options().add("nk",1).add("ifiles",-1);
}

void CLoptions::show_usage() const
{
  std::cerr
    << "usage: " << progname << "[options] nk ifile [ifile ....]\n\n"
    << "Computes C(k) for nk wavevectors kn (Delta k computed automatically).\n"
    << "Directions must be specified with option -d.\n"
    << "This computes\n\n"
    << "              C(k) = (1/N) \\sum_{ij} v_i v_i exp(i k r_{ij}) / v0^2\n\n"
    << "or, if -C (connect) is given, the connected C(k), Roman style:\n\n"
    << "      C(k) = (1/N) \\sum_{ij} \\delta \\hat v_i \\delta \\hat v_i exp(i k r_{ij})\n\n"
    << "where \\delta v_i = v_i/v_0 - (1/N) \\sum_i v_i/v_0 and \\delta\\hat v_i is\n"
    << "normalized.  The La Plata style \n"
    << "connection can be obtained by computing C(k) and setting C(k=0) to 0.\n"
    << "\n"
    << " Options:\n";
  show_command_line_options(std::cerr);
}


/*****************************************************************************
 *
 * main and normalize
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
  fn= fn==0 ? 1 : sqrt(fn/conf.N);
  for (int i=0; i<conf.N; ++i) {
    conf.v[i][0]/=fn;
    conf.v[i][1]/=fn;
    conf.v[i][2]/=fn;
  }
}

void wmain(int argc,char *argv[])
{
  CLoptions o;
  o.parse_command_line(argc,argv);
  
  glsim::OLconfiguration conf,conft;
  glsim::OLconfig_file   cfile(&conf);
  glsim::H5_multi_file   ifs(options.ifiles,cfile);

  if (options.kdir<0 || options.kdir>3)
    throw glsim::Runtime_error("kdir must be 0, 1, 2, or 3");

  Ck C(conf.box_length,options.nk,options.kdir==3 ? 0 : options.kdir,options.connect);
  if (options.kdir==3) {
    while (ifs.read()) {
      conft=conf;
      if (options.connect) normalize_vel(conft);
      for (int n=0; n<47; ++n) {
	glsim::apply_cubic_operation(conf,conft.r,n);
	C.push_config(conft);
      }
    }
  } else 
    while (ifs.read()) {
      if (options.connect) normalize_vel(conf);
      C.push_config(conf);
    }
  std::cout << C;
}

int main(int argc, char *argv[])
{
  return glsim::StandardEC(argc,argv,wmain);
}
