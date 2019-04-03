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


int find_k(std::string& find_k_file)
{
  FILE *f=fopen(find_k_file.c_str(),"r");
  if (f==0) throw glsim::Clib_error(HERE);
  char buf[5000];
  double Ckmax=-1;
  int     n=0,nkmax;
  while (ungetc(fgetc(f),f)!=EOF) {
    double k,Ci,Cr,C;
    fgets(buf,5000,f);
    if (*buf=='#') continue;
    errno=0;
    if (sscanf(buf,"%lg %lg %lg",&k,&Cr,&Ci)!=3) {
      if (errno!=0) throw glsim::Clib_error(HERE);
      else throw glsim::Runtime_error("Cannot read, bad format");
    }
    C=Cr*Cr+Ci*Ci;
    if (C>Ckmax) {Ckmax=C; nkmax=n;}
    ++n;
  }
  return nkmax;
}

/*****************************************************************************
 *
 *  The following class will compute the current j(k,t) from which the
 *  correlation function can be computed.
 *
 */

class Current_jk {
public:
  Current_jk(double box_length[],int kn,int kdir);

  Current_jk& push_config(glsim::OLconfiguration &conf);
  const double* wavevector() {return k;}
  void clear() {jkx_.clear(); jky_.clear(); jkz_.clear();}

  const vcomplex& jkx() {return jkx_;}
  const vcomplex& jky() {return jky_;}
  const vcomplex& jkz() {return jkz_;}

private:
  int  kn,kdir;
  double k[3],deltak_[3];

  vcomplex jkx_,jky_,jkz_;
} ;

Current_jk::Current_jk(double box_length[],int kn_,int kdir_) :
  kn(kn_), kdir(kdir_)
{
  // choose deltak
  deltak_[0]=2*M_PI/box_length[0];
  deltak_[1]=2*M_PI/box_length[1];
  deltak_[2]=2*M_PI/box_length[2];

  memset(k,0,3*sizeof(double));
  k[kdir]=deltak_[kdir]*kn;
}

Current_jk& Current_jk::push_config(glsim::OLconfiguration &conf)
{
  vcomplex jjk(3,dcomplex(0.));

  for (int i=0; i<conf.N; i++) {
    double kr= k[0]*conf.r[i][0] + k[1]*conf.r[i][1] + k[2]*conf.r[i][2];
    jjk[0]+=conf.v[i][0]*exp(dcomplex(0,-1)*kr);
    jjk[1]+=conf.v[i][1]*exp(dcomplex(0,-1)*kr);
    jjk[2]+=conf.v[i][2]*exp(dcomplex(0,-1)*kr);
  }

  jkx_.push_back(jjk[0]/sqrt(conf.N));
  jky_.push_back(jjk[1]/sqrt(conf.N));
  jkz_.push_back(jjk[2]/sqrt(conf.N));

  return *this;
}

/*****************************************************************************
 *
 *  Class for C(k,t).  I'll choose one direction of k, parallel to the
 *  cooordinate axes (set on construction), this way PBCs are
 *  correctly enforced.
 *
 */

class Ckt {
public:
  Ckt(double deltat_,Current_jk *jk=0);
  Ckt& compute_Ckt();
  
  vcomplex     Cktx,Ckty,Cktz;

private:
  double      deltat;
  Current_jk* jk;

  friend std::ostream& operator<<(std::ostream&,const Ckt&); 
} ;

Ckt::Ckt(double deltat_,Current_jk *jk_) :
  deltat(deltat_),
  jk(jk_)
{}

Ckt& Ckt::compute_Ckt()
{
  int clen=jk->jkx().size()/2;
  cFFT FF(glsim::FFT::in_place);

  correlation_1d_tti_fft(jk->jkx(),FF,clen);
  Cktx=FF.tdata();
  correlation_1d_tti_fft(jk->jky(),FF,clen);
  Ckty=FF.tdata();
  correlation_1d_tti_fft(jk->jkz(),FF,clen);
  Cktz=FF.tdata();

  return *this;
}

struct optlst {
public:
  int         kn,kdir;
  bool        normalize,connect;
  bool        find_k;
  std::string find_k_file;
  bool        write_jk;
  std::string jk_file;
  std::vector<std::string> ifiles;

  optlst() : kn(-1),find_k(false) {}
} options;

std::ostream& operator<<(std::ostream& o,const Ckt& Ckt_)
{
  o << "# time   Ck'   Ck''  Cktx (etc)\n";
  dcomplex fac=1.;
  if (options.normalize) {
    dcomplex C0=Ckt_.Cktx[0]+Ckt_.Ckty[0]+Ckt_.Cktz[0];
    fac=1./C0;
  }
  for (int i=0; i<Ckt_.Cktx.size(); i++) {
    dcomplex cc=fac*(Ckt_.Cktx[i]+Ckt_.Ckty[i]+Ckt_.Cktz[i]);
    o << i*Ckt_.deltat << "  " << cc.real() << "  " << cc.imag() << "  "
      << Ckt_.Cktx[i].real() << ' ' << Ckt_.Cktx[i].imag() << "   "
      << Ckt_.Ckty[i].real() << ' ' << Ckt_.Ckty[i].imag() << "   "
      << Ckt_.Cktz[i].real() << ' ' << Ckt_.Cktz[i].imag() << "   "
      << '\n';
  }
  return o;
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
    ("ifiles",po::value<std::vector<std::string> >(&options.ifiles)->required(),"input files")
    ;
  command_line_options().add_options()
    ("wavevector,k",po::value<int>(&options.kn),
     "nk (integer): specify wavevector as multiple of Deltak")
    ("find-nk-in-file,f",po::value<std::string>(&options.find_k_file),
     "find nk as that which maximizes C(k), to be read from the given file")
    ("print-jk,j",po::value<std::string>(&options.jk_file),"write current j(k) to given file")
    ("kdirection,d",po::value<int>(&options.kdir)->required(),
     "direction of k: 0=x, 1=y, 2=z, 3=average of 0,1,2 (not good for noncubic boxes)")
    ("connect,C",po::bool_switch(&options.connect)->default_value(false),
     "Computed connected C(k,t) (Roman style)")
    ("normalize,N",po::bool_switch(&options.normalize)->default_value(false),
     "Normalize correlation to 1 at t=0")
     ;

  positional_options().add("ifiles",-1);
}

void CLoptions::show_usage() const
{
  std::cerr
    << "usage: " << progname << "[options] ifile [ifile ....]\n\n"
    << "Computes C(k,t) (def2) at wave vector kn (the desired multiple of Delta k).\n"
    << "Directions must be specified with option -d.\n"
    << "This computes\n\n"
    << "         C(k,t) = (1/N) \\sum_{ij} v_i(0) v_i(t) exp(i k r_{ij}(t)) / v0^2\n\n"
    << "or, if -C (connect) is given, the connected C(k), Roman style:\n\n"
    << "      C(k) = (1/N) \\sum_{ij} \\delta \\hat v_i \\delta \\hat v_i exp(i k r_{ij})\n\n"
    << "where \\delta v_i = v_i/v_0 - (1/N) \\sum_i v_i/v_0 and \\delta\\hat v_i is\n"
    << "normalized.  The La Plata style \n"
    << "connection can be obtained by setting C(k=0,t) to 0.\n"
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
  
  glsim::OLconfiguration conf;
  glsim::OLconfig_file   cfile(&conf);
  glsim::H5_multi_file   ifs(options.ifiles,cfile);

  if (options.kdir<0 || options.kdir>3)
    throw glsim::Runtime_error("kdir must be 0,1,2,3");

  // Check options and find k
  if (options.kn>=0 && options.find_k_file.size()>0)
    throw glsim::Runtime_error("Must give either -k or -f");
  if (options.find_k_file.size()>0) {
    options.kn=find_k(options.find_k_file);
    options.find_k=true;
  }
  if (options.kn<0)
    throw glsim::Runtime_error("Must give one of -k or -f");
  options.write_jk = options.jk_file.size()>0;

  double deltat=get_deltat(ifs,conf);
  Current_jk jk(conf.box_length,options.kn,options.kdir>2 ? 0 : options.kdir);

  double k=sqrt(modsq(jk.wavevector()));
  std::cout << "# C(k,t) at |k|=" << k
	    << options.find_k ? " (found from file " + options.find_k_file + ")\n" :
    "(given in command line)\n";

  if (options.kdir<=2)  { // Single direction

    while (ifs.read()) {
      conf.unfold_coordinates();
      substract_cm(conf);
      if (options.connect) normalize_vel(conf);
      jk.push_config(conf);
    }

    if (options.write_jk) {
      std::ofstream rfile(options.jk_file);
      rfile << "# t    jkx'  jkx''  jky'  jky''  jkz'  jkz''\n";
      for (int i=0; i<jk.jkx().size(); ++i)
	rfile << i*deltat << "  "
	      << jk.jkx()[i].real() << ' ' << jk.jkx()[i].imag() << "  "
	      << jk.jky()[i].real() << ' ' << jk.jky()[i].imag() << "  "
	      << jk.jkz()[i].real() << ' ' << jk.jkz()[i].imag() << "\n";
    }

    Ckt C(deltat,&jk);
    C.compute_Ckt();
    std::cout << "#\n";
    std::cout << C;

  } else { // average over rotations

    std::cout << "# AVERAGED over cubic lattice transformations\n";
    std::cout << "#\n";

    Ckt C(deltat);
    Ckt C2(deltat,&jk);

    for (int nop=0; nop<46; ++nop) {
      ifs.rewind();
      jk.clear();
      while (ifs.read()) {
	conf.unfold_coordinates();
	substract_cm(conf);
	if (options.connect) normalize_vel(conf);
	glsim::apply_cubic_operation(conf,conf.r,nop);
	jk.push_config(conf);
      }
      C2.compute_Ckt();

      if (C.Cktx.empty()) {
	C.Cktx.resize(C2.Cktx.size());
	for (int j=0; j<C.Cktx.size(); ++j) C.Cktx[j]=C2.Cktx[j]/dcomplex(46,0);
	C.Ckty.resize(C2.Ckty.size());
	for (int j=0; j<C.Ckty.size(); ++j) C.Ckty[j]=C2.Ckty[j]/dcomplex(46,0);
	C.Cktz.resize(C2.Cktz.size());
	for (int j=0; j<C.Cktz.size(); ++j) C.Cktz[j]=C2.Cktz[j]/dcomplex(46,0);
      } else {
	for (int j=0; j<C.Cktx.size(); ++j) C.Cktx[j]+=C2.Cktx[j]/dcomplex(46,0);
	for (int j=0; j<C.Ckty.size(); ++j) C.Ckty[j]+=C2.Ckty[j]/dcomplex(46,0);
	for (int j=0; j<C.Cktz.size(); ++j) C.Cktz[j]+=C2.Cktz[j]/dcomplex(46,0);
      }

    }

    std::cout << C;
  }

}

int main(int argc, char *argv[])
{
  return glsim::StandardEC(argc,argv,wmain);
}
