/*
 * ckiso.cc -- Compute the isotropic (averaged over k) C(k)
 *
 */

#include <valarray>
#include <gsl/gsl_sf_trig.h>

#include "glsim/blib.hh"
#include "glsim/offlattice.hh"
#include "glsim/timecorr.hh"
//#include "glsim/binvec.hh"
#include "glsim/Grk.hh"

#include "3dvecs.hh"


/******************************************************************************
 *
 * This computes C(r)
 *
 */
// class Cr {
// public:
//   Cr(const glsim::OLconfiguration& c,double rnn);
//   ~Cr();

//   /// Number of bins
//   int    size() const {return esize;}
//   /// Number of bins including all particles (up to distances larger that box/2)
//   int    extended_size() const {return vec->nbins();}
//   /// \f$\Delta r\f$
//   double deltar() const {return vec->delta();}
//   /// Radius corresponding to ith bin
//   double r(int i) const {return vec->binc(i);}
//   /// 4\pi\rho r^2 \Delta r C(r) at ith bin, excluding the self (i=j) contribution
//   double corr_nonorm(std::size_t i) const
//   {return (*vec)[i]/nconf;}

//   void push(const glsim::OLconfiguration&);
//   void push(const Cr&);

// private:
//   double Cr0;
//   glsim::Binned_vector<double> *vec;
//   double rnn;
//   int    esize;   // The size as seen from the outside.  Internally we
// 		  // keep a few more bins for convenience, but we must
// 		  // not report values of r larger than half the box
// 		  // side
//   int    nconf;
//   double rho;

//   friend std::ostream& operator<<(std::ostream& o,const Cr& corr);
// } ;

// Cr::Cr(const glsim::OLconfiguration &c,double rnn_) :
//   rnn(rnn_), nconf(0), Cr0(0.)
// {
//   double rmax=sqrt(c.box_length[0]*c.box_length[0] + c.box_length[1]*c.box_length[1] +
// 		   c.box_length[2]*c.box_length[2])/1.99;
//   double lmax=std::min(c.box_length[0],c.box_length[1]);
//   lmax=std::min(lmax,c.box_length[2]);
//   int nbins=(int) ceil(rmax/(0.1*rnn));
//   vec=new glsim::Binned_vector<double>(nbins,0,rmax);
//   esize=vec->binn(lmax/2.);
//   for (int i=0; i<nbins; (*vec)[i++]=0.) ;
//   rho=c.number_density();
// }

// Cr::~Cr()
// {
//   delete vec;
// }

// void Cr::push(const glsim::OLconfiguration& conf)
// {
//   for (int i=0; i<conf.N; ++i) {
//     Cr0+=dotp(conf.v[i],conf.v[i])/conf.N;
//   }
//   for (int i=0; i<conf.N-1; ++i) {
//     for (int j=i+1; j<conf.N; ++j) {
//       double r=sqrt(conf.distancesq(i,j));
//       (*vec)[r]+=2*dotp(conf.v[i],conf.v[j])/conf.N;
//     }
//   }
//   nconf++;
// }

// void Cr::push(const Cr& C)
// {
//   Cr0+=C.Cr0;
//   for (int i=0; i<vec->nbins(); ++i)
//     (*vec)[i]+=(*(C.vec))[i];
//   nconf+=C.nconf;
// }
 
// std::ostream& operator<<(std::ostream& o,const Cr& corr)
// {
//   double f=1./(4*M_PI*corr.rho*corr.nconf*corr.vec->delta());
//   o << 0. << " " << corr.Cr0/corr.nconf << '\n';
//   for (int i=0; i < corr.size(); ++i) {
//     double r=corr.vec->binc(i);
//     o << r << "  " << f*(*corr.vec)[i]/(r*r) << '\n';
//   }
// }

// /*****************************************************************************
//  *
//  * C(k)
//  *
//  */

// class Ckiso {
// public:
//   Ckiso(const glsim::OLconfiguration&,const Cr&,int);
//   std::vector<double> Ck;
//   double              deltak;
// } ;

// Ckiso::Ckiso(const glsim::OLconfiguration &conf,const Cr& C,int nk)
// {
//   deltak=0.2*M_PI/conf.box_length[0];

//   Ck.resize(nk);
//   for (int ik=0; ik<nk; ++ik) {
//     double k=ik*deltak/M_PI; // Divide by PI because of GSL's definition 
//                              // of the sinc function
//     // We start with 1 instead of 0 because corr_nonorm, which is used
//     // in the loop below, does note include the self (i=j) contribution
//     // also note that we don't multiply by \Delta r because we use corr_nonorm
//     // which return value is already multiplied by \Delta r
//     Ck[ik]=1;
//     for (int i=0; i<C.extended_size(); ++i) {  // extended_size ensures all pairs are included in sum
//       double kr=k*C.r(i);
//       Ck[ik]+=C.corr_nonorm(i)*gsl_sf_sinc(kr);
//     }
//   }
// }

// std::ostream& operator<<(std::ostream& o,const Ckiso& ck)
// {
//   for (int i=0; i < ck.Ck.size(); ++i) {
//     double k=i*ck.deltak;
//     o << k << "  " << ck.Ck[i] << '\n';
//   }
// }


/******************************************************************************
 *
 * Options and main
 *
 */

struct optlst {
public:
  bool   multithreaded;
  bool   phase_average;
  bool   space_ave_fields;
  int    nr,nk;
  // double rnn;
  std::vector<std::string> ifiles;
  std::string c_of_r_file;

} options;

double aveV[3];

class CLoptions : public glsim::UtilityCL {
public:
  CLoptions();
  void show_usage() const;
} ;

CLoptions::CLoptions() : glsim::UtilityCL("ckiso")
{
  hidden_command_line_options().add_options()
    // ("rnn",po::value<double>(&options.rnn)->required(),"mean NN distance")
    ("nr",po::value<int>(&options.nr)->required(),"number of bins for C(r)")
    ("nk",po::value<int>(&options.nk)->required(),"number of wavevectors to compute")
    ("ifiles",po::value<std::vector<std::string> >(&options.ifiles)->required(),"input files")
    ;
  command_line_options().add_options()
    ("multithreaded,m",po::bool_switch(&options.multithreaded),
     "process several configurations in parallel")
    ("c-of-r,r",po::value<std::string>(&options.c_of_r_file),
     "also write C(r) to given file")
    ("phase-average,P",po::bool_switch(&options.phase_average),
     "connect using phase average (default is space, or single-configuration, average)")
    ("space-average-fields,S",po::bool_switch(&options.space_ave_fields),
     "use the definition for space-averaged fields")
     ;

  positional_options().add("nr",1).add("nk",1).add("ifiles",-1);
}

void CLoptions::show_usage() const
{
  std::cerr
    << "usage: " << progname << "[options] nr nk ifile [ifile ....]\n\n"
    << "This program computes the (connected) velocity correlation in real\n"
    << "and Fourier space (the Fourier space correlation is analytically\n"
    << "averaged over directions of the wavevector).\n"
    << "nr is the number of space bins, nk is the number of k values.\n"
    // << "nk is the number of k values, rnn is the mean nearest-neighbour distance.\n"
    // << "This computes the isotropic C(k) for nk wavevectors kn (Delta k computed automatically).\n"
    << "The definitions are\n\n"
    << "      C(k) = (1/N) \\sum_{ij} \\delta \\hat v_i \\delta \\hat v_i \\sinc( k r_{ij})\n\n"
    << "where \\delta v_i = v_i/v_0 - (1/N) \\sum_i v_i/v_0 and \\delta\\hat v_i is\n"
    << "normalized,\n"
    << "         \\sum_{ij} \\delta \\hat v_i \\delta \\hat v_i \\delta(r-r_{ij})\n"
    << "  C(r) = ---------------------------------------------------------------\n"
    << "                     \\sum_ij \\delta(r-r_{ij}).\n\n"
    << "The alternative definitions for the real space correlation (option -S below)\n"
    << "uses the space-average fields, so that\n"
    << "  CS(r) = (1/Nrho) \\sum_{ij} \\delta \\hat v_i \\delta \\hat v_i \\delta(r-r_{ij)).\n\n"
    << "\n"
    << " Options:\n";
  show_command_line_options(std::cerr);
}


/*****************************************************************************
 *
 * main and prepare velocities
 *
 */
void compute_phase_ave(glsim::H5_multi_file &ifs,glsim::OLconfiguration& conf)
{
  aveV[0]=aveV[1]=aveV[2]=0.;
  long n=0;
  while (ifs.read()) {
    for (int i=0; i<conf.N; ++i) {
      ++n;
      aveV[0]+=conf.v[i][0]/conf.N;
      aveV[1]+=conf.v[i][1]/conf.N;
      aveV[2]+=conf.v[i][2]/conf.N;
    }
  }
  ifs.rewind();
  aveV[0]/=n;
  aveV[1]/=n;
  aveV[2]/=n;
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

void prepare_vel(glsim::OLconfiguration& conf)
{
  if (options.phase_average) { // substract phase-space average of V
    for (int i=0; i<conf.N; ++i) {
      conf.v[i][0]-=aveV[0];
      conf.v[i][1]-=aveV[1];
      conf.v[i][2]-=aveV[2];
    }
  }
  else normalize_vel(conf);
}

// void wmain(int argc,char *argv[])
// {
//   CLoptions o;
//   o.parse_command_line(argc,argv);

//   glsim::OLconfiguration conf;
//   glsim::OLconfig_file   cfile(&conf);
//   glsim::H5_multi_file   ifs(options.ifiles,cfile);

//   ifs.read();
//   Cr C(conf,options.rnn);
//   ifs.rewind();

//   if (options.multithreaded) {
//      #pragma omp parallel
//      {
//        glsim::OLconfiguration confloc=conf;
//        Cr Cloc(confloc,options.rnn);

//        #pragma omp for schedule(static) nowait
//        for (hsize_t rec=0; rec<ifs.size(); ++rec) {
// 	 #pragma omp critical
// 	 {
// 	   ifs.read();
// 	   confloc=conf;
// 	 }
// 	 normalize_vel(confloc);
// 	 Cloc.push(confloc);
//        }

//        #pragma omp critical
//        C.push(Cloc);
//      }
//   } else {
//     while (ifs.read()) {
//       normalize_vel(conf);
//       C.push(conf);
//     }
//   }

//   if (options.c_of_r_file.length()>0) {
//     std::ofstream rfile(options.c_of_r_file);
//     int ixi=0;
//     while (C.corr_nonorm(ixi)>0 && ixi <C.size()) ++ixi;
//     rfile << "# Velocity correlation (real space)\n";
//     rfile << "#\n# r0 (crossing zero) = " << C.r(ixi) << '\n';
//     rfile << "#\n# r    C(r)\n";
//     rfile << C;
//   }
//   Ckiso Ck(conf,C,options.nk);
//   std::cout << "# k   C(k)\n";
//   std::cout << Ck;
// }

void wmain(int argc,char *argv[])
{
  CLoptions o;
  o.parse_command_line(argc,argv);

  glsim::OLconfiguration conf;
  glsim::OLconfig_file   cfile(&conf);
  glsim::H5_multi_file   ifs(options.ifiles,cfile);

  ifs.read();
  glsim::Grk C(conf.box_length,options.nr,options.nk,options.space_ave_fields);
  ifs.rewind();

  if (options.phase_average) compute_phase_ave(ifs,conf);

  if (options.multithreaded) {
     // #pragma omp parallel
     // {
     //   glsim::OLconfiguration confloc=conf;
     //   Cr Cloc(confloc,options.rnn);

     //   #pragma omp for schedule(static) nowait
     //   for (hsize_t rec=0; rec<ifs.size(); ++rec) {
     // 	 #pragma omp critical
     // 	 {
     // 	   ifs.read();
     // 	   confloc=conf;
     // 	 }
     // 	 normalize_vel(confloc);
     // 	 Cloc.push(confloc);
     //   }

     //   #pragma omp critical
     //   C.push(Cloc);
     // }
  } else {
    while (ifs.read()) {
      prepare_vel(conf);
      C.push(conf,conf.v);
    }
  }

  if (options.c_of_r_file.length()>0) {
    std::ofstream rfile(options.c_of_r_file);
    rfile << "# Velocity correlation (real space)\n";
    if (!options.phase_average) {
      int ixi=0;
      while (C.Gr(ixi)>0 && ixi <C.sizer()) ++ixi;
      rfile << "# (space average)\n# r0 (crossing zero) = " << C.r(ixi) << '\n';
    } else rfile << "# (phase average, <v> = " << aveV[0] << ' ' << aveV[1] << ' ' << aveV[2] << '\n';
    rfile << "#\n# r    C(r)\n";
    rfile << C.pGr();
  }
  std::cout << "# k   C(k)\n";
  std::cout << C.pGk();
}

int main(int argc, char *argv[])
{
  return glsim::StandardEC(argc,argv,wmain);
}
