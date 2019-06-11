
#include "glsim/blib.hh"
#include "glsim/offlattice.hh"

using namespace glsim;

static struct ooptions{
	std::vector<std::string> ifiles;
}options;


class CLoptions : public glsim::UtilityCL {
public:
  CLoptions();
  void show_usage() const;
} ;


CLoptions::CLoptions() : UtilityCL("ct_polarization")
{
  hidden_command_line_options().add_options()
    ("ifiles",po::value<std::vector<std::string> >(&options.ifiles)->required(),"input files")
    ;
  positional_options().add("ifiles",-1);
}



void CLoptions::show_usage() const
{
  std::cerr
    << "\nusage: " << progname << " [options] ifile [..]\n\n"
    << "Reads the given glsim trajectory/configuration files and writes them in ascii\n"
    << "COBBS 12-column format and compute the time correlation of polarization.\n"
    << "\nOptions:\n";
  show_command_line_options(std::cerr);
  std::cerr  << "\n";
  show_parameters(std::cerr);
  std::cerr  << "\n";
}



double polarization(glsim::OLconfiguration& conf,double V[])
{
  memset(V,0,2*sizeof(double));
  for (int i=0; i<conf.N; ++i) {

    double v0=sqrt(conf.v[i][0]*conf.v[i][0] + conf.v[i][1]*conf.v[i][1]);
    V[0]+=conf.v[i][0]/v0;
    V[1]+=conf.v[i][1]/v0;
  
  }

  V[0]/=conf.N;
  V[1]/=conf.N;

  return sqrt(V[0]*V[0]+V[1]*V[1]);
}



void wmain(int argc,char *argv[])
{
  CLoptions    opt;
  opt.parse_command_line(argc,argv);

  OLconfiguration conf;
  OLconfig_file cfile(&conf);
  H5_multi_file ifs(options.ifiles,cfile);
  glsim::logs.set_stream(std::cout,glsim::error);

  std::vector<double> Phi;
  double Phiav=0;
  double phi;
  

  	do{
    		double V[2];
		phi = polarization(conf,V);
		//std::cout << phi << '\n';
    		Phi.push_back(phi);
		
 	 }while(ifs.read());

  int Tmax = Phi.size();
  


for(int t=1; t < Tmax; t++){
	//  std::cout << Phi[t] << '\n';
	  Phiav += Phi[t];
}

  Phiav = Phiav/Tmax;
 // std::cout << Phiav << '\n';
 
  double Cnorm;

  for(int t=0; t < Tmax; t++){
	
	double C=0;

	for(int to = 1 ; to < (Tmax-t) ; to++){
 		C += (Phi[to]-Phiav)*(Phi[to+t]-Phiav);
        }
    	
	if(t==0){ Cnorm=C/(Tmax-t);
	//std::cout << Cnorm << '\n';
	}
	
   std::cout << t << "  " << C/(Cnorm*(Tmax-t)) <<  '\n';
  }

}


int main(int argc, char *argv[])
{
  return glsim::StandardEC(argc,argv,wmain);
}
