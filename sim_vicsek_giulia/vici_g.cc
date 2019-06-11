#include <cmath>
#include "vici_g.hh"


//costruttore di VicPar


VicPar::VicPar(const char* scope):
	Parameters(scope){
//questa e' solo una definizione dei parametri
	parameter_file_options().add_options()
    	("Vic.time_step",po::value<double>()->required(),"time step for Vic integrator")
        ("Vic.steps",po::value<long>()->required(),"number of steps to run")
        ("Vic.eta",po::value<double>()->required(),"dissipation coefficient")
         ("Vic.v0",po::value<double>()->required(),"fixed modulus of velocity")
("Vic.rc",po::value<double>()->required(),"metric cutoff")
	;
}


//costruttore di VicEnvironment, prima crea un environment di default
VicEnvironment::VicEnvironment(const char* scope):
	SimEnvironment(scope),
	Dt(0),
	Vicsteps(0),
	eta(0),
	v0(0),
	rc(0),
	polarization(0),
	SE(scope),
//sto chiamando il costruttore dei parametri, sa che deve leggere quella lista sopra
	par(scope){}

// quando la chiama simulation ha capito che bisogna costruirsi un env da capo
void VicEnvironment::init_local()
{
//init local della base class
  SimEnvironment::init_local();
  Dt=par.value("Vic.time_step").as<double>();
  Vicsteps=par.value("Vic.steps").as<long>();
  eta=par.value("Vic.eta").as<double>();
  v0= par.value("Vic.v0").as<double>();
  rc= par.value("Vic.rc").as<double>();
}

//warm = quando il run e' gia finito devo aggiornare l'environment se e' cambiato rispetto all'altro run
void VicEnvironment::warm_init_local()
{
  SimEnvironment::warm_init_local();
  Dt=par.value("Vic.time_step").as<double>();
  Vicsteps=par.value("Vic.steps").as<long>();
  eta=par.value("Vic.eta").as<double>();
  v0= par.value("Vic.v0").as<double>();
  rc= par.value("Vic.rc").as<double>();
}

//SPIEGAZIONI SU QUESTO:interfaccia con boost serialize, deve salvare tutto lo stato interno della classe e rileggerlo. Questa e' la classe che serve per leggere e salvare

template <typename Archive>
inline void VicEnvironment::serialize(Archive &ar,const unsigned int version)
{
  if (version!=class_version)
    throw glsim::Environment_wrong_version("VicEnvironment",version,class_version);
  ar & boost::serialization::base_object<SimEnvironment>(*this);
  ar & Dt;
  ar & Vicsteps;
  ar & eta;
  ar & v0;
  ar & rc;
}
BOOST_CLASS_VERSION(VicEnvironment,VicEnvironment::class_version);


//costruttore di VicSimulation


VicSimulation::VicSimulation(VicEnvironment& e, glsim::OLconfiguration &c):
//inizializzare Simulation (base class) sempre con gli argomenti e e c
Simulation(e,c),
env(e),
conf(c),
//inizializzazione dell'oggetto NN
NN(e.rc),
//inizializzare ran
ran(-e.eta*M_PI,+ e.eta*M_PI)
{
//inizializzo vmedio, array doppio= alloco memoria (new=malloc)
    vmedio= new double[conf.N][2];
//dichiaro rcsq
    rcsq=env.rc*env.rc;

//rebuild: crea tutte le coppie di primi vicini con il cutoff rc e con un deltar del 30% di rc
    conf.box_length[2]=20000;
    NN.rebuild(c,env.rc);

//PER LE INIZIALIZZAZIONI DI VELOCITA E POSIZIONI: il programma ism_create_conf puo' essere usato per creare qualsiasi tipo di configurazione, inizializza sia le velocita' che le posizioni iniziali, nel caso di sistema planare c'e' l'opzione -P
}

//distruttore della classe: tutte le variabili standard vengono automaticamente distrutte mentre devi liberare memoria per gli array, in questo caso vmedio
VicSimulation::~VicSimulation(){
    //libera memoria alla fine del programma
    delete[] vmedio;
}



//dichiarazione della classe step
void VicSimulation::step(){

    double distsq;
    double tetanew;
    double modvmedio;
  
    /*
    for(int i=0; i < conf.N; i++){
        //un vmedio per ogni particella
        vmedio[i][0] = conf.v[i][0];
        vmedio[i][1] = conf.v[i][1];
        
        //calcolo i vicini della particella i
        for(int j=0; j < conf.N; j++){
            if(j != i){
	

//sto calcolando la distanza senza PBC, sbagliato---> vedere ciclo sulle coppie
            distsq=(conf.r[i][0] - conf.r[j][0])*(conf.r[i][0] - conf.r[j][0])+ (conf.r[i][1] - conf.r[j][1])*(conf.r[i][1] - conf.r[j][1]);
            
            if(distsq > rcsq) continue;
            //se è un vicino lo metto nel vettore vmedio
            vmedio[i][0]+= conf.v[j][0];
            vmedio[i][1]+= conf.v[j][1];
           }
        } //fine for j
*/ 

for(int i=0; i < conf.N; i++){
	vmedio[i][0] = conf.v[i][0];
        vmedio[i][1] = conf.v[i][1];
}



//CICLO SULLE COPPIE, uso la classe nearestneighbours subcell. Questo e' un ciclo su tutte le coppie che si sono create nel costruttore con il rebuild
for(auto p = NN.pairs_begin(); p != NN.pairs_end(); p++){
	
   double rxmn=conf.ddiff(conf.r[p->first][0],conf.r[p->second][0],conf.box_length[0]);
   double rymn=conf.ddiff(conf.r[p->first][1],conf.r[p->second][1],conf.box_length[1]);
   double distsq=rxmn*rxmn+rymn*rymn;

   if(distsq > rcsq) continue;

   vmedio[p->first][0] += vmedio[p->second][0];
   vmedio[p->first][1] += vmedio[p->second][1];
   vmedio[p->second][0] += vmedio[p->first][0];
   vmedio[p->second][1] += vmedio[p->first][1];

}


//adesso ogni particella ha il suo vmedio-> ora devo ruotarlo di un numero random: per ogni i calcolo l'angolo del vettore e gli aggiungo un numero random generato da ran
for(int i=0; i < conf.N; i++){
	tetanew = atan2(vmedio[i][1],vmedio[i][0]) + ran();

	//calcolo il vmedio nuovo normalizzato
	vmedio[i][0] = cos(tetanew);
	vmedio[i][1] = sin(tetanew);

	// aggiornamento temporale delle velocita'
	conf.v[i][0] = env.v0*vmedio[i][0];
	conf.v[i][1] = env.v0*vmedio[i][1];

	// aggiornamento forward delle posizioni
	conf.r[i][0] += env.Dt*conf.v[i][0];
        conf.r[i][1] += env.Dt*conf.v[i][1];

    } 

	conf.fold_coordinates();
	NN.update(env.v0*env.Dt);


  env.time_completed+=env.Dt;
  env.time_in_run+=env.Dt;
  env.run_completed = env.steps_in_run>=env.Vicsteps;
  conf.step=env.steps_completed;
  conf.time=env.time_completed;
  update_observables();
//  glsim::logs(glsim::debug) << "Steps " << conf.step << " time " <<   env.time_completed << '\n';

}

void VicSimulation::update_observables(){

 	//qui calcolo ad esempio la polarizzazione
	
	double V[2];
	memset(V,0,2*sizeof(double));

	for(int i=0; i < conf.N; i++){
		V[0]+=conf.v[i][0];
    		V[1]+=conf.v[i][1];
	}

	env.polarization = sqrt(V[0]*V[0]+V[1]*V[1])/(conf.N*env.v0);
	env.Vcm[0]=V[0]/conf.N;
  	env.Vcm[1]=V[1]/conf.N;

}


//LOG: scrittura su terminale sullo stato della simulazione

void VicSimulation::log_start_sim(){
 char buff[300];
  
  Simulation::log_start_sim();
  glsim::logs(glsim::info) << "    Step       Time        Phi   \n";

  sprintf(buff," Initial            %10.3e \n",
	  env.polarization);
  glsim::logs(glsim::info) << buff;

}


void VicSimulation::log(){
  update_observables();
  static char buff[300];
  sprintf(buff,"%8ld %10.3e %10.3e\n", env.steps_completed, env.time_completed, env.polarization);
glsim::logs(glsim::info) << buff;
}



//OBSERVABLES


VicObservable_parameters::VicObservable_parameters(const char* scope) :
  glsim::Parameters(scope)
{
//acquisizione dei parametri dal file.ini
  parameter_file_options().add_options()
    ("Vic.obs_interval",po::value<int>()->default_value(0),
     "Interval for standard observation, 0=skip")
    ("Vic.obs_file_prefix",po::value<std::string>(),"Observation file prefix")
    ;
}

//assegnazione di questi parametri alle variabili e apertura file	
void VicObservable::interval_and_file()
{
  obs_interval=par.value("Vic.obs_interval").as<int>();
  obs_file_prefix=par.value("Vic.obs_file_prefix").as<std::string>();
}

//impostazione del file obsXXX
void VicObservable::write_header()
{
/*
  fprintf(of,"#   (1)| |     (2)| |     (3)| |     (4)| |     (5)| |     (6)| |     (7)| |     (8)| |     (9)| |     (10)| |    (11)| |    (12)| |    (13)| |    (14)| |    (15)| |    (16)| |    (17)|\n");
  fprintf(of,"#- Step and time -| |------- Social energy --------| | Av v^2 | |--- Center of mass velocity --| | Polariz.| |---------- Total spin --------|  |---------- Spin (single conf) ----------|\n");
  fprintf(of,"#   Step       Time  Potential    Kinetic      Total  <|v_i|^2>       VCMx       VCMy       VXMz         Phi         Sx         Sy         Sz    Average   Variance        Min        Max\n"); */

  fprintf(of,"#     (1)| |     (2)| |     (3)| | 	 (4)| |  	(5)| |  \n"); 
  fprintf(of,"#- Step and time -| | Polariz.|  Vcmx | |  Vcmy| |  \n");
  fprintf(of, "#   Step       Time      Phi	Vcmx	Vcmy \n");

}

//effettiva funzione di stampa
void VicObservable::observe() {
	update();
	fprintf(of, "%8ld %10.3e %10.3e %10.3e  %10.3e\n", env.steps_completed, env.time_completed, env.polarization, env.Vcm[0], env.Vcm[1]);
}


//funzione che rifa l'update "per sicurezza", e' uguale alla update_observables chiamata nella SIM
void VicObservable::update(){
	
	double V[2];
	memset(V,0,2*sizeof(double));

	for(int i=0; i < conf.N; i++){
		V[0]+=conf.v[i][0];
    		V[1]+=conf.v[i][1];
	}

	env.polarization = sqrt(V[0]*V[0]+V[1]*V[1])/(conf.N*env.v0);
	env.Vcm[0]=V[0]/conf.N;
  	env.Vcm[1]=V[1]/conf.N;
}







