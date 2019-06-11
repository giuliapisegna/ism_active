
#include "glsim/random.hh"
#include "glsim/stochastic.hh"
#include "glsim/simulation.hh"
#include "glsim/olconfiguration.hh"
#include "glsim/observable.hh"
#include "glsim/md.hh"
#include "glsim/avevar.hh"
#include "glsim/nneighbours.hh"

//ENVIRONMENT E PARAMETERS

class VicPar: public glsim::Parameters{
	//classe che contiene i parametri che prendiamo dal file.ini
	public:
//il costruttore fa solo lettura dei parametri, questa classe ha lo scopo di inizializzarsi con un .ini
	 VicPar(const char *scope);
};


class VicEnvironment : public glsim::SimEnvironment {
	public:
	 VicEnvironment(const char* scope=glsim::Parameters::default_scope);

         double Dt;
         double Vicsteps;
	 double eta;
	 double v0;
	 double rc;

//Observables
	 
 	 double polarization;
	 double Vcm[2];
	

	protected:
 	 void init_local(), warm_init_local();
	// void update_observables();

	private:
	VicPar par;
	glsim::StochasticEnvironment SE;

//---> dichiarazione della classe di boost serialize
	void vserial(oarchive_t &ar) {ar << *this;}
	void vserial(iarchive_t &ar) {ar >> *this;}
 	template <typename Archive>
  	void serialize(Archive &ar,const unsigned int version);
//serve a serializzare anche la parte privata della classe, altrimenti non puo' accedere
  	friend class boost::serialization::access;
	

	public:
  	static const int class_version=0;
};



class VicSimulation : public glsim::Simulation{
    public:
    VicSimulation(VicEnvironment& e, glsim::OLconfiguration& c);
    ~VicSimulation();
    
    const char* name() const {return "Vicsek 2d";}
    void step();
    void log();
    void log_start_sim();

    protected:

    void update_observables();
    VicEnvironment&         env;
    glsim::OLconfiguration& conf;
//creo l'oggetto dei primi vicini nella classe Sim
    glsim::NeighbourList_subcells NN;

    private:
    //al momento della dichiarazione del rumore devi dare solo il nome, nel costruttore gli dici chi deve usare
    glsim::Uniform_real ran;
 
    //vmedio= puntatore a un array di due double, lo dichiaro dentro la classe perchè non ha senso ridichiararlo ad ogni step
    double (*vmedio)[2];
//variabile interna che serve solo per fare i calcoli, la metto privato e la inizializzo nel costruttore
    double rcsq;
};



//OBSERVABLES

//classe che serve per acquisire i parametri per la stampa degli osservabili, tipo con che intervallo stampare gli osservabili
class VicObservable_parameters : public glsim::Parameters {
public:
  VicObservable_parameters(const char* scope);
} ;


//classe osservabili, sottoclasse di SBObservable nella base class
class VicObservable :  public glsim::SBObservable {
public:
  VicObservable(VicEnvironment&,glsim::OLconfiguration&);
//queste sono le funzioni da scrivere: 1) apertura file, 2) scrittura dell'intestazione e 3) scrittura del file
  void interval_and_file();
  void write_header();
  void observe();

protected:
  VicEnvironment  &env;
  glsim::OLconfiguration &conf;
  VicObservable_parameters par;

//questa e'una funzione che mi rifa l'update perche' questa classe e' meglio non si parli con quella della simulazioni. e' piu' una questione di sicurezza che altro, perche' tutto dovrebbe essere gia' aggiornato dalle chiamate dentro la sim.

  void update();
} ;


//perche' sta qui il costruttore? Le inline possono andare solo nel .hh, non nel .cc
inline VicObservable::VicObservable(VicEnvironment& e,glsim::OLconfiguration &c) :
  SBObservable(e),
  env(e),
  conf(c),
  par(e.scope())
{}

