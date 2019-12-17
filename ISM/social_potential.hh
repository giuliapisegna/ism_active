/*
 * social.hh -- Header for social interactions
 *
 * Tomás S. Grigera <tgrigera@iflysib.unlp.edu.ar>
 *
 * La Plata, July 2015
 *
 */

#ifndef SOCIAL_HH
#define SOCIAL_HH

#include "glsim/parameters.hh"
#include "glsim/olconfiguration.hh"
#include "glsim/nneighbours.hh"
#include "3dvecs.hh"
#include <cmath> 

/*****************************************************************************
 * 
 * SocialInteracions
 * 
 */

class SocialInteractions {
public:
  SocialInteractions() {}
  virtual ~SocialInteractions() {}

  virtual double social_mass(short type) const=0;
  virtual double social_potential_energy_and_force(glsim::OLconfiguration&,double b[][3], double rl[3], double drl_x)=0;
  virtual double social_potential_energy_and_acceleration(glsim::OLconfiguration&,double b[][3], double rl[3], double drl_x)=0;
} ;


/*****************************************************************************
 *
 * Vicsek with external potential
 * 
 */

class VicsekParameters : public glsim::Parameters {
public:
  VicsekParameters(const char *scope=glsim::Parameters::default_scope);
} ;

class VicsekInteraction : public SocialInteractions {
public:
  VicsekInteraction(VicsekParameters &par,glsim::OLconfiguration& c);

  double social_mass(short type) const {return chi/v0sq;}
  double speed(short type=0) const {return v0;}
  double coupling() const {return J;}
  double cutoff() const {return rc;}

protected:
   bool              metric;
   bool		     withfield;
   bool		     withpotential;
   double            chi,J,v0,v0sq,rc,rcsq;
//stiffness del potenziale armonico
   double  	     k;


private:
  VicsekParameters& par;
  
} ;

/*
 * Metric Vicsek with external potential
 *
 */

template <typename NeighboursT=glsim::NeighbourList_subcells>
class MetricVicsekInteractionwPotential : public VicsekInteraction {
public:
  MetricVicsekInteractionwPotential(VicsekParameters &par,glsim::OLconfiguration& c,
			  NeighboursT *NN=0);
  double social_potential_energy_and_force(glsim::OLconfiguration&, double b[][3], double rl[3], double drl_x);
  double social_potential_energy_and_acceleration(glsim::OLconfiguration&, double b[][3], double rl[3], double drl_x);

private:
  double implement_social_interactions_wpotential(glsim::OLconfiguration&, double b[][3], double rl[3], double drl_x, double ffac, double ffpp);

  bool         own_NN;
  NeighboursT *NN;
} ;

template <typename NeighboursT>
MetricVicsekInteractionwPotential<NeighboursT>::MetricVicsekInteractionwPotential(VicsekParameters &p,
							      glsim::OLconfiguration &c,
							      NeighboursT *n) :
  VicsekInteraction(p,c)
{
  if (!metric) throw glsim::Runtime_error("You asked for topological interactions but created the metric object");
  rcsq=rc*rc;

  if (n) {
    NN=n;
    own_NN=false;
  } else {
    NN=new NeighboursT(rc);
    own_NN=true;
  }
  NN->rebuild(c,rc);
}

/*
   This computes the social interactions, but only for the metric case and in the case without magnetic field.

   These routine assumes that the input velocities have modulus v0 (the
   modulus of the velocity is fixed in the Vicsek model).

   The different between computing accelerations and force is just a prefactor, so we precompute the prefactor and call a common implementation routine.

   The acceleration can be computed as F/m or as v0sq*F/chi.

*/
template <typename NeighboursT>
double MetricVicsekInteractionwPotential<NeighboursT>::
social_potential_energy_and_force(glsim::OLconfiguration &conf, double b[][3],double rl[3], double drl_x)
{
  double ffac=J/v0sq;
  double ffpp= (k*chi)/v0sq;
  return implement_social_interactions_wpotential(conf, b,rl,drl_x, ffac,ffpp);
}

template <typename NeighboursT>
double MetricVicsekInteractionwPotential<NeighboursT>::
social_potential_energy_and_acceleration(glsim::OLconfiguration &conf, double b[][3],double rl[3], double drl_x)
{
  double ffac=J/chi;
  double ffpp = k;
  return implement_social_interactions_wpotential(conf, b, rl, drl_x, ffac, ffpp);
}

/* ffac should be J/chi to compute acceleration, J/v0sq for force; ffpp should be k to compute acceleration and k*chi/v0sq to compute force  */
template <typename NeighboursT>
double MetricVicsekInteractionwPotential<NeighboursT>::
implement_social_interactions_wpotential(glsim::OLconfiguration &conf, double b[][3], double rl[3], double drl_x, double ffac, double ffpp)
{
  double E=0;
  memset(b,0,conf.N*3*sizeof(double));

  double efac=-J/v0sq;

  for (auto p = NN->pairs_begin(); p!=NN->pairs_end(); ++p) {

//cambiare tutta questa parte di calcolo delle distanze, distanze vere senza pbc
    double rxmn= conf.r[p->first][0] - conf.r[p->second][0];
    double rymn= conf.r[p->first][1] - conf.r[p->second][1];
    double rzmn= conf.r[p->first][2] - conf.r[p->second][2];
    double dsq=rxmn*rxmn+rymn*rymn+rzmn*rzmn;

    if (dsq>rcsq) continue;

    E += efac*dotp(conf.v[p->first],conf.v[p->second]);

    b[p->first][0] += ffac * conf.v[p->second][0];
    b[p->first][1] += ffac * conf.v[p->second][1];
    b[p->first][2] += ffac * conf.v[p->second][2];
    b[p->second][0] += ffac * conf.v[p->first][0];
    b[p->second][1] += ffac * conf.v[p->first][1];
    b[p->second][2] += ffac * conf.v[p->first][2];

  }
 
      for(int i=0; i < conf.N; i++){

	 E += 0.5*k*((conf.r[i][0] - rl[0])*(conf.r[i][0] - rl[0])+ (conf.r[i][1] - rl[1])*(conf.r[i][1] - rl[1])+(conf.r[i][2] - rl[2])*(conf.r[i][2] - rl[2])); 


         b[i][0] += ffpp*(-conf.v[i][0]+drl_x);
	 b[i][1] += ffpp*(-conf.v[i][1]);
	 b[i][2] += ffpp*(-conf.v[i][2]); 
	 
       }

 return E;

}




#endif /* SOCIAL_HH */
