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
  virtual double social_potential_energy_and_force(glsim::OLconfiguration&,double b[][3])=0;
  virtual double social_potential_energy_and_acceleration(glsim::OLconfiguration&,double b[][3])=0;
  virtual void fold_coordinates(glsim::OLconfiguration&,double maxdisp=-1);
} ;

inline void SocialInteractions::fold_coordinates(glsim::OLconfiguration& conf,double maxdisp)
{
  conf.fold_coordinates();
}

/*****************************************************************************
 *
 * Vicsek
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

protected:
  bool              metric;
  double            chi,J,v0,v0sq,rc,rcsq;

private:
  VicsekParameters& par;
} ;

/*
 * Metric Vicsek 
 *
 */

template <typename NeighboursT=glsim::NeighbourList_subcells>
class MetricVicsekInteraction : public VicsekInteraction {
public:
  MetricVicsekInteraction(VicsekParameters &par,glsim::OLconfiguration& c,
			  NeighboursT *NN=0);
  double social_potential_energy_and_force(glsim::OLconfiguration&,double b[][3]) {}
  double social_potential_energy_and_acceleration(glsim::OLconfiguration&,double b[][3]);
  void   fold_coordinates(glsim::OLconfiguration&,double maxdisp=-1);

private:
  bool         own_NN;
  NeighboursT *NN;
} ;

template <typename NeighboursT>
MetricVicsekInteraction<NeighboursT>::MetricVicsekInteraction(VicsekParameters &p,
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

template <typename NeighboursT>
inline void MetricVicsekInteraction<NeighboursT>::
fold_coordinates(glsim::OLconfiguration& conf,double maxdisp)
{
  conf.fold_coordinates();
  if (maxdisp<0)
    NN->rebuild(conf,rc);
  else
    NN->update(maxdisp);
}

/*
   This is computes the social interactions, but only for the metric case.

   This routine assumes that the input velocities have modulus v0 (the
   modulus of the velocity is fixed in the Vicsek model).

   Note that the acceleration can be computed as F/m or as v0sq*F/chi.

*/
template <typename NeighboursT>
double MetricVicsekInteraction<NeighboursT>::
social_potential_energy_and_acceleration(glsim::OLconfiguration &conf,double b[][3])
{
  double E=0;
  memset(b,0,conf.N*3*sizeof(double));

  double efac=-J/v0sq;
  double ffac=J/chi;

  for (auto p = NN->pairs_begin(); p!=NN->pairs_end(); ++p) {

    double rxmn=conf.ddiff(conf.r[p->first][0],conf.r[p->second][0],conf.box_length[0]);
    double rymn=conf.ddiff(conf.r[p->first][1],conf.r[p->second][1],conf.box_length[1]);
    double rzmn=conf.ddiff(conf.r[p->first][2],conf.r[p->second][2],conf.box_length[2]);
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
  return E;
}

/*
 * Topological Vicsek 
 *
 */

class TopologicalVicsekInteraction : public VicsekInteraction {
public:
  TopologicalVicsekInteraction(VicsekParameters &par,glsim::OLconfiguration& c,
		    glsim::TopologicalNearestNeighbours *NN=0);
  double social_potential_energy_and_force(glsim::OLconfiguration&,double b[][3]) {}
  double social_potential_energy_and_acceleration(glsim::OLconfiguration&,double b[][3]);
  void   fold_coordinates(glsim::OLconfiguration&,double maxdisp=-1);

private:
  bool                                 own_NN;
  int                                  NNeighbours;
  glsim::TopologicalNearestNeighbours *NN;
} ;

inline void TopologicalVicsekInteraction::fold_coordinates(glsim::OLconfiguration& conf,double maxdisp)
{
  conf.fold_coordinates();
  if (maxdisp<0)
    NN->rebuild(conf,rc);
  else
    NN->update(maxdisp);
}

#endif /* SOCIAL_HH */
