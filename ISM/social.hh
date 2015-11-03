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

class MetricVicsekInteraction : public VicsekInteraction {
public:
  MetricVicsekInteraction(VicsekParameters &par,glsim::OLconfiguration& c,
			  glsim::MetricNearestNeighbours *NN=0);
  double social_potential_energy_and_force(glsim::OLconfiguration&,double b[][3]) {}
  double social_potential_energy_and_acceleration(glsim::OLconfiguration&,double b[][3]);
  void   fold_coordinates(glsim::OLconfiguration&,double maxdisp=-1);

private:
  bool                            own_NN;
  glsim::MetricNearestNeighbours *NN;
} ;

inline void MetricVicsekInteraction::fold_coordinates(glsim::OLconfiguration& conf,double maxdisp)
{
  conf.fold_coordinates();
  if (maxdisp<0)
    NN->rebuild(conf,rc);
  else
    NN->update(conf,maxdisp);
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
    NN->update(conf,maxdisp);
}

#endif /* SOCIAL_HH */
