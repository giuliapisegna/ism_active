#include <math.h>
#include <iostream>

#ifndef T3DVECS_HH
#define T3DVECS_HH

inline double dotp(double x[3],double y[3])
{
  return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}

inline double modsq(const double x[3])
{
  return x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
}

inline void vprod(double *prod,double *a,double *b)
{
  prod[0]=a[1]*b[2] - a[2]*b[1];
  prod[1]=a[2]*b[0] - a[0]*b[2];
  prod[2]=a[0]*b[1] - a[1]*b[0];
}

inline void normalize(double *v)
{
  double mod=sqrt(modsq(v));
  v[0]/=mod;
  v[1]/=mod;
  v[2]/=mod;
}

inline std::ostream& operator<<(std::ostream& o,double *v)
{
  o << '(' << v[0] << ',' << v[1] << ',' << v[2] << ')';
  return o;
}

#endif /* T3DVECS_HH */

