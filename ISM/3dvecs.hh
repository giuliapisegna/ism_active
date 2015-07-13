

#ifndef T3DVECS_HH
#define T3DVECS_HH

inline double dotp(double x[3],double y[3])
{
  return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}

inline double modsq(double x[3])
{
  return x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
}

inline void vprod(double *prod,double *a,double *b)
{
  prod[0]=a[1]*b[2] - a[2]*b[1];
  prod[1]=a[2]*b[0] - a[0]*b[2];
  prod[2]=a[0]*b[1] - a[1]*b[0];
}


#endif /* T3DVECS_HH */

