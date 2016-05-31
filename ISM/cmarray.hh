/*
 * cmarray.hh -- column-major order arrays
 *
 */

// Column_major means members of one column are contiguous i.e.
// M(0,0) M(1,0) M(2,0) are contiguous
//
// In other words, the first dimension is contiguous

template <typename T=int>
class column_major_array {
public:
  column_major_array(int ROWS,int COLS) :
    nrows(ROWS), ncols(COLS)
  {arr=new T[nrows*ncols];}
  ~column_major_array() {delete[] arr;}
  T& operator()(int row,int col)
  {return *(arr+lpos(row,col));}
  T* data() {return arr;}

  const int nx() {return nrows;}
  const int ny() {return ncols;}

private:
  T   *arr;
  int nrows,ncols;

  int lpos(int row,int col)
  {return row+nrows*col;}
} ;
