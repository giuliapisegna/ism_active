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
    nrows_(ROWS), ncols_(COLS)
  {arr=new T[nrows_*ncols_];}
  ~column_major_array() {delete[] arr;}
  T& operator()(int row,int col)
  {return *(arr+lpos(row,col));}
  T* data() {return arr;}

  const int nrows() const {return nrows_;}
  const int ncols() const {return ncols_;}

private:
  T   *arr;
  int nrows_,ncols_;

  int lpos(int row,int col)
  {return row+nrows_*col;}
} ;
