//////////////////////////////////////////////////////
// Matrix TCL Pro 2.2
// Copyright (c) 2000-2011 Techsoft . All Rights Reserved.
// File name: Matrix.h - Header file for matrix class (do NOT include this file)
// Web: http://www.techsoftpl.com/matrix/
//
#include"iostream"

#ifndef _MATRIX_H_12a08ba0_8940_11d3_910c_002018382069_INCLUDED
#define _MATRIX_H_12a08ba0_8940_11d3_910c_002018382069_INCLUDED

namespace techsoft {

using std::istream;
using std::ostream;
using std::valarray;
using std::complex;
using std::numeric_limits;
using std::slice;
using std::slice_array;
using std::abs;
using std::sqrt;
using std::swap;
using std::setw;

#ifdef RANGE_CHECK_
using std::out_of_range;
using std::invalid_argument;
using std::runtime_error;
using std::exception;
#endif

#if defined(_MSC_VER)
using ::ceil;
using ::floor;
using ::srand;
using ::rand;
using ::sqrt;
using ::sqrtl;
#else
using std::ceil;
using std::floor;
using std::srand;
using std::rand;
using std::size_t;
#endif

//////////////////
// Temporary helping classes
//
class mslice;
class gmslice;
template <class T> class Mref;
template <class T> class Cmat_iter;
template <class T> class Mat_iter;
template <class T> class mslice_matrix;
template <class T> class gmslice_matrix;

enum MatArrayType { C_ARRAY, FORTRAN_ARRAY };

template <class T>
class matrix
{
public:
   typedef T value_type;
   typedef valarray<T> array_type;

   // Constructors
   matrix ();
   matrix (size_t nRow, size_t nCol);
   matrix (size_t nRow, size_t nCol, const T& e);
   matrix (size_t nRow, size_t nCol, const T* Array, MatArrayType aType = C_ARRAY);
   matrix (size_t nRow, size_t nCol, const valarray<T>& vec);
   ~matrix ();

   // Copy constructors
   matrix (const matrix<T>& m);
   matrix (const mslice_matrix<T>& sm);
   matrix (const gmslice_matrix<T>& sm);
#ifndef _TS_NO_MEMBER_TEMPLATES
   template <class X> matrix (const matrix<X>& m);
#endif

   // Assignment operators
   matrix<T>& operator= (const matrix<T>& m);
   matrix<T>& operator= (const mslice_matrix<T>& sm);
   matrix<T>& operator= (const gmslice_matrix<T>& sm);
   matrix<T>& operator= (const T& e);
#ifndef _TS_NO_MEMBER_TEMPLATES
   template <class X> matrix<T>& operator= (const matrix<X>& m);
#endif

   size_t size () const throw() { return mPtr->nrow * mPtr->ncol; }
   size_t typesize () const throw() { return sizeof(T); }

   size_t colsize () const throw() { return mPtr->nrow; }	// No. of elements in a column, i.e., it's row number
   size_t rowsize () const throw() { return mPtr->ncol; }	// No. of elements in a row, i.e., it's column number

   size_t rowno () const throw() { return mPtr->nrow; }		// No. of rows
   size_t colno () const throw() { return mPtr->ncol; }		// No. of columns

   void resize (size_t nRow, size_t nCol, const T& dval = T(0));
   void free ();

   Mat_iter<T> row (size_t i);
   Cmat_iter<T> row (size_t i) const;
   Mat_iter<T> column (size_t i);
   Cmat_iter<T> column (size_t i) const;
   Mat_iter<T> diag (int i=0);
   Cmat_iter<T> diag (int i=0) const;
   T trace (int i=0) const;

   Mat_iter<T> operator[] (size_t i) { return row(i); }
   Cmat_iter<T> operator[] (size_t i) const { return row(i); }
   Mat_iter<T> operator[] (int i) { return row(i); }
   Cmat_iter<T> operator[] (int i) const { return row(i); }

   Mat_iter<T> operator() (size_t i) { return column(i); }
   Cmat_iter<T> operator() (size_t i) const { return column(i); }
   Mat_iter<T> operator() (int i) { return column(i); }
   Cmat_iter<T> operator() (int i) const { return column(i); }

   Mref<T> operator() (size_t i, size_t j);
   const T& operator() (size_t i, size_t j) const;

   mslice_matrix<T> operator[] (mslice ms);
   gmslice_matrix<T> operator[] (gmslice gm);

   const matrix<T> operator[] (mslice ms) const;
   const matrix<T> operator[] (gmslice gm) const;

   // Unary operators
   matrix<T> operator+ () const { return *this; }
   matrix<T> operator- () const;
   matrix<T> operator~ () const;
   matrix<T> operator! () const;

   // Computed assignment
   matrix<T>& operator+= (const matrix<T>& m);
   matrix<T>& operator-= (const matrix<T>& m);
   matrix<T>& operator*= (const matrix<T>& m);
   matrix<T>& operator/= (const matrix<T>& m);
   matrix<T>& operator*= (const T& x);
   matrix<T>& operator/= (const T& x);

   // Miscellaneous methods
   void null ();
   void unit ();
   void rand (int rmin=-1, int rmax=1, int rseed=0);

   T sum () const { return mPtr->val.sum(); }
   T min () const { return mPtr->val.min(); }
   T max () const { return mPtr->val.max(); }

   matrix<T> apply (T (*fn)(T)) const;
   matrix<T> apply (T (*fn)(const T&)) const;

   matrix<T> apply (T (*fn)(size_t,size_t,T)) const;
   matrix<T> apply (T (*fn)(size_t,size_t,const T&)) const;

   // Utility methods
   bool solve (const valarray<T>& v, valarray<T>& s) const;
   bool solve (const matrix<T>& v, matrix<T>& s) const;

   matrix<T> adj () const;
   T cofact (size_t row, size_t col) const;

   T det () const;
   T cond () const;
   size_t rank () const;

   T norm1 () const;
   T norm2 () const;
   T normI () const;
   T normF () const;
   
   bool lud (valarray<size_t>& ri, T* pDet = NULL);
   void lubksb (const valarray<size_t>& ri, const valarray<T>& v, valarray<T>& s) const;
   void lumpove (const matrix<T>& ludm, const valarray<size_t>& ri, const valarray<T>& v, valarray<T>& s) const;
   bool solve_lu (const valarray<T>& v, valarray<T>& s) const { return solve( v, s); }
   
   bool inv ();
   bool inv_lu ();
   bool inv_sv ();
   bool inv_qr ();

   bool svd (matrix<T>& vc, valarray<T>& w);
   void svbksb (const matrix<T>& vc, const valarray<T>& w, const valarray<T>& b, valarray<T>& x) const;
   bool solve_sv (const valarray<T>& v, valarray<T>& s) const;

   void qrd (matrix<T>& r);
   void qrbksb (const matrix<T>& r, const valarray<T>& v, valarray<T>& s) const;
   bool solve_qr (const valarray<T>& v, valarray<T>& s) const;

   bool chold ();
   void cholbksb (const valarray<T>& v, valarray<T>& s) const;
   bool solve_chol (const valarray<T>& v, valarray<T>& s) const;

   bool eigen (valarray<T>& eival) const;
   bool eigen (valarray<T>& eival, matrix<T>& eivec) const;

   bool eigen (valarray<T>& rev, valarray<T>& iev) const;
   bool eigen (valarray<T>& rev, valarray<T>& iev, matrix<T>& eivec) const;

   // Type of matrices
   bool isSquare () const throw() { return (mPtr->nrow == mPtr->ncol); }
   bool isSingular () const;
   bool isDiagonal () const;
   bool isScalar () const;
   bool isUnit () const throw();
   bool isNull () const throw();
   bool isSymmetric () const;
   bool isSkewSymmetric () const;
   bool isUpperTriangular () const;
   bool isLowerTriangular () const;
   bool isRowOrthogonal () const;
   bool isColOrthogonal () const;

private:

   struct base_mat
   {
      valarray<T> val;
      size_t nrow, ncol;
      int	 refcnt;

      base_mat (size_t row, size_t col)
            : nrow( row), ncol( col), val( row * col) { refcnt = 1; }

      base_mat (const T& v, size_t row, size_t col)
            : nrow( row), ncol( col), val( v, row * col) { refcnt = 1; }

      base_mat (const T* v, size_t row, size_t col)
            : nrow( row), ncol( col), val( v, row * col) { refcnt = 1; }

      base_mat (size_t row, size_t col, const valarray<T>& v)
            : nrow( row), ncol( col), val( v) { refcnt = 1; }

      ~base_mat () {;}
   };

   base_mat *mPtr;
   enum AllocType { Allocate, Free, Record };

   bool shared () { return (mPtr->refcnt > 1); }
   void clone ();
   base_mat *allocator (AllocType bOpt, size_t nrow = 0, size_t ncol = 0);

   void tred2 (valarray<T>& d, valarray<T>& e, bool eivec);
   bool tql2 (valarray<T>& d, valarray<T>& e, bool eivec);
   void balanc (matrix<T>& v, bool eivec);
   bool hqr2 (valarray<T>& d, valarray<T>& e, matrix<T>& v, bool eivec);

   friend class Mat_iter<T>;
   friend class Cmat_iter<T>;
   friend class Mref<T>;
};


/////////////
// Non-member binary operators
//
template <class T> matrix<T> operator+ (const matrix<T>& m1, const matrix<T>& m2);
template <class T> matrix<T> operator- (const matrix<T>& m1, const matrix<T>& m2);

template <class T> matrix<T> operator* (const matrix<T>& m1, const matrix<T>& m2);
template <class T> matrix<T> operator* (const matrix<T>& m, const T& el);
template <class T> matrix<T> operator* (const T& el, const matrix<T>& m);

template <class T> matrix<T> operator/ (const matrix<T>& m1, const matrix<T>& m2);
template <class T> matrix<T> operator/ (const matrix<T>& m, const T& el);
template <class T> matrix<T> operator/ (const T& el, const matrix<T>& m);

template <class T> valarray<T> operator* (const matrix<T>& m, const valarray<T>& v);
template <class T> valarray<T> operator* (const valarray<T>& v, const matrix<T>& m);
template <class T> valarray<T> operator/ (const matrix<T>& m, const valarray<T>& v);
template <class T> valarray<T> operator/ (const valarray<T>& v, const matrix<T>& m);

template <class T> bool operator== (const matrix<T>& m1, const matrix<T>& m2);
template <class T> bool operator!= (const matrix<T>& m1, const matrix<T>& m2);


/////////////
// Non-member functions
//
template <class T> matrix<T> pow (const matrix<T>& m, size_t n);
template <class T> matrix<T> abs (const matrix<T>& m);
template <class T> matrix<T> floor (const matrix<T>& m);
template <class T> matrix<T> ceil (const matrix<T>& m);

template <class T> void mswap (const Mref<T>& x, const Mref<T>& y);
template <class T> void mswap (const Mat_iter<T>& x, const Mat_iter<T>& y);
template <class T> void mswap (matrix<T>& x, matrix<T>& y);

template <class T> istream& operator>> (istream& is, Mref<T> el);
template <class T> istream& operator>> (istream& is, valarray<T>& v);
template <class T> istream& operator>> (istream& is, Mat_iter<T> v);
template <class T> istream& operator>> (istream& is, matrix<T>& m);

template <class T> ostream& operator<< (ostream& os, const Mref<T>& el);
template <class T> ostream& operator<< (ostream &os, const valarray<T>& v);
template <class T> ostream& operator<< (ostream& os, const Mat_iter<T>& v);
template <class T> ostream& operator<< (ostream& os, const Cmat_iter<T>& v);
template <class T> ostream& operator<< (ostream &os, const matrix<T>& m);

float epsilon (complex<float> v);
double epsilon (complex<double> v);
long double epsilon (complex<long double> v);

template <class T> T epsilon (const T& v);
template <class T> bool isVecEq (const valarray<T>& v1, const valarray<T>& v2);


//////////////////
// Temporary helping classes
//

template <class T>
class Mref
{
public:
   Mref (matrix<T>& mm, size_t ii) : m(mm), i(ii) {;}

   operator T () const { return m.mPtr->val[i]; }
   T* operator& () const { if (m.shared()) m.clone(); return &m.mPtr->val[i]; }

   T& operator= (const Mref<T>& e) const { if (m.shared()) m.clone(); T te = e; return m.mPtr->val[i] = te; }
   T& operator= (const T& e) const { if (m.shared()) m.clone(); return m.mPtr->val[i] = e; }
   T& operator+= (const T& e) const { if (m.shared()) m.clone(); return m.mPtr->val[i] += e; }
   T& operator-= (const T& e) const { if (m.shared()) m.clone(); return m.mPtr->val[i] -= e; }
   T& operator*= (const T& e) const { if (m.shared()) m.clone(); return m.mPtr->val[i] *= e; }
   T& operator/= (const T& e) const { if (m.shared()) m.clone(); return m.mPtr->val[i] /= e; }

#ifndef _TS_NO_MEMBER_TEMPLATES
   template <class X> T& operator= (const X& e) const { if (m.shared()) m.clone(); return m.mPtr->val[i] = e; }
   template <class X> T& operator+= (const X& e) const { if (m.shared()) m.clone(); return m.mPtr->val[i] += e; }
   template <class X> T& operator-= (const X& e) const { if (m.shared()) m.clone(); return m.mPtr->val[i] -= e; }
   template <class X> T& operator*= (const X& e) const { if (m.shared()) m.clone(); return m.mPtr->val[i] *= e; }
   template <class X> T& operator/= (const X& e) const { if (m.shared()) m.clone(); return m.mPtr->val[i] /= e; }

   template <class X> T operator+ (const X& e) { T y = *this; return y += e;  }
   template <class X> T operator- (const X& e) { T y = *this; return y -= e;  }
   template <class X> T operator* (const X& e) { T y = *this; return y *= e;  }
   template <class X> T operator/ (const X& e) { T y = *this; return y /= e;  }
#endif

   T operator++ () { if (m.shared()) m.clone(); return ++m.mPtr->val[i]; }
   T operator++ (int) { if (m.shared()) m.clone(); return m.mPtr->val[i]++; }
   T operator-- () { if (m.shared()) m.clone(); return --m.mPtr->val[i]; }
   T operator-- (int) { if (m.shared()) m.clone(); return m.mPtr->val[i]--; }

   T operator+ () { return m.mPtr->val[i]; }
   T operator- () { return -m.mPtr->val[i]; }

   Mref (const Mref& mr) : m(mr.m), i(mr.i) {;}
private:
   Mref ();
   matrix<T>& m;
   size_t i;
};


template <class T>
class Cmat_iter
{
public:
   Cmat_iter (const matrix<T>& mm, const slice& ss) : m(mm), s(ss) {;}

   const T& operator[] (size_t i) const;
   operator valarray<T> () const { return m.mPtr->val[s]; }
   size_t size () const { return s.size(); }

   Cmat_iter(const Cmat_iter& mi) : m(mi.m), s(mi.s) {;}
private:
   Cmat_iter();
   Cmat_iter& operator= (const Cmat_iter&);

   const matrix<T>& m;
   slice s;
};


template <class T>
class Mat_iter
{
public:
   Mat_iter (matrix<T>& mm, const slice& ss) : m(mm), s(ss) {;}

   Mref<T> operator[] (size_t i) const;
   operator valarray<T>() const { return m.mPtr->val[s]; }

   void operator= (const valarray<T>& v) const;
   void operator= (const T& e) const;
   void operator*= (const valarray<T>& v) const;
   void operator/= (const valarray<T>& v) const;
   void operator+= (const valarray<T>& v) const;
   void operator-= (const valarray<T>& v) const;
   size_t size () const { return s.size(); }

   Mat_iter(const Mat_iter& mi) : m(mi.m), s(mi.s) {;}
private:
   Mat_iter();
   Mat_iter<T>& operator= (const Mat_iter<T>& mat);

   matrix<T>& m;
   slice s;
};

class mslice
{
public:
   mslice (size_t sRow, size_t sCol, size_t nRow, size_t nCol)
      : srow_(sRow), scol_(sCol), nrow_(nRow), ncol_(nCol) {;}
   mslice (const mslice& ms)
      : srow_(ms.srow_), scol_(ms.scol_),
        nrow_(ms.nrow_), ncol_(ms.ncol_) {;}

   size_t size () const { return nrow_ * ncol_; }
   size_t start (size_t colno) const { return srow_*colno+scol_; }
   size_t end (size_t colno) const { return (srow_+nrow_-1)*colno+(scol_+ncol_); }

   size_t colsize () const { return nrow_; }
   size_t rowsize () const { return ncol_; }

   size_t rowno () const { return nrow_; }
   size_t colno () const { return ncol_; }

private:
   mslice ();
   size_t srow_,scol_, nrow_, ncol_;
};

template <class T>
class mslice_matrix
{
public:
   typedef T value_type;
   mslice_matrix (matrix<T>& mm, const mslice& ss)
      : s(ss), m(mm) {;}

   const mslice& get_slice() const { return s; }
   matrix<T>& get_ref_mat () const { return m; }

   void operator= (const matrix<T>& m2) const;
   void operator= (const T& e) const;
   void operator+= (const matrix<T>& m2) const;
   void operator-= (const matrix<T>& m2) const;
   void operator*= (const T& e) const;
   void operator/= (const T& e) const;

   mslice_matrix (const mslice_matrix& ms) : s(ms.s), m(ms.m) {;}
private:
   mslice_matrix ();
   mslice_matrix& operator= (const mslice_matrix&);

   mslice s;
   matrix<T>& m;
};


enum MATIRX_TYPE { GENERAL, DIAGONAL, TRIDIAGONAL, UTRIANG, LTRIANG };

class gmslice
{
public:
   gmslice (MATIRX_TYPE mtype) 
      : MType_(mtype), ncol_(0) {;}
   gmslice (size_t nCol, const valarray<size_t>& sRow, const valarray<size_t>& nRow)
      : MType_(GENERAL), ncol_(nCol), srow_(sRow), rsize_(nRow) {;}
   gmslice (const gmslice& gm)
      : MType_(gm.MType_), ncol_(gm.ncol_),
        srow_(gm.srow_), rsize_(gm.rsize_) {;}

   size_t size () const { return srow_.size() * ncol_; }
   size_t start (size_t rowno) const { return srow_[rowno]; }
   size_t end (size_t rowno) const { return srow_[rowno]+rsize_[rowno]; }

   size_t colsize () const { return srow_.size(); }
   size_t rowsize () const { return ncol_; }

   size_t rowno () const { return srow_.size(); }
   size_t colno () const { return ncol_; }

   MATIRX_TYPE getype() { return MType_; }

   size_t ncol_;
   valarray<size_t> srow_,rsize_;

private:
   gmslice();
   MATIRX_TYPE MType_;
};


template <class T>
class gmslice_matrix
{
public:
   typedef T value_type;
   gmslice_matrix (matrix<T>& mm, const gmslice& ss)
      : s(ss), m(mm) {;}

   const gmslice& get_slice() const { return s; }
   matrix<T>& get_ref_mat () const { return m; }

   void operator= (const matrix<T>& m2) const;
   void operator= (const T& e) const;
   void operator+= (const matrix<T>& m2) const;
   void operator-= (const matrix<T>& m2) const;
   void operator*= (const T& e) const;
   void operator/= (const T& e) const;

   gmslice_matrix (const gmslice_matrix& ms) : s(ms.s), m(ms.m) {;}
private:
   gmslice_matrix ();
   gmslice_matrix& operator= (const gmslice_matrix&);

   gmslice s;
   matrix<T>& m;
};


}   // namespace techsoft

#include "matrix.cc"

#endif // _MATRIX_H
