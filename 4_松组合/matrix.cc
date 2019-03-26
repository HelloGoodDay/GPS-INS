//////////////////////////////////////////////////////
// Matrix TCL Pro 2.2
// Copyright (c) 2000-2011 Techsoft. All Rights Reserved.
// File name: Matrix.cc - Implementation of matrix class (do NOT include this file)
// Web: http://www.techsoftpl.com/matrix/
//
#ifndef _MATRIX_CC_12a08ba0_8940_11d3_910c_002018382069_INCLUDED
#define _MATRIX_CC_12a08ba0_8940_11d3_910c_002018382069_INCLUDED

#include"cmatrix"

namespace techsoft {

/////////////////////////
// Overloaded versions of "abs" functions
//
#if (defined(__BORLANDC__) && (__BORLANDC__ < 0x0560)) || (defined(_MSC_VER) && (_MSC_VER < 1310))

inline float abs (float v) { return (float)fabs( v); }
inline double abs (double v) { return fabs( v); }
inline long double abs (long double v) { return fabsl( v); }

inline float sqrt (float v) { return (float)sqrt( (double)v); }
inline long double sqrt (long double v) { return sqrtl( v); }

#endif //

template<class T> inline
const T& max (const T& x, const T& y)
{
   return (x < y) ? y : x;
}

template<class T> inline
const T& min (const T& x, const T& y)
{
   return (y < x) ? y : x;
}

template <class T> inline 
void mswap (const Mref<T>& x, const Mref<T>& y)
{
   T z = x;
   x = y;
   y = z;
}

template <class T> inline 
void mswap (const Mat_iter<T>& x, const Mat_iter<T>& y)
{
   size_t n = min( x.size(), y.size());
   for (size_t i=0; i < n; i++)
   {
      T z = x[i];
      x[i] = y[i];
      y[i] = z;
   }
}

template <class T> inline 
void mswap (matrix<T>& x, matrix<T>& y)
{
   matrix<T> z = x;
   x = y;
   y = z;
}

inline float epsilon (complex<float>) 
{
   return std::numeric_limits<float>::epsilon() * 100.0f;
}

inline double epsilon (complex<double>) 
{
   return std::numeric_limits<double>::epsilon() * 1e4;
}

inline long double epsilon (complex<long double>) 
{
   return std::numeric_limits<long double>::epsilon() * 1e4;
}

inline double epsilon (const double& v) 
{ 
   double ep = std::numeric_limits<double>::epsilon() * 1e4;
   return v > 1.0 ? v * ep : ep;
}

inline long double epsilon (const long double& v) 
{ 
   long double ep = std::numeric_limits<long double>::epsilon() * 1e4;
   return v > 1.0 ? v * ep : ep;
}

template <class T> inline 
T epsilon (const T& v) 
{ 
   T e = std::numeric_limits<T>::epsilon() * 100;
   return v > T(1) ? v * e : e;
}

//////////////////////////////////////////////////////////////////
// Implementation of Cmat_iter class
//

template <class T> inline const T&
Cmat_iter<T>::operator[] (size_t i) const
{
   #ifdef RANGE_CHECK_
   if (i >= s.size())
      throw out_of_range( "matrix<T>::operator[][]: Subscript out of range!");
   #endif
   return m.mPtr->val[s.start()+i*s.stride()];
}


//////////////////////////////////////////////////////////////////
// Implementation of Mat_iter class
//

template <class T> inline Mref<T>
Mat_iter<T>::operator[] (size_t i) const
{
   #ifdef RANGE_CHECK_
   if (i >= s.size())
      throw out_of_range( "matrix<T>::operator[][]: Subscript out of range!");
   #endif
   return Mref<T>( m, s.start()+i*s.stride());
}

template <class T> inline void
Mat_iter<T>::operator= (const valarray<T>& v) const
{
   #ifdef RANGE_CHECK_
   if (size() != v.size())
      throw invalid_argument( "matrix<T>::operator[]: size mismatch in assignment operation!");
   #endif

   if (m.shared())
      m.clone();
   m.mPtr->val[s] = v;
}

template <class T> inline void
Mat_iter<T>::operator= (const T& e) const
{
   if (m.shared())
      m.clone();
   m.mPtr->val[s] = e;
}

template <class T> inline void
Mat_iter<T>::operator*= (const valarray<T>& v) const
{
   #ifdef RANGE_CHECK_
   if (size() != v.size())
      throw invalid_argument( "matrix<T>::operator*=: size mismatch!");
   #endif

   if (m.shared())
      m.clone();
   m.mPtr->val[s] *= v;
}

template <class T> inline void
Mat_iter<T>::operator/= (const valarray<T>& v) const
{
   #ifdef RANGE_CHECK_
   if (size() != v.size())
      throw invalid_argument( "matrix<T>::operator/=: size mismatch!");
   #endif

   if (m.shared())
      m.clone();
   m.mPtr->val[s] /= v;
}

template <class T> inline void
Mat_iter<T>::operator+= (const valarray<T>& v) const
{
   #ifdef RANGE_CHECK_
   if (size() != v.size())
      throw invalid_argument( "matrix<T>::operator+=: size mismatch!");
   #endif

   if (m.shared())
      m.clone();
   m.mPtr->val[s] += v;
}

template <class T> inline void
Mat_iter<T>::operator-= (const valarray<T>& v) const
{
   #ifdef RANGE_CHECK_
   if (size() != v.size())
      throw invalid_argument( "matrix<T>::operator-=: size mismatch!");
   #endif

   if (m.shared())
      m.clone();
   m.mPtr->val[s] -= v;
}


//////////////////////////////////////////////////////////////////
// Implementation of mslice_matrix class
//

template <class T>
void mslice_matrix<T>::operator= (const matrix<T>& m2) const
{
   #ifdef RANGE_CHECK_
   if (s.rowno() != m2.rowno() || s.colno() != m2.colno())
      throw out_of_range( "mslice_matrix<T>operator= Non-conforming matrix size!");
   #endif

   T *pv = &m(0,0);
   size_t start = s.start( m.colno());

   for (size_t i=0; i < s.rowno(); i++)
   {
      size_t ipos = m.colno() * i;
      for (size_t j=0; j < s.colno(); j++)
         pv[start+ipos+j] = m2(i,j);
   }
}

template <class T>
void mslice_matrix<T>::operator+= (const matrix<T>& m2) const
{
   #ifdef RANGE_CHECK_
   if (s.rowno() != m2.rowno() || s.colno() != m2.colno())
      throw out_of_range( "mslice_matrix<T>::operator+= Non-conforming matrix size!");
   #endif

   T *pv = &m(0,0);
   size_t start = s.start( m.colno());

   for (size_t i=0; i < s.rowno(); i++)
   {
      size_t ipos = m.colno() * i;
      for (size_t j=0; j < s.colno(); j++)
         pv[start+ipos+j] += m2(i,j);
   }
}

template <class T>
void mslice_matrix<T>::operator-= (const matrix<T>& m2) const
{
   #ifdef RANGE_CHECK_
   if (s.rowno() != m2.rowno() || s.colno() != m2.colno())
      throw out_of_range( "mslice_matrix<T>::operator+= Non-conforming matrix size!");
   #endif

   T *pv = &m(0,0);
   size_t start = s.start( m.colno());

   for (size_t i=0; i < s.rowno(); i++)
   {
      size_t ipos = m.colno() * i;
      for (size_t j=0; j < s.colno(); j++)
         pv[start+ipos+j] -= m2(i,j);
   }
}

template <class T>
void mslice_matrix<T>::operator= (const T& e) const
{
   T *pv = &m(0,0);
   size_t start = s.start( m.colno());

   for (size_t i=0; i < s.rowno(); i++)
   {
      size_t ipos = m.colno() * i;
      for (size_t j=0; j < s.colno(); j++)
         pv[start+ipos+j] = e;
   }
}

template <class T>
void mslice_matrix<T>::operator*= (const T& e) const
{
   T *pv = &m(0,0);
   size_t start = s.start( m.colno());

   for (size_t i=0; i < s.rowno(); i++)
   {
      size_t ipos = m.colno() * i;
      for (size_t j=0; j < s.colno(); j++)
         pv[start+ipos+j] *= e;
   }
}

template <class T>
void mslice_matrix<T>::operator/= (const T& e) const
{
   T *pv = &m(0,0);
   size_t start = s.start( m.colno());

   for (size_t i=0; i < s.rowno(); i++)
   {
      size_t ipos = m.colno() * i;
      for (size_t j=0; j < s.colno(); j++)
         pv[start+ipos+j] /= e;
   }
}


//////////////////////////////////////////////////////////////////
// Implementation of mslice_matrix class
//

template <class T>
void gmslice_matrix<T>::operator= (const matrix<T>& m2) const
{
   #ifdef RANGE_CHECK_
   if (s.rowno() != m2.rowno() || s.colno() != m2.colno())
      throw out_of_range( "mslice_matrix<T>operator= Non-conforming matrix size!");
   #endif

   T *pv = &m(0,0);
   const T *pm = &m2(0,0);

   for (size_t i=0,ipos=0; i < s.rowno(); ++i,ipos+=s.colno())
      for (size_t j=s.start(i); j < s.end(i); j++)
      {
         size_t jpos = ipos+j;
         pv[jpos] = pm[jpos];
      }
}

template <class T>
void gmslice_matrix<T>::operator+= (const matrix<T>& m2) const
{
   #ifdef RANGE_CHECK_
   if (s.rowno() != m2.rowno() || s.colno() != m2.colno())
      throw out_of_range( "mslice_matrix<T>::operator+= Non-conforming matrix size!");
   #endif

   T *pv = &m(0,0);
   const T *pm = &m2(0,0);

   for (size_t i=0,ipos=0; i < s.rowno(); ++i,ipos+=s.colno())
      for (size_t j=s.start(i); j < s.end(i); j++)
      {
         size_t jpos = ipos+j;
         pv[jpos] += pm[jpos];
      }
}

template <class T>
void gmslice_matrix<T>::operator-= (const matrix<T>& m2) const
{
   #ifdef RANGE_CHECK_
   if (s.rowno() != m2.rowno() || s.colno() != m2.colno())
      throw out_of_range( "mslice_matrix<T>::operator+= Non-conforming matrix size!");
   #endif

   T *pv = &m(0,0);
   const T *pm = &m2(0,0);

   for (size_t i=0,ipos=0; i < s.rowno(); ++i,ipos+=s.colno())
      for (size_t j=s.start(i); j < s.end(i); j++)
      {
         size_t jpos = ipos+j;
         pv[jpos] -= pm[jpos];
      }
}

template <class T>
void gmslice_matrix<T>::operator= (const T& e) const
{
   T *pv = &m(0,0);

   for (size_t i=0,ipos=0; i < s.rowno(); ++i,ipos+=s.colno())
      for (size_t j=s.start(i); j < s.end(i); j++)
         pv[ipos+j] = e;
}

template <class T>
void gmslice_matrix<T>::operator*= (const T& e) const
{
   T *pv = &m(0,0);

   for (size_t i=0,ipos=0; i < s.rowno(); ++i,ipos+=s.colno())
      for (size_t j=s.start(i); j < s.end(i); j++)
         pv[ipos+j] *= e;
}

template <class T>
void gmslice_matrix<T>::operator/= (const T& e) const
{
   T *pv = &m(0,0);

   for (size_t i=0,ipos=0; i < s.rowno(); ++i,ipos+=s.colno())
      for (size_t j=s.start(i); j < s.end(i); j++)
         pv[ipos+j] /= e;
}


//////////////////////////////////////////////////////////////////
// Implementation of matrix class
//

template <class T> inline
matrix<T>::matrix ()
{
    mPtr = allocator( Allocate, 0, 0);
}

template <class T> inline
matrix<T>::matrix (size_t nRow, size_t nCol)
{
    mPtr = allocator( Allocate, nRow, nCol);
}

template <class T>
matrix<T>::matrix (size_t nRow, size_t nCol, const T& e)
{
    mPtr = new base_mat( e, nRow, nCol);
    allocator( Record);
}

template <class T>
matrix<T>::matrix (size_t nRow, size_t nCol, const valarray<T>& vec)
{
   mPtr = new base_mat( nRow, nCol, vec);
   allocator( Record);
}

template <class T>
matrix<T>::matrix (size_t row, size_t col, const T* val, MatArrayType aType)
{
    if (aType == C_ARRAY)
    {
        mPtr = new base_mat( val, row, col);
        allocator( Record);
    }
    else
    {
        mPtr = allocator( Allocate, row, col);
        T *pv = &mPtr->val[0];

        for (size_t k=0, j=0; j < col; j++)
            for (size_t i=0; i < row; i++)
                pv[i*col+j] = val[k++];
    }
}

template <class T> inline
matrix<T>::matrix (const matrix<T>& m)
{
   mPtr = m.mPtr;
   mPtr->refcnt++;
}

template <class T>
matrix<T>::matrix (const mslice_matrix<T>& sm)
{
   const mslice& s = sm.get_slice();
   matrix<T>& m =  sm.get_ref_mat();
   mPtr = allocator( Allocate, s.rowno(), s.colno());

   for (size_t i=0; i < s.rowno(); i++)
      for (size_t j=0; j < s.colno(); j++)
         mPtr->val[i*s.colno()+j] = m.mPtr->val[s.start(m.colno())+i*m.colno()+j];
}

template <class T>
matrix<T>::matrix (const gmslice_matrix<T>& sm)
{
   const gmslice& s = sm.get_slice();
   const matrix<T>& m =  sm.get_ref_mat();
   mPtr = allocator( Allocate, s.rowno(), s.colno());
   mPtr->val = T(0);
   T *pv = &mPtr->val[0];
   const T *pm = &m.mPtr->val[0];

   for (size_t i=0,ipos=0; i < s.rowno(); ++i,ipos+=s.colno())
      for (size_t j=s.start(i); j < s.end(i); j++)
      {
         size_t jpos = ipos+j;
         pv[jpos] = pm[jpos];
      }
}

#ifndef _TS_NO_MEMBER_TEMPLATES
template <class T> template <class X>
matrix<T>::matrix (const matrix<X>& m)
{
   mPtr = allocator( Allocate, m.rowno(), m.colno());
   for (size_t i=0; i < rowno(); i++)
      for (size_t j=0; j < colno(); j++)
      {
         X x = m(i,j);
         mPtr->val[i*colno()+j] = x;
      }
}
#endif

// Internal copy constructor
template <class T>
void matrix<T>::clone ()
{
   if (mPtr->refcnt == 1)
      return;
   base_mat *mp = allocator( Allocate, mPtr->nrow, mPtr->ncol);
   mp->val = mPtr->val;
   mPtr->refcnt--;
   mPtr = mp;
}

// Internal matrix allocator/deallocator
template <class T>
#if defined(__MWERKS__) || defined(_MSC_VER) || defined(__GNUC__)
typename
#endif
matrix<T>::base_mat*
matrix<T>::allocator (AllocType bOpt, size_t nrow, size_t ncol)
{
   static base_mat *mlist[MAX_MATRIX_CACHE+1];
   static int instCount = 0;
   static volatile int semCount = 0;
   static bool first = true;

   if (first)
   {
      first = false;
      for (int k=0; k <= MAX_MATRIX_CACHE; k++)
         mlist[k] = NULL;   
   }

   switch (bOpt)
   {
   case Record:
      instCount++;
      return NULL;

   case Allocate:
      instCount++;
      base_mat *mp;
      if (++semCount == 1)
         for (int i=0; i < MAX_MATRIX_CACHE; i++)
            if (mlist[i] != NULL && mlist[i]->nrow == nrow && mlist[i]->ncol == ncol)
            {
               mp = mlist[i];
               mlist[i] = NULL;
               --semCount;
               mp->refcnt = 1;
               return mp;
            }
      --semCount;
      mp = new base_mat( nrow, ncol);
      return mp;

   case Free:
      if (--instCount)
      {
         if (++semCount == 1)
            for (int i=0; i < MAX_MATRIX_CACHE; i++)
            {
               swap( mlist[i], mPtr);
               if (mPtr == NULL)
               {
                  --semCount;
                  return NULL;
               }
            }
         --semCount;
      }
      else                               // clear the cached matrices
      {
         ++semCount;
         for (int i=0; i < MAX_MATRIX_CACHE; i++)
            if (mlist[i] != NULL)
            {
               delete mlist[i];
               mlist[i] = NULL;
            }
         --semCount;
      }
      delete mPtr;
      mPtr = NULL;
      break;
   }
   return NULL;
}

// destructor
template <class T> inline
matrix<T>::~matrix ()
{
   if (--mPtr->refcnt == 0)
     allocator( Free);
}

template <class T> void
matrix<T>::resize (size_t nRow, size_t nCol, const T& dval)
{
   if (mPtr->nrow == nRow && mPtr->ncol == nCol)
      return;

   base_mat *mp = allocator( Allocate, nRow, nCol);

   for (size_t i=0; i < nRow; i++)
      for (size_t j=0; j < nCol; j++)
         mp->val[i*nCol+j] = (i >= mPtr->nrow || j >= mPtr->ncol) ?
                                  dval : mPtr->val[i*mPtr->ncol+j];
   if (--mPtr->refcnt == 0)
      allocator( Free);
   mPtr = mp;
}

template <class T> void
matrix<T>::free ()
{
   if (--mPtr->refcnt == 0)
      allocator( Free);
   mPtr = allocator( Allocate, 0, 0);
}

template <class T> inline mslice_matrix<T>
matrix<T>::operator[] (mslice ms)
{
   #ifdef RANGE_CHECK_
   if (ms.end( colno()-1) > size())
       throw invalid_argument( "matrix<T>::operator[]: sub-matrix is out of range!");
   #endif
   return mslice_matrix<T>( *this, ms);
}

template <class T> inline gmslice_matrix<T>
matrix<T>::operator[] (gmslice gm)
{
   #ifdef RANGE_CHECK_
   if (gm.size() > size())
       throw invalid_argument( "matrix<T>::operator[]: sub-matrix is out of range!");
   #endif
   
   if (gm.getype() != GENERAL)
   {
      gm.ncol_ = mPtr->ncol;
      gm.srow_.resize( mPtr->nrow);
      gm.rsize_.resize( mPtr->nrow);
   }
   size_t i;
   switch(gm.getype())
   {
   case DIAGONAL:
      gm.rsize_ = 1;
      for (i=0; i < gm.ncol_; ++i)
         gm.srow_[i] = i;
      break;

   case TRIDIAGONAL:
      gm.rsize_ = 3;
      gm.rsize_[0] = 2;
      gm.rsize_[gm.ncol_-1] = 2;
      gm.srow_[0] = 0;
      for (i=1; i < gm.ncol_; ++i)
         gm.srow_[i] = i-1;
      break;
      
   case UTRIANG:
      for (i=0; i < gm.ncol_; ++i)
      {
         gm.srow_[i] = i;
         gm.rsize_[i] = gm.ncol_-i;
      }
      break;
      
   case LTRIANG:
      gm.srow_ = 0;
      for (i=0; i < gm.ncol_; ++i)
         gm.rsize_[i] = i+1;
      break;
   }
   return gmslice_matrix<T>( *this, gm);
}

template <class T> inline matrix<T>&
matrix<T>::operator= (const matrix<T>& m)
{
   m.mPtr->refcnt++;
   if (--mPtr->refcnt == 0)
      allocator( Free);
   mPtr = m.mPtr;
   return *this;
}

template <class T> matrix<T>&
matrix<T>::operator= (const mslice_matrix<T>& sm)
{
   const mslice& s = sm.get_slice();
   matrix<T>& m =  sm.get_ref_mat();
   base_mat *mp = allocator( Allocate, s.rowno(), s.colno());
   T *pv = &mp->val[0];
   const T *pm = &m.mPtr->val[0];

   for (size_t i=0; i < s.rowno(); i++)
      for (size_t j=0; j < s.colno(); j++)
         pv[i*s.colno()+j] = pm[s.start(m.colno())+i*m.colno()+j];

   if (--mPtr->refcnt == 0)
      allocator( Free);
   mPtr = mp;

   return *this;
}

template <class T> matrix<T>&
matrix<T>::operator= (const gmslice_matrix<T>& sm)
{
   const gmslice& s = sm.get_slice();
   const matrix<T>& m =  sm.get_ref_mat();
   base_mat *mp = allocator( Allocate, s.rowno(), s.colno());
   mp->val = T(0);
   T *pv = &mp->val[0];
   const T *pm = &m(0,0);

   for (size_t i=0,ipos=0; i < s.rowno(); ++i,ipos+=s.colno())
      for (size_t j=s.start(i); j < s.end(i); j++)
      {
         size_t jpos = ipos+j;
         pv[jpos] = pm[jpos];
      }

   if (--mPtr->refcnt == 0)
      allocator( Free);
   mPtr = mp;

   return *this;
}

template <class T>
matrix<T>& matrix<T>::operator= (const T& el)
{
   if (shared())
   {
      --mPtr->refcnt;
      mPtr = allocator( Allocate, mPtr->nrow, mPtr->ncol);
   }
   mPtr->val = el;

   return *this;
}

#ifndef _TS_NO_MEMBER_TEMPLATES
template <class T> template <class X>
matrix<T>& matrix<T>::operator= (const matrix<X>& m)
{
   if (--mPtr->refcnt == 0)
      allocator( Free);

   mPtr = allocator( Allocate, m.rowno(), m.colno());
   const X *pm = &m(0,0);
   T *pv = &mPtr->val[0];
   size_t n = size();

   for (size_t i=0; i < n; i++)
      pv[i] = pm[i];

   return *this;
}
#endif

template <class T> const matrix<T>
matrix<T>::operator[] (mslice ms) const
{
   #ifdef RANGE_CHECK_
   if (ms.end( colno()-1) > size())
       throw invalid_argument( "matrix<T>::operator[]: sub-matrix is out of range!");
   #endif

   matrix<T> tm(ms.rowno(),ms.colno());
   for (size_t i=0; i < tm.rowno(); i++)
      for (size_t j=0; j < tm.colno(); j++)
         tm.mPtr->val[i*tm.colno()+j] = mPtr->val[ms.start(colno())+i*colno()+j];
   return tm;
}

template <class T> const matrix<T>
matrix<T>::operator[] (gmslice gm) const
{
   #ifdef RANGE_CHECK_
   if (gm.size() > size())
       throw invalid_argument( "matrix<T>::operator[]: sub-matrix is out of range!");
   #endif
   
   if (gm.getype() != GENERAL)
   {
      gm.ncol_ = mPtr->ncol;
      gm.srow_.resize( mPtr->nrow);
      gm.rsize_.resize( mPtr->nrow);
   }
   size_t k;
   switch(gm.getype())
   {
   case DIAGONAL:
      gm.rsize_ = 1;
      for (k=0; k < gm.ncol_; ++k)
         gm.srow_[k] = k;
      break;

   case TRIDIAGONAL:
      gm.rsize_ = 3;
      gm.rsize_[0] = 2;
      gm.rsize_[gm.ncol_-1] = 2;
      gm.srow_[0] = 0;
      for (k=1; k < gm.ncol_; ++k)
         gm.srow_[k] = k-1;
      break;
      
   case UTRIANG:
      for (k=0; k < gm.ncol_; ++k)
      {
         gm.srow_[k] = k;
         gm.rsize_[k] = gm.ncol_-k;
      }
      break;
      
   case LTRIANG:
      gm.srow_ = 0;
      for (k=0; k < gm.ncol_; ++k)
         gm.rsize_[k] = k+1;
      break;
   }

   matrix<T> tm(gm.rowno(), gm.colno(), T(0));
   const T *pv = &mPtr->val[0];
   T *pm = &tm.mPtr->val[0];

   for (size_t i=0,ipos=0; i < gm.rowno(); ++i,ipos+=gm.colno())
      for (size_t j=gm.start(i); j < gm.end(i); j++)
      {
         size_t jpos = ipos+j;
         pm[jpos] = pv[jpos];
      }
   return tm;
}

template <class T> inline Mat_iter<T>
matrix<T>::row (size_t i)
{
   #ifdef RANGE_CHECK_
   if (i >= mPtr->nrow)
      throw out_of_range( "matrix<T>::row: subscript is out of range!");
   #endif
   return Mat_iter<T>( *this, slice( i * mPtr->ncol, mPtr->ncol, 1));
}

template <class T> inline Cmat_iter<T>
matrix<T>::row (size_t i) const
{
   #ifdef RANGE_CHECK_
   if (i >= mPtr->nrow)
      throw out_of_range( "matrix<T>::row: subscript is out of range!");
   #endif
   return Cmat_iter<T>( *this, slice( i * mPtr->ncol, mPtr->ncol, 1));
}

template <class T> inline Mat_iter<T>
matrix<T>::column (size_t i)
{
   #ifdef RANGE_CHECK_
   if (i >= mPtr->ncol)
      throw out_of_range( "matrix<T>::column: subscript is out of range!");
   #endif
   return Mat_iter<T>( *this, slice( i, mPtr->nrow, mPtr->ncol));
}

template <class T> inline Cmat_iter<T>
matrix<T>::column (size_t i) const
{
   #ifdef RANGE_CHECK_
   if (i >= mPtr->ncol)
      throw out_of_range( "matrix<T>::column: subscript is out of range!");
   #endif
   return Cmat_iter<T>( *this, slice( i, mPtr->nrow, mPtr->ncol));
}

template <class T> inline Mat_iter<T>
matrix<T>::diag (int i)
{
   #ifdef RANGE_CHECK_
   if (i >= int(mPtr->ncol) || int(mPtr->nrow)+i < 1)
      throw out_of_range( "matrix<T>::diag: subscript is out of range!");
   #endif
   size_t st,sz;
   if (i >= 0)
   {
      st = i;
      sz = mPtr->ncol - i;
      if (sz > mPtr->nrow)
        sz = mPtr->nrow;
   }
   else
   {
      st = -i * int(mPtr->ncol);
      sz = int(mPtr->nrow) + i;
      if (sz > mPtr->ncol)
         sz = mPtr->ncol;
   }
   return Mat_iter<T>( *this, slice( st, sz, mPtr->ncol+1));
}

template <class T> inline Cmat_iter<T>
matrix<T>::diag (int i) const
{
   #ifdef RANGE_CHECK_
   if (i >= int(mPtr->ncol) || int(mPtr->nrow)+i < 1)
      throw out_of_range( "matrix<T>::diag: subscript is out of range!");
   #endif

   size_t st,sz;
   if (i >= 0)
   {
      st = i;
      sz = mPtr->ncol - i;
      if (sz > mPtr->nrow)
        sz = mPtr->nrow;
   }
   else
   {
      st = -i * int(mPtr->ncol);
      sz = int(mPtr->nrow) + i;
      if (sz > mPtr->ncol)
         sz = mPtr->ncol;
   }
   return Cmat_iter<T>( *this, slice( st, sz, mPtr->ncol+1));
}

template <class T> inline Mref<T>
matrix<T>::operator() (size_t i, size_t j)
{
   #ifdef RANGE_CHECK_
   if (i >= mPtr->nrow || j >= mPtr->ncol)
      throw out_of_range( "matrix<T>::operator(): subscript is out of range!");
   #endif
   return Mref<T>( *this, i * mPtr->ncol + j);
}

template <class T> inline const T&
matrix<T>::operator() (size_t i, size_t j) const
{
   #ifdef RANGE_CHECK_
   if (i >= mPtr->nrow || j >= mPtr->ncol)
      throw out_of_range( "matrix<T>::operator(): subscript is out of range!");
   #endif
   return mPtr->val[i*mPtr->ncol+j];
}

template <class T> matrix<T>&
matrix<T>::operator+= (const matrix<T>& m)
{
   #ifdef RANGE_CHECK_
   if (rowno() != m.rowno() || colno() != m.colno())
      throw invalid_argument( "matrix<T>::operator+=: Inconsistent matrix size in addition!");
   #endif

   base_mat *mp;
   if (shared())
      mp = allocator( Allocate, mPtr->nrow, mPtr->ncol);
   else
      mp = mPtr;

   const T *pm = &m(0,0);
   const T *pv = &mPtr->val[0];
   T *pv2 = &mp->val[0];
   size_t n = size();

   for (size_t i=0; i < n; i++)
      pv2[i] = pv[i] + pm[i];

   if (shared())
   {
      --mPtr->refcnt;
      mPtr = mp;
   }
   return *this;
}

template <class T> matrix<T>&
matrix<T>::operator-= (const matrix<T>& m)
{
   #ifdef RANGE_CHECK_
   if (rowno() != m.rowno() || colno() != m.colno())
      throw invalid_argument( "matrix<T>::operator-=: Inconsistent matrix size in addition!");
   #endif

   base_mat *mp;
   if (shared())
      mp = allocator( Allocate, mPtr->nrow, mPtr->ncol);
   else
      mp = mPtr;

   const T *pm = &m(0,0);
   const T *pv = &mPtr->val[0];
   T *pv2 = &mp->val[0];
   size_t n = size();

   for (size_t i=0; i < n; i++)
      pv2[i] = pv[i] - pm[i];

   if (shared())
   {
      --mPtr->refcnt;
      mPtr = mp;
   }
   return *this;
}

template <class T> matrix<T>&
matrix<T>::operator*= (const T& x)
{
   base_mat *mp;

   if (shared())
      mp = allocator( Allocate, mPtr->nrow, mPtr->ncol);
   else
      mp = mPtr;

   const T *pv = &mPtr->val[0];
   T *pv2 = &mp->val[0];
   size_t n = size();

   for (size_t i=0; i < n; i++)
      pv2[i] = pv[i] * x;

   if (shared())
   {
      --mPtr->refcnt;
      mPtr = mp;
   }
   return *this;
}

template <class T> matrix<T>&
matrix<T>::operator*= (const matrix<T>& m)
{
   size_t nr1 = mPtr->nrow;
   size_t nc1 = mPtr->ncol;
   size_t nc2 = m.mPtr->ncol;

   #ifdef RANGE_CHECK_
   size_t nr2 = m.mPtr->nrow;
   if (nc1 != nr2)
      throw invalid_argument( "matrix<T>::operator*=: Inconsistent matrix size in multiplication!");
   #endif

   matrix<T> tm( nr1, nc2);
   const T *pm1 = &mPtr->val[0];
   const T *pm2 = &m(0,0);
   T *pm = &tm.mPtr->val[0];

   for (size_t i=0; i < nr1; i++)
      for (size_t j=0; j < nc2; j++)
      {
         pm[i*nc2+j] = T(0);
         for (size_t k=0; k < nc1; k++)
            pm[i*nc2+j] += pm1[i*nc1+k] * pm2[k*nc2+j];
      }
   *this = tm;
   return *this;
}

template <class T> inline matrix<T>&
matrix<T>::operator/= (const T& x)
{
   *this *= 1/x;
   return *this;
}

template <class T> matrix<T>&
matrix<T>::operator/= (const matrix<T>& m)
{
   *this *= !m;
   return *this;
}

template <class T> matrix<T> 
matrix<T>::operator- () const
{
   matrix<T> m( mPtr->nrow, mPtr->ncol);
   T *pv2 = &m.mPtr->val[0];
   const T *pv = &mPtr->val[0];
   size_t n = size();

   for (size_t i=0; i < n; i++)
      pv2[i] = -pv[i];

   return m;
}

template <class T>
matrix<T> matrix<T>::operator~ () const
{
   size_t m = mPtr->nrow;
   size_t n = mPtr->ncol;
   matrix<T> tm( n, m);
   
   const T *pm = &mPtr->val[0];
   T *ptm = &tm.mPtr->val[0];

   for (size_t i=0; i < m; i++)
      for (size_t j=0; j < n; j++)
         ptm[j*m+i] = pm[i*n+j];
         
   return tm;
}

template <class T>
inline matrix<T> matrix<T>::operator! () const
{
   matrix<T> m = *this;
   if (!m.inv())
      m.null();
   return m;
}

template <class T> inline matrix<T>
operator+ (const matrix<T>& m1, const matrix<T>& m2)
{
   matrix<T> m = m1;
   return m += m2;
}

template <class T> inline matrix<T>
operator- (const matrix<T>& m1, const matrix<T>& m2)
{
   matrix<T> m = m1;
   return m -= m2;
}

template <class T> inline matrix<T>
operator* (const matrix<T>& m1, const matrix<T>& m2)
{
   matrix<T> m = m1;
   return m *= m2;
}

template <class T> inline matrix<T>
operator* (const matrix<T>& m, const T& el)
{
   matrix<T> m2 = m;
   return m2 *= el;
}

template <class T> inline matrix<T>
operator* (const T& el, const matrix<T>& m)
{
   return m * el;
}

template <class T> valarray<T>
operator* (const matrix<T>& m, const valarray<T>& v)
{
   #ifdef RANGE_CHECK_
   if (m.rowsize() != v.size())
      throw invalid_argument( "matrix<T>::operator*: Inconsistent matrix size!");
   #endif

   valarray<T> vt( m.rowno());
   for (size_t i=0; i < m.rowno(); i++)
   {
      vt[i] = T(0);
      for (size_t j=0; j < m.colno(); j++)
      {
         T el = m(i,j);
         vt[i] += el * v[j];
      }
   }
   return vt;
}

template <class T> valarray<T>
operator* (const valarray<T>& v, const matrix<T>& m)
{
   #ifdef RANGE_CHECK_
   if (v.size() != m.rowno())
      throw invalid_argument( "matrix<T>::operator*: Inconsistent matrix size!");
   #endif

   valarray<T> vt( m.colno());
   for (size_t i=0; i < m.colno(); i++)
   {
      vt[i] = T(0);
      for (size_t j=0; j < v.size(); ++j)
      {
         T el = m(j,i);
         vt[i] += el * v[j];
      }
   }
   return vt;
}

template <class T> inline matrix<T>
operator/ (const matrix<T>& m1, const matrix<T>& m2)
{
   matrix<T> m = m1;
   return m *= !m2;
}

template <class T> inline matrix<T>
operator/ (const matrix<T>& m, const T& el)
{
   matrix<T> m2 = m;
   return m2 *= 1/el;
}

template <class T> inline matrix<T>
operator/ (const T& el, const matrix<T>& m)
{
   matrix<T> m2 = !m;
   return m2 *= el;
}

template <class T> inline valarray<T>
operator/ (const matrix<T>& m, const valarray<T>& v)
{
   valarray<T> v2 = T(1) / v;
   return m * v2;
}

template <class T> inline valarray<T>
operator/ (const valarray<T>& v, const matrix<T>& m)
{
   return v * !m;
}

template <class T> bool
operator== (const matrix<T>& m1, const matrix<T>& m2)
{
   if (m1.rowno() != m2.rowno() || m1.colno() != m2.colno())
      return false;

   size_t n = m1.size();
   const T *pm1 = &m1(0,0);
   const T *pm2 = &m2(0,0);

   for (size_t i=0; i < n; i++)
      if (abs( pm1[i] - pm2[i]) > epsilon( pm1[i]))
         return false;

   return true;
}

template <class T> inline bool
operator!= (const matrix<T>& m1, const matrix<T>& m2)
{
   return !(m1 == m2);
}

template <class T> bool 
isVecEq (const valarray<T>& v1, const valarray<T>& v2)
{
   if (v1.size() != v2.size())
      return false;

   for (size_t i=0; i < v1.size(); i++)
   {
      T d = v1[i] - v2[i];
      if (abs( d) > epsilon( v1[i]))
         return false;
   }
   return true;
}

template <class T> inline void 
matrix<T>::unit ()
{
   null();
   diag(0) = T(1);
}

template <class T> inline void 
matrix<T>::null ()
{
   if (shared())
   {
      --mPtr->refcnt;
      mPtr = allocator( Allocate, mPtr->nrow, mPtr->ncol);
   }
   mPtr->val = T(0);
}

inline int rand_ () { return rand(); }

template <class T> void
matrix<T>::rand (int rmin, int rmax, int rseed)
{
   if (!rmin && !rmax)
   {
      rmin = -1;
      rmax = 1;
   }
   if (rmin > rmax)
   {
      int tmp = rmin;
      rmin = rmax;
      rmax = tmp;
   }
   int rdiff = rmax - rmin;

   if (rseed)
      srand( rseed);

   if (shared())
   {
      --mPtr->refcnt;
      mPtr = allocator( Allocate, mPtr->nrow, mPtr->ncol);
   }
   
   T *pv = &mPtr->val[0];
   size_t n = size();

   for (size_t i=0; i < n; i++)
   {
      int rn = rand_();
      if (rdiff > 1)
      {
         rn %= rdiff;
         pv[i] = T(rn) - T(rdiff) / T(2);
         if (rdiff > 1)
         {
            T f = ((T)rand_()) / T(RAND_MAX);
            if (rn > rdiff)
               pv[i] -= f;
            else
               pv[i] += f;
         }
      }
      else
         pv[i] = T(rn) / T(RAND_MAX);
   }
}

template <class T> matrix<T>
pow (const matrix<T>& m, size_t n)
{
   #ifdef RANGE_CHECK_
   if (m.colno() != m.rowno() || m.rowno() == 0)
      throw invalid_argument( "pow(): Non-square matrix!");
   #endif

   matrix<T> x(m.rowno(), m.colno());
   if (n == 0)
   {
      x.unit();
      return x;
   }
   matrix<T> m1 = m;
   bool first = true;

   while (1)
   {
     if (n & 1)
     {
       if (first)
       {
         x = m1;
         first = false;
       }
       else
         x *= m1;
     }
     n >>= 1;
     if (n == 0)
        break;

     m1 *= m1;
   }
   return x;
}

template <class T> matrix<T>
abs (const matrix<T>& m)
{
   matrix<T> a(m.rowno(), m.colno());
   T *pa = &a(0,0);
   const T *pm = &m(0,0);
   size_t n = m.size();
   for (size_t i=0; i < n; i++)
      pa[i] = T(abs( pm[i]));
   return a;
}

template <class T> matrix<T>
floor (const matrix<T>& m)
{
   matrix<T> a(m.rowno(), m.colno());
   T *pa = &a(0,0);
   const T *pm = &m(0,0);
   size_t n = m.size();
   for (size_t i=0; i < n; i++)
      pa[i] = T(floor( pm[i]));
   return a;
}

template <class T> matrix<T>
ceil (const matrix<T>& m)
{
   matrix<T> a(m.rowno(), m.colno());
   T *pa = &a(0,0);
   const T *pm = &m(0,0);
   size_t n = m.size();
   for (size_t i=0; i < n; i++)
      pa[i] = T(ceil( pm[i]));
   return a;
}

template <class T> inline T
matrix<T>::trace (int i) const
{
   valarray<T> d = diag( i);
   return d.sum();
}

///////////
// Computes (a^2 + b^2 )^1/2 without destructive underflow or overflow.
//
template <class T> inline 
T hypot (T a, T b)
{
   a = abs( a);
   b = abs( b);
   if (a > b)
   {
      T c = b/a;
      return (a * sqrt( T(1) + c * c));
   }
   if (b < epsilon(b))
      return T(0);
   T c = a/b;
   return (b * sqrt( T(1) + c * c));
}

template <class T> inline
T sign (T a, T b)
{
   return (b >= T(0) ? abs(a) : -abs(a));
}

template <class T>
bool matrix<T>::lud (valarray<size_t>& ri, T* pDet)
{
   size_t i,j,k;
   double ta,tb;

   #ifdef RANGE_CHECK_
   if (mPtr->nrow != mPtr->ncol)
      throw invalid_argument( "matrix<T>::lud: Non-square matrix!");
   #endif

   if (mPtr->nrow != ri.size())
      ri.resize( mPtr->nrow);

   if (shared())
      clone();

   if (pDet != NULL)
      *pDet = T(1);

   size_t n = mPtr->nrow;
   T *pv = &mPtr->val[0];

   for (i=0; i < n; i++)
      ri[i] = i;

   for (k=0; k < n-1; k++)
   {
      j = k;
      ta = abs( pv[ri[k]*n+k]);
      for (i=k+1; i < n; i++)
         if ((tb = abs( pv[ri[i]*n+k])) > ta)
         {
            ta = tb;
            j = i;
         }
      if (j != k)
      {
         swap( ri[j], ri[k]);
         if (pDet != NULL)
            *pDet = - *pDet;
      }
      size_t kpos = ri[k] * n;

      if (abs( pv[kpos+k]) < epsilon( pv[kpos+k]))
         return false;

      if (pDet != NULL)
         *pDet *= pv[kpos+k];

      for (i=k+1; i < n; i++)
      {
         size_t ipos = ri[i] * n;
         T a = pv[ipos+k] /= pv[kpos+k];

         for (j=k+1; j < n; j++)
            pv[ipos+j] -= a * pv[kpos+j];
      }
   }
   if (pDet != NULL)
      *pDet *= pv[ri[k]*n+k];
   
   return true;
}

template <class T> void
matrix<T>::lubksb (const valarray<size_t>& ri, const valarray<T>& v, valarray<T>& s) const
{
   size_t i,j,k,ip;
   bool nonzero = false;
   size_t n = mPtr->nrow;
   const T *pv = &mPtr->val[0];

   #ifdef RANGE_CHECK_
   if (n != v.size())
      throw invalid_argument( "matrix<T>::lubksb: Invalid vector size!");
   #endif

   if (n != s.size())
      s.resize( n);

   for (k=0,i=0; i < n; i++)
   {
      ip = ri[i];
      s[i] = v[ip];

      if (nonzero)
      {
         ip *= n;
         for (j=k; j < i; j++)
            s[i] -= pv[ip+j] * s[j];
      }
      else if (s[i] != T(0))
      {
         k = i;
         nonzero = true;
      }
   }
   for (i=n-1; ; i--)
   {
      ip = ri[i] * n;
      for (j=i+1; j < n; j++)
         s[i] -= pv[ip+j] * s[j];
      s[i] /= pv[ip+i];
      if (i == 0) break;
   }
}

template <class T> void 
matrix<T>::lumpove (const matrix<T>& ludm, const valarray<size_t>& ri, const valarray<T>& v, valarray<T>& s) const
{
   size_t i,j;
   size_t n = mPtr->ncol;
   const T *pa = &mPtr->val[0];

   valarray<T> rv(n), rs(n);
   for (i=0; i < n; i++)
   {
      T tmp = -v[i];
      for (j=0; j < n; j++)
         tmp += pa[i*n+j] * s[j];
      rv[i] = tmp;
   }
   ludm.lubksb( ri, rv, rs);
   for (i=0; i < n; i++)
      s[i] -= rs[i];
}

template <class T> bool 
matrix<T>::solve (const valarray<T>& v, valarray<T>& s) const
{
   matrix<T> tm = *this;
   valarray<size_t> ri( mPtr->nrow);

   if (tm.lud( ri))
   {
      tm.lubksb( ri, v, s);
      return true;
   }
   return false;
}

template <class T> bool 
matrix<T>::solve_chol (const valarray<T>& v, valarray<T>& s) const
{
   matrix<T> tm = *this;

   if (tm.chold())
   {
      tm.cholbksb( v, s);
      return true;
   }
   return false;
}

template <class T> bool 
matrix<T>::solve_sv (const valarray<T>& v, valarray<T>& s) const
{
   #ifdef RANGE_CHECK_
   if (mPtr->nrow < mPtr->ncol)
      throw invalid_argument( "matrix<T>::solve_sv(): No. of rows must be equal or greater than column!");
   #endif

   matrix<T> vc( mPtr->ncol, mPtr->ncol);
   matrix<T> tm = *this;
   valarray<T> w( mPtr->ncol);

   if (tm.svd( vc, w))
   {
      tm.svbksb( vc, w, v, s);
      return true;
   }
   return false;
}

template <class T> bool 
matrix<T>::solve_qr (const valarray<T>& v, valarray<T>& s) const
{
   #ifdef RANGE_CHECK_
   if (mPtr->nrow < mPtr->ncol)
      throw invalid_argument( "matrix<T>::solve_qr(): No. of rows must be equal or greater than column!");
   #endif

   matrix<T> r( mPtr->ncol, mPtr->ncol);
   matrix<T> tm = *this;

   tm.qrd( r);
   tm.qrbksb( r, v, s);
   
   return true;
}

template <class T> bool
matrix<T>::solve (const matrix<T>& v, matrix<T>& s) const
{
   size_t n = mPtr->nrow;
   #ifdef RANGE_CHECK_
   if (n != v.rowno())
      throw invalid_argument( "matrix<T>::operator*: Inconsistent number of vectors in solve!");
   #endif

   if (v.size() != s.size())
      s.resize( v.rowno(), v.colno());

   matrix<T> tm = *this;
   valarray<size_t> ri(n);
   valarray<T> vt(n), st(n);

   if (!tm.lud( ri))
      return false;

   for (size_t i=0; i < v.colno(); i++)
   {
      vt = v.column(i);
      tm.lubksb( ri, vt, st);
      s.column(i) = st;
   }
   return true;
}

template <class T> bool 
matrix<T>::chold ()
{
   #ifdef RANGE_CHECK_
   if (!isSymmetric())
      throw invalid_argument( "matrix<T>::chold: Non-symmetric matrix!");
   #endif

   int n = int(mPtr->nrow);
   if (shared())
      clone();
   T *pa = &mPtr->val[0];
   int i,j;

   for (i=0; i < n; i++)
      for (j=i; j < n; j++)
      {
         T s = pa[i*n+j];
         for (int k=i-1; k >= 0; k--)
            s -= pa[i*n+k] * pa[j*n+k];

         if (i == j)
         {
            if (s  <= epsilon(s))
               return false;
            pa[i*n+i] = sqrt( s);
         }
         else
            pa[j*n+i] = s / pa[i*n+i];
      }
	for (j=1; j < n; j++)
		for (i=0; i < j; i++)
			pa[i*n+j] = T(0);
   
   return true;
}

template <class T> void
matrix<T>::cholbksb (const valarray<T>& v, valarray<T>& s) const
{
   int i,k,n;
   T a;

   const T *pm = &mPtr->val[0];
   n = int(mPtr->nrow);
   i = int(s.size());

   if (n != i)
      s.resize( n);

   for (i=0; i < n; i++)
   {
      a = v[i];
      for (k=i-1; k >= 0; k--)
         a -= pm[i*n+k] * s[k];
      s[i] = a / pm[i*n+i];
   }
   for (i=n-1; i >= 0; i--)
   {
      a = s[i];
      for (k=i+1; k < n; k++)
         a -= pm[k*n+i] * s[k];
      s[i] = a / pm[i*n+i];
   }
}

template <class T> void        
matrix<T>::qrd (matrix<T>& r)
{
   size_t i,j,k;
   size_t m = mPtr->nrow;
   size_t n = mPtr->ncol;

   #ifdef RANGE_CHECK_
   if (m < n)
      throw invalid_argument( "matrix<T>::qrd: Inconsistent matrix size!");
   #endif

   if (shared())
      clone();

   if (n != r.rowno() || r.rowno() != r.colno())
      r.resize( n, n);

   T *pa = &mPtr->val[0];
   T *pr = &r(0,0);

   for (k=0; k < n; k++)
   {
      T nrm(0);
      for (i = k; i < m; i++)
         nrm = hypot( nrm, pa[i*n+k]);

      if (nrm > epsilon(nrm))
      {
         if (pa[k*n+k] < T(0))
            nrm = -nrm;

         for (i = k; i < m; i++)
            pa[i*n+k] /= nrm;
         pa[k*n+k] += T(1);

         for (j=k+1; j < n; j++)
         {
            T s(0);
            for (i=k; i < m; i++)
               s += pa[i*n+k] * pa[i*n+j];

            s /= -pa[k*n+k];
            for (i=k; i < m; i++)
               pa[i*n+j] += s * pa[i*n+k];
         }
      }
      pr[k*n+k] = -nrm;
   }

	for (j=1; j < n; j++)
		for (size_t i=0; i < j; i++)
      {
         pr[i*n+j] = pa[i*n+j];
         pr[j*n+i] = T(0);
      }

   matrix<T> q(m,n);
   T *pq = &q(0,0);

   for (k=n-1; ; k--)
   {
      for (i = 0; i < m; i++)
         pq[i*n+k] = T(0);

      pq[k*n+k] = T(1);
      for (j=k; j < n; j++)
      {
         if (pa[k*n+k] != T(0))
         {
            T s(0);
            for (i=k; i < m; i++)
               s += pa[i*n+k] * pq[i*n+j];
            s /= -pa[k*n+k];
            for (i = k; i < m; i++)
               pq[i*n+j] += s * pa[i*n+k];
         }
      }
      if (k == 0)
         break;
   }
   *this = q;
}

template <class T> void
matrix<T>::qrbksb (const matrix<T>& r, const valarray<T>& v, valarray<T>& s) const
{
   size_t m = mPtr->nrow;
   size_t n = mPtr->ncol;
   size_t i,j;

   #ifdef RANGE_CHECK_
   if (m != v.size())
      throw std::invalid_argument( "matrix<T>::qrsolve: Inconsistent vector size!");
   #endif

   if (n != s.size())
      s.resize( n);

   for (i=0; i < n; i++)
   {
      s[i] = T(0);
      for (j=0; j < m; j++)
         s[i] += mPtr->val[j*n+i] * v[j];
   }

   const T *pr = &r(0,0);
   for (i=n-1; ; i--)
   {
      size_t ip = i * n;
      for (j=i+1; j < n; j++)
         s[i] -= pr[ip+j] * s[j];
      s[i] /= pr[ip+i];
      if (i == 0) break;
   }
}

template <class T> bool 
matrix<T>::svd (matrix<T>& vc, valarray<T>& w)
{
   size_t flag,i,its,j,jj,k,l,nm;
   T c,f,h,s,x,y,z,tmp;

   if (shared())
      clone();

   size_t m = mPtr->nrow;
   size_t n = mPtr->ncol;

   if (vc.rowno() != n || vc.colno() != n)
      vc.resize( n, n);
   if (w.size() != n)
      w.resize( n);

   T *a = &mPtr->val[0];
   T *v = &vc(0,0);
   valarray<T> rv1(n);

   T g(0), scale(0), anorm(0);
   for (i=0; i < n; i++)
   {
      l = i + 1;
      rv1[i] = scale * g;
      g = s = scale = T(0);
      if (i < m)
      {
         for (k=i; k < m; k++)
            scale += abs( a[k*n+i]);

         if (scale > epsilon( scale))
         {
            for (k=i; k < m; k++)
            {
               tmp = a[k*n+i] /= scale;
               s += tmp * tmp;
            }
            f = a[i*n+i];
            g = -sign( sqrt(s), f);
            h = f * g - s;
            a[i*n+i] = f - g;

            for (j=l; j < n; j++)
            {
               for (s=T(0),k=i; k < m; k++)
                  s += a[k*n+i] * a[k*n+j];
               f = s / h;
               for (k=i; k < m; k++)
                  a[k*n+j] += f * a[k*n+i];
            }
            for (k=i; k < m; k++)
               a[k*n+i] *= scale;
         }
      }
      w[i] = scale * g;
      g = s = scale = T(0);
      if (i < m && i != n-1)
      {
         for (k=l; k < n; k++)
            scale += abs(a[i*n+k]);

         if (scale > epsilon(scale))
         {
            for (k=l; k < n; k++)
            {
               tmp = a[i*n+k] /= scale;
               s += tmp * tmp;
            }
            f = a[i*n+l];
            g = -sign( sqrt(s), f);
            h = f * g - s;
            a[i*n+l] = f - g;

            for (k=l; k < n; k++)
               rv1[k] = a[i*n+k] / h;

            for (j=l; j < m; j++)
            {
               for (s=T(0),k=l; k < n; k++)
                  s += a[j*n+k] * a[i*n+k];
               for (k=l; k < n; k++)
                  a[j*n+k] += s * rv1[k];
            }
            for (k=l; k < n; k++)
               a[i*n+k] *= scale;
         }
      }
      anorm = techsoft::max( anorm, abs( w[i]) + abs( rv1[i]) );
   }

   for (i=n-1; ; i--)
   {
      if (i < n-1)
      {
         if (abs(g) > epsilon(g))
         {
            for (j=l; j < n; j++)
                v[j*n+i] = (a[i*n+j] / a[i*n+l]) / g;
            for (j=l; j < n; j++)
            {
               for (s=T(0),k=l; k < n; k++)
                  s += a[i*n+k] * v[k*n+j];

               for (k=l; k < n; k++)
                  v[k*n+j] += s * v[k*n+i];
            }
         }
         for (j=l; j < n; j++)
            v[i*n+j] = v[j*n+i] = T(0);
      }
      v[i*n+i] = T(1);
      g = rv1[i];
      l = i;
      if (i == 0)
         break;
   }

   for (i=techsoft::min(m,n)-1; ; i--)
   {
      l = i + 1;
      g = w[i];
      for (j=l; j < n; j++)
         a[i*n+j] = T(0);
      if (abs(g) > epsilon(g))
      {
         g = T(1) / g;
         for (j=l; j < n; j++)
         {
            for (s=T(0),k=l; k < m; k++)
               s += a[k*n+i] * a[k*n+j];

            f = (s / a[i*n+i]) * g;

            for (k=i; k < m; k++)
               a[k*n+j] += f * a[k*n+i];
         }
         for (j=i; j < m; j++)
            a[j*n+i] *= g;

      }
      else
      {
          for (j=i; j < m; j++)
              a[j*n+i] = T(0);
      }
      ++a[i*n+i];
      if (i == 0)
         break;
   }

   for (k=n-1; ; k--)
   {
      for (its=1; its <= 30; its++)
      {
         flag = 1;
         for (l=k; ; l--)
         {
            nm = l-1;
            if (abs( rv1[l]) < epsilon(rv1[l]))
            {
               flag = 0;
               break;
            }
            if (abs( w[nm]) < epsilon(w[nm]))
               break;
            if (l == 0)
               break;
         }
         if (flag)
         {
            c = T(0);
            s = T(1);
            for (i=l; i <= k; i++)
            {
               f = s * rv1[i];
               rv1[i] = c * rv1[i];
               if (abs(f) < epsilon(f))
                  break;
               g = w[i];
               h = hypot( f, g);
               w[i] = h;
               h = T(1) / h;
               c = g * h;
               s = -f * h;
               for (j=0; j < m; j++)
               {
                  y = a[j*n+nm];
                  z = a[j*n+i];
                  a[j*n+nm] = y * c + z * s;
                  a[j*n+i] = z * c - y * s;
               }
            }
         }
         z = w[k];
         if (l == k)                 
         {
            if (z < T(0))
            {
               w[k] = -z;
               for (j=0; j < n; j++)
                  v[j*n+k] = -v[j*n+k];
            }
            break;
         }
         
         if (its == 30)
            return false;

         x = w[l];
         nm = k - 1;
         y = w[nm];
         g = rv1[nm];
         h = rv1[k];
         f = ((y-z) * (y+z) + (g-h) * (g+h)) / (T(2) * h * y);
         g = hypot( f, T(1));
         f = ((x-z) * (x+z) + h * ((y / (f+sign( g, f))) - h)) / x;
         c = s = T(1);                    
         for (j=l; j <= nm; j++)
         {
            i = j + 1;
            g = rv1[i];
            y = w[i];
            h = s * g;
            g = c * g;
            z = hypot( f, h);
            rv1[j] = z;
            c = f / z;
            s = h / z;
            f = x * c + g * s;
            g = g*c - x*s;
            h = y * s;
            y *= c;
            for (jj=0; jj < n; jj++)
            {
               x = v[jj*n+j];
               z = v[jj*n+i];
               v[jj*n+j] = x * c + z * s;
               v[jj*n+i] = z * c - x * s;
            }
            z = hypot( f, h);
            w[j] = z;
            if (abs(z) > epsilon(z))
            {
               z = 1.0 / z;
               c = f * z;
               s = h * z;
            }
            f = c * g + s * y;
            x = c * y - s * g;
            for (jj=0; jj < m; jj++)
            {
               y = a[jj*n+j];
               z = a[jj*n+i];
               a[jj*n+j] = y * c + z * s;
               a[jj*n+i] = z * c - y * s;
            }
         }
         rv1[l] = 0.0;
         rv1[k] = f;
         w[k] = x;
      }
      if (k == 0)
         break;
   }
   return true;
}

template <class T> void
matrix<T>::svbksb (const matrix<T>& vc, const valarray<T>& w, const valarray<T>& b, valarray<T>& s) const
{
   size_t i,j,k;
   const T *pm;
   size_t m = mPtr->nrow;
   size_t n = mPtr->ncol;
   valarray<T> tmp(n);
   
   pm = &mPtr->val[0];
   for (j=0; j < n; j++)                  
   { 
      T ta(0);
      if (abs(w[j]) > epsilon(w[j]))
      { 
         for (i=0; i < m; i++) 
            ta += pm[i*n+j] * b[i];
         ta /= w[j];                       
      }
      tmp[j] = ta;
   }

   pm = &vc.mPtr->val[0];
   for (j=0; j < n; j++)                  
   { 
      T ta(0);
      for (k=0; k < n; k++) 
         ta += pm[j*n+k] * tmp[k];
      s[j] = ta;
   }
}

template <class T> T
matrix<T>::norm1 () const
{
   T nr(0);
   size_t m = mPtr->nrow;
   size_t n = mPtr->ncol;
   const T *pm = &mPtr->val[0];

   for (size_t j=0; j < n; j++)
   {
      T s(0);
      for (size_t i=0; i < m; i++)
         s += abs( pm[i*n+j]);
      nr = techsoft::max( nr, s);
   }
   return nr;
}

template <class T> T
matrix<T>::norm2 () const
{
   matrix<T> m(*this);
   matrix<T> v(mPtr->ncol,mPtr->ncol);
   valarray<T> w(mPtr->ncol);
   if (m.svd( v, w))
      return w.max();
   return T(0);
}

template <class T> T
matrix<T>::normI () const
{
   T nr(0);
   size_t m = mPtr->nrow;
   size_t n = mPtr->ncol;
   const T *pm = &mPtr->val[0];

   for (size_t i=0; i < m; i++)
   {
      T s(0);
      for (size_t j=0; j < n; j++)
         s += abs( pm[i*n+j]);
      nr = techsoft::max( nr, s);
   }
   return nr;
}

template <class T> T
matrix<T>::normF () const
{
   T nr(0);
   size_t n = size();
   const T *pm = &mPtr->val[0];

   for (size_t i=0; i < n; i++)
      nr += pm[i] * pm[i];
   nr = sqrt( nr);

   return nr;
}

template <class T> T
matrix<T>::cond () const
{
   matrix<T> m(*this);
   matrix<T> v(mPtr->ncol,mPtr->ncol);
   valarray<T> w(mPtr->ncol);
   if (m.svd( v, w))
      return w.max()/w.min();
   return T(0);
}

template <class T> size_t 
matrix<T>::rank () const
{
   matrix<T> m(*this);
   matrix<T> v(mPtr->ncol,mPtr->ncol);
   valarray<T> w(mPtr->ncol);

   if (!m.svd( v, w))
      return 0;

   size_t r = 0;
   for (size_t i=0; i < mPtr->ncol; i++)
      if (abs(w[i]) > epsilon( w[i]))
         ++r;
   return r;
}

template <class T> T
matrix<T>::det () const
{
   #ifdef RANGE_CHECK_
   if (mPtr->nrow != mPtr->ncol)
      throw invalid_argument( "matrix<T>::det: Determinant of a non-square matrix!");
   #endif

   T d;

   matrix<T> m(*this);
   valarray<size_t> ri( mPtr->nrow);
   if (!m.lud( ri, &d))
      d = T(0);

   return d;
}

template <class T> T
matrix<T>::cofact (size_t row, size_t col) const
{
   size_t i,i1,j,j1;

   #ifdef RANGE_CHECK_
   if (mPtr->nrow != mPtr->ncol)
      throw invalid_argument( "matrix<T>::Cofact(): Cofactor of a non-square matrix!");
   
   if (row >= mPtr->nrow || col >= mPtr->ncol)
      throw out_of_range( "matrix<T>::Cofact(): Index out of range!");
   #endif

   size_t n = mPtr->nrow;
   matrix<T> tm (n-1,n-1);
   const T *pm = &mPtr->val[0];
   T *ptm = &tm.mPtr->val[0];

   for (i=0,i1=0; i < n; i++)
   {
      if (i == row)
        continue;
      for (j=0,j1=0; j < n; j++)
      {
      	 if (j == col)
            continue;
      	 ptm[i1*(n-1)+j1] = pm[i*n+j];
         j1++;
      }
      i1++;
   }
   T cof = tm.det();
   if ((row+col)%2 == 1)
      cof = -cof;

   return cof;
}

template <class T> matrix<T>
matrix<T>::adj () const
{
   #ifdef RANGE_CHECK_
   if (mPtr->nrow != mPtr->ncol)
      throw invalid_argument( "matrix<T>::adj(): Adjoin of a non-square matrix!");
   #endif

   size_t n = mPtr->nrow;
   matrix<T> tm(n,n);
   T *pm = &tm.mPtr->val[0];

   for (size_t i=0; i < n; i++)
      for (size_t j=0; j < n; j++)
         pm[j*n+i] = cofact(i,j);
   return tm;
}

template <class T>
bool matrix<T>::inv ()
{
   #ifdef RANGE_CHECK_
   if (mPtr->nrow != mPtr->ncol)
      throw invalid_argument( "matrix<T>::inv(): Inversion of non-square matrix!");
   #endif

   if (shared())
      clone();

   size_t n = mPtr->nrow;
   T *pv = &mPtr->val[0];
   size_t i,j,k,ipos,kpos;

   valarray<size_t> ri(n);
   for (i=0; i < n; i++)
      ri[i] = i;

   for (k=0; k < n; k++)
   {
      double ta,tb;
      T a(0);

      kpos = k * n;
      i = k;
      ta = abs( pv[kpos+k]);
      for (j=k+1; j < n; j++)
         if ((tb = abs(pv[j*n+k])) > ta)
         {
            ta = tb;
            i = j;
         }

      if (ta < epsilon( a))
         return false;

      if (i != k)
      {
         swap( ri[k], ri[i]);
         for (ipos=i*n,j=0; j < n; j++)
            swap( pv[kpos+j], pv[ipos+j]);
      }

      a = T(1) / pv[kpos+k];
      pv[kpos+k] = T(1);

      for (j=0; j < n; j++)
         pv[kpos+j] *= a;

      for (i=0; i < n; i++)
      {
         if (i != k)
         {
            ipos = i * n;
            a = pv[ipos+k];
            pv[ipos+k] = T(0);
            for (j=0; j < n; j++)
               pv[ipos+j] -= a * pv[kpos+j];
         }
      }
   }
   for (j=0; j < n; j++)
   {
      if (j != ri[j])         // Column is out of order
      {
         k = j+1;
         while (j != ri[k])
             k++;
         for (i=0; i < n; i++)
            swap( pv[i*n+j], pv[i*n+k]);
         swap( ri[j], ri[k]);
      }
   }
   return true;
}

template <class T>
bool matrix<T>::inv_lu ()
{
   #ifdef RANGE_CHECK_
   if (mPtr->nrow != mPtr->ncol)
      throw invalid_argument( "matrix<T>::inv_lu(): Inversion of non-square matrix!");
   #endif

   size_t n = mPtr->nrow;
   size_t i,j,k;

   valarray<size_t> ri(n);
   valarray<T> v(n), s(n);

   if (!lud( ri))
      return false;

   matrix<T> tm(n,n);
   T *pv = &tm.mPtr->val[0];
   for (j=0; j < n; j++)
   {
      for (i=0; i < n; i++)
         v[i] = T(0);
      v[j] = T(1);
      lubksb( ri, v, s);
      for (k=0,i=0; i < n; i++,k+=n)
         pv[k+j] = s[i];
   }
   *this = tm;
   return true;
}

template <class T>
bool matrix<T>::inv_qr ()
{
   #ifdef RANGE_CHECK_
   if (mPtr->nrow != mPtr->ncol)
      throw invalid_argument( "matrix<T>::inv_qr(): Inversion of non-square matrix!");
   #endif

   size_t n = mPtr->nrow;
   size_t i,j,k;

   valarray<T> v(n), s(n);
   matrix<T> tm(n,n), r(n,n);

   qrd( r);
   T *pv = &tm.mPtr->val[0];
   for (j=0; j < n; j++)
   {
      for (i=0; i < n; i++)
         v[i] = T(0);
      v[j] = T(1);
      qrbksb( r, v, s);
      for (k=0,i=0; i < n; i++,k+=n)
         pv[k+j] = s[i];
   }
   *this = tm;
   return true;
}

template <class T>
bool matrix<T>::inv_sv ()
{
   #ifdef RANGE_CHECK_
   if (mPtr->nrow != mPtr->ncol)
      throw invalid_argument( "matrix<T>::inv2(): Inversion of non-square matrix!");
   #endif

   size_t n = mPtr->ncol;
   size_t i,j,k;

   valarray<T> w(n), v(n), s(n);
   matrix<T> vc(n,n);

   if (!svd( vc, w))
      return false;

   matrix<T> tm(n,n);
   T *pv = &tm.mPtr->val[0];
   for (j=0; j < n; j++)
   {
      for (i=0; i < n; i++)
         v[i] = T(0);
      v[j] = T(1);
      svbksb( vc, w, v, s);
      for (k=0,i=0; i < n; i++,k+=n)
         pv[k+j] = s[i];
   }
   *this = tm;
   return true;
}

template <class T> inline bool
matrix<T>::isSingular () const
{
   if (mPtr->nrow != mPtr->ncol)
      return false;
   return (abs( det()) < epsilon(T(0)));
}

template <class T> bool
matrix<T>::isDiagonal () const
{
   if (mPtr->nrow != mPtr->ncol)
      return false;

   size_t n = mPtr->nrow;
   const T *pm = &mPtr->val[0];

   for (size_t i=0; i < n; i++)
      for (size_t j=0; j < n; j++)
         if (i != j && abs( pm[i*n+j]) > epsilon( pm[i]))
            return false;
   return true;
}

template <class T> bool
matrix<T>::isScalar () const
{
   if (!isDiagonal())
     return false;

   T v = mPtr->val[0];
   size_t n = mPtr->nrow;
   const T *pm = &mPtr->val[0];

   for (size_t i=1; i < n; i++)
     if (abs( pm[i*n+i] - v) > epsilon( v))
        return false;
        
   return true;
}

template <class T> bool
matrix<T>::isUnit () const throw()
{
   if (mPtr->nrow != mPtr->ncol)
      return false;

   size_t n = mPtr->nrow;
   const T *pm = &mPtr->val[0];

   for (size_t k=0; k < n; k++)
     if (abs( pm[k*n+k] - T(1)) > epsilon( T(1)))
       return false;

	for (size_t j=1; j < n; j++)
	  for (size_t i=0; i < j; i++)
       if (abs( pm[i*n+j]) > epsilon( T(1)) || abs( pm[j*n+i]) > epsilon( T(1)))
          return false;

   return true;
}

template <class T> bool
matrix<T>::isNull () const throw()
{
   size_t n = size();
   const T *pm = &mPtr->val[0];

   for (size_t i=0; i < n; i++)
      if (abs( pm[i]) > epsilon( pm[i]))
         return false;

   return true;
}

template <class T> bool
matrix<T>::isSymmetric () const
{
   if (mPtr->nrow != mPtr->ncol)
      return false;

   size_t n = mPtr->nrow;
   const T *pm = &mPtr->val[0];

	for (size_t j=1; j < n; j++)
		for (size_t i=0; i < j; i++)
			if (abs( pm[i*n+j] - pm[j*n+i]) > epsilon( pm[0]))
            return false;
   return true;
}

template <class T> bool
matrix<T>::isSkewSymmetric () const
{
   if (mPtr->nrow != mPtr->ncol)
      return false;

   size_t n = mPtr->nrow;
   const T *pm = &mPtr->val[0];

	for (size_t j=1; j < n; j++)
		for (size_t i=0; i < j; i++)
			if (abs( pm[i*n+j] + pm[j*n+i]) > epsilon( pm[0]))
            return false;
   return true;
}

template <class T> bool
matrix<T>::isUpperTriangular () const
{
   if (mPtr->nrow != mPtr->ncol)
      return false;

   size_t n = mPtr->nrow;
   const T *pm = &mPtr->val[0];

   for (size_t i=1; i < n; i++)
      for (size_t j=0; j < i-1; j++)
         if (abs( pm[i*n+j]) > epsilon( pm[0]))
            return false;
   return true;
}

template <class T> bool
matrix<T>::isLowerTriangular () const
{
   if (mPtr->nrow != mPtr->ncol)
      return false;

   size_t n = mPtr->nrow;
   const T *pm = &mPtr->val[0];

   for (size_t j=1; j < n; j++)
      for (size_t i=0; i < j-1; i++)
         if (abs( pm[i*n+j]) > epsilon( pm[0]))
            return false;
   return true;
}

template <class T> inline bool
matrix<T>::isRowOrthogonal () const
{
   matrix<T> tm = *this * ~(*this);
   return tm.isUnit();
}

template <class T> inline bool
matrix<T>::isColOrthogonal () const
{
   matrix<T> tm = ~(*this) * (*this);
   return tm.isUnit();
}

template <class T> matrix<T> 
matrix<T>::apply (T (*fn)(T)) const
{
   matrix<T> a( mPtr->nrow, mPtr->ncol);
   size_t n = a.size();

   const T *pm = &mPtr->val[0];
   T *pa = &a.mPtr->val[0];

   for (size_t i=0; i < n; i++)
      pa[i] = fn( pm[i]);

   return a;
}

template <class T> matrix<T> 
matrix<T>::apply (T (*fn)(size_t,size_t,T)) const
{
   size_t m = mPtr->nrow;
   size_t n = mPtr->ncol;
   matrix<T> a( m, n);

   const T *pm = &mPtr->val[0];
   T *pa = &a.mPtr->val[0];
   size_t i,j,im,jn;

   for (i=0,im=0; i < m; im+=n,i++)
      for (j=0,jn=im; j < n; jn++,j++)
         pa[jn] = fn( i, j, pm[jn]);

   return a;
}

template <class T> matrix<T> 
matrix<T>::apply (T (*fn)(size_t,size_t,const T&)) const
{
   size_t m = mPtr->nrow;
   size_t n = mPtr->ncol;
   matrix<T> a( m, n);

   const T *pm = &mPtr->val[0];
   T *pa = &a.mPtr->val[0];
   size_t i,j,im,jn;

   for (i=0,im=0; i < m; im+=n,i++)
      for (j=0,jn=im; j < n; jn++,j++)
         pa[jn] = fn( i, j, pm[jn]);

   return a;
}

template <class T> inline matrix<T> 
matrix<T>::apply (T (*fn)(const T&)) const
{
   matrix<T> a( mPtr->nrow, mPtr->ncol);
   size_t n = a.size();

   const T *pm = &mPtr->val[0];
   T *pa = &a.mPtr->val[0];

   for (size_t i=0; i < n; i++)
      pa[i] = fn( pm[i]);

   return a;
}

template <class T> bool 
matrix<T>::eigen (valarray<T>& eival) const
{
   #ifdef RANGE_CHECK_
   if (!isSymmetric())
      return false;
   #endif
   
   if (eival.size() != mPtr->nrow)
      eival.resize( mPtr->nrow);

   matrix<T> a(*this);
   valarray<T> e(mPtr->nrow);

   a.tred2( eival, e, false);
   return a.tql2( eival, e, false);
}

template <class T> bool 
matrix<T>::eigen (valarray<T>& eival, matrix<T>& eivec) const
{
   #ifdef RANGE_CHECK_
   if (!isSymmetric())
      return false;
   #endif
   
   valarray<T> e(mPtr->nrow);

   if (eival.size() != mPtr->nrow)
      eival.resize( mPtr->nrow);

   eivec = *this;
   eivec.tred2( eival, e, true);

   return eivec.tql2( eival, e, true);
}

template <class T> bool 
matrix<T>::eigen (valarray<T>& rev, valarray<T>& iev) const
{
   #ifdef RANGE_CHECK_
   if (mPtr->nrow != mPtr->ncol)
      return false;
   #endif
   
   if (rev.size() != mPtr->nrow)
      rev.resize( mPtr->nrow);
   if (iev.size() != mPtr->nrow)
      iev.resize( mPtr->nrow);

   matrix<T> eivec, m(*this);
   m.balanc( eivec, false);

   return m.hqr2( rev, iev, eivec, false);
}

template <class T> bool 
matrix<T>::eigen (valarray<T>& rev, valarray<T>& iev, matrix<T>& eivec) const
{
   #ifdef RANGE_CHECK_
   if (mPtr->nrow != mPtr->ncol)
      return false;
   #endif
   
   if (rev.size() != mPtr->nrow)
      rev.resize( mPtr->nrow);
   if (iev.size() != mPtr->nrow)
      iev.resize( mPtr->nrow);
   if (eivec.rowno() != mPtr->nrow || eivec.colno() != mPtr->nrow)
      eivec.resize( mPtr->nrow, mPtr->nrow);

   matrix<T> m(*this);
   m.balanc( eivec, true);

   return m.hqr2( rev, iev, eivec, true);
}

template <class T> inline
void cdiv (const T& xr, const T& xi, const T& yr, const T& yi, T& cdivr, T& cdivi) 
{
   T r,d;

   if (abs( yr) > abs( yi)) 
   {
      r = yi/yr;
      d = yr + r*yi;
      cdivr = (xr + r*xi)/d;
      cdivi = (xi - r*xr)/d;
   } 
   else 
   {
      r = yr/yi;
      d = yi + r*yr;
      cdivr = (r*xr + xi)/d;
      cdivi = (r*xi - xr)/d;
   }
}

////////////////////////
//  This is derived from the Algol procedure hqr2, by Martin 
//  and Wilkinson, Handbook for Auto. Comp.,Vol.ii-Linear Algebra, 
//  and the corresponding Fortran subroutine in EISPACK.
//
template <class T> bool
matrix<T>::hqr2 (valarray<T>& d, valarray<T>& e, matrix<T>& v, bool eivec)
{
   int i,j,k,l;

   int nn = int(mPtr->ncol);
   int n = nn-1;
   int low = 0;
   int high = nn-1;
   T exshift(0);
   T p(0),q(0),r(0),s(0),z(0),t,w,x,y;
   T cdivr, cdivi;

   if (shared())
      clone();

   T *pv, *pm = &mPtr->val[0];
   if (eivec)
      pv = &v(0,0);

   T norm(0);
   for (i=0; i < nn; i++) 
   {
      if (i < low || i > high) 
      {
         d[i] = pm[i*nn+i];
         e[i] = T(0);
      }
      for (j = techsoft::max(i-1,0); j < nn; j++)
         norm += abs( pm[i*nn+j]);
   }

   size_t iter = 0;
   while (n >= low) 
   {
      l = n;
      while (l > low) 
      {
         s = abs( pm[(l-1)*nn+l-1]) + abs( pm[l*nn+l]);
         if (s < epsilon(s))
            s = norm;
        
         if (abs( pm[l*nn+l-1]) < epsilon(pm[l*nn+l-1]) * s)
            break;
        
         l--;     
      }
    
      if (l == n) 
      {
         pm[n*nn+n] += exshift;
         d[n] = pm[n*nn+n];
         e[n] = T(0);
         n--;
         iter = 0;
      } 
      else if (l == n-1) 
      {
         w = pm[n*nn+n-1] * pm[(n-1)*nn+n];
         p = (pm[(n-1)*nn+n-1] - pm[n*nn+n]) / T(2);
         q = p * p + w;
         z = sqrt( abs( q));
         pm[n*nn+n] += exshift;
         pm[(n-1)*nn+n-1] += exshift;
         x = pm[n*nn+n];

         if (q >= T(0)) 
         {
            if (p >= T(0))
               z = p + z;
            else
               z = p - z;
           
            d[n-1] = x + z;
            d[n] = d[n-1];
            if (z != T(0))
               d[n] = x - w / z;
           
            e[n-1] = T(0);
            e[n] = T(0);

            x = pm[n*nn+n-1];
            s = abs( x) + abs( z);
            p = x / s;
            q = z / s;
            r = sqrt( p*p + q*q);
            p = p / r;
            q = q / r;

            for (j = n-1; j < nn; j++) 
            {
               z = pm[(n-1)*nn+j];
               pm[(n-1)*nn+j] = q * z + p * pm[n*nn+j];
               pm[n*nn+j] = q * pm[n*nn+j] - p * z;
            }

            for (i = 0; i <= n; i++) 
            {
               z = pm[i*nn+n-1];
               pm[i*nn+n-1] = q * z + p * pm[i*nn+n];
               pm[i*nn+n] = q * pm[i*nn+n] - p * z;
            }

            if (eivec)
               for (i = low; i <= high; i++) 
               {
                  z = pv[i*nn+n-1];
                  pv[i*nn+n-1] = q * z + p * pv[i*nn+n];
                  pv[i*nn+n] = q * pv[i*nn+n] - p * z;
               }
         } 
         else 
         {
            d[n-1] = x + p;
            d[n] = x + p;
            e[n-1] = z;
            e[n] = -z;
         }
         n = n - 2;
         iter = 0;
      } 
      else 
      {
         x = pm[n*nn+n];
         y = T(0);
         w = T(0);
         if (l < n) 
         {
            y = pm[(n-1)*nn+n-1];
            w = pm[n*nn+n-1] * pm[(n-1)*nn+n];
         }

         if (iter == 10) 
         {
            exshift += x;
            for (i = low; i <= n; i++)
               pm[i*nn+i] -= x;
           
            s = abs( pm[n*nn+n-1]) + abs( pm[(n-1)*nn+n-2]);
            x = y = 0.75 * s;
            w = -0.4375 * s * s;
         }

         if (iter == 30) 
         {
             s = (y - x) / T(2);
             s = s * s + w;
             if (s > T(0)) 
             {
                 s = sqrt( s);
                 if (y < x)
                    s = -s;
                
                 s = x - w / ((y - x) / T(2) + s);
                 for (i = low; i <= n; i++)
                    pm[i*nn+i] -= s;
                 
                 exshift += s;
                 x = y = w = 0.964;
             }
         }

         if (++iter > 250)
            return false;

         int m = n-2;
         while (m >= l) 
         {
            z = pm[m*nn+m];
            r = x - z;
            s = y - z;
            p = (r * s - w) / pm[(m+1)*nn+m] + pm[m*nn+m+1];
            q = pm[(m+1)*nn+m+1] - z - r - s;
            r = pm[(m+2)*nn+m+1];
            s = abs( p) + abs( q) + abs( r);
            p = p / s;
            q = q / s;
            r = r / s;
            if (m == l)
               break;
            
            if (abs( pm[m*nn+m-1]) * (abs( q) + abs( r)) <
               epsilon(r) * (abs( p) * (abs( pm[(m-1)*nn+m-1]) + abs( z) +
               abs( pm[(m+1)*nn+m+1])))) {
                  break;
            }
            m--;
         }

         for (i = m+2; i <= n; i++) 
         {
            pm[i*nn+i-2] = T(0);
            if (i > m+2)
               pm[i*nn+i-3] = T(0);
         }

         for (k = m; k <= n-1; k++) 
         {
            bool notlast = (k != n-1);
            if (k != m) 
            {
               p = pm[k*nn+k-1];
               q = pm[(k+1)*nn+k-1];
               r = notlast ? pm[(k+2)*nn+k-1] : T(0);
               x = abs( p) + abs( q) + abs( r);
               if (x > epsilon(x)) 
               {
                  p = p / x;
                  q = q / x;
                  r = r / x;
               }
            }
            if (x < epsilon(x))
               break;
            
            s = sqrt( p * p + q * q + r * r);
            if (p < T(0))
               s = -s;
            
            if (s != T(0)) 
            {
               if (k != m)
                  pm[k*nn+k-1] = -s * x;
               else if (l != m)
                  pm[k*nn+k-1] = -pm[k*nn+k-1];
               
               p = p + s;
               x = p / s;
               y = q / s;
               z = r / s;
               q = q / p;
               r = r / p;

               for (j = k; j < nn; j++) 
               {
                  p = pm[k*nn+j] + q * pm[(k+1)*nn+j];
                  if (notlast) 
                  {
                     p = p + r * pm[(k+2)*nn+j];
                     pm[(k+2)*nn+j] = pm[(k+2)*nn+j] - p * z;
                  }
                  pm[k*nn+j] = pm[k*nn+j] - p * x;
                  pm[(k+1)*nn+j] = pm[(k+1)*nn+j] - p * y;
               }

               for (i = 0; i <= techsoft::min( n, k+3); i++) 
               {
                  p = x * pm[i*nn+k] + y * pm[i*nn+k+1];
                  if (notlast) 
                  {
                     p = p + z * pm[i*nn+k+2];
                     pm[i*nn+k+2] = pm[i*nn+k+2] - p * r;
                  }
                  pm[i*nn+k] = pm[i*nn+k] - p;
                  pm[i*nn+k+1] = pm[i*nn+k+1] - p * q;
               }

               if (eivec)
                  for (i = low; i <= high; i++) 
                  {
                     p = x * pv[i*nn+k] + y * pv[i*nn+k+1];
                     if (notlast) 
                     {
                        p += z * pv[i*nn+k+2];
                        pv[i*nn+k+2] -= p * r;
                     }
                     pv[i*nn+k] -= p;
                     pv[i*nn+k+1] -= p * q;
                  }
            }
         }
      }
   }

   if (norm < epsilon(norm))
      return true;

   for (n = nn-1; n >= 0; n--) 
   {
      p = d[n];
      q = e[n];

      if (q == T(0)) 
      {
         int l = n;
         pm[n*nn+n] = T(1);
         for (i = n-1; i >= 0; i--) 
         {
            w = pm[i*nn+i] - p;
            r = T(0);
            for (j = l; j <= n; j++)
               r += pm[i*nn+j] * pm[j*nn+n];
            
            if (e[i] < T(0)) 
            {
               z = w;
               s = r;
            } 
            else 
            {
               l = i;
               if (e[i] == T(0)) 
               {
                  if (w != T(0))
                     pm[i*nn+n] = -r / w;
                  else
                     pm[i*nn+n] = -r / (epsilon(norm) * norm);
               } 
               else 
               {
                  x = pm[i*nn+i+1];
                  y = pm[(i+1)*nn+i];
                  q = (d[i] - p) * (d[i] - p) + e[i] * e[i];
                  t = (x * s - z * r) / q;
                  pm[i*nn+n] = t;
                  if (abs( x) > abs( z))
                     pm[(i+1)*nn+n] = (-r - w * t) / x;
                  else
                     pm[(i+1)*nn+n] = (-s - y * t) / z;
               }

               t = abs( pm[i*nn+n]);
               if ((epsilon(t) * t) * t > T(1)) 
                  for (j = i; j <= n; j++)
                     pm[j*nn+n] /= t;
            }
         }
      } 
      else if (q < T(0)) 
      {
         int l = n-1;

         if (abs( pm[n*nn+n-1]) > abs( pm[(n-1)*nn+n])) 
         {
            pm[(n-1)*nn+n-1] = q / pm[n*nn+n-1];
            pm[(n-1)*nn+n] = -(pm[n*nn+n] - p) / pm[n*nn+n-1];
         } 
         else 
         {
            cdiv( 0.0, -pm[(n-1)*nn+n], pm[(n-1)*nn+n-1]-p, q, cdivr, cdivi);
            pm[(n-1)*nn+n-1] = cdivr;
            pm[(n-1)*nn+n] = cdivi;
         }
         pm[n*nn+n-1] = T(0);
         pm[n*nn+n] = T(1);
         for (i = n-2; i >= 0; i--) 
         {
            T ra(0),sa(0),vr,vi;

            for (j = l; j <= n; j++) 
            {
               ra += pm[i*nn+j] * pm[j*nn+n-1];
               sa += pm[i*nn+j] * pm[j*nn+n];
            }
            w = pm[i*nn+i] - p;

            if (e[i] < T(0)) 
            {
               z = w;
               r = ra;
               s = sa;
            } 
            else 
            {
               l = i;
               if (e[i] == T(0))
               {
                  cdiv( -ra, -sa, w, q, pm[i*nn+n-1], pm[i*nn+n]);
               } 
               else 
               {
                  x = pm[i*nn+i+1];
                  y = pm[(i+1)*nn+i];
                  vr = (d[i] - p) * (d[i] - p) + e[i] * e[i] - q * q;
                  vi = (d[i] - p) * 2.0 * q;
                  if (vr == T(0) && vi == T(0)) 
                     vr = epsilon(norm) * norm * (abs( w) + abs(q) + abs(x) + abs(y) + abs(z));

                  cdiv( x*r - z*ra + q*sa, x*s - z*sa - q*ra, vr, vi, cdivr, cdivi);
                  pm[i*nn+n-1] = cdivr;
                  pm[i*nn+n] = cdivi;
                  if (abs(x) > (abs(z) + abs(q))) 
                  {
                     pm[(i+1)*nn+n-1] = (-ra - w * pm[i*nn+n-1] + q * pm[i*nn+n]) / x;
                     pm[(i+1)*nn+n] = (-sa - w * pm[i*nn+n] - q * pm[i*nn+n-1]) / x;
                  } 
                  else 
                  {
                     cdiv( -r - y * pm[i*nn+n-1], -s - y * pm[i*nn+n], z, q, cdivr, cdivi);
                     pm[(i+1)*nn+n-1] = cdivr;
                     pm[(i+1)*nn+n] = cdivi;
                  }
               }

               t = techsoft::max( abs( pm[i*nn+n-1]), abs( pm[i*nn+n]));
               if ((epsilon(t) * t) * t > T(1)) 
                  for (j = i; j <= n; j++) 
                  {
                     pm[j*nn+n-1] /= t;
                     pm[j*nn+n] /= t;
                  }
            }
         }
      }
   }

   if (eivec)
   {
      for (i = 0; i < nn; i++) 
         if (i < low || i > high) 
            for (j = i; j < nn; j++)
               pv[i*nn+j] = pm[i*nn+j];
         
      for (j = nn-1; j >= low; j--) 
         for (i = low; i <= high; i++) 
         {
            z = T(0);
            for (k = low; k <= techsoft::min( j, high); k++)
               z += pv[i*nn+k] * pm[k*nn+j];
         
            pv[i*nn+j] = z;
         }
   }
   return true;
}

///////////////////////////
//  This is derived from the Algol procedures orthes and ortran,
//  by Martin and Wilkinson, Handbook for Auto. Comp.,Vol.ii-Linear 
//  Algebra, and the corresponding Fortran subroutines in EISPACK.
//
template <class T> void 
matrix<T>::balanc (matrix<T>& v, bool eivec) 
{
   size_t i,j,lo,hi,m,n;

   if (shared())
      clone();

   n = mPtr->ncol;
   lo = 0;
   hi = n-1;
   valarray<T> ort(n);
   T *pv, *pm = &mPtr->val[0];
   if (eivec)
      pv = &v(0,0);

   for (m=lo+1; m <= hi-1; m++) 
   {
      T scale(0);
      for (i=m; i <= hi; i++)
         scale += abs( pm[i*n+m-1]);

      if (scale > epsilon(scale)) 
      {
         T h(0);
         for (i=hi; i >= m; i--) 
         {
            ort[i] = pm[i*n+m-1] / scale;
            h += ort[i] * ort[i];
         }
         T g = sqrt( h);
         if (ort[m] > T(0)) 
            g = -g;
         
         h -= ort[m] * g;
         ort[m] -= g;

         for (j=m; j < n; j++) 
         {
            T f(0);
            for (i = hi; i >= m; i--) 
               f += ort[i] * pm[i*n+j];
            f /= h;
            for (i=m; i <= hi; i++)
               pm[i*n+j] -= f * ort[i];
        }

        for (i=0; i <= hi; i++) 
        {
            T f(0);
            for (j=hi; j >= m; j--)
               f += ort[j] * pm[i*n+j];
            f /= h;
            for (j=m; j <= hi; j++)
               pm[i*n+j] -= f * ort[j];
         }
         ort[m] = scale * ort[m];
         pm[m*n+m-1] = scale * g;
      }
   }

   if (eivec)
      for (i=0; i < n; i++)
         for (j = 0; j < n; j++)
            pv[i*n+j] = (i == j ? T(1) : T(0));

   for (m=hi-1; m >= lo+1; m--) 
   {
      if (abs(pm[m*n+m-1]) > T(0)) 
      {
         for (i=m+1; i <= hi; i++)
            ort[i] = pm[i*n+m-1];

         if (eivec)
            for (j=m; j <= hi; j++) 
            {
               T g(0);
               for (i=m; i <= hi; i++)
                  g += ort[i] * pv[i*n+j];

               g = (g / ort[m]) / pm[m*n+m-1];
               for (i=m; i <= hi; i++)
                  pv[i*n+j] += g * ort[i];
            }
      }
   }
}

template <class T> bool
matrix<T>::tql2 (valarray<T>& d, valarray<T>& e, bool eivec)
{
   size_t m,l,iter,i,k;
   T s,r,p,g,f,c,b;
   
   size_t n = mPtr->nrow;
   T *pm = &mPtr->val[0];

   for (i=1; i < n; i++)
      e[i-1] = e[i]; 
   e[n-1] = T(0);

   for (l=0; l < n; l++) 
   {
      iter=0;
      do 
      {
         for (m=l; m < n-1; m++)
         { 
            T dd = abs( d[m]) + abs( d[m+1]);
            if (abs( e[m]) <= epsilon(e[m])*dd) 
               break;
         }
         if (m != l) 
         {
            if (iter++ == 30) 
               return false;

            g = (d[l+1] - d[l]) / (T(2) * e[l]);
            r = hypot( g, T(1));
            g = d[m] - d[l] + e[l] / (g + sign( r, g));
            s = c = T(1);
            p = T(0);
            for (i=m-1; ; i--)
            {
               f = s * e[i];
               b = c * e[i];
               r = hypot( f, g);
               e[i+1] = r;
               if (r < epsilon(r))
               { 
                  d[i+1] -= p;
                  e[m] = T(0);
                  break;
               }
               s = f / r;
               c = g / r;
               g = d[i+1] - p;
               r = (d[i] - g) * s + T(2) * c * b;
               p = s * r;
               d[i+1] = g + p;
               g = c * r - b;

               if (eivec)
                  for (k=0; k < n; k++)
                  { 
                     f = pm[k*n+i+1];
                     pm[k*n+i+1] = s * pm[k*n+i] + c * f;
                     pm[k*n+i] = c * pm[k*n+i] - s * f;
                  }
               if (i == l)
                  break;
            }
            if (abs(r) < epsilon(r) && i >= l) 
               continue;
            d[l] -= p;
            e[l] = g;
            e[m] = T(0);
         }
      } 
      while (m != l);

      if (abs(d[l]) < epsilon(d[l]))
         d[l] = T(0);
   }
   return true;
}

template <class T> void
matrix<T>::tred2 (valarray<T>& d, valarray<T>& e, bool eivec)
{
   size_t i,j,k;
   size_t n = mPtr->nrow;

   #ifdef RANGE_CHECK_
   if (n != mPtr->ncol)
      throw std::invalid_argument( "matrix<T>::tred2: non-square matrix!");
   #endif

   if (shared())
      clone();
   T *pm = &mPtr->val[0];

   for (i=n-1; i > 0; i--) 
   {
      T h(0), scale(0);

      if (i > 1) 
      {
         for (k=0; k < i; k++)
            scale += abs( pm[i*n+k]);

         if (scale < epsilon(scale))
            e[i] = pm[i*n+i-1];
         else 
         {
            for (k=0; k < i; k++) 
            {
               pm[i*n+k] /= scale;
               T tmp = pm[i*n+k]; 
               h += tmp * tmp;
            }
            T f = pm[i*n+i-1];
            T g = sqrt( h);
            if (f >= T(0))
               g = -g;

            e[i] = scale * g;
            h -= f * g;
            pm[i*n+i-1] = f - g;
            f = T(0);
            for (j=0; j < i; j++) 
            {
               if (eivec)
                  pm[j*n+i] = pm[i*n+j] / h;
               g = T(0);
               for (k=0; k <= j; k++)
                  g += pm[j*n+k] * pm[i*n+k];
               for (k=j+1; k < i; k++)
                  g += pm[k*n+j] * pm[i*n+k];
               e[j] = g / h;
               f += e[j] * pm[i*n+j];
            }
            T hh = f / (h+h);
            for (j=0; j < i; j++)
            { 
               f = pm[i*n+j];
               g = e[j] - hh * f;
               e[j] = g;
               for (k=0; k <= j; k++)
                  pm[j*n+k] -= (f * e[k] + g * pm[i*n+k]);
            }
         }
      } 
      else
         e[i] = pm[i*n+i-1];
      d[i] = h;
   }

   if (eivec)
      d[0] = T(0);
   e[0] = T(0);

   for (i=0; i < n; i++)
   { 
      if (eivec)
      {
         if (d[i] != T(0))
         { 
            for (j=0; j < i; j++) 
            {
               T g(0);
               for (k=0; k < i; k++)
                  g += pm[i*n+k] * pm[k*n+j];
               for (k=0; k < i; k++)
                  pm[k*n+j] -= g * pm[k*n+i];
            }
         }
         d[i] = pm[i*n+i];
         pm[i*n+i] = T(1);
         for (j=0; j < i; j++) 
            pm[j*n+i] = pm[i*n+j] = T(0);
      }
      else
         d[i] = pm[i*n+i];
   }
}


/////////////////////////////////
//  I/O functions
//

#ifndef MATRIX_ROW_SEP
#define MATRIX_ROW_SEP   ('\n')
#endif

#ifndef MATRIX_COL_SEP
#define MATRIX_COL_SEP   ('\t')
#endif

// Input stream function for a single element
template <class T> inline istream&
operator>> (istream& is, Mref<T> el)
{
    T x; is >> x; el = x; return is;
}

// Output stream function for a single element
template <class T> inline ostream&
operator<< (ostream& os, const Mref<T>& el)
{
    const T x = el;
    os << x;
    return os;
}

// Output stream function for matrix row/column
template <class T> ostream&
operator<< (ostream& os, const Mat_iter<T>& v)
{
    int w = os.width();
    for (size_t i=0; i < v.size(); i++)
        os << setw( w) << v[i] << MATRIX_COL_SEP;
    return os;
}

// Iutput stream function for matrix row/column
template <class T> istream&
operator>> (istream& is, Mat_iter<T> v)
{
    for (size_t i=0; i < v.size(); i++)
        is >> v[i];
    return is;
}


// Output stream function for const matrix row/column
template <class T> ostream&
operator<< (ostream& os, const Cmat_iter<T>& v)
{
    int w = os.width();
    for (size_t i=0; i < v.size(); i++)
        os << setw( w) << v[i] << MATRIX_COL_SEP;
    return os;
}

// Input stream function
template <class T> istream&
operator>> (istream& is, matrix<T>& m)
{
    T *pm = &m(0,0);
    size_t n = m.size();
    for (size_t i=0; i < n; i++)
       is >> pm[i];
    return is;
}

// Output stream function for matrix
template <class T> ostream&
operator<< (ostream &os, const matrix<T>& m)
{
    int w = os.width();
    for (size_t i=0; i < m.rowno(); i++)
       for (size_t j=0; j < m.colno(); j++)
          os << setw( w) << m(i,j) << (j == m.colno()-1 ? MATRIX_ROW_SEP : MATRIX_COL_SEP);
   return os;
}


// Iutput stream function for valarray
template <class T> istream&
operator>> (istream& is, valarray<T>& v)
{
    for (size_t i=0; i < v.size(); i++)
        is >> v[i];
   return is;
}

// Output stream function for valarray
template <class T> ostream&
operator<< (ostream& os, const valarray<T>& v)
{
    int w = os.width();
    for (size_t i=0; i < v.size(); i++)
        os << setw( w) << v[i] << MATRIX_COL_SEP;
    return os;
}


//
// Specializations
//

template <> inline 
void matrix<complex<float> >::rand (int rmin, int rmax, int rseed)
{
   typedef value_type T;
   typedef float VT;

   if (!rmin && !rmax)
   {
      rmin = -1;
      rmax = 1;
   }
   if (rmin > rmax)
   {
      int tmp = rmin;
      rmin = rmax;
      rmax = tmp;
   }
   int rdiff = rmax - rmin;

   if (rseed)
      srand( rseed);

   if (shared())
   {
      --mPtr->refcnt;
      mPtr = allocator( Allocate, mPtr->nrow, mPtr->ncol);
   }
   
   T *pv = &mPtr->val[0];
   size_t n = size();

   for (size_t i=0; i < n; ++i)
   {
      VT rd[2];

      for (size_t j=0; j < 2; ++j)
      {
         int rn = rand_();
         if (rdiff > 1)
         {
            rn %= rdiff;
            rd[j] = VT(rn) - VT(rdiff) / VT(2);
            if (rdiff > 1)
            {
               VT f = ((VT)rand_()) / VT(RAND_MAX);
               if (rn > rdiff)
                  rd[j] -= f;
               else
                  rd[j] += f;
            }
         }
         else
            rd[j] = VT(rn) / VT(RAND_MAX);
      }
      pv[i] = T( rd[0], rd[1]);
   }
}

template <> inline 
void matrix<complex<double> >::rand (int rmin, int rmax, int rseed)
{
   typedef value_type T;
   typedef double VT;

   value_type y = value_type(1.0);

   if (!rmin && !rmax)
   {
      rmin = -1;
      rmax = 1;
   }
   if (rmin > rmax)
   {
      int tmp = rmin;
      rmin = rmax;
      rmax = tmp;
   }
   int rdiff = rmax - rmin;

   if (rseed)
      srand( rseed);

   if (shared())
   {
      --mPtr->refcnt;
      mPtr = allocator( Allocate, mPtr->nrow, mPtr->ncol);
   }
   
   T *pv = &mPtr->val[0];
   size_t n = size();

   for (size_t i=0; i < n; ++i)
   {
      VT rd[2];

      for (size_t j=0; j < 2; ++j)
      {
         int rn = rand_();
         if (rdiff > 1)
         {
            rn %= rdiff;
            rd[j] = VT(rn) - VT(rdiff) / VT(2);
            if (rdiff > 1)
            {
               VT f = ((VT)rand_()) / VT(RAND_MAX);
               if (rn > rdiff)
                  rd[j] -= f;
               else
                  rd[j] += f;
            }
         }
         else
            rd[j] = VT(rn) / VT(RAND_MAX);
      }
      pv[i] = T( rd[0], rd[1]);
   }
}

template <> inline 
void matrix<complex<long double> >::rand (int rmin, int rmax, int rseed)
{
   typedef value_type T;
   typedef long double VT;

   if (!rmin && !rmax)
   {
      rmin = -1;
      rmax = 1;
   }
   if (rmin > rmax)
   {
      int tmp = rmin;
      rmin = rmax;
      rmax = tmp;
   }
   int rdiff = rmax - rmin;

   if (rseed)
      srand( rseed);

   if (shared())
   {
      --mPtr->refcnt;
      mPtr = allocator( Allocate, mPtr->nrow, mPtr->ncol);
   }
   
   T *pv = &mPtr->val[0];
   size_t n = size();

   for (size_t i=0; i < n; ++i)
   {
      VT rd[2];

      for (size_t j=0; j < 2; ++j)
      {
         int rn = rand_();
         if (rdiff > 1)
         {
            rn %= rdiff;
            rd[j] = VT(rn) - VT(rdiff) / VT(2);
            if (rdiff > 1)
            {
               VT f = ((VT)rand_()) / VT(RAND_MAX);
               if (rn > rdiff)
                  rd[j] -= f;
               else
                  rd[j] += f;
            }
         }
         else
            rd[j] = VT(rn) / VT(RAND_MAX);
      }
      pv[i] = T( rd[0], rd[1]);
   }
}

template <> inline
complex<float> sign (complex<float> a, complex<float> b)
{
   if (abs( a-b) < epsilon( a))
      return a;

   return (b.real() < 0.0f || b.imag() < 0.0f)? -a : a;
}

template <> inline
complex<double> sign (complex<double> a, complex<double> b)
{
   if (abs( a-b) < epsilon( a))
      return a;

   return (b.real() < double(0) || b.imag() < double(0))? -a : a;
}

template <> inline
complex<long double> sign (complex<long double> a, complex<long double> b)
{
   if (abs( a-b) < epsilon( a))
      return a;

   return (b.real() < 0.0 || b.imag() < 0.0)? -a : a;
}

template <> inline 
bool matrix<complex<float> >::svd (matrix<complex<float> >& vc, valarray<complex<float> >& w)
{
   typedef value_type T;
   typedef float VT;

   size_t flag,i,its,j,jj,k,l,nm;
   T c,f,g(0),h,s,x,y,z,tmp;

   if (shared())
      clone();

   size_t m = mPtr->nrow;
   size_t n = mPtr->ncol;

   if (vc.rowno() != n || vc.colno() != n)
      vc.resize( n, n);
   if (w.size() != n)
      w.resize( n);

   T *a = &mPtr->val[0];
   T *v = &vc(0,0);
   valarray<T> rv1(n);

   VT scale(0), anorm(0);
   for (i=0; i < n; i++)
   {
      l = i + 1;
      rv1[i] = scale * g;
      scale = VT( 0);
      s = T( 0);
      g = T( 0);

      if (i < m)
      {
         for (k=i; k < m; k++)
            scale += abs( a[k*n+i]);

         if (scale > epsilon( scale))
         {
            for (k=i; k < m; k++)
            {
               a[k*n+i] /= scale;
               tmp = a[k*n+i];
               s += tmp * tmp;
            }
            f = a[i*n+i];
            g = -sign( sqrt(s), f);
            h = f * g - s;
            a[i*n+i] = f - g;

            for (j=l; j < n; j++)
            {
               for (s=T(0),k=i; k < m; k++)
                  s += a[k*n+i] * a[k*n+j];
               f = s / h;
               for (k=i; k < m; k++)
                  a[k*n+j] += f * a[k*n+i];
            }
            for (k=i; k < m; k++)
               a[k*n+i] *= scale;
         }
      }
      
      w[i] = scale * g;
      scale = VT( 0);
      s = T(0);
      g = T(0);

      if (i < m && i != n-1)
      {
         for (k=l; k < n; k++)
            scale += abs(a[i*n+k]);

         if (scale > epsilon(scale))
         {
            for (k=l; k < n; k++)
            {
               a[i*n+k] /= scale;
               tmp = a[i*n+k];
               s += tmp * tmp;
            }
            f = a[i*n+l];
            g = -sign( sqrt(s), f);
            h = f * g - s;
            a[i*n+l] = f - g;

            for (k=l; k < n; k++)
               rv1[k] = a[i*n+k] / h;

            for (j=l; j < m; j++)
            {
               for (s=T(0),k=l; k < n; k++)
                  s += a[j*n+k] * a[i*n+k];
               for (k=l; k < n; k++)
                  a[j*n+k] += s * rv1[k];
            }
            for (k=l; k < n; k++)
               a[i*n+k] *= scale;
         }
      }
      anorm = techsoft::max( anorm, abs( w[i]) + abs( rv1[i]) );
   }

   for (i=n-1; ; i--)
   {
      if (i < n-1)
      {
         if (abs(g) > epsilon(g))
         {
            for (j=l; j < n; j++)
                v[j*n+i] = (a[i*n+j] / a[i*n+l]) / g;
            for (j=l; j < n; j++)
            {
               for (s=T(0),k=l; k < n; k++)
                  s += a[i*n+k] * v[k*n+j];

               for (k=l; k < n; k++)
                  v[k*n+j] += s * v[k*n+i];
            }
         }
         for (j=l; j < n; j++)
            v[i*n+j] = v[j*n+i] = T(0);
      }

      v[i*n+i] = T(1);
      g = rv1[i];
      l = i;
      if (i == 0)
         break;
   }

   for (i=techsoft::min(m,n)-1; ; i--)
   {
      l = i + 1;
      g = w[i];
      for (j=l; j < n; j++)
         a[i*n+j] = T(0);
      if (abs(g) > epsilon(g))
      {
         g = T(1) / g;
         for (j=l; j < n; j++)
         {
            for (s=T(0),k=l; k < m; k++)
               s += a[k*n+i] * a[k*n+j];

            f = (s / a[i*n+i]) * g;

            for (k=i; k < m; k++)
               a[k*n+j] += f * a[k*n+i];
         }
         for (j=i; j < m; j++)
            a[j*n+i] *= g;

      }
      else
      {
          for (j=i; j < m; j++)
              a[j*n+i] = T(0);
      }
      a[i*n+i] = a[i*n+i] + T(1);
      if (i == 0)
         break;
   }

   for (k=n-1; ; k--)
   {
      for (its=1; its <= 30; its++)
      {
         flag = 1;
         for (l=k; ; l--)
         {
            nm = l-1;
            if (abs( rv1[l]) < epsilon(rv1[l]))
            {
               flag = 0;
               break;
            }
            if (abs( w[nm]) < epsilon(w[nm]))
               break;
            if (l == 0)
               break;
         }
         if (flag)
         {
            c = T(0);
            s = T(1);
            for (i=l; i <= k; i++)
            {
               f = s * rv1[i];
               rv1[i] = c * rv1[i];
               if (abs(f) < epsilon(f))
                  break;
               g = w[i];
               h = sqrt( f * f + g * g);
               w[i] = h;
               h = VT(1) / h;
               c = g * h;
               s = -f * h;
               for (j=0; j < m; j++)
               {
                  y = a[j*n+nm];
                  z = a[j*n+i];
                  a[j*n+nm] = y * c + z * s;
                  a[j*n+i] = z * c - y * s;
               }
            }
         }
         z = w[k];
         if (l == k)                 
         {
            if (z.real() < VT(0))
            {
               w[k] = -z;
               for (j=0; j < n; j++)
                  v[j*n+k] = -v[j*n+k];
            }
            break;
         }

         if (its == 30)
            return false;

         x = w[l];
         nm = k - 1;
         y = w[nm];
         g = rv1[nm];
         h = rv1[k];
         f = ((y-z) * (y+z) + (g-h) * (g+h)) / (VT(2) * h * y);
         g = sqrt( f * f + T(1) * T(1));
         f = ((x-z) * (x+z) + h * ((y / (f+sign( g, f))) - h)) / x;
         c = s = T(1);                    
         for (j=l; j <= nm; j++)
         {
            i = j + 1;
            g = rv1[i];
            y = w[i];
            h = s * g;
            g = c * g;
            z = sqrt( f * f + h * h);
            rv1[j] = z;
            c = f / z;
            s = h / z;
            f = x * c + g * s;
            g = g*c - x*s;
            h = y * s;
            y *= c;
            for (jj=0; jj < n; jj++)
            {
               x = v[jj*n+j];
               z = v[jj*n+i];
               v[jj*n+j] = x * c + z * s;
               v[jj*n+i] = z * c - x * s;
            }
            z = sqrt( f * f + h * h);
            w[j] = z;
            if (abs(z) > epsilon(z))
            {
               z = VT(1) / z;
               c = f * z;
               s = h * z;
            }
            f = c * g + s * y;
            x = c * y - s * g;
            for (jj=0; jj < m; jj++)
            {
               y = a[jj*n+j];
               z = a[jj*n+i];
               a[jj*n+j] = y * c + z * s;
               a[jj*n+i] = z * c - y * s;
            }
         }
         rv1[l] = 0.0;
         rv1[k] = f;
         w[k] = x;
      }
      if (k == 0)
         break;
   }
   return true;
}

template <> inline 
bool matrix<complex<double> >::svd (matrix<complex<double> >& vc, valarray<complex<double> >& w)
{
   typedef value_type T;
   typedef double VT;

   size_t flag,i,its,j,jj,k,l,nm;
   T c,f,g(0),h,s,x,y,z,tmp;

   if (shared())
      clone();

   size_t m = mPtr->nrow;
   size_t n = mPtr->ncol;

   if (vc.rowno() != n || vc.colno() != n)
      vc.resize( n, n);
   if (w.size() != n)
      w.resize( n);

   T *a = &mPtr->val[0];
   T *v = &vc(0,0);
   valarray<T> rv1(n);

   VT scale(0), anorm(0);
   for (i=0; i < n; i++)
   {
      l = i + 1;
      rv1[i] = scale * g;
      scale = VT( 0);
      s = T( 0);
      g = T( 0);

      if (i < m)
      {
         for (k=i; k < m; k++)
            scale += abs( a[k*n+i]);

         if (scale > epsilon( scale))
         {
            for (k=i; k < m; k++)
            {
               a[k*n+i] /= scale;
               tmp = a[k*n+i];
               s += tmp * tmp;
            }
            f = a[i*n+i];
            g = -sign( sqrt(s), f);
            h = f * g - s;
            a[i*n+i] = f - g;

            for (j=l; j < n; j++)
            {
               for (s=T(0),k=i; k < m; k++)
                  s += a[k*n+i] * a[k*n+j];
               f = s / h;
               for (k=i; k < m; k++)
                  a[k*n+j] += f * a[k*n+i];
            }
            for (k=i; k < m; k++)
               a[k*n+i] *= scale;
         }
      }
      
      w[i] = scale * g;
      scale = VT( 0);
      s = T(0);
      g = T(0);

      if (i < m && i != n-1)
      {
         for (k=l; k < n; k++)
            scale += abs(a[i*n+k]);

         if (scale > epsilon(scale))
         {
            for (k=l; k < n; k++)
            {
               a[i*n+k] /= scale;
               tmp = a[i*n+k];
               s += tmp * tmp;
            }
            f = a[i*n+l];
            g = -sign( sqrt(s), f);
            h = f * g - s;
            a[i*n+l] = f - g;

            for (k=l; k < n; k++)
               rv1[k] = a[i*n+k] / h;

            for (j=l; j < m; j++)
            {
               for (s=T(0),k=l; k < n; k++)
                  s += a[j*n+k] * a[i*n+k];
               for (k=l; k < n; k++)
                  a[j*n+k] += s * rv1[k];
            }
            for (k=l; k < n; k++)
               a[i*n+k] *= scale;
         }
      }
      anorm = techsoft::max( anorm, abs( w[i]) + abs( rv1[i]) );
   }

   for (i=n-1; ; i--)
   {
      if (i < n-1)
      {
         if (abs(g) > epsilon(g))
         {
            for (j=l; j < n; j++)
                v[j*n+i] = (a[i*n+j] / a[i*n+l]) / g;
            for (j=l; j < n; j++)
            {
               for (s=T(0),k=l; k < n; k++)
                  s += a[i*n+k] * v[k*n+j];

               for (k=l; k < n; k++)
                  v[k*n+j] += s * v[k*n+i];
            }
         }
         for (j=l; j < n; j++)
            v[i*n+j] = v[j*n+i] = T(0);
      }

      v[i*n+i] = T(1);
      g = rv1[i];
      l = i;
      if (i == 0)
         break;
   }

   for (i=techsoft::min(m,n)-1; ; i--)
   {
      l = i + 1;
      g = w[i];
      for (j=l; j < n; j++)
         a[i*n+j] = T(0);
      if (abs(g) > epsilon(g))
      {
         g = T(1) / g;
         for (j=l; j < n; j++)
         {
            for (s=T(0),k=l; k < m; k++)
               s += a[k*n+i] * a[k*n+j];

            f = (s / a[i*n+i]) * g;

            for (k=i; k < m; k++)
               a[k*n+j] += f * a[k*n+i];
         }
         for (j=i; j < m; j++)
            a[j*n+i] *= g;

      }
      else
      {
          for (j=i; j < m; j++)
              a[j*n+i] = T(0);
      }
      a[i*n+i] = a[i*n+i] + T(1);
      if (i == 0)
         break;
   }

   for (k=n-1; ; k--)
   {
      for (its=1; its <= 30; its++)
      {
         flag = 1;
         for (l=k; ; l--)
         {
            nm = l-1;
            if (abs( rv1[l]) < epsilon(rv1[l]))
            {
               flag = 0;
               break;
            }
            if (abs( w[nm]) < epsilon(w[nm]))
               break;
            if (l == 0)
               break;
         }
         if (flag)
         {
            c = T(0);
            s = T(1);
            for (i=l; i <= k; i++)
            {
               f = s * rv1[i];
               rv1[i] = c * rv1[i];
               if (abs(f) < epsilon(f))
                  break;
               g = w[i];
               h = sqrt( f * f + g * g);
               w[i] = h;
               h = VT(1) / h;
               c = g * h;
               s = -f * h;
               for (j=0; j < m; j++)
               {
                  y = a[j*n+nm];
                  z = a[j*n+i];
                  a[j*n+nm] = y * c + z * s;
                  a[j*n+i] = z * c - y * s;
               }
            }
         }
         z = w[k];
         if (l == k)                 
         {
            if (z.real() < VT(0))
            {
               w[k] = -z;
               for (j=0; j < n; j++)
                  v[j*n+k] = -v[j*n+k];
            }
            break;
         }
         
         if (its == 30)
            return false;

         x = w[l];
         nm = k - 1;
         y = w[nm];
         g = rv1[nm];
         h = rv1[k];
         f = ((y-z) * (y+z) + (g-h) * (g+h)) / (VT(2) * h * y);
         g = sqrt( f * f + T(1) * T(1));
         f = ((x-z) * (x+z) + h * ((y / (f+sign( g, f))) - h)) / x;
         c = s = T(1);                    
         for (j=l; j <= nm; j++)
         {
            i = j + 1;
            g = rv1[i];
            y = w[i];
            h = s * g;
            g = c * g;
            z = sqrt( f * f + h * h);
            rv1[j] = z;
            c = f / z;
            s = h / z;
            f = x * c + g * s;
            g = g*c - x*s;
            h = y * s;
            y *= c;
            for (jj=0; jj < n; jj++)
            {
               x = v[jj*n+j];
               z = v[jj*n+i];
               v[jj*n+j] = x * c + z * s;
               v[jj*n+i] = z * c - x * s;
            }
            z = sqrt( f * f + h * h);
            w[j] = z;
            if (abs(z) > epsilon(z))
            {
               z = VT(1) / z;
               c = f * z;
               s = h * z;
            }
            f = c * g + s * y;
            x = c * y - s * g;
            for (jj=0; jj < m; jj++)
            {
               y = a[jj*n+j];
               z = a[jj*n+i];
               a[jj*n+j] = y * c + z * s;
               a[jj*n+i] = z * c - y * s;
            }
         }
         rv1[l] = 0.0;
         rv1[k] = f;
         w[k] = x;
      }
      if (k == 0)
         break;
   }
   return true;
}

template <> inline 
bool matrix<complex<long double> >::svd (matrix<complex<long double> >& vc, valarray<complex<long double> >& w)
{
   typedef value_type T;
   typedef long double VT;

   size_t flag,i,its,j,jj,k,l,nm;
   T c,f,g(0),h,s,x,y,z,tmp;

   if (shared())
      clone();

   size_t m = mPtr->nrow;
   size_t n = mPtr->ncol;

   if (vc.rowno() != n || vc.colno() != n)
      vc.resize( n, n);
   if (w.size() != n)
      w.resize( n);

   T *a = &mPtr->val[0];
   T *v = &vc(0,0);
   valarray<T> rv1(n);

   VT scale(0), anorm(0);
   for (i=0; i < n; i++)
   {
      l = i + 1;
      rv1[i] = scale * g;
      scale = VT( 0);
      s = T( 0);
      g = T( 0);

      if (i < m)
      {
         for (k=i; k < m; k++)
            scale += abs( a[k*n+i]);

         if (scale > epsilon( scale))
         {
            for (k=i; k < m; k++)
            {
               a[k*n+i] /= scale;
               tmp = a[k*n+i];
               s += tmp * tmp;
            }
            f = a[i*n+i];
            g = -sign( sqrt(s), f);
            h = f * g - s;
            a[i*n+i] = f - g;

            for (j=l; j < n; j++)
            {
               for (s=T(0),k=i; k < m; k++)
                  s += a[k*n+i] * a[k*n+j];
               f = s / h;
               for (k=i; k < m; k++)
                  a[k*n+j] += f * a[k*n+i];
            }
            for (k=i; k < m; k++)
               a[k*n+i] *= scale;
         }
      }
      
      w[i] = scale * g;
      scale = VT( 0);
      s = T(0);
      g = T(0);

      if (i < m && i != n-1)
      {
         for (k=l; k < n; k++)
            scale += abs(a[i*n+k]);

         if (scale > epsilon(scale))
         {
            for (k=l; k < n; k++)
            {
               a[i*n+k] /= scale;
               tmp = a[i*n+k];
               s += tmp * tmp;
            }
            f = a[i*n+l];
            g = -sign( sqrt(s), f);
            h = f * g - s;
            a[i*n+l] = f - g;

            for (k=l; k < n; k++)
               rv1[k] = a[i*n+k] / h;

            for (j=l; j < m; j++)
            {
               for (s=T(0),k=l; k < n; k++)
                  s += a[j*n+k] * a[i*n+k];
               for (k=l; k < n; k++)
                  a[j*n+k] += s * rv1[k];
            }
            for (k=l; k < n; k++)
               a[i*n+k] *= scale;
         }
      }
      anorm = techsoft::max( anorm, abs( w[i]) + abs( rv1[i]) );
   }

   for (i=n-1; ; i--)
   {
      if (i < n-1)
      {
         if (abs(g) > epsilon(g))
         {
            for (j=l; j < n; j++)
                v[j*n+i] = (a[i*n+j] / a[i*n+l]) / g;
            for (j=l; j < n; j++)
            {
               for (s=T(0),k=l; k < n; k++)
                  s += a[i*n+k] * v[k*n+j];

               for (k=l; k < n; k++)
                  v[k*n+j] += s * v[k*n+i];
            }
         }
         for (j=l; j < n; j++)
            v[i*n+j] = v[j*n+i] = T(0);
      }

      v[i*n+i] = T(1);
      g = rv1[i];
      l = i;
      if (i == 0)
         break;
   }

   for (i=techsoft::min(m,n)-1; ; i--)
   {
      l = i + 1;
      g = w[i];
      for (j=l; j < n; j++)
         a[i*n+j] = T(0);
      if (abs(g) > epsilon(g))
      {
         g = T(1) / g;
         for (j=l; j < n; j++)
         {
            for (s=T(0),k=l; k < m; k++)
               s += a[k*n+i] * a[k*n+j];

            f = (s / a[i*n+i]) * g;

            for (k=i; k < m; k++)
               a[k*n+j] += f * a[k*n+i];
         }
         for (j=i; j < m; j++)
            a[j*n+i] *= g;

      }
      else
      {
          for (j=i; j < m; j++)
              a[j*n+i] = T(0);
      }
      a[i*n+i] = a[i*n+i] + T(1);
      if (i == 0)
         break;
   }

   for (k=n-1; ; k--)
   {
      for (its=1; its <= 30; its++)
      {
         flag = 1;
         for (l=k; ; l--)
         {
            nm = l-1;
            if (abs( rv1[l]) < epsilon(rv1[l]))
            {
               flag = 0;
               break;
            }
            if (abs( w[nm]) < epsilon(w[nm]))
               break;
            if (l == 0)
               break;
         }
         if (flag)
         {
            c = T(0);
            s = T(1);
            for (i=l; i <= k; i++)
            {
               f = s * rv1[i];
               rv1[i] = c * rv1[i];
               if (abs(f) < epsilon(f))
                  break;
               g = w[i];
               h = sqrt( f * f + g * g);
               w[i] = h;
               h = VT(1) / h;
               c = g * h;
               s = -f * h;
               for (j=0; j < m; j++)
               {
                  y = a[j*n+nm];
                  z = a[j*n+i];
                  a[j*n+nm] = y * c + z * s;
                  a[j*n+i] = z * c - y * s;
               }
            }
         }
         z = w[k];
         if (l == k)                 
         {
            if (z.real() < VT(0))
            {
               w[k] = -z;
               for (j=0; j < n; j++)
                  v[j*n+k] = -v[j*n+k];
            }
            break;
         }

         if (its == 30)
            return false;

         x = w[l];
         nm = k - 1;
         y = w[nm];
         g = rv1[nm];
         h = rv1[k];
         f = ((y-z) * (y+z) + (g-h) * (g+h)) / (VT(2) * h * y);
         g = sqrt( f * f + T(1) * T(1));
         f = ((x-z) * (x+z) + h * ((y / (f+sign( g, f))) - h)) / x;
         c = s = T(1);                    
         for (j=l; j <= nm; j++)
         {
            i = j + 1;
            g = rv1[i];
            y = w[i];
            h = s * g;
            g = c * g;
            z = sqrt( f * f + h * h);
            rv1[j] = z;
            c = f / z;
            s = h / z;
            f = x * c + g * s;
            g = g*c - x*s;
            h = y * s;
            y *= c;
            for (jj=0; jj < n; jj++)
            {
               x = v[jj*n+j];
               z = v[jj*n+i];
               v[jj*n+j] = x * c + z * s;
               v[jj*n+i] = z * c - x * s;
            }
            z = sqrt( f * f + h * h);
            w[j] = z;
            if (abs(z) > epsilon(z))
            {
               z = VT(1) / z;
               c = f * z;
               s = h * z;
            }
            f = c * g + s * y;
            x = c * y - s * g;
            for (jj=0; jj < m; jj++)
            {
               y = a[jj*n+j];
               z = a[jj*n+i];
               a[jj*n+j] = y * c + z * s;
               a[jj*n+i] = z * c - y * s;
            }
         }
         rv1[l] = 0.0;
         rv1[k] = f;
         w[k] = x;
      }
      if (k == 0)
         break;
   }
   return true;
}

template <> inline
complex<float> matrix<complex<float> >::normF () const
{
   float nr = 0.0f;
   size_t n = size();

   for (size_t i=0; i < n; i++)
   {
      float e = abs( mPtr->val[i]);
      nr += e * e;
   }
   nr = sqrt( nr);

   return complex<float>( nr, 0.0f);
}

template <> inline
complex<double> matrix<complex<double> >::normF () const
{
   double nr = 0.0;
   size_t n = size();

   for (size_t i=0; i < n; i++)
   {
      double e = abs( mPtr->val[i]);
      nr += e * e;
   }
   nr = sqrt( nr);

   return complex<double>( nr, 0.0);
}

template <> inline
complex<long double> matrix<complex<long double> >::normF () const
{
   long double nr = 0.0;
   size_t n = size();

   for (size_t i=0; i < n; i++)
   {
      long double e = abs( mPtr->val[i]);
      nr += e * e;
   }
   nr = sqrt( nr);

   return complex<long double>( nr, 0.0);
}

template <> inline
complex<float> matrix<complex<float> >::norm1 () const
{
   typedef value_type T;
   float nr = 0.0;
   size_t m = mPtr->nrow;
   size_t n = mPtr->ncol;
   const T *pm = &mPtr->val[0];

   for (size_t j=0; j < n; j++)
   {
      float s = 0.0f;
      for (size_t i=0; i < m; i++)
         s += abs( pm[i*n+j]);
      nr = techsoft::max( nr, s);
   }
   return T( nr, 0.0f);
}

template <> inline
complex<double> matrix<complex<double> >::norm1 () const
{
   typedef value_type T;
   double nr = 0.0;
   size_t m = mPtr->nrow;
   size_t n = mPtr->ncol;
   const T *pm = &mPtr->val[0];

   for (size_t j=0; j < n; j++)
   {
      double s = 0.0;
      for (size_t i=0; i < m; i++)
         s += abs( pm[i*n+j]);
      nr = techsoft::max( nr, s);
   }
   return T( nr, 0.0);
}

template <> inline
complex<long double> matrix<complex<long double> >::norm1 () const
{
   typedef value_type T;
   long double nr = 0.0;
   size_t m = mPtr->nrow;
   size_t n = mPtr->ncol;
   const T *pm = &mPtr->val[0];

   for (size_t j=0; j < n; j++)
   {
      long double s = 0.0;
      for (size_t i=0; i < m; i++)
         s += abs( pm[i*n+j]);
      nr = techsoft::max( nr, s);
   }
   return T( nr, 0.0);
}

template <>  inline
complex<float> matrix<complex<float> >::norm2 () const
{
   matrix<complex<float> > m(*this);
   matrix<complex<float> > v(mPtr->ncol,mPtr->ncol);

   size_t i, n = mPtr->ncol;
   float nr = 0.0f;
   valarray<complex<float> > w( n);

   if (m.svd( v, w))
      for (i=0; i < n; i++)
         nr = techsoft::max( nr, abs( w[i]));

   return complex<float>( nr, 0.0f);
}

template <>  inline
complex<double> matrix<complex<double> >::norm2 () const
{
   matrix<complex<double> > m(*this);
   matrix<complex<double> > v(mPtr->ncol,mPtr->ncol);

   size_t i, n = mPtr->ncol;
   double nr = 0.0;
   valarray<complex<double> > w( n);

   if (m.svd( v, w))
      for (i=0; i < n; i++)
         nr = techsoft::max( nr, abs( w[i]));

   return complex<double>( nr, 0.0);
}

template <>  inline
complex<long double> matrix<complex<long double> >::norm2 () const
{
   matrix<complex<long double> > m(*this);
   matrix<complex<long double> > v(mPtr->ncol,mPtr->ncol);

   size_t i, n = mPtr->ncol;
   long double nr = 0.0;
   valarray<complex<long double> > w( n);

   if (m.svd( v, w))
      for (i=0; i < n; i++)
         nr = techsoft::max( nr, abs( w[i]));

   return complex<long double>( nr, 0.0);
}

template <> inline
complex<float> matrix<complex<float> >::normI () const
{
   typedef value_type T;
   float nr = 0.0f;
   size_t m = mPtr->nrow;
   size_t n = mPtr->ncol;
   const T *pm = &mPtr->val[0];

   for (size_t i=0; i < m; i++)
   {
      float s = 0.0f;
      for (size_t j=0; j < n; j++)
         s += abs( pm[i*n+j]);
      nr = techsoft::max( nr, s);
   }
   return T( nr, 0.0f);
}

template <> inline
complex<double> matrix<complex<double> >::normI () const
{
   typedef value_type Ty;
   double nr = 0.0;
   size_t m = mPtr->nrow;
   size_t n = mPtr->ncol;
   const Ty *pm = &mPtr->val[0];

   for (size_t i=0; i < m; i++)
   {
      double s = 0.0;
      for (size_t j=0; j < n; j++)
         s += abs( pm[i*n+j]);
      nr = techsoft::max( nr, s);
   }
   return Ty( nr, 0.0);
}

template <> inline
complex<long double> matrix<complex<long double> >::normI () const
{
   typedef value_type Ty;
   long double nr = 0.0;
   size_t m = mPtr->nrow;
   size_t n = mPtr->ncol;
   const Ty *pm = &mPtr->val[0];

   for (size_t i=0; i < m; i++)
   {
      long double s = 0.0;
      for (size_t j=0; j < n; j++)
         s += abs( pm[i*n+j]);
      nr = techsoft::max( nr, s);
   }
   return Ty( nr, 0.0);
}


} // namespace techsoft

#endif // _MATRIX_CC
