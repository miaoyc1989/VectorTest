/*
 * Copyright (c) 2008-2011 Zhang Ming (M. Zhang), zmjerry@163.com
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 2 or any later version.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details. A copy of the GNU General Public License is available at:
 * http://www.fsf.org/licensing/licenses
 */


/*****************************************************************************
 *                               matrix-impl.h
 *
 * Implementation for Matrix class.
 *
 * Zhang Ming, 2010-01 (revised 2010-12), Xi'an Jiaotong University.
 *****************************************************************************/


/**
 * initialize
 */
template <typename Type>
void Matrix<Type>::init( int rows, int columns )
{
	_nRow = rows;
	_nColumn = columns;
	_nTotal = _nRow * _nColumn;

	_pMatrixData = new Type[_nTotal];
	prow0 = new Type*[_nRow];
	prow1 = new Type*[_nRow];

	assert( _pMatrixData != NULL );
	assert( prow0 != NULL );
	assert( prow1 != NULL );

	Type *p = _pMatrixData;
	pv1 = _pMatrixData - 1;

	for( int i=0; i<_nRow; ++i )
	{
		prow0[i] = p;
		prow1[i] = p-1;
		p += _nColumn;
	}
	prow1--;

	memset(_pMatrixData, 0, GetSize() * (sizeof (Type)));
}


/**
 * copy matrix from normal array
 */
template <typename Type>
inline void Matrix<Type>::copyFromArray( const Type *v )
{
	memcpy(_pMatrixData, v, _nTotal*sizeof(Type));
}


/**
 * set matrix by a scalar
 */
template <typename Type>
inline void Matrix<Type>::setByScalar( const Type &x )
{
	for( long i=0; i<_nTotal; ++i )
	{
		_pMatrixData[i] = x;
	}
}


/**
 * destroy the matrix
 */
template <typename Type>
void Matrix<Type>::destroy()
{
	if( _pMatrixData == NULL )
		return ;
	else
		delete []_pMatrixData;

	if( prow0 != NULL )
		delete []prow0;

	prow1++;
	if( prow1 != NULL )
		delete []prow1;
}


/**
 * constructors and destructor
 */
template <typename Type>
Matrix<Type>::Matrix()
: _pMatrixData(0), pv1(0), prow0(0), prow1(0), _nRow(0), _nColumn(0), _nTotal(0)
{
}

template <typename Type>
Matrix<Type>::Matrix( const Matrix<Type> &A )
{
	init( A._nRow, A._nColumn );
	copyFromArray( A._pMatrixData );
}

template <typename Type>
Matrix<Type>::Matrix( int rows, int columns )
{
	init( rows,columns );
}

template <typename Type>
Matrix<Type>::Matrix( int rows, int columns, const Type &x )
{
	init( rows,columns );
//	setByScalar(x);
}

template <typename Type>
Matrix<Type>::Matrix( int rows, int columns, const Type *arrays )
{
	init( rows,columns );
	copyFromArray( arrays );
}

template <typename Type>
Matrix<Type>::~Matrix()
{
	destroy();
}


/**
 * overload evaluate operator = from matrix to matrix
 */
template <typename Type>
Matrix<Type>& Matrix<Type>::operator=( const Matrix<Type> &A )
{
	if( _pMatrixData == A._pMatrixData )
		return *this;

	if( _nRow == A._nRow && _nColumn == A._nColumn )
		copyFromArray( A._pMatrixData );
	else
	{
		destroy();
		init( A._nRow, A._nColumn );
		copyFromArray( A._pMatrixData );
	}

	return *this;
}


/**
 * overload evaluate operator = from scalar to matrix
 */
template <typename Type>
inline Matrix<Type>& Matrix<Type>::operator=( const Type &x )
{
	setByScalar( x );

	return *this;
}


/**
 * overload operator [] for 0-offset access
 */
template <typename Type>
inline Type* Matrix<Type>::operator[]( int i )
{
#ifdef BOUNDS_CHECK
	assert( 0 <= i );
	assert( i < _nRow );
#endif

	return prow0[i];
}

template <typename Type>
inline const Type* Matrix<Type>::operator[]( int i ) const
{
#ifdef BOUNDS_CHECK
	assert( 0 <= i );
	assert( i < _nRow );
#endif

	return prow0[i];
}


/**
 * overload operator () for 1-offset access
 */
template <typename Type>
inline Type& Matrix<Type>::operator()( int row, int column )
{
#ifdef BOUNDS_CHECK
	assert( 1 <= row );
	assert( row <= _nRow ) ;
	assert( 1 <= column);
	assert( column <= _nColumn );
#endif

	return  prow1[row][column];
}

template <typename Type>
inline const Type& Matrix<Type>::operator()( int row, int column ) const
{
#ifdef BOUNDS_CHECK
	assert( 1 <= row );
	assert( row <= _nRow ) ;
	assert( 1 <= column);
	assert( column <= _nColumn );
#endif

	return  prow1[row][column];
}


/**
 * type conversion functions
 */
template <typename Type>
inline Matrix<Type>::operator Type*()
{
	return _pMatrixData;
}

template <typename Type>
inline Matrix<Type>::operator const Type*() const
{
	return _pMatrixData;
}


/**
 * get the matrix's size
 */
template <typename Type>
inline long Matrix<Type>::GetSize() const
{
	return _nTotal;
}


/**
 * get the matrix's dimension
 */
template <typename Type>
int Matrix<Type>::dim( int dimension ) const
{
#ifdef BOUNDS_CHECK
	assert( dimension >= 1);
	assert( dimension <= 2);
#endif

	if( dimension == 1 )
		return _nRow;
	else if( dimension == 2 )
		return _nColumn;
	else
		return 0;
}

template <typename Type>
inline int Matrix<Type>::GetRows() const
{
    return _nRow;
}

template <typename Type>
inline int Matrix<Type>::GetCols() const
{
    return _nColumn;
}


/**
 * reallocate matrix's size
 */
template <typename Type>
Matrix<Type>& Matrix<Type>::resize( int rows, int columns )
{
	if(  rows == _nRow && columns == _nColumn )
		return *this;

	destroy();
	init( rows, columns );
	return *this;
}

/**
 * set the matrix's row vector
 */
template <typename Type>
void Matrix<Type>::setRow( const Vector<Type> &v, int row )
{
#ifdef BOUNDS_CHECK
	assert( row >= 0 );
	assert( row < _nRow );
	assert( v.dim() == _nColumn );
#endif

	for( int j=0; j<_nColumn; ++j )
		prow0[row][j] = v[j];
}


/**
 * set the matrix's column vector
 */
template <typename Type>
void Matrix<Type>::setColumn( const Vector<Type> &v, int column )
{
#ifdef BOUNDS_CHECK
	assert( column >= 0 );
	assert( column < _nColumn );
	assert( v.dim() == _nRow );
#endif

	for( int i=0; i<_nRow; ++i )
		prow0[i][column] = v[i];
}


/**
 * compound assignment operators +=
 */
template <typename Type>
Matrix<Type>& Matrix<Type>::operator+=( const Type &x )
{
    Type **rowPtr = prow0;
    Type *colPtr = 0;

    for( int i=0; i<_nRow; ++i )
    {
        colPtr = *rowPtr++;
        for( int j=0; j<_nColumn; ++j )
            *colPtr++ += x;
    }

	return *this;
}

template <typename Type>
Matrix<Type>& Matrix<Type>::operator+=( const Matrix<Type> &rhs )
{
    assert( _nRow == rhs.GetRows() );
    assert( _nColumn == rhs.GetCols() );

    Type **rowPtrL = prow0;
    Type *colPtrL = 0;
    Type **rowPtrR = rhs.prow0;
    const Type *colPtrR = 0;

    for( int i=0; i<_nRow; ++i )
    {
        colPtrL = *rowPtrL++;
        colPtrR = *rowPtrR++;
        for( int j=0; j<_nColumn; ++j )
            *colPtrL++ += *colPtrR++;
    }

	return *this;
}


/**
 * compound assignment operators -=
 */
template <typename Type>
Matrix<Type>& Matrix<Type>::operator-=( const Type &x )
{
    Type **rowPtr = prow0;
    Type *colPtr = 0;

    for( int i=0; i<_nRow; ++i )
    {
        colPtr = *rowPtr++;
        for( int j=0; j<_nColumn; ++j )
            *colPtr++ -= x;
    }

	return *this;
}

template <typename Type>
Matrix<Type>& Matrix<Type>::operator-=( const Matrix<Type> &rhs )
{
    assert( _nRow == rhs.GetRows() );
    assert( _nColumn == rhs.GetCols() );

    Type **rowPtrL = prow0;
    Type *colPtrL = 0;
    Type **rowPtrR = rhs.prow0;
    const Type *colPtrR = 0;

    for( int i=0; i<_nRow; ++i )
    {
        colPtrL = *rowPtrL++;
        colPtrR = *rowPtrR++;
        for( int j=0; j<_nColumn; ++j )
            *colPtrL++ -= *colPtrR++;
    }

	return *this;
}


/**
 * compound assignment operators *=
 */
template <typename Type>
Matrix<Type>& Matrix<Type>::operator*=( const Type &x )
{
    Type **rowPtr = prow0;
    Type *colPtr = 0;

    for( int i=0; i<_nRow; ++i )
    {
        colPtr = *rowPtr++;
        for( int j=0; j<_nColumn; ++j )
            *colPtr++ *= x;
    }

	return *this;
}

// WARNING: this is element-by-element multiplication
template <typename Type>
Matrix<Type>& Matrix<Type>::operator*=( const Matrix<Type> &rhs )
{
    assert( _nRow == rhs.GetRows() );
    assert( _nColumn == rhs.GetCols() );

    Type **rowPtrL = prow0;
    Type *colPtrL = 0;
    Type **rowPtrR = rhs.prow0;
    const Type *colPtrR = 0;

    for( int i=0; i<_nRow; ++i )
    {
        colPtrL = *rowPtrL++;
        colPtrR = *rowPtrR++;
        for( int j=0; j<_nColumn; ++j )
            *colPtrL++ *= *colPtrR++;
    }

	return *this;
}


/**
 * compound assignment operators /=
 */
template <typename Type>
Matrix<Type>& Matrix<Type>::operator/=( const Type &x )
{
    Type **rowPtr = prow0;
    Type *colPtr = 0;

    for( int i=0; i<_nRow; ++i )
    {
        colPtr = *rowPtr++;
        for( int j=0; j<_nColumn; ++j )
            *colPtr++ /= x;
    }

	return *this;
}

// WARNING: this is element-by-element division
template <typename Type>
Matrix<Type>& Matrix<Type>::operator/=( const Matrix<Type> &rhs )
{
    assert( _nRow == rhs.GetRows() );
    assert( _nColumn == rhs.GetCols() );

    Type **rowPtrL = prow0;
    Type *colPtrL = 0;
    Type **rowPtrR = rhs.prow0;
    const Type *colPtrR = 0;

    for( int i=0; i<_nRow; ++i )
    {
        colPtrL = *rowPtrL++;
        colPtrR = *rowPtrR++;
        for( int j=0; j<_nColumn; ++j )
            *colPtrL++ /= *colPtrR++;
    }

	return *this;
}


/**
 * Overload the output stream function.
 */
template <typename Type>
ostream& operator<<( ostream &out, const Matrix<Type> &A )
{
	int rows = A.GetRows();
	int columns = A.GetCols();

	out << "size: " << rows << " by " << columns << "\n";
	for( int i=0; i<rows; ++i )
	{
		for( int j=0; j<columns; ++j )
			out << A[i][j] << "\t";
		out << "\n";
	}

	return out;
}


/**
 * Overload the intput stream function.
 */
template <typename Type>
istream& operator>>( istream &in, Matrix<Type> &A )
{
	int rows, columns;
	in >> rows >> columns;

	if( !( rows == A.GetRows() && columns == A.GetCols() ) )
		A.resize( rows, columns );

	for( int i=0; i<rows; ++i )
		for( int j=0; j<columns; ++j )
			in >> A[i][j];

	return in;
}


/**
 * get negative matrix
 */
template<typename Type>
Matrix<Type> operator-( const Matrix<Type> &A )
{
	int rows = A.GetRows();
	int columns = A.GetCols();

	Matrix<Type> tmp( rows, columns );
	for( int i=0; i<rows; ++i )
		for( int j=0; j<columns; ++j )
			tmp[i][j] = -A[i][j];

	return tmp;
}


/**
 * matrix-scalar addition
 */
template<typename Type>
inline Matrix<Type> operator+( const Matrix<Type> &A, const Type &x )
{
	Matrix<Type> tmp( A );
	return tmp += x;
}

template<typename Type>
inline Matrix<Type> operator+( const Type &x, const Matrix<Type> &A )
{
	return A + x;
}


/**
 * matrix-matrix addition
 */
template<typename Type>
inline Matrix<Type> operator+( const Matrix<Type> &A1, const Matrix<Type> &A2 )
{
	Matrix<Type> tmp( A1 );
	return tmp += A2;
}


/**
 * matrix-scalar subtraction
 */
template<typename Type>
inline Matrix<Type> operator-( const Matrix<Type> &A, const Type &x )
{
	Matrix<Type> tmp( A );
	return tmp -= x;
}

template<typename Type>
inline Matrix<Type> operator-( const Type &x, const Matrix<Type> &A )
{
	Matrix<Type> tmp( A );
	return -tmp += x;
}


/**
 * matrix-matrix subtraction
 */
template<typename Type>
inline Matrix<Type> operator-( const Matrix<Type> &A1, const Matrix<Type> &A2 )
{
	Matrix<Type> tmp( A1 );
	return tmp -= A2;
}


/**
 * matrix-scaling multiplication
 */
template <typename Type>
inline Matrix<Type> operator*( const Matrix<Type> &A, const Type &x )
{
	Matrix<Type> tmp( A );
	return tmp *= x;
}

template <typename Type>
inline Matrix<Type> operator*( const Type &x, const Matrix<Type> &A )
{
	return A * x;
}


/**
 * matrix-matrix multiplication
 */
template <typename Type>
Matrix<Type> operator*( const Matrix<Type> &A1, const Matrix<Type> &A2 )
{
	assert( A1.GetCols() == A2.GetRows() );

	int rows = A1.GetRows();
	int columns = A2.GetCols();

	Matrix<Type> tmp( rows, columns );

    mult( A1, A2, tmp );

	return tmp;
}


/**
 * matrix-vector multiplication
 */
template <typename Type>
Vector<Type> operator*( const Matrix<Type> &A, const Vector<Type> &b )
{
	assert( A.GetCols() == b.dim() );

	int rows = A.GetRows();

	Vector<Type> tmp(rows);

    mult( A, b, tmp );

	return tmp;
}


/**
 * matrix-scalar division
 */
template <typename Type>
inline Matrix<Type> operator/( const Matrix<Type> &A, const Type &x )
{
	Matrix<Type> tmp( A );
	return tmp /= x;
}

template <typename Type>
Matrix<Type> operator/( const Type &x, const Matrix<Type> &A )
{
	int rows = A.GetRows();
	int clumns = A.GetCols();

	Matrix<Type> tmp( rows,clumns );
	for( int i=0; i<rows; ++i )
		for( int j=0; j<clumns; ++j )
			tmp[i][j] = x / A[i][j];

	return tmp;
}


/**
 * This is an optimized version of matrix multiplication,
 * where the destination matrix has already been allocated.
 */
template <typename Type>
Matrix<Type>& mult( const Matrix<Type> &A, const Matrix<Type> &B,
                    Matrix<Type> &C )
{
    int M = A.GetRows();
    int N = B.GetCols();
    int K = A.GetCols();

    assert( B.GetRows() == K );

    C.resize( M, N );
    Type        sum;
    const Type  *pRow,
                *pCol;

    for( int i=0; i<M; i++ )
        for( int j=0; j<N; ++j )
        {
            pRow  = &A[i][0];
            pCol  = &B[0][j];
            sum = 0;

            for( int k=0; k<K; ++k )
            {
                sum += (*pRow) * (*pCol);
                pRow++;
                pCol += N;
            }
            C[i][j] = sum;
        }
    return C;
}


/**
 * This is an optimized version of matrix and vector multiplication,
 * where the destination vector has already been allocated.
 */
template <typename Type>
Vector<Type>& mult( const Matrix<Type> &A, const Vector<Type> &b,
                    Vector<Type> &c )
{
    int M = A.GetRows();
    int N = A.GetCols();

    assert( b.GetSize() == N );

    c.resize( M );
    Type        sum;
    const Type  *pRow,
                *pCol;

    for( int i=0; i<M; i++ )
    {
        pRow  = &A[i][0];
        pCol  = &b[0];
        sum = 0;

        for( int j=0; j<N; ++j )
        {
            sum += (*pRow) * (*pCol);
            pRow++;
            pCol++;
        }
        c[i] = sum;
    }
    return c;
}


/**
 * matrix-matrix elementwise multiplication
 */
template<typename Type>
inline Matrix<Type> elemMult( const Matrix<Type> &A1, const Matrix<Type> &A2 )
{
	Matrix<Type> tmp( A1 );
	return tmp *= A2;
}

template <typename Type>
inline Matrix<Type>& elemMultEq( Matrix<Type> &A1, const Matrix<Type> &A2 )
{
    return A1 *= A2;
}


/**
 * matrix-matrix elementwise division
 */
template <typename Type>
inline Matrix<Type> elemDivd( const Matrix<Type> &A1, const Matrix<Type> &A2 )
{
	Matrix<Type> tmp( A1 );
	return tmp /= A2;
}

template <typename Type>
inline Matrix<Type>& elemDivdEq( Matrix<Type> &A1, const Matrix<Type> &A2 )
{
    return A1 /= A2;
}


/**
 * matrix tranpose
 */
template <typename Type>
Matrix<Type> trT( const Matrix<Type> &A )
{
	int rows = A.GetCols();
	int clumns = A.GetRows();

	Matrix<Type> tmp( rows, clumns );
	for( int i=0; i<rows; ++i )
		for( int j=0; j<clumns; ++j )
			tmp[i][j] = A[j][i];

	return tmp;
}


/**
 * matrix conjugate tranpose
 */
template <typename Type>
Matrix<Type> trH( const Matrix<Type> &A )
{
	int rows = A.GetCols();
	int clumns = A.GetRows();

	Matrix<Type> tmp( rows, clumns );
	for( int i=0; i<rows; ++i )
		for( int j=0; j<clumns; ++j )
			tmp[i][j] = conj(A[j][i]);

	return tmp;
}


/**
 * matrix-matrix tranpose multiplication: A^T * B.
 */
template <typename Type>
Matrix<Type> trMult( const Matrix<Type> &A1, const Matrix<Type> &A2 )
{
	assert( A1.GetRows() == A2.GetRows() );

	int rows = A1.GetCols();
	int columns = A2.GetCols();
	int K = A1.GetRows();

	Matrix<Type> tmp( rows, columns );
    for( int i=0; i<rows; ++i )
		for( int j=0; j<columns; ++j )
			for( int k=0; k<K; ++k )
			   tmp[i][j] += A1[k][i] * A2[k][j];

	return tmp;
}


/**
 * matrix-vector tranpose multiplication: A^T * b.
 */
template <typename Type>
Vector<Type> trMult( const Matrix<Type> &A, const Vector<Type> &v )
{
	assert( A.GetRows() == v.dim() );

	int rows = A.GetRows();
	int columns = A.GetCols();

	Vector<Type> tmp( columns );
    for( int i=0; i<columns; ++i )
		for( int j=0; j<rows; ++j )
			tmp[i] += A[j][i] * v[j];

	return tmp;
}


/**
 * matrix-matrix tranpose multiplication: A * B^T.
 */
template <typename Type>
Matrix<Type> multTr( const Matrix<Type> &A1, const Matrix<Type> &A2 )
{
	assert( A1.GetCols() == A2.GetCols() );

	int rows = A1.GetRows();
	int columns = A2.GetRows();
	int K = A1.GetCols();

	Matrix<Type> tmp( rows, columns );
    for( int i=0; i<rows; ++i )
		for( int j=0; j<columns; ++j )
			for( int k=0; k<K; ++k )
			   tmp[i][j] += A1[i][k] * A2[j][k];

	return tmp;
}


/**
 * vector-vector tranpose multiplication: a * b^T.
 */
template <typename Type>
Matrix<Type> multTr( const Vector<Type> &a, const Vector<Type> &b )
{
	int rows = a.dim();
	int columns = b.dim();

	Matrix<Type> tmp( rows, columns );
	for( int i=0; i<rows; ++i )
		for( int j=0; j<columns; ++j )
			tmp[i][j] = a[i]*b[j];

	return tmp;
}


/**
 * matrix-matrix tranpose multiplication: A^H * B.
 */
template <typename Type>
Matrix<complex<Type> > trMult( const Matrix<complex<Type> > &A1,
                               const Matrix<complex<Type> > &A2 )
{
	assert( A1.GetRows() == A2.GetRows() );

	int rows = A1.GetCols();
	int columns = A2.GetCols();
	int K = A1.GetRows();

	Matrix<complex<Type> > tmp( rows, columns );
    for( int i=0; i<rows; ++i )
		for( int j=0; j<columns; ++j )
			for( int k=0; k<K; ++k )
			   tmp[i][j] += conj(A1[k][i]) * A2[k][j];

	return tmp;
}


/**
 * matrix-vector tranpose multiplication: A^H * b.
 */
template <typename Type>
Vector<complex<Type> > trMult( const Matrix<complex<Type> > &A,
                               const Vector<complex<Type> > &v )
{
	assert( A.GetRows() == v.dim() );

	int rows = A.GetRows();
	int columns = A.GetCols();

	Vector<complex<Type> > tmp( columns );
    for( int i=0; i<columns; ++i )
		for( int j=0; j<rows; ++j )
			tmp[i] += conj(A[j][i]) * v[j];

	return tmp;
}


/**
 * matrix-matrix tranpose multiplication: A * B^H.
 */
template <typename Type>
Matrix<complex<Type> > multTr( const Matrix<complex<Type> > &A1,
                               const Matrix<complex<Type> > &A2 )
{
	assert( A1.GetCols() == A2.GetCols() );

	int rows = A1.GetRows();
	int columns = A2.GetRows();
	int K = A1.GetCols();

	Matrix<complex<Type> > tmp( rows, columns );
    for( int i=0; i<rows; ++i )
		for( int j=0; j<columns; ++j )
			for( int k=0; k<K; ++k )
			   tmp[i][j] += A1[i][k] * conj(A2[j][k]);

	return tmp;
}


/**
 * vector-vector tranpose multiplication: a * b^H.
 */
template <typename Type>
Matrix<complex<Type> > multTr( const Vector<complex<Type> > &a,
                               const Vector<complex<Type> > &b )
{
	int rows = a.dim();
	int columns = b.dim();

	Matrix<complex<Type> > tmp( rows, columns );
	for( int i=0; i<rows; ++i )
		for( int j=0; j<columns; ++j )
			tmp[i][j] = a[i]*conj(b[j]);

	return tmp;
}


/**
 * Generate the identity matrix.
 */
template <typename Type>
Matrix<Type> eye( int N, const Type &x )
{
    Matrix<Type> tmp( N, N );
	for( int i=0; i<N; ++i )
		tmp[i][i] = x;

	return tmp;
}


/**
 * Get the diagonal entries of matrix.
 */
template <typename Type>
Vector<Type> diag( const Matrix<Type> &A )
{
	int _nColumn = A.GetRows();
	if( _nColumn > A.GetCols() )
		_nColumn = A.GetCols();

	Vector<Type> tmp( _nColumn );
	for( int i=0; i<_nColumn; ++i )
		tmp[i] = A[i][i];

	return tmp;
}


/**
 * Generate the diagonal of matrix by given its diagonal elements.
 */
template <typename Type>
Matrix<Type> diag( const Vector<Type> &d )
{
	int N = d.GetSize();

	Matrix<Type> tmp( N, N );
	for( int i=0; i<N; ++i )
		tmp[i][i] = d[i];

	return tmp;
}


/**
 * Compute Frobenius norm of matrix.
 */
template <typename Type>
Type norm( const Matrix<Type> &A )
{
	int m = A.GetRows();
	int n = A.GetCols();

	Type sum = 0;
	for( int i=1; i<=m; ++i )
		for( int j=1; j<=n; ++j )
            sum += A(i,j) * A(i,j);

	return sqrt(sum);
}

template <typename Type>
Type norm( const Matrix<complex<Type> > &A )
{
	int m = A.GetRows();
	int n = A.GetCols();

	Type sum = 0;
	for( int i=1; i<=m; ++i )
		for( int j=1; j<=n; ++j )
            sum += norm(A(i,j));

	return sqrt(sum);
}


/**
 * Swap two matrixes.
 */
template <typename Type> void swap( Matrix<Type> &lhs, Matrix<Type> &rhs )
{
    int m = lhs.GetRows();
	int n = lhs.GetCols();

	assert( m == rhs.GetRows() );
	assert( n == rhs.GetCols() );

	for( int i=1; i<=m; ++i )
		for( int j=1; j<=n; ++j )
            swap( lhs(i,j), rhs(i,j) );
}


/**
 * Matrix's column vecotrs sum.
 */
template <typename Type>
Vector<Type> sum( const Matrix<Type> &A )
{
	int m = A.GetRows();
	int n = A.GetCols();
	Vector<Type> sum(n);

	for( int j=1; j<=n; ++j )
		for( int i=1; i<=m; ++i )
            sum(j) += A(i,j);

	return sum;
}


/**
 * Minimum of matrix's column vecotrs.
 */
template <typename Type>
Vector<Type> min( const Matrix<Type> &A )
{
	int m = A.GetRows();
	int n = A.GetCols();
	Vector<Type> sum(n);

	for( int j=1; j<=n; ++j )
	{
	    Type tmp = A(1,j);
        for( int i=2; i<m; ++i )
            if( tmp > A(i,j) )
                tmp = A(i,j);
        sum(j) = tmp;
	}

	return sum;
}


/**
 * Maximum of matrix's column vecotrs.
 */
template <typename Type>
Vector<Type> max( const Matrix<Type> &A )
{
	int m = A.GetRows();
	int n = A.GetCols();
	Vector<Type> sum(n);

	for( int j=1; j<=n; ++j )
	{
	    Type tmp = A(1,j);
        for( int i=2; i<m; ++i )
            if( tmp < A(i,j) )
                tmp = A(i,j);
        sum(j) = tmp;
	}

	return sum;
}


/**
 * Matrix's column vecotrs mean.
 */
template <typename Type>
inline Vector<Type> mean( const Matrix<Type> &A )
{
	return sum(A) / Type(A.GetRows());
}


/**
 * Convert real matrix to complex matrix.
 */
template <typename Type>
Matrix<complex<Type> > complexMatrix( const Matrix<Type> &rA )
{
	int rows = rA.GetRows();
	int columns = rA.GetCols();

    Matrix<complex<Type> > cA( rows, columns );
    for( int i=0; i<rows; ++i )
		for( int j=0; j<columns; ++j )
			cA[i][j] = rA[i][j];

    return cA;
}

template <typename Type>
Matrix<complex<Type> > complexMatrix( const Matrix<Type> &mR,
                                      const Matrix<Type> &mI )
{
	int rows = mR.GetRows();
	int columns = mR.GetCols();

	assert( rows == mI.GetRows() );
	assert( columns == mI.GetCols() );

    Matrix<complex<Type> > cA( rows, columns );
    for( int i=0; i<rows; ++i )
		for( int j=0; j<columns; ++j )
			cA[i][j] = complex<Type>( mR[i][j], mI[i][j] );

    return cA;
}


/**
 * Get magnitude of a complex matrix.
 */
template <typename Type>
Matrix<Type> abs( const Matrix<complex<Type> > &A )
{
    int m = A.GetRows(),
        n = A.GetCols();
    Matrix<Type> tmp( m, n );

    for( int i=0; i<m; ++i )
        for( int j=0; j<n; ++j )
            tmp[i][j] = abs( A[i][j] );

    return tmp;
}


/**
 * Get angle of a complex matrix.
 */
template <typename Type>
Matrix<Type> arg( const Matrix<complex<Type> > &A )
{
    int m = A.GetRows(),
        n = A.GetCols();
    Matrix<Type> tmp( m, n );

    for( int i=0; i<m; ++i )
        for( int j=0; j<n; ++j )
            tmp[i][j] = arg( A[i][j] );

    return tmp;
}


/**
 * Get real part of a complex matrix.
 */
template <typename Type>
Matrix<Type> real( const Matrix<complex<Type> > &A )
{
    int m = A.GetRows(),
        n = A.GetCols();
    Matrix<Type> tmp( m, n );

    for( int i=0; i<m; ++i )
        for( int j=0; j<n; ++j )
            tmp[i][j] = A[i][j].real();

    return tmp;
}


/**
 * Get imaginary part of a complex matrix.
 */
template <typename Type>
Matrix<Type> imag( const Matrix<complex<Type> > &A )
{
    int m = A.GetRows(),
        n = A.GetCols();
    Matrix<Type> tmp( m, n );

    for( int i=0; i<m; ++i )
        for( int j=0; j<n; ++j )
            tmp[i][j] = A[i][j].imag();

    return tmp;
}


//==============================================================================
//
//  Matrix::SetSize
//
//  Description:
//      Sets the size of the matrix.
//
//  Parameters:
//      source  - object to match size of
//      newSize - size of square matrix
//      newRows - number of rows
//      newCols - number of columns
//
//==============================================================================
 template <class Type>
 void Matrix<Type>::SetSize (const Matrix<Type>& source) 
 {
 	SetSize(source.GetRows(), source.GetCols());
 }
 
 
 template <class Type>
 void Matrix<Type>::SetSize (unsigned int newSize) 
 {
 	SetSize(newSize, newSize);
 }
 
 
 template <class Type>
 void Matrix<Type>::SetSize (unsigned int newRows, unsigned int newCols) 
 {
 	unsigned int oldRows = GetRows();
 	unsigned int oldCols = GetCols();
 	//--------------------------------------------------------------------------
 	//  No size change?
 	//--------------------------------------------------------------------------
 	if ((newRows == oldRows) && (newCols == oldCols)) 
 	{
		memset(_pMatrixData, 0, GetSize() * (sizeof (Type)));
 		return;
 	}
 	//--------------------------------------------------------------------------
 	//  Delete the previous matrix data, then allocate and initialize the new 
 	//  matrix data.
 	//--------------------------------------------------------------------------
 	delete [] _pMatrixData;
 	_pMatrixData = NULL;
 	_nRow = newRows;
 	_nColumn = newCols;
 	_nTotal = newRows * newCols;
 
 	init(_nRow, _nColumn);
 }
 
 //==============================================================================
 //
 //  Matrix::Resize
 //
 //  Description:
 //      Resizes the matrix, preserving the previous data.
 //
 //  Parameters:
 //      newRows - number of rows
 //      newCols - number of columns
 //
 //==============================================================================
template <class Type>
void Matrix<Type>::Resize(unsigned int newRows, unsigned int newCols) 
{
 	unsigned int oldRows = GetRows();
 	unsigned int oldCols = GetCols();
 	unsigned int minRows = oldRows < newRows ? oldRows : newRows;
 	//--------------------------------------------------------------------------
 	//  No size change?
 	//--------------------------------------------------------------------------
 	if ((newRows == oldRows) && (newCols == oldCols))
	{
		return;
	}
 	//--------------------------------------------------------------------------
 	//  Update the matrix extents and allocate the new matrix data.
 	//--------------------------------------------------------------------------
 	Type* oldData = _pMatrixData;
 	try
 	{
 		_pMatrixData = NULL;
		
		//set row number, column number and total number
		_nRow    = newRows;
		_nColumn = newCols;
		_nTotal  = newRows * newCols;

 		unsigned int elementCount = _nTotal;
 		if (elementCount != 0) 
		{
 			_pMatrixData = new Type[elementCount];
 			//----------------------------------------------------------------------
 			//  Compose the new matrix from the old matrix data
 			//----------------------------------------------------------------------
 			if(newCols > oldCols) 
			{
 				Type *ptr = _pMatrixData;
 				for(unsigned int ind = 0; ind < minRows; ++ind) 
				{
 					memmove(ptr,(oldData+ind*oldCols),oldCols*sizeof(Type));
 					memset((ptr+oldCols), 0, (newCols-oldCols) * (sizeof (Type)));
					ptr+=newCols;
 				}   
 				// Set remaining cells to zero
 				if(elementCount > (unsigned)(((Type*)ptr) - _pMatrixData)) 
				{
 					unsigned int cells_left = elementCount - (((Type*)ptr) - _pMatrixData);
 					memset(ptr, 0, cells_left * (sizeof (Type)));
 				}   
 		
 			}
 			else 
			{
 				Type *ptr = _pMatrixData;
 				for(unsigned int ind = 0; ind < minRows; ++ind) 
				{
 					memmove(ptr,(oldData+ind*oldCols),newCols*sizeof(Type));
 					ptr+=newCols;
 				}
 				// Set remaining cells to zero
 				if(elementCount > (unsigned)(((Type*)ptr) - _pMatrixData)) 
				{
 					unsigned int cells_left = elementCount - (((Type*)ptr) - _pMatrixData);
 					memset(ptr, 0, cells_left * (sizeof (Type)));
 				}
 			}
 		}
 		else 
		{
 			_pMatrixData = NULL;
 		}
 		delete [] oldData;
 	}
 	catch (...)
 	{
 		delete _pMatrixData;
 		_pMatrixData = oldData;
 		SetSize(oldRows,oldCols);
 		throw;
	}
}

//==============================================================================
//
//  Matrix::Scale
//
//  Description:
//      Scales the matrix.
//
//  Parameters:
//      scale  - scaling value
//      offset - offset value
//
//==============================================================================
template <class Type>
void Matrix<Type>::Scale (const Type& scale) 
{
	const unsigned int elementCount = GetSize();
	for (unsigned int i=0; i < elementCount; i++) 
	{
		_pMatrixData[i] *= scale;
	}
}

template <class Type>
void Matrix<Type>::Scale (const Type& scale, const Type& offset) 
{
	const unsigned int elementCount = GetSize();
	for (unsigned int i=0; i < elementCount; i++) 
	{
		_pMatrixData[i] *= scale;
		_pMatrixData[i] += offset;
	}
}

//==============================================================================
//
//  Matrix::CopyRow
//  Matrix::CopyColumn
//
//  Description:
//      Copies the specified row/column of the matrix into the given vector.
//      The vector must be of the same type as the matrix (e.g. 
//      a row/column from a CNiComplexReal64Matrix can be copied to a 
//      CNiComplexReal64Vector).
//
//  Parameters:
//      rowIndex - index of the row to copy
//      columnIndex - index of the column to copy
//      destVector - vector which will hold the copy
//
//  Return Value:
//      None.
//
//==============================================================================
template <class Type>
void Matrix<Type>::CopyRow (unsigned int rowIndex, Vector<Type>& destVector) const
{
#ifdef BOUNDS_CHECK
	assert( rowIndex >= 0 );
	assert( rowIndex < _nRow );
#endif

	destVector.SetSize(GetCols());
	Type *vec = destVector; // pointer to internal vector data
	memmove(vec,(_pMatrixData+rowIndex*GetCols()),GetCols()*sizeof(Type));
}


template <class Type>
void Matrix<Type>::CopyColumn (unsigned int columnIndex, Vector<Type>& destVector) const
{
#ifdef BOUNDS_CHECK
	assert( columnIndex >= 0 );
	assert( columnIndex < _nColumn );
#endif

	destVector.SetSize(GetRows());
	Type *vec = destVector; // pointer to internal vector data
	for (int rowIndex=0; rowIndex < GetRows(); ++rowIndex) 
	{
		vec[rowIndex] = _pMatrixData[rowIndex*GetCols()+columnIndex];
	}
}

//==============================================================================
//
//  Matrix::AssignRow
//  Matrix::AssignColumn
//
//  Description:
//      Assigns the given vector to the specified row/column of the matrix.  The
//      vector must be of the same type as the matrix and have at least the 
//      same number of elements as the vector has columns/rows.  For example, 
//      a row/column from a 4-by-5 CNiComplexReal64Matrix must be copied from a 
//      CNiComplexReal64Vector with at least 5/4 elements.
//
//  Parameters:
//      rowIndex - index of the row to assign
//      columnIndex - index of the column to assign
//      sourceVector - vector to assign from
//
//==============================================================================
template <class Type>
void Matrix<Type>::AssignRow (unsigned int rowIndex, const Vector<Type>& sourceVector)
{
#ifdef BOUNDS_CHECK
	assert( rowIndex >= 0 );
	assert( rowIndex < _nRow );
	assert( v.dim() == _nColumn );
#endif

	const Type *vec = sourceVector; // pointer to internal vector data
	memmove((_pMatrixData+rowIndex*GetCols()),vec,GetCols()*sizeof(Type));
}


template <class Type>
void Matrix<Type>::AssignColumn (unsigned int columnIndex, const Vector<Type>& sourceVector)
{
#ifdef BOUNDS_CHECK
	assert( columnIndex >= 0 );
	assert( columnIndex < _nColumn );
	assert( v.dim() == _nRow );
#endif

	const Type *vec = sourceVector; // pointer to internal vector data
	for (unsigned int rowIndex=0; rowIndex < GetRows(); ++rowIndex) 
	{
		_pMatrixData[rowIndex*GetCols()+columnIndex] = vec[rowIndex];
	}
}

//==============================================================================
//
//  Matrix::MinMax
//
//  Description:
//      Finds the minimum and maximum elements of the matrix.
//
//  Parameters:
//      minData     - minimum value
//      minRowIndex - row index of minimum value
//      minColIndex - column index of minimum value
//      maxData     - maximum value         
//      maxRowIndex - row index of maximum value
//      maxColIndex - column index of maximum value
//
//==============================================================================
template <class Type>
void Matrix<Type>::MinMax (Type& minData, unsigned int& minRowIndex, unsigned int& minColIndex, Type& maxData, unsigned int& maxRowIndex, unsigned int& maxColIndex) const 
{
	if (GetRows() * GetCols() > 0) 
	{
		const unsigned int rows = GetRows();
		const unsigned int cols = GetCols();
		//----------------------------------------------------------------------
		//  Set initial min and max extents.
		//----------------------------------------------------------------------
		minData     = _pMatrixData[0];
		minRowIndex = 0;
		minColIndex = 0;
		maxData     = _pMatrixData[0];
		maxRowIndex = 0;
		maxColIndex = 0;
		//----------------------------------------------------------------------
		//  Traverse matrix comparing data values.
		//----------------------------------------------------------------------
		for (unsigned int i=0; i < rows; i++) 
		{
			for (unsigned int j=0; j < cols; j++) 
			{
				Type value = _pMatrixData[i*cols+j];
				if (value < minData) 
				{
					   minData  = value;
					minRowIndex = i;
					minColIndex = j;
				}
				if (value > maxData) 
				{
					   maxData  = value;
					maxRowIndex = i;
					maxColIndex = j;
				}
			}
		}
	}
	else 
	{
		//----------------------------------------------------------------------
		//  Indicate absent data.
		//----------------------------------------------------------------------
		minRowIndex = 0xffffffff;  /* maximum unsigned int value */
		maxRowIndex = 0xffffffff;
		minColIndex = 0xffffffff;
		maxColIndex = 0xffffffff;
	}
}
