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
 *                               vector-impl.h
 *
 * Implementation for Vector class.
 *
 * Zhang Ming, 2010-01 (revised 2010-12), Xi'an Jiaotong University.
 *****************************************************************************/


/**
 * initialize
 */
template <typename Type>
void Vector<Type>::init( int length )
{
	assert( _pVectorData == NULL );
	_pVectorData = new Type[length];

	assert( _pVectorData != NULL );
	pv1 = _pVectorData - 1;
	_nSize = length;

	memset(_pVectorData, 0, GetSize() * (sizeof (Type)));
}


/**
 * copy vector from normal array
 */
template <typename Type>
inline void Vector<Type>::copyFromArray( const Type *v )
{
	for( int i=0; i<_nSize; ++i )
		_pVectorData[i] = v[i];
}

/**
 * destroy the vector
 */
template <typename Type>
void Vector<Type>::destroy()
{
	if( _pVectorData == NULL )
		return;

	delete []_pVectorData;

	_pVectorData = NULL;
	pv1 = NULL;
}


/**
 * constructors and destructor
 */
template <typename Type>
Vector<Type>::Vector()
: _pVectorData(0), pv1(0), _nSize(0)
{
}

template <typename Type>
Vector<Type>::Vector( const Vector<Type> &v )
: _pVectorData(0), pv1(0), _nSize(0)
{
	init( v._nSize );
	copyFromArray( v._pVectorData );
}

template <typename Type>
Vector<Type>::Vector( int length )
:  _pVectorData(0), pv1(0), _nSize(0)
{
	init( length );
}

template <typename Type>
Vector<Type>::Vector( int length, const Type &x )
:  _pVectorData(0), pv1(0), _nSize(0)
{
	init( length );
	Set( x );
}

template <typename Type>
Vector<Type>::Vector( int length, const Type *array )
:  _pVectorData(0), pv1(0), _nSize(0)
{
	init( length );
	copyFromArray( array );
}

template <typename Type>
Vector<Type>::~Vector()
{
	destroy();
}


/**
 * overload evaluate operator= from vector to vector
 */
template <typename Type>
Vector<Type>& Vector<Type>::operator=( const Vector<Type> &v )
{
	if( _pVectorData == v._pVectorData )
		return *this;

	if( _nSize == v._nSize )
	{
		copyFromArray( v._pVectorData );
	}
	
	else
	{
		destroy();
		init( v._nSize );
		copyFromArray( v._pVectorData );
	}

	return *this;
}


/**
 * overload evaluate operator= from scalar to vector
 */
template <typename Type>
inline Vector<Type>& Vector<Type>::operator=( const Type &x )
{
	Set( x );

	return *this;
}


/**
 * overload operator [] for 0-offset access
 */
template <typename Type>
inline Type& Vector<Type>::operator[]( int i )
{
#ifdef BOUNDS_CHECK
	assert( 0 <= i );
	assert( i < _nSize );
#endif

	return _pVectorData[i];
}

template <typename Type>
inline const Type& Vector<Type>::operator[]( int i ) const
{
#ifdef BOUNDS_CHECK
	assert( 0 <= i );
	assert( i < _nSize );
#endif

	return _pVectorData[i];
}


/**
 * overload operator () for 1-offset access
 */
template <typename Type>
inline Type& Vector<Type>::operator()( int i )
{
#ifdef BOUNDS_CHECK
	assert( 1 <= i );
	assert( i <= _nSize );
#endif

	return pv1[i];
}

template <typename Type>
inline const Type& Vector<Type>::operator()( int i ) const
{
#ifdef BOUNDS_CHECK
	assert( 1 <= i );
	assert( i <= _nSize );
#endif

	return pv1[i];
}


/**
 * iterators
 */
template <typename Type>
inline typename Vector<Type>::iterator Vector<Type>::begin()
{
    return _pVectorData;
}

template <typename Type>
inline typename Vector<Type>::const_iterator Vector<Type>::begin() const
{
    return _pVectorData;
}

template <typename Type>
inline typename Vector<Type>::iterator Vector<Type>::end()
{
    return _pVectorData + _nSize;
}

template <typename Type>
inline typename Vector<Type>::const_iterator Vector<Type>::end() const
{
    return _pVectorData + _nSize;
}


/**
 * type conversion functions
 */
template <typename Type>
inline Vector<Type>::operator Type*()
{
	return _pVectorData;
}

template <typename Type>
inline Vector<Type>::operator const Type*() const
{
	return _pVectorData;
}


/**
 * get the vector's total size
 */
template <typename Type>
inline int Vector<Type>::GetSize() const
{
	return  _nSize;
}

/**
 * get the vector's dimension
 */
template <typename Type>
inline int Vector<Type>::dim() const
{
	return  _nSize;
}

/**
 * compound assignment operators +=
 */
template <typename Type>
Vector<Type>& Vector<Type>::operator+=( const Type &x )
{
    iterator itr = (*this).begin();
    while( itr != (*this).end() )
        *itr++ += x;

	return *this;
}

template <typename Type>
Vector<Type>& Vector<Type>::operator+=( const Vector<Type> &rhs )
{
    assert( _nSize == rhs.dim() );

    iterator itrL = (*this).begin();
    const_iterator itrR = rhs.begin();
    while( itrL != (*this).end() )
        *itrL++ += *itrR++;

	return *this;
}


/**
 * compound assignment operators -=
 */
template <typename Type>
Vector<Type>& Vector<Type>::operator-=( const Type &x )
{
    iterator itr = (*this).begin();
    while( itr != (*this).end() )
        *itr++ -= x;

	return *this;
}

template <typename Type>
Vector<Type>& Vector<Type>::operator-=( const Vector<Type> &rhs )
{
    assert( _nSize == rhs.dim() );

    iterator itrL = (*this).begin();
    const_iterator itrR = rhs.begin();
    while( itrL != (*this).end() )
        *itrL++ -= *itrR++;

	return *this;
}


/**
 * compound assignment operators *=
 */
template <typename Type>
Vector<Type>& Vector<Type>::operator*=( const Type &x )
{
    iterator itr = (*this).begin();
    while( itr != (*this).end() )
        *itr++ *= x;

	return *this;
}

template <typename Type>
Vector<Type>& Vector<Type>::operator*=( const Vector<Type> &rhs )
{
    assert( _nSize == rhs.dim() );

    iterator itrL = (*this).begin();
    const_iterator itrR = rhs.begin();
    while( itrL != (*this).end() )
        *itrL++ *= *itrR++;

	return *this;
}


/**
 * compound assignment operators /=
 */
template <typename Type>
Vector<Type>& Vector<Type>::operator/=( const Type &x )
{
    iterator itr = (*this).begin();
    while( itr != (*this).end() )
        *itr++ /= x;

	return *this;
}

template <typename Type>
Vector<Type>& Vector<Type>::operator/=( const Vector<Type> &rhs )
{
    assert( _nSize == rhs.dim() );

    iterator itrL = (*this).begin();
    const_iterator itrR = rhs.begin();
    while( itrL != (*this).end() )
        *itrL++ /= *itrR++;

	return *this;
}


/**
 * Overload the output stream function.
 */
template <typename Type>
ostream& operator<<( ostream &out, const Vector<Type> &v )
{
	int N = v.dim();
	out << "size: " << N << " by 1" << "\n";

	for( int i=0; i<N; ++i )
		out << v[i] << " " << "\n";

	return out;
}


/**
 * Overload the input stream function.
 */
template <typename Type>
istream& operator>>( istream &in, Vector<Type> &v )
{
	int N;
	in >> N;

	if( !( N == v.dim() ) )
		v.Resize( N );

	for( int i=0; i<N; ++i )
		in >> v[i];

	return in;
}


/**
 * get negative vector
 */
template <typename Type>
Vector<Type> operator-( const Vector<Type> &v )
{
	Vector<Type> tmp( v.dim() );
    typename Vector<Type>::iterator itrL = tmp.begin();
    typename Vector<Type>::const_iterator itrR = v.begin();

    while( itrL != tmp.end() )
        *itrL++ = -(*itrR++);

    return tmp;
}


/**
 * vector-scalar addition.
 */
template <typename Type>
inline Vector<Type> operator+( const Vector<Type> &v, const Type &x )
{
	Vector<Type> tmp( v );
	return tmp += x;
}

template <typename Type>
inline Vector<Type> operator+( const Type &x, const Vector<Type> &v )
{
	return v+x;
}


/**
 * vector-scalar substraction.
 */
template <typename Type>
inline Vector<Type> operator-( const Vector<Type> &v, const Type &x )
{
	Vector<Type> tmp( v );
	return tmp -= x;
}

template <typename Type>
inline Vector<Type> operator-( const Type &x, const Vector<Type> &v )
{
	Vector<Type> tmp( v );
	return -tmp += x;
}


/**
 * vector-scalar multiplication.
 */
template <typename Type>
inline Vector<Type> operator*( const Vector<Type> &v, const Type &x )
{
	Vector<Type> tmp( v );
	return tmp *= x;
}

template <typename Type>
inline Vector<Type> operator*( const Type &x, const Vector<Type> &v )
{
	return v*x;
}


/**
 * vector-scalar division.
 */
template <typename Type>
inline Vector<Type> operator/( const Vector<Type> &v, const Type &x )
{
	Vector<Type> tmp( v );
	return tmp /= x;
}

template <typename Type>
inline Vector<Type> operator/( const Type &x, const Vector<Type> &v )
{
	int N = v.dim();
	Vector<Type> tmp( N );

	for( int i=0; i<N; ++i )
		tmp[i] = x / v[i];

	return tmp;
}


/**
 * vector-vector addition.
 */
template <typename Type>
inline Vector<Type> operator+( const Vector<Type> &v1, const Vector<Type> &v2 )
{
    Vector<Type> tmp( v1 );
	return tmp += v2;
}


/**
 * vector-vector substraction.
 */
template <typename Type>
inline Vector<Type> operator-( const Vector<Type> &v1, const Vector<Type> &v2 )
{
    Vector<Type> tmp( v1 );
	return tmp -= v2;
}


/**
 * vector-vector multiplication.
 */
template <typename Type>
inline Vector<Type> operator*( const Vector<Type> &v1, const Vector<Type> &v2 )
{
    Vector<Type> tmp( v1 );
	return tmp *= v2;
}


/**
 * vector-vector division.
 */
template <typename Type>
inline Vector<Type> operator/( const Vector<Type> &v1, const Vector<Type> &v2 )
{
    Vector<Type> tmp( v1 );
	return tmp /= v2;
}


/**
 * Inner product for vectors.
 */
template <typename Type>
Type dotProd( const Vector<Type> &v1, const Vector<Type> &v2 )
{
	assert( v1.dim() == v2.dim() );

    Type sum = 0;
    typename Vector<Type>::const_iterator itr1 = v1.begin();
    typename Vector<Type>::const_iterator itr2 = v2.begin();

    while( itr1 != v1.end() )
		sum += (*itr1++) * (*itr2++);

	return sum;
}


/**
 * Inner product for vectors.
 */
template <typename Type>
complex<Type> dotProd( const Vector<complex<Type> > &v1,
                       const Vector<complex<Type> > &v2 )
{
	assert( v1.dim() == v2.dim() );

    complex<Type> sum = 0;
    typename Vector<complex<Type> >::const_iterator itr1 = v1.begin();
    typename Vector<complex<Type> >::const_iterator itr2 = v2.begin();

    while( itr1 != v1.end() )
		sum += (*itr1++) * conj(*itr2++);

	return sum;
}


/**
 * Vector's sum.
 */
template <typename Type>
Type sum( const Vector<Type> &v )
{
    Type sum = 0;
    typename Vector<Type>::const_iterator itr = v.begin();

    while( itr != v.end() )
		sum += *itr++;

	return sum;
}


/**
 * Minimum value of vector.
 */
template <typename Type>
Type min( const Vector<Type> &v )
{
    Type m = v[0];
    for( int i=1; i<v.GetSize(); ++i )
        if( m > v[i] )
            m = v[i];

    return m;
}


/**
 * Maximum value of vector.
 */
template <typename Type>
Type max( const Vector<Type> &v )
{
    Type M = v[0];
    for( int i=1; i<v.GetSize(); ++i )
        if( M < v[i] )
            M = v[i];

    return M;
}


/**
 * Vector's norm in Euclidean space.
 */
template <typename Type>
Type norm( const Vector<Type> &v )
{
	Type sum = 0;
	typename Vector<Type>::const_iterator itr = v.begin();

	while( itr != v.end() )
	{
	    sum += (*itr) * (*itr);
	    itr++;
	}

	return Type(sqrt(1.0*sum));
}


/**
 * Vector's norm in Euclidean space.
 */
template <typename Type>
Type norm( const Vector<complex<Type> > &v )
{
	Type sum = 0;
	typename Vector<complex<Type> >::const_iterator itr = v.begin();

	while( itr != v.end() )
	    sum += norm(*itr++);

	return Type(sqrt(1.0*sum));
}


/**
 * return vector's reversion
 */
template <typename Type>
void swap( Vector<Type> &lhs, Vector<Type> &rhs )
{
    typename Vector<Type>::iterator itrL = lhs.begin(),
                                    itrR = rhs.begin();

    while( itrL != lhs.end() )
        std::swap( *itrL++, *itrR++ );
}


/**
 * Generates a vector of n points linearly spaced between and
 * including a and b.
 */
template <typename Type>
Vector<Type> linspace( Type a, Type b, int n )
{
    if( n < 1 )
        return Vector<Type>();
    else if( n == 1 )
        return Vector<Type>( 1, a );
    else
    {
        Type dx = (b-a) / (n-1);

        Vector<Type> tmp(n);
        for( int i=0; i<n; ++i )
            tmp[i] = a + i*dx;

        return tmp;
    }
}


/**
 * Get magnitude of a complex vector.
 */
template <typename Type>
Vector<Type> abs( const Vector< complex<Type> > &v )
{
    Vector<Type> tmp( v.dim() );
    typename Vector<Type>::iterator itrL = tmp.begin();
    typename Vector< complex<Type> >::const_iterator itrR = v.begin();

    while( itrL != tmp.end() )
        *itrL++ = abs(*itrR++);

    return tmp;
}


/**
 * Get angle of a complex vector.
 */
template <typename Type>
Vector<Type> arg( const Vector< complex<Type> > &v )
{
    Vector<Type> tmp( v.dim() );
    typename Vector<Type>::iterator itrL = tmp.begin();
    typename Vector< complex<Type> >::const_iterator itrR = v.begin();

    while( itrL != tmp.end() )
        *itrL++ = arg(*itrR++);

    return tmp;
}


/**
 * Get real part of a complex vector.
 */
template <typename Type>
Vector<Type> real( const Vector< complex<Type> > &v )
{
    Vector<Type> tmp( v.dim() );
    typename Vector<Type>::iterator itrL = tmp.begin();
    typename Vector< complex<Type> >::const_iterator itrR = v.begin();

    while( itrL != tmp.end() )
        *itrL++ = (*itrR++).real();

    return tmp;
}


/**
 * Get imaginary part of a complex vector.
 */
template <typename Type>
Vector<Type> imag( const Vector< complex<Type> > &v )
{
    Vector<Type> tmp( v.dim() );
    typename Vector<Type>::iterator itrL = tmp.begin();
    typename Vector< complex<Type> >::const_iterator itrR = v.begin();

    while( itrL != tmp.end() )
        *itrL++ = (*itrR++).imag();

    return tmp;
}


/**
 * Convert real vector to complex vector.
 */
template <typename Type>
Vector<complex<Type> > complexVector( const Vector<Type> &rv )
{
	int N = rv.dim();

    Vector<complex<Type> > cv( N );
    typename Vector<complex<Type> >::iterator itrL = cv.begin();
    typename Vector<Type>::const_iterator itrR = rv.begin();

    while( itrR != rv.end() )
        *itrL++ = *itrR++;

    return cv;
}

template <typename Type>
Vector<complex<Type> > complexVector( const Vector<Type> &vR,
                                      const Vector<Type> &vI )
{
	int N = vR.dim();

	assert( N == vI.dim() );

    Vector<complex<Type> > cv( N );
    typename Vector<complex<Type> >::iterator itrC = cv.begin();
    typename Vector<Type>::const_iterator itrR = vR.begin(),
                                          itrI = vI.begin();

    while( itrC != cv.end() )
        *itrC++ = complex<Type>( *itrR++, *itrI++ );

    return cv;
}

//==============================================================================
//
//  Vector::SetSize
//
//  Description:
//      Sets the size of the vector
//
//  Parameters:
//      source  - vector to match size of
//      newSize - new vector size
//
//==============================================================================
template <typename Type>
inline void Vector<Type>::SetSize(const Vector<Type>& source)
{
	SetSize(source.GetSize());
}

template <class Type>
inline void Vector<Type>::SetSize (const unsigned int newSize) 
{
	const unsigned int oldSize = GetSize();
	//--------------------------------------------------------------------------
	//  No size change?
	//--------------------------------------------------------------------------
	if (newSize == oldSize) 
	{
		init(newSize);
		return;
	}
	//--------------------------------------------------------------------------
	//  Delete the previous vector data, then allocate and initialize the new 
	//  vector data.
	//--------------------------------------------------------------------------
	delete [] _pVectorData;
	_pVectorData = NULL;
	_nSize = newSize;
	if (newSize != 0) 
	{
		init(newSize);
	}
}

//==============================================================================
//
//  Vector::Set
//
//  Description:
//      Sets every element in the vector to the specfied value.
//
//  Parameters:
//      value - value to set
//
//==============================================================================
template <typename Type>
inline void Vector<Type>::Set( const Type &x )
{
	for( int i=0; i<_nSize; ++i )
	{
		_pVectorData[i] = x;
	}
}

//==============================================================================
//
//  Vector::Scale
//
//  Description:
//      Scales the vector.
//
//  Parameters:
//      scale  - scaling value
//      offset - offset value
//
//==============================================================================
template <typename Type>
void Vector<Type>::Scale (const Type& scale) 
{
	for (int i=0; i < _nSize; i++)
	{
		_pVectorData[i] *= scale;
	}
}

template <class Type>
void Vector<Type>::Scale (const Type& scale, const Type& offset) 
{
	for (int i=0; i < _nSize; i++) 
	{
		_pVectorData[i] *= scale;
		_pVectorData[i] += offset;
	}
}

//==============================================================================
//
//  Vector::Resize
//
//  Description:
//      Resizes the vector, preserving the previous data.
//
//  Parameters:
//      newSize - new vector size
//
//==============================================================================
template <class Type>
void Vector<Type>::Resize (unsigned int newSize) 
{
	const unsigned int oldSize = GetSize();
	//--------------------------------------------------------------------------
	//  No size change?
	//--------------------------------------------------------------------------
	if (newSize == oldSize) 
	{
		return;
	}
	//--------------------------------------------------------------------------
	//  Update the vector extents and allocate the new vector data.
	//--------------------------------------------------------------------------
	Type* oldData = _pVectorData;
	try
	{
		_pVectorData = NULL;
		_nSize = newSize;
		if (newSize != 0) 
		{
			_pVectorData = new Type[newSize];
			//----------------------------------------------------------------------
			//  Compose the new vector from the old vector data, padding any extra
			//  cells with zero.
			//----------------------------------------------------------------------
			unsigned int size = newSize < oldSize? newSize : oldSize;
			memmove(_pVectorData,oldData,size*sizeof(Type));
			if(newSize > oldSize)
			{
				memset((_pVectorData + oldSize), 0, (newSize - oldSize) * (sizeof (Type)));
			}
		}
		delete [] oldData;
	}
	catch (...)
	{
		delete _pVectorData;
		_pVectorData = oldData;
		_nSize = oldSize;
		throw;
	}
}

//==============================================================================
//
//  Vector::Append
//
//  Description:
//      Appends a vector or single element to this vector.
//
//  Parameters:
//      other   - vector to append
//      element - element to append
//
//  Return Value:
//      Reference to this object
//
//==============================================================================
template <class Type>
Vector<Type>& Vector<Type>::Append (const Vector<Type>& other) 
{
	const unsigned int sizeOfThis  = GetSize();
	const unsigned int sizeOfOther = other.GetSize();
	Resize(sizeOfThis + sizeOfOther);
	memmove(_pVectorData+sizeOfThis,other._pVectorData,sizeOfOther * sizeof(Type));
	return *this;
}

template <class Type>
Vector<Type>& Vector<Type>::Append (const Type& element) 
{
	const unsigned int sizeOfThis = GetSize();
	Resize(sizeOfThis + 1);
	_pVectorData[sizeOfThis] = element;
	return *this;
}

//==============================================================================
//
//  Vector::MinMax
//
//  Description:
//      Finds the minimum and maximum elements of the vector.
//
//  Parameters:
//      minData  - minimum value in vector
//      minIndex - index of minimum value
//      maxData  - maximum value in vector
//      maxIndex - index of maximum value
//
//==============================================================================
template <class Type>
void Vector<Type>::MinMax (Type& minData, unsigned int& minIndex, Type& maxData, unsigned int& maxIndex) const 
{
	const unsigned int size = GetSize();
	if (size > 0)
	{
		//----------------------------------------------------------------------
		//  Set initial min and max extents.
		//----------------------------------------------------------------------
		minData  = _pVectorData[0];
		minIndex = 0;
		maxData  = _pVectorData[0];
		maxIndex = 0;
		//----------------------------------------------------------------------
		//  Traverse vector interrogating data values.
		//----------------------------------------------------------------------
		for (unsigned int i=1; i < size; i++) 
		{
			if (_pVectorData[i] < minData) 
			{
				minData  = _pVectorData[i];
				minIndex = i;
			}
			if (_pVectorData[i] > maxData) 
			{
				maxData  = _pVectorData[i];
				maxIndex = i;
			}
		}
	}
	else {
		//----------------------------------------------------------------------
		//  Indicate absent data with negative indices.
		//----------------------------------------------------------------------
		minIndex = 0xffffffff;
		maxIndex = 0xffffffff;
	}
}

//==============================================================================
//
//  Vector::Abs
//
//  Description:
//      Calculates the absolute value of each element in the vector.
//
//==============================================================================
template <class Type>
void Vector<Type>::Abs() 
{
	const unsigned int size = GetSize();
	for (unsigned int i=0; i < size; i++)
	{
		if (_pVectorData[i] < 0)
		{
			_pVectorData[i] = -_pVectorData[i];
		}
	}
}

//==============================================================================
//
//  Vector::Subset
//
//  Description:
//      Creates a subset from the vector object.
//
//  Parameters:
//      source - source vector
//      subset - on return, contains the vector subset
//      offset - starting index of subset
//      size   - size of subset (-1 for all elements)
//
//==============================================================================

template <class Type>
Vector<Type> Vector<Type>::Subset (unsigned int offset, int size) const 
{
	Vector<Type> subset;

	const unsigned int sourceSize = GetSize();
	//--------------------------------------------------------------------------
	//  Constrain subset to the extents of the vector.
	//--------------------------------------------------------------------------
	const unsigned int sourceOffset = offset < sourceSize ? offset : sourceSize;
	const unsigned int subsetSize   = (size < 0) ? (sourceSize - sourceOffset) : ((unsigned int)size < (sourceSize-sourceOffset) ? size : (sourceSize-sourceOffset));
	//--------------------------------------------------------------------------
	//  Size the subset and extract the elements.
	//--------------------------------------------------------------------------
	subset.SetSize(subsetSize);
	memmove(subset._pVectorData,(_pVectorData+sourceOffset),subsetSize * sizeof(Type));

	return subset;
}

template <class Type>
void Vector<Type>::Subset ( Vector<Type>& source, Vector<Type>& subset, unsigned int offset, int size) 
{
	const unsigned int sourceSize = source.GetSize();
	//--------------------------------------------------------------------------
	//  Constrain subset to the extents of the vector.
	//--------------------------------------------------------------------------
	const unsigned int sourceOffset = offset < sourceSize ? offset : sourceSize;
	const unsigned int subsetSize   = (size < 0) ? (sourceSize - sourceOffset) : ((unsigned int)size < (sourceSize-sourceOffset) ? size : (sourceSize-sourceOffset));
	//--------------------------------------------------------------------------
	//  Size the subset and extract the elements.
	//--------------------------------------------------------------------------
	subset.SetSize(subsetSize);
 	memmove(subset._pVectorData,(source._pVectorData+sourceOffset),subsetSize * sizeof(Type));
}

//==============================================================================
//
//  Vector::Sum
//
//  Description:
//      Calculates the sum of every element in the vector.
//
//  Return Value:
//      Result
//
//==============================================================================
template <class Type>
Type Vector<Type>::Sum () const 
{
	Type result(0);
	const unsigned int size = GetSize();
	for (unsigned int i=0; i < size; i++) 
	{
		result += _pVectorData[i];
	}
	return result;
}