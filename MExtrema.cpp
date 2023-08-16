#include <DataFilters2/include/MExtrema.h>
#include <Log/include/Log.h>
#include <MDR/include/MVoxelData.h>
#include <ScientificCore/include/MMathematicalUtils.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <iterator>
#include <numeric>
#include <set>
#include <tuple>
#include <vector>

namespace msc
{

LOG_INIT( "MSC_SDK.MExtrema" );

MExtrema::ReturnType MExtrema::calculate( const mcp::MVoxelData& aData )
{
	const std::uint64_t M = aData.sizeX();
	const std::uint64_t N = aData.sizeY();

	if ( !M || !N )
	{
		LOG_ERROR( "calculate: x or y dimension was 0" );
		return {};
	}

	// Search peaks through columns:
	Indices2DType smaxcolM, smincolM;
	std::tie( smaxcolM, smincolM ) = extremos( aData );

	if ( smaxcolM.empty() || smincolM.empty() )
	{
		LOG_ERROR( "calculate: No extrema values found in columns" );
		return {};
	}

	// Search peaks through rows, on columns with extrema points:
	std::set<std::uint64_t> columnExtremaRowIndices; // Rows with column extrema
	for ( auto& maxCol : smaxcolM )
	{
		columnExtremaRowIndices.emplace( maxCol.x() );
	}

	for ( auto& minCol : smincolM )
	{
		columnExtremaRowIndices.emplace( minCol.x() );
	}

	const std::size_t imSize = columnExtremaRowIndices.size();
	mcp::MVoxelData	  filteredData( imSize, N, 1 );

	const float* dataPtr		 = aData.slice( 0 );
	float* const filteredDataPtr = filteredData.slice( 0 );
	int			 idx			 = 0;
	for ( const auto& i : columnExtremaRowIndices )
	{
		for ( std::size_t j = 0; j < N; ++j )
		{
			filteredDataPtr[idx + j * imSize] = dataPtr[j + i * M];
		}
		++idx;
	}

	Indices2DType smaxfilM, sminfilM;
	std::tie( smaxfilM, sminfilM ) = extremos( filteredData );

	if ( smaxfilM.empty() || sminfilM.empty() )
	{
		LOG_ERROR( "calculate: no extrema values found in rows where it was "
				   "already found "
				   "in the column search " ); // unlikely
		return {};
	}

	// Convertion from 2 to 1 index:
	auto firstAndSecondColOf = []( const Indices2DType& aCont ) {
		const std::size_t rowSize = aCont.size();

		IndicesType firstCol;
		IndicesType secondCOl;

		firstCol.reserve( rowSize );
		secondCOl.reserve( rowSize );

		for ( auto& row : aCont )
		{
			firstCol.emplace_back( row.x() );
			secondCOl.emplace_back( row.y() );
		}

		return std::make_pair( firstCol, secondCOl );
	};

	auto temp = firstAndSecondColOf( smaxcolM );

	IndicesType smaxcolV{
		MMathematicalUtils::sub2ind( { M, N }, temp.first, temp.second ) };

	temp = firstAndSecondColOf( smincolM );

	IndicesType smincolV{
		MMathematicalUtils::sub2ind( { M, N }, temp.first, temp.second ) };

	auto filterFrom = [&columnExtremaRowIndices]( const Indices2DType& aCont ) {
		const std::size_t rowSize = aCont.size();

		IndicesType secondCol;
		IndicesType firstCol;

		secondCol.reserve( rowSize );
		firstCol.reserve( rowSize );

		for ( auto& row : aCont )
		{
			firstCol.emplace_back( row.x() );
			secondCol.emplace_back(
				*std::next( columnExtremaRowIndices.begin(), row.y() ) );
		}

		return std::make_pair( secondCol, firstCol );
	};

	temp = filterFrom( smaxfilM );

	IndicesType smaxfilV{
		MMathematicalUtils::sub2ind( { M, N }, temp.first, temp.second ) };

	temp = filterFrom( sminfilM );

	IndicesType sminfilV{
		MMathematicalUtils::sub2ind( { M, N }, temp.first, temp.second ) };

	// Peaks in rows and in columns:

	IndicesType smax = intersect( smaxcolV, smaxfilV ); // sorted unique values
	IndicesType smin = intersect( smincolV, sminfilV );

	// Search peaks through diagonals?
	std::set<std::uint64_t> extremaIndices;

	for ( auto& v : smax )
	{
		extremaIndices.emplace( v );
	}

	for ( auto& v : smin )
	{
		extremaIndices.emplace( v );
	}

	// Check peaks on down-up diagonal:
	IndicesType iext;
	IndicesType jext;
	std::tie( iext, jext ) = MMathematicalUtils::ind2sub( { M, N }, extremaIndices );

	IndicesType sextmax;
	IndicesType sextmin;
	std::tie( sextmax, sextmin ) = extremosDiag( iext, jext, aData, 1 );

	// Check peaks on up-down diagonal:
	smax = intersect( smax, sextmax, { M - 1, N * M - M - 1 } );
	smin = intersect( smin, sextmin, { M - 1, N * M - M - 1 } );

	// Peaks on up-down diagonals:
	extremaIndices.clear();
	for ( auto& v : smax )
	{
		extremaIndices.emplace( v );
	}

	for ( auto& v : smin )
	{
		extremaIndices.emplace( v );
	}

	std::tie( iext, jext ) = MMathematicalUtils::ind2sub( { M, N }, extremaIndices );
	std::tie( sextmax, sextmin ) = extremosDiag( iext, jext, aData, -1 );

	// Peaks on columns, rows and diagonals:
	smax = intersect( smax, sextmax, { 0, ( N * M ) - 1 } );
	smin = intersect( smin, sextmin, { 0, ( N * M ) - 1 } );

	auto filter = [&dataPtr, &M]( const IndicesType& indices ) {
		ContType temp;
		temp.reserve( indices.size() );
		int i, j;
		for ( auto& idx : indices )
		{
			std::tie( i, j ) = MMathematicalUtils::linerIndexToSubsrcipt( (int)idx, M );
			temp.emplace_back( dataPtr[j + i * M] );
		}

		return temp;
	};

	// MExtrema points:
	ContType xymax( filter( smax ) );
	ContType xymin( filter( smin ) );
	std::tie( xymax, smax ) =
		sortDataAndTheirIndices<std::greater<float>>( xymax, smax );
	std::tie( xymin, smin ) = sortDataAndTheirIndices( xymin, smin );

	return std::make_tuple( xymax, smax, xymin, smin );
}

std::pair<MExtrema::Indices2DType, MExtrema::Indices2DType>
	MExtrema::extremos( const mcp::MVoxelData& aData )
{
	const int	  M = aData.sizeX();
	const int	  N = aData.sizeY();
	Indices2DType smax, smin;

	const float* dataPtr = aData.slice( 0 );
#pragma omp for
	for ( int n = 0; n < M; ++n )
	{
		ContType	column;
		IndicesType imaxfil, iminfil;
		column.reserve( N );
		for ( int i = 0; i < N; ++i )
		{
			column.emplace_back( dataPtr[n + i * M] );
		}

		std::tie( std::ignore, imaxfil, std::ignore, iminfil ) =
			extrema( std::move( column ) );

		if ( imaxfil.size() ) // Maxima indexes
		{
			for ( const auto& idx : imaxfil )
			{
				smax.emplace_back(
					mcp::MIndex2D<std::uint64_t>{ idx, (std::uint64_t)n } );
			}
		}

		if ( iminfil.size() ) // Minima indexes
		{
			for ( const auto& idx : iminfil )
			{
				smin.emplace_back(
					mcp::MIndex2D<std::uint64_t>{ idx, (std::uint64_t)n } );
			}
		}
	}

	return std::make_pair( smax, smin );
}

MExtrema::ReturnType MExtrema::extrema( ContType&& aData )
{
	IndicesType imax, imin;

	std::size_t Nt = aData.size();
	if ( !Nt )
	{
		LOG_ERROR( "extrema: dimension were zero" );
		return {};
	}

	IndicesType indx;
	indx.reserve( Nt );
	for ( int i = 0; i < Nt; ++i )
	{
		indx.emplace_back( i );
	}

	auto isNanIt = std::find_if( aData.begin(), aData.end(), []( const auto& aValue ) {
		return std::isnan( aValue );
	} );
	bool hasNat = false;
	while ( isNanIt != aData.end() )
	{
		if ( !hasNat )
		{
			hasNat = true;
		}

		--Nt;
		int dist = std::distance( aData.begin(), isNanIt );
		for ( std::size_t i = dist; i < indx.size() - 1; ++i )
		{
			std::swap( indx[i], indx[i + 1] );
		}

		indx.pop_back();

		auto it = aData.erase( isNanIt );
		isNanIt = std::find_if( it, aData.end(), []( const auto& aValue ) {
			return std::isnan( aValue );
		} );
	}

	ContType dx;
	dx.reserve( Nt - 1 );

	bool horizontalLine = true;

	IndicesType a;

	for ( std::uint64_t i = 0; i < Nt - 1; ++i )
	{
		const float diff = aData[i + 1] - aData[i];
		dx.emplace_back( diff );
		if ( horizontalLine && diff != 0 )
		{
			horizontalLine = false;
		}

		if ( diff != 0 )
		{
			a.emplace_back( i ); // Indexes where x changes
		}
	}

	if ( horizontalLine )
	{
		LOG_ERROR( "extrema: data was a horizontalLine" );
		return {};
	}

	// Flat peaks? Put the middle element:
	IndicesType lm;
	for ( std::uint64_t i = 0; i < a.size() - 1; ++i )
	{
		if ( a[i + 1] - a[i] != 1 ) // Indexes where a do not changes
		{
			lm.emplace_back( i + 1 );
		}
	}

	for ( auto& idx : lm )
	{
		const std::uint64_t d = a[idx] - a[idx - 1]; // Number of elements in the flat peak
		a[idx] -= std::floor( d / 2.0 );   // Save middle elements
	}

	a.emplace_back( Nt - 1 );

	// Peaks?
	std::vector<bool> b;
	b.reserve( a.size() - 1 );
	for ( std::size_t i = 0; i < a.size() - 1; ++i )
	{
		const float xa = aData[a[i + 1]] - aData[a[i]]; // Serie without flat peaks
		b.emplace_back( xa > 0 ); // 1  =>  positive slopes (minima begin)
								  // 0  =>  negative slopes (maxima begin)
	}

	for ( std::size_t i = 0; i < b.size() - 1; ++i )
	{
		const int xb = b[i + 1] - b[i]; // -1 =>  maxima indexes (but one)
										// +1 =>  minima indexes (but one)

		if ( xb == -1 ) // maxima indexes
		{
			imax.emplace_back( a[i + 1] );
		}
		else if ( xb == 1 ) // minima indexes
		{
			imin.emplace_back( a[i + 1] );
		}
	}

	std::size_t nmaxi = imax.size();
	std::size_t nmini = imin.size();

	ContType xmax, xmin;

	// Maximum or minumim on a flat peak at the ends?
	if ( !nmaxi && !nmini )
	{
		xmax.reserve( 2 );
		xmin.reserve( 2 );
		imax.reserve( 2 );
		imin.reserve( 2 );

		if ( aData.front() > aData.back() )
		{
			xmax.emplace_back( aData.front() );
			imax.emplace_back( indx.front() );
			xmin.emplace_back( aData.back() );
			imin.emplace_back( indx.back() );
		}
		else if ( aData.front() < aData.back() )
		{
			xmax.emplace_back( aData.back() );
			imax.emplace_back( indx.back() );
			xmin.emplace_back( aData.front() );
			imin.emplace_back( indx.front() );
		}

		return std::make_tuple( xmax, imax, xmin, imin );
	}

	// Maximum or minimum at the ends?
	if ( !nmaxi )
	{
		imax.emplace_back( 0 );
		imax.emplace_back( Nt - 1 );
	}
	else if ( !nmini )
	{
		imin.emplace_back( 0 );
		imin.emplace_back( Nt - 1 );
	}
	else
	{
		auto shiftRight = []( IndicesType& cont ) {
			std::rotate( cont.begin(), std::prev( cont.end() ), cont.end() );

			// last element currently at front
			cont.emplace_back( cont.front() );
			cont.front() = 0;
		};

		if ( imax.front() < imin.front() )
		{
			shiftRight( imin );
		}
		else
		{
			shiftRight( imax );
		}

		if ( imax.back() > imin.back() )
		{
			imin.emplace_back( Nt - 1 );
		}
		else
		{
			imax.emplace_back( Nt - 1 );
		}
	}

	xmax.reserve( imax.size() );
	for ( auto& idx : imax )
	{
		xmax.emplace_back( aData[idx] );
	}

	xmin.reserve( imin.size() );
	for ( auto& idx : imin )
	{
		xmin.emplace_back( aData[idx] );
	}

	// NaN's:
	if ( hasNat )
	{
		auto filterIndicies = [&indx]( const IndicesType& indicies ) {
			IndicesType temp;
			temp.reserve( indicies.size() );
			for ( auto& i : indicies )
			{
				temp.emplace_back( indx[i] );
			}

			return temp;
		};

		imax = filterIndicies( imax );
		imin = filterIndicies( imin );
	}

	// Descending order:
	std::tie( xmax, imax ) =
		sortDataAndTheirIndices<std::greater<float>>( xmax, imax );
	std::tie( xmin, imin ) = sortDataAndTheirIndices( xmin, imin );

	return std::make_tuple( xmax, imax, xmin, imin );
}

std::pair<MExtrema::IndicesType, MExtrema::IndicesType> MExtrema::extremosDiag(
	IndicesType& aIext, const IndicesType& aJext, const mcp::MVoxelData& aData, int aA )
{
	// Peaks through diagonals (down-up aA=-1)

	const std::uint64_t M = aData.sizeX();
	const std::uint64_t N = aData.sizeY();

	if ( M < 2 || N < 1 )
	{
		LOG_INFO( "extremosDiag: dimensions were smaller than expected" );
		return {};
	}

	if ( aA == -1 )
	{
		for ( auto& i : aIext )
		{
			i = ( M - 2 ) - i + 1;
		}
	}

	IndicesType iini, jini;
	std::tie( iini, jini ) = cruce( aIext, aJext, 0, 0 );

	const IndicesType temp = MMathematicalUtils::sub2ind( { M, N }, iini, jini );

	const std::set<std::uint64_t> uniqueIndices( temp.begin(), temp.end() );

	std::tie( iini, jini ) = MMathematicalUtils::ind2sub( { M, N }, uniqueIndices );

	IndicesType ifin, jfin;
	std::tie( ifin, jfin ) = cruce( iini, jini, M - 1, N - 1 );

	IndicesType sextmax;
	IndicesType sextmin;

	const float* dataPtr = aData.slice( 0 );

#pragma omp for
	for ( int n = 0; n < (int)iini.size(); ++n )
	{
		IndicesType ises;
		IndicesType jses;

		for ( std::uint64_t i = iini[n]; i < ifin[n] + 1; ++i )
		{
			ises.emplace_back( i );
		}

		for ( std::uint64_t i = jini[n]; i < jfin[n] + 1; ++i )
		{
			jses.emplace_back( i );
		}

		if ( aA == -1 )
		{
			std::transform( ises.begin(), ises.end(), ises.begin(), [&M]( mint v ) {
				return ( M - 2 ) - v + 1;
			} );
		}

		ContType extremaData;
		extremaData.reserve( ises.size() );

		const IndicesType s = MMathematicalUtils::sub2ind( { M, N }, ises, jses );

		for ( std::uint64_t i = 0; i < ises.size(); ++i )
		{
			extremaData.emplace_back( dataPtr[jses[i] + ises[i] * M] );
		}

		IndicesType imax, imin;
		std::tie( std::ignore, imax, std::ignore, imin ) =
			extrema( std::move( extremaData ) );

		for ( auto& i : imax )
		{
			sextmax.emplace_back( s[i] );
		}

		for ( auto& i : imin )
		{
			sextmin.emplace_back( s[i] );
		}
	}

	return { sextmax, sextmin };
}

std::pair<MExtrema::IndicesType, MExtrema::IndicesType>
	MExtrema::cruce( const IndicesType& aI0, const IndicesType& aJ0, mint aI, mint aJ )
{
	// Indexes where the diagonal of the element io,jo crosses the left/superior
	// (aI=0,aJ=0) or right/inferior (aI=M - 1,aJ=N - 1) side of an MxN matrix.

	const mint arriba = 2 * ( aI * aJ == 0 ) - 1;

	const std::size_t i0Size = aI0.size();
	IndicesType		  i;
	IndicesType		  j;
	i.reserve( i0Size );
	j.reserve( i0Size );

	for ( std::size_t n = 0; n < i0Size; ++n )
	{
		const bool si = ( arriba * ( aJ0[n] - aJ ) ) > ( arriba * ( aI0[n] - aI ) );
		i.emplace_back( ( aI - ( aJ + aI0[n] - aJ0[n] ) ) * si + aJ + aI0[n] - aJ0[n] );
		j.emplace_back( ( aI + aJ0[n] - aI0[n] - aJ ) * si + aJ );
	}

	return std::make_pair( i, j );
}

MExtrema::IndicesType
	MExtrema::intersect( IndicesType& aLhs, IndicesType& aRhs, const IndicesType& aExtra )
{
	std::sort( aLhs.begin(), aLhs.end() );
	std::copy( aExtra.begin(), aExtra.end(), std::back_inserter( aRhs ) );
	std::sort( aRhs.begin(), aRhs.end() );

	IndicesType rv;
	rv.reserve( aLhs.size() );

	std::set_intersection(
		aLhs.begin(),
		aLhs.end(),
		aRhs.begin(),
		aRhs.end(),
		std::back_inserter( rv ) );

	return rv;
}
} // namespace msc