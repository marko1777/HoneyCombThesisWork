/*!
 * \file
 * MExtrema class definition. This file is part of DataFilters module of MSC SDK.
 *
 * \remarks
 * Copyright (C) 2018 Mediso. All rights reserved.
 *
 * \authors
 * pmark
 */

#pragma once

#include <DataFilters2/include/Export.h>
#include <MDR/include/MIndex.h>
#include <MPlatform/include/IntegralTypes.h>
#include <tuple>
#include <vector>

namespace mcp
{
class MVoxelData;
}

namespace msc
{

/*!	\brief Gets the extrema points from a surface.
 * \details The extrema points are searched only through the column, the row and
 * the diagonals crossing each matrix element
 *
 * This is a port of the matlap implementation originally
 * Written by
 * Lic. on Physics Carlos Adriďż˝n Vargas Aguilera
 * Physical Oceanography MS candidate
 * UNIVERSIDAD DE GUADALAJARA
 * Mexico, 2005
 *
 * You can get see the original one
 * From       : http://www.mathworks.com/matlabcentral/fileexchange
 * File ID    : 12275
 * \todo fix usings not showing up in doxygen
 */
class DataFilters2_API MExtrema
{
public:
	using ContType = std::vector<float>;

	using MatrixType = std::vector<ContType>;

	using IndicesType = std::vector<std::uint64_t>;

	using Indices2DType = std::vector<mcp::MIndex2D<std::uint64_t>>;

	/*! \typedef
	 * \brief Return type
	 *
	 * \tparam extrema max elements
	 * \tparam extrema max indices
	 * \tparam extrema min elements
	 * \tparam extrema min indices
	 */
	using ReturnType = std::tuple<ContType, IndicesType, ContType, IndicesType>;

	/*!
	 * \brief Calculate extrema points
	 *
	 * \param [in] aData: row aligned vector<vector> matrix data
	 * \param [out] [max, maxIndices, min, minIndices] extrema max/min elements
	 * in descendant order and max/min linear indices
	 */
	static ReturnType calculate( const mcp::MVoxelData& aData );

private:
	/*!
	 * \brief Sort data and it's linear in index in the matrix in respect of the data
	 *
	 * \param [in] aData: the data to be sorted from a matrix
	 * \param [in] aDataIndices: the data linear indices in the original matrix
	 * \param [out] [data, indicies] sorted data and linear indices in respect of the data
	 */
	template <typename Tcomp = std::less<double>>
	static std::pair<ContType, IndicesType>
		sortDataAndTheirIndices( const ContType& aData, const IndicesType& aDataIndices );

	/*!
	 * \brief Gets the global extrema points from a time series.
	 *
	 * \param [in] aData: the data to be used for searching the extrema points
	 * \param [out] [max, min] extrema max/min elements in descendant order and max/min linear indices
	 */
	static ReturnType extrema( ContType&& aData );

	/*!
	 * \brief Peaks through columns or rows.
	 *
	 * \param [in] aData: the data to be used for searching the extrema points
	 * \param [out] [maxIndices, minIndices] extrema max/min linear indices
	 */
	static std::pair<Indices2DType, Indices2DType>
		extremos( const mcp::MVoxelData& aData );

	/*!
	 * \brief Indexes where the diagonal of the element io,jo crosses the left/superior.
	 * (aI=0,aJ=0) or right/inferior (aI=M-1,aJ=N-1) side of an MxN matrix.
	 *
	 * \param [in] aI0: the first dimension elements of the subscripts
	 * \param [in] aJ0: the second dimension elements of the subscripts
	 * \param [in] aI: specifies which inferior should be checked
	 * \param [in] aJ: specifies which inferior should be checked
	 * \param [out] Indeces where the diagonal of the element io,jo crosses the left/superior
	 */
	static std::pair<IndicesType, IndicesType>
		cruce( const IndicesType& aI0, const IndicesType& aJ0, mint aI, mint aJ );

	/*!
	 * \brief Peaks through diagonals
	 *
	 * \param [in] aI0: the first dimension elements of the extrema subscripts
	 * \param [in] aJ0: the second dimension elements of the extrema subscripts
	 * \param [in] aData: data which the extrema points were optained
	 * \param [in] aA: direction which the diagonal should be processed.
	 * (up-down aA=1) (down-up aA=-1) \param [out] Indeces where the diagonal of
	 * the element io,jo crosses the left/superior
	 */
	static std::pair<IndicesType, IndicesType> extremosDiag(
		IndicesType& aIext, const IndicesType& aJext,
		const mcp::MVoxelData& aData, int aA );

	/*!
	 * \brief Unique intersection of linear indices
	 *
	 * \param [in] aLhs: left hand side of the intersection
	 * \param [in] aRhs: right hand side of the intersection
	 * \param [in] aExtra: extra indices to be considered
	 * \param [out] [indices] Unique intersection of indices
	 */
	static IndicesType intersect(
		IndicesType& aLhs, IndicesType& aRhs, const IndicesType& aExtra = {} );
};

template <typename Tcomp>
std::pair<MExtrema::ContType, MExtrema::IndicesType>
	MExtrema::sortDataAndTheirIndices( const ContType& aData, const IndicesType& aDataIndices )
{
	std::vector<std::pair<double, size_t>> dataWithIndex;
	dataWithIndex.reserve( aData.size() );

	for ( size_t i = 0; i < aData.size(); ++i )
	{
		dataWithIndex.emplace_back( std::make_pair( aData[i], i ) );
	}

	std::sort(
		dataWithIndex.begin(),
		dataWithIndex.end(),
		[]( const std::pair<double, size_t>& lhs, const std::pair<double, size_t>& rhs ) {
			return Tcomp()( lhs.first, rhs.first );
		} );

	ContType sortedData;
	sortedData.reserve( aData.size() );

	for ( auto& pair : dataWithIndex )
	{
		sortedData.emplace_back( pair.first );
	}

	IndicesType sortedIdata;
	sortedData.reserve( aDataIndices.size() );

	for ( auto& pair : dataWithIndex )
	{
		sortedIdata.emplace_back( aDataIndices[pair.second] );
	}

	return { sortedData, sortedIdata };
}

} // namespace msc