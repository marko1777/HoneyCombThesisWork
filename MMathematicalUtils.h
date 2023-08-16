#pragma once

#include <ScientificCore/include/Export.h>
#include <array>
#include <cstdint>
#include <set>
#include <tuple>
#include <type_traits>
#include <vector>

namespace msc
{

struct ScientificCore_API MMathematicalUtils
{
	/*!
	 * \brief Convert linear index to subsript.
	 *
	 * \param [in] aIdx: linear index to be converted
	 * \param [in] aFirstDIm: the matrix first dimension
	 * \param [out] [x, y] the subscript of the linear index
	 */
	constexpr static std::pair<int, int> linerIndexToSubsrcipt( int aIdx, int aFirstDim )
	{
		int i = aIdx % aFirstDim;
		int j = aIdx / aFirstDim;
		return { i, j };
	}

	/*!
	 * \brief Convert linear indices to subsripts.
	 * https://www.mathworks.com/help/matlab/ref/ind2sub.html
	 *
	 * \param [in] aDim: dimension of the matrix which the indices were obtained
	 * \param [in] aIndices: linear indices to be converted
	 * \param [out] [row, col] the arrays row and col containing the equivalent
	 * row and column subscripts corresponding to the linear indices aIndices
	 * for a matrix of size aDim \todo make cont T param uniform
	 */
	template <typename T>
	static std::enable_if_t<
		std::is_same<T, std::set<std::uint64_t>>::value ||
			std::is_same<T, std::vector<std::uint64_t>>::value,
		std::pair<std::vector<std::uint64_t>, std::vector<std::uint64_t>>>
		ind2sub( const std::pair<int, int>& aDim, const T& aIndices );
	/*!
	 * \brief Convert subscripts to linear indices. https://www.mathworks.com/help/matlab/ref/sub2ind.html
	 *
	 * \param [in] aDim: the matrix dimensions where the subscripts were optained
	 * \param [in] aRow: the first dimension subscripts
	 * \param [in] aCol: the second dimension subscripts
	 * \param [in] aDim: the matrix dimensions where the subscripts were optained
	 * \param [out] [indices] linear indices corresponding to the row and column
	 * subscripts in aRow and aCol for a matrix of size aDim.
	 */
	static std::vector<std::uint64_t> sub2ind(
		const std::pair<int, int>& aDim, const std::vector<std::uint64_t>& aRow,
		const std::vector<std::uint64_t>& aCol );

	/*!
	 * \brief 2D grid https://www.mathworks.com/help/matlab/ref/meshgrid.html
	 * \param [in] aX: coordinates to be used
	 * \param [in] aY: coordinates to be used
	 * \param [out] [X, Y]: X is a matrix where each row is a copy of aX; Y is a
	 * matrix where each column is a copy of aY
	 */
	template <size_t tLhsSize, size_t tRhsSize>
	constexpr static std::pair<
		std::array<std::array<double, tLhsSize>, tRhsSize>,
		std::array<std::array<double, tLhsSize>, tRhsSize>>
		meshgrid( const std::array<double, tLhsSize>& aX, const std::array<double, tRhsSize>& aY );

	/*!
	 * \brief Shift zero-frequency component to center of spectrum
	 * https://www.mathworks.com/help/matlab/ref/fftshift.html
	 *
	 * \param [in] aData: data
	 * \param [in] aDImX: first dimension
	 * \param [in] aDImY: second dimension
	 * \param [out] Zero-frequency components in the middle
	 * \todo refactor these so they use the same container/ uniform T and make a
	 * constraint on it to have a subscript operator
	 */
	template <typename T>
	constexpr static std::vector<T>
		fftshift( const std::vector<T>& aData, int aDimX, int aDimY );

	/*!
	 * \brief Inverse zero-frequency shift
	 * https://www.mathworks.com/help/matlab/ref/ifftshift.html
	 *
	 * \param [in] aData: data
	 * \param [in] aDImX: first dimension
	 * \param [in] aDImY: second dimension
	 * \param [out] Inverse zero-frequency shift
	 * \todo refactor these so they use the same container/ uniform T and make a
	 * constraint on it to have a subscript operator
	 */
	template <typename T>
	constexpr static std::vector<T>
		ifftshift( const std::vector<T>& aData, int aDimX, int aDimY );

	/*!
	 * \brief Generate linearly spaced vector https://www.mathworks.com/help/matlab/ref/linspace.html
	 *
	 * \param [in] aBegin: start of the interval
	 * \param [in] aEnd: end of the interval
	 * \param [in] aLength: length of the interval
	 * \param [out] aLength lenght interval from aBegin to aEnd
	 */
	template <typename Trv, typename Tbegin, typename Tend>
	static std::vector<Trv> linspace( Tbegin aBegin, Tend aEnd, std::uint64_t aLength );
};

template <size_t tLhsSize, size_t tRhsSize>
constexpr std::pair<std::array<std::array<double, tLhsSize>, tRhsSize>, std::array<std::array<double, tLhsSize>, tRhsSize>>
	MMathematicalUtils::meshgrid(
		const std::array<double, tLhsSize>& aX, const std::array<double, tRhsSize>& aY )
{
	if ( !tLhsSize || !tRhsSize )
	{
		return {};
	}

	std::array<std::array<double, tLhsSize>, tRhsSize> rx;
	std::array<std::array<double, tLhsSize>, tRhsSize> ry;

	for ( size_t i = 0; i < tRhsSize; ++i )
	{
		rx[i] = aX;
	}

	for ( size_t i = 0; i < tRhsSize; ++i )
	{
		for ( size_t j = 0; j < tLhsSize; ++j )
		{
			ry[i][j] = aY[i];
		}
	}
	return { rx, ry };
}

template <typename T>
std::enable_if_t<
	std::is_same<T, std::set<std::uint64_t>>::value ||
		std::is_same<T, std::vector<std::uint64_t>>::value,
	std::pair<std::vector<std::uint64_t>, std::vector<std::uint64_t>>>
	MMathematicalUtils::ind2sub( const std::pair<int, int>& aDim, const T& aIndices )
{
	std::vector<std::uint64_t> rows;
	rows.reserve( aIndices.size() );

	std::vector<std::uint64_t> cols;
	cols.reserve( aIndices.size() );

	int i, j;

	for ( auto& idx : aIndices )
	{
		std::tie( i, j ) = linerIndexToSubsrcipt( idx, aDim.first );
		rows.emplace_back( i );
		cols.emplace_back( j );
	}

	return { rows, cols };
}

// TODO refactor these so they use the same container/ uniform T and make a constraint on it to has subscript operator
template <typename T>
constexpr std::vector<T>
	MMathematicalUtils::fftshift( const std::vector<T>& aData, int aDimX, int aDimY )
{
	if ( aData.empty() )
	{
		return {};
	}

	std::vector<T> rv( aDimX * aDimY );

	for ( int k = 0; k < ( aDimX * aDimY ); ++k )
	{
		int i;
		int j;
		std::tie( i, j ) = linerIndexToSubsrcipt( k, aDimX );
		i				 = ( i + aDimX / 2 ) % aDimX;
		j				 = ( j + aDimY / 2 ) % aDimY;
		int idx			 = i + ( j * aDimX );

		rv[idx] = aData[k];
	}

	return rv;
}

// TODO refactor as above
template <typename T>
constexpr std::vector<T>
	MMathematicalUtils::ifftshift( const std::vector<T>& aData, int aDimX, int aDimY )
{
	if ( aData.empty() )
	{
		return {};
	}

	std::vector<T> rv( aDimX * aDimY );

	for ( int k = 0; k < ( aDimX * aDimY ); ++k )
	{
		int i;
		int j;
		std::tie( i, j ) = linerIndexToSubsrcipt( k, aDimX );
		i				 = ( aDimX + ( ( i - aDimX / 2 ) % aDimX ) ) % aDimX;
		j				 = ( aDimY + ( ( j - aDimY / 2 ) % aDimY ) ) % aDimY;
		int idx			 = i + ( j * aDimX );

		rv[idx] = aData[k];
	}

	return rv;
}

template <typename Trv, typename Tbegin, typename Tend>
std::vector<Trv> MMathematicalUtils::linspace( Tbegin aBegin, Tend aEnd, std::uint64_t aLength )
{
	if ( !aLength )
	{
		return {};
	}
	else if ( aLength == 1 )
	{
		return std::vector<Trv>{ static_cast<Trv>( aEnd ) };
	}

	std::vector<Trv> rv;
	rv.reserve( aLength );

	double h = ( aEnd - aBegin ) / double( aLength - 1 );

	for ( uint64_t i = 0; i < aLength; ++i )
	{
		rv.emplace_back( aBegin + i * h );
	}

	return rv;
}

} // namespace msc