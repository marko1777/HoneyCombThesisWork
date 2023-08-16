/*!
 * \file
 * Member function definitions for MHoneyCombFilter class. This file is part of
 * DataFilters module of MSC SDK.
 *
 * \remarks
 * Copyright (C) 2018 Mediso. All rights reserved.
 *
 * \authors
 * pmark
 */

#include <DataFilters2/include/MExtrema.h>
#include <DataFilters2/include/MFourierTransform.h>
#include <DataFilters2/include/MFspecial.h>
#include <DataFilters2/include/MHoneyCombFilter.h>
#include <Log/include/Log.h>
#include <MPlatform/include/IntegralTypes.h>
#include <MPlatform/include/VectorTypes.h>
#include <QString>
#include <ScientificCore/include/MMathematicalUtils.h>
#include <ScientificCore/include/MRoipoly.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>

namespace msc
{

LOG_INIT( "MSC_SDK.MHoneyCombFilter" );

// MxN-> from (N,N) take a MxM subMatrix
// http://www.songho.ca/dsp/convolution/convolution.html#convolution_2d
void conv2d( const mcp::MVoxelData& aImage, const mcp::MVoxelData& aMask, mcp::MVoxelData& aResult )
{
	int m = aImage.sizeX(); // == offset[0]
	int n = aImage.sizeY();

	int kCols = aMask.sizeX();
	int kRows = aMask.sizeY();

	int kCenterX = kCols / 2;
	int kCenterY = kRows / 2;

	aResult.resize( aImage.sizeX(), aImage.sizeX(), 1 );
	aResult.fill( 0.f );

	float*			   resultPtr = aResult.slice( 0 );
	const float* const imagePtr	 = aImage.slice( 0 );
	const float* const maskPtr	 = aMask.slice( 0 );

	for ( int i = 0; i < m; ++i ) // rows
	{
		for ( int j = 0; j < n; ++j ) // columns
		{
			for ( int k = 0; k < kRows; ++k ) // kernel rows
			{
				int mm = kRows - 1 - k; // row index of flipped kernel

				for ( int l = 0; l < kCols; ++l ) // kernel columns
				{
					int nn = kCols - 1 - l; // column index of flipped kernel

					// index of input signal, used for checking boundary
					int ii = i + ( kCenterY - mm );
					int jj = j + ( kCenterX - nn );
					// ignore input samples which are out of bound
					if ( ii >= 0 && ii < m && jj >= 0 && jj < n )
					{
						resultPtr[i + ( j * m )] +=
							imagePtr[ii + ( jj * m )] * maskPtr[mm + nn * kRows];
					}
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------

MHoneyCombFilter::MHoneyCombFilter()
		: MFilter( false, true, false )
{
}

std::vector<std::vector<float>> MHoneyCombFilter::voxelToMatrix( const mcp::MVoxelData& aData )
{
	const mint dimX = aData.sizeX();
	const mint dimY = aData.sizeY();

	std::vector<std::vector<float>> matrixData( dimX, std::vector<float>( dimY ) );

	for ( mint j = 0; j < dimY; ++j )
	{
		for ( mint i = 0; i < dimX; ++i )
		{
			matrixData[i][j] = aData.at( i, j );
		}
	}

	return matrixData;
}

mcp::MVoxelData
	MHoneyCombFilter::matrixToVoxel( const std::vector<std::vector<float>>& aData )
{
	const mint dimX = aData[0].size();
	const mint dimY = aData.size();

	mcp::MVoxelData voxelData( dimX, dimY, 1 );

	for ( mint j = 0; j < dimY; ++j )
	{
		for ( mint i = 0; i < dimX; ++i )
		{
			voxelData.at( i, j ) = aData[i][j];
		}
	}

	return voxelData;
}

void MHoneyCombFilter::transformVoxelData( const std::function<void( mint, mint )>& aFunction )
{
	for ( mint j = 0; j < m_DimY; ++j )
	{
		for ( mint i = 0; i < m_DimX; ++i )
		{
			aFunction( i, j );
		}
	}
}

bool MHoneyCombFilter::applyFilter( mcp::MVoxelData& aData )
{
	LOG_INFO( "applyFilter" );
	m_DimX = aData.sizeX();
	m_DimY = aData.sizeY();

	if ( !m_DimX || !m_DimY )
	{
		LOG_ERROR( "x or y dimension was 0" );
		return false;
	}

	std::unique_ptr<mcp::MVoxelData> smoothedFreqDomainData = preprocessData( aData );

	if ( !smoothedFreqDomainData )
	{
		LOG_ERROR( "preprocessData failed" );
		return false;
	}

	std::vector<float>		   maxValues;
	std::vector<std::uint64_t> maxLinearIndices;
	std::tie( maxValues, maxLinearIndices, std::ignore, std::ignore ) =
		MExtrema::calculate( *smoothedFreqDomainData );

	if ( maxValues.empty() || maxLinearIndices.empty() )
	{
		LOG_ERROR( "no extrema values were found" );
		return false;
	}

	std::tie( m_ExtremaIndicesX, m_ExtremaIndicesY ) = MMathematicalUtils::ind2sub(
		{ smoothedFreqDomainData->sizeX(), smoothedFreqDomainData->sizeY() },
		maxLinearIndices );

	std::unique_ptr<std::vector<mcp::MIndex2D<std::uint64_t>>> honeyCombIndices =
		collectHighFrequencyRegionsCenter( maxValues );
	if ( !honeyCombIndices )
	{
		LOG_ERROR( "collectHighFrequencyRegionsCenter failed" );
		return false;
	}

	std::unique_ptr<std::vector<std::vector<bool>>> mask =
		drawAroundHighFrequencyRegions( *honeyCombIndices, *smoothedFreqDomainData );

	if ( !mask )
	{
		LOG_ERROR( "drawAroundHighFrequencyRegions failed" );
		return false;
	}

	mcp::MVoxelData finalizedFilter =
		createAndConvoluteFilter( *smoothedFreqDomainData, *mask );

	filterImage( aData, finalizedFilter );

	return true;
}

std::unique_ptr<mcp::MVoxelData>
	MHoneyCombFilter::preprocessData( const mcp::MVoxelData& aData )
{
	LOG_INFO( "prerpocessData" );

	transformDataToFrequencyDomain( aData );

	m_FrequencyDomainData =
		MMathematicalUtils::fftshift( m_FrequencyDomainData, m_DimX, m_DimY );

	std::unique_ptr<mcp::MVoxelData> smoothedFreqDomainData(
		new mcp::MVoxelData( dataSmoothing() ) );

	const bool validRegionFound =
		findAndMaskLowFrequencyRegion( *smoothedFreqDomainData );

	if ( validRegionFound )
	{
		LOG_INFO( "valid low frequency region found" );
		return smoothedFreqDomainData;
	}
	else
	{
		LOG_ERROR( "No valid low frequency region found" );
		return nullptr;
	}
}

void MHoneyCombFilter::transformDataToFrequencyDomain( const mcp::MVoxelData& aData )
{
	LOG_INFO( "transformDataToFrequencyDomain" );

	m_FrequencyDomainData.reserve( m_DimX * m_DimY );

	auto inserter = std::back_inserter( m_FrequencyDomainData );

	transformVoxelData( [&aData, &inserter]( mint i, mint j ) {
		inserter++ = aData.at( i, j );
	} );

	MFourierTransform fftObject;
	fftObject.forward2DInPlace( m_FrequencyDomainData, std::log2( m_DimX ), std::log2( m_DimY ) );
}

mcp::MVoxelData MHoneyCombFilter::dataSmoothing()
{
	LOG_INFO( "dataSmoothing" );

	mcp::MVoxelData highFreqRegionsLifted( m_DimX, m_DimY, 1 );

	auto fftDataIt = m_FrequencyDomainData.begin();
	transformVoxelData( [&highFreqRegionsLifted, fftDataIt]( mint i, mint j ) mutable {
		highFreqRegionsLifted.at( i, j ) = (float)std::log( std::abs( *fftDataIt++ ) );
	} );

	m_Dimc = m_DimX / 8;

	mcp::MVoxelData smoothedFreqDomainData( m_DimX, m_DimY, 1 );

	conv2d( highFreqRegionsLifted, arrayToVoxel( MFspecial::disk<2>() ), smoothedFreqDomainData );

	return smoothedFreqDomainData;
}

bool MHoneyCombFilter::findAndMaskLowFrequencyRegion( mcp::MVoxelData& aData )
{
	LOG_INFO( "findAndMaskLowFrequencyRegion" );
	const bool validRegionFound = findLowFrequencyRegionIn( aData );

	if ( validRegionFound )
	{
		maskLowFrequencyRegionIn( aData );
	}

	return validRegionFound;
}

bool MHoneyCombFilter::findLowFrequencyRegionIn( const mcp::MVoxelData& aData )
{
	LOG_INFO( "findLowFrequencyRegionIn" );
	m_Avg = averageOfVoxelData( aData );

	const mint dimcPlus = m_DimX / 32;

	const float m0 = aData.at( m_DimX / 2 - 1, m_DimY / 2 - 1 );
	float m = m0;

	mint pos = m_DimX / 2 - 1 + dimcPlus;

	while ( m > ( m_Avg + ( m0 - m_Avg ) / 4 ) )
	{
		m =	( aData.at( m_DimX / 2 - 1, pos ) + aData.at( pos, m_DimX / 2 - 1 ) ) / 2.f;
		++pos;
	}

	m_Dimc = pos - m_DimX / 2;

	const bool validRegionFound = pos < m_DimX;

	return validRegionFound;
}

void MHoneyCombFilter::maskLowFrequencyRegionIn( mcp::MVoxelData& aData )
{
	LOG_INFO( "maskLowFrequencyRegionIn" );

	const mint lb = std::max( (long long)0, m_DimY / 2 - ( m_Dimc + 2 ) );
	const mint ub = std::min( m_DimY, m_DimY / 2 + ( m_Dimc + 1 ) );

	for ( mint j = lb; j < ub; ++j )
	{
		for ( mint i = lb; i < ub; ++i )
		{
			aData.at( i, j ) = m_Avg;
		}
	}
}

float MHoneyCombFilter::averageOfVoxelData( const mcp::MVoxelData& aData )
{
	float sum = 0.0;

	for ( mint i = 0; i < m_Dimc; ++i )
	{
		for ( mint j = 0; j < m_Dimc; ++j )
		{
			sum += aData.at( i, j );
		}
	}

	return sum / float( m_Dimc * m_Dimc );
}

std::unique_ptr<std::vector<mcp::MIndex2D<std::uint64_t>>>
	MHoneyCombFilter::collectHighFrequencyRegionsCenter( const std::vector<float>& aExtremaMaxValues )
{
	LOG_INFO( "collectHighFrequencyRegionsCenter" );
	const float extremaMax =
		*std::max_element( aExtremaMaxValues.begin(), aExtremaMaxValues.end() );
	int highFreqRegionCount = 0;

	std::unique_ptr<std::vector<mcp::MIndex2D<std::uint64_t>>> centerIndices(
		new std::vector<mcp::MIndex2D<std::uint64_t>> );

	const std::uint64_t iterations = 40;

	std::vector<float> centerValueLowerBounds;
	centerValueLowerBounds.reserve( iterations );

	for ( int it = 1; it < iterations && highFreqRegionCount != 6; ++it )
	{
		const float centerValueLowerBound =
			m_Avg + ( extremaMax - m_Avg ) * ( 1.f - 0.025f * (float)it );

		highFreqRegionCount = findHighFrequencyRegionsCenter(
			*centerIndices,
			aExtremaMaxValues,
			centerValueLowerBound );

		centerValueLowerBounds.emplace_back( centerValueLowerBound );
	}

	if ( highFreqRegionCount != 6 )
	{
		highFreqRegionCount = 0;

		for ( int it = 1; it < iterations && highFreqRegionCount != 4; ++it )
		{
			highFreqRegionCount = findHighFrequencyRegionsCenter(
				*centerIndices,
				aExtremaMaxValues,
				centerValueLowerBounds[it - 1] );
		}
	}

	if ( highFreqRegionCount == 6 || highFreqRegionCount == 4 )
	{
		return centerIndices;
	}
	else
	{
		LOG_ERROR( QString( "%1 piece of high frequency region was found" ).arg( highFreqRegionCount ) );
		return nullptr;
	}
}

int MHoneyCombFilter::findHighFrequencyRegionsCenter(
	std::vector<mcp::MIndex2D<std::uint64_t>>& aCenterIndices,
	const std::vector<float>& aExtremaMaxValues, float aCenterValueLowerBound )
{
	aCenterIndices.clear();

	for ( std::size_t i = 0; i < m_ExtremaIndicesX.size(); ++i )
	{
		if ( ( (int)m_ExtremaIndicesY[i] < m_DimY - m_Dimc / 4 ) &&
			 ( (int)m_ExtremaIndicesY[i] > m_Dimc / 4 ) &&
			 ( (int)m_ExtremaIndicesX[i] < m_DimX - m_Dimc / 4 ) &&
			 ( (int)m_ExtremaIndicesX[i] > m_Dimc / 4 ) &&
			 ( aExtremaMaxValues[i] > aCenterValueLowerBound ) )
		{
			aCenterIndices.emplace_back( m_ExtremaIndicesX[i], m_ExtremaIndicesY[i] );
		}
	}

	std::uint64_t pars = 0;

	for ( std::uint64_t i = 0; i < aCenterIndices.size(); ++i )
	{
		for ( std::uint64_t j = 0; j < aCenterIndices.size(); ++j )
		{
			if ( i != j )
			{
				const std::uint64_t diffX =
					aCenterIndices[i].x() - aCenterIndices[j].x();
				const std::uint64_t diffY =
					aCenterIndices[i].y() - aCenterIndices[j].y();

				const float r = (float)std::sqrt( diffX * diffX + diffY * diffY );

				if ( r < m_Dimc / 2.0 )
				{
					++pars;
				}
			}
		}
	}

	return (int)aCenterIndices.size() - (int)pars / 2;
}

std::unique_ptr<std::vector<std::vector<bool>>> MHoneyCombFilter::drawAroundHighFrequencyRegions(
	const std::vector<mcp::MIndex2D<std::uint64_t>>& aHoneyCombIndices,
	const mcp::MVoxelData&							 aData )
{
	const std::uint64_t		  degree = 10;
	const std::vector<float> phi =
		MMathematicalUtils::linspace<float>( 0, 2 * PIValue, degree );

	std::vector<float> cosX;
	std::vector<float> sinY;
	cosX.reserve( degree );
	sinY.reserve( degree );

	std::for_each( phi.begin(), phi.end(), [&cosX, &sinY]( float aV ) {
		cosX.emplace_back( std::cos( aV ) );
		sinY.emplace_back( std::sin( aV ) );
	} );

	std::vector<mcp::MIndex2D<std::uint64_t>> circleIndicesAroundHoneyCombIndices;
	circleIndicesAroundHoneyCombIndices.resize( degree );

	std::vector<mcp::MIndex2D<std::uint64_t>>		circlePoints;
	std::unique_ptr<std::vector<std::vector<bool>>> mask(
		new std::vector<std::vector<bool>>( m_DimX, std::vector<bool>( m_DimY ) ) );

	std::vector<mcp::MIndex2D<std::uint64_t>> roiIndices;
	roiIndices.reserve( circlePoints.size() );

	MRoipoly roipoly( m_DimX, m_DimY );

	for ( std::size_t j = 0; j < aHoneyCombIndices.size(); ++j )
	{
		circlePoints.clear();
		circlePoints.resize( degree );

		const float regionHalfHeight =
			( m_Avg +
			  ( aData.at( aHoneyCombIndices[j].y(), aHoneyCombIndices[j].x() ) - m_Avg ) *
				  0.5f );

		for ( mint R = 0; R < m_Dimc; ++R )
		{
			for ( std::size_t i = 0; i < degree; ++i )
			{
				const float xValue = std::round( float( R + 1 ) * cosX[i] ) + (float)aHoneyCombIndices[j].x();
				const float yValue = std::round( float( R + 1 ) * sinY[i] ) + (float)aHoneyCombIndices[j].y();
				circleIndicesAroundHoneyCombIndices[i] = {
					std::uint64_t( xValue < 0 ? 0 : xValue ),
					std::uint64_t( yValue < 0 ? 0 : yValue ) };
			}

			for ( std::uint64_t i = 0; i < degree; ++i )
			{
				if ( circlePoints[i].y() )
				{
					continue;
				}

				const std::uint64_t xCircleIndex =
					circleIndicesAroundHoneyCombIndices[i].x();
				const std::uint64_t yCircleIndex =
					circleIndicesAroundHoneyCombIndices[i].y();

				float circleIdxValue = 0.0;

				if ( xCircleIndex < m_DimX && yCircleIndex < m_DimY )
				{
					circleIdxValue = aData.at( yCircleIndex, xCircleIndex );
				}

				if ( circleIdxValue < regionHalfHeight )
				{
					circlePoints[i] = { xCircleIndex, yCircleIndex };
				}
			}
		}

		for ( int i = 0; i < degree; ++i )
		{
			if ( !circlePoints[i].x() || !circlePoints[i].y() )
			{
				LOG_ERROR(
					QString( "When checking (%1, %2) honeyComb coordinates in "
							 "circle point: (%3, %4) values doesn't "
							 "satifisfied the conditions" )
						.arg( aHoneyCombIndices[j].x() )
						.arg( aHoneyCombIndices[j].y() )
						.arg( circleIndicesAroundHoneyCombIndices[i].x() )
						.arg( circleIndicesAroundHoneyCombIndices[i].y() ) );
				return nullptr;
			}
		}

		roiIndices.clear();

		for ( std::size_t i = 0; i < circlePoints.size(); ++i )
		{
			roiIndices.emplace_back( circlePoints[i].x(), circlePoints[i].y() );
		}

		const std::vector<mcp::MIndex2D<std::uint64_t>> indicesInRoi =
			roipoly.mask( roiIndices );

		for ( const auto& point : indicesInRoi )
		{
			( *mask )[point.x()][point.y()] = true;
		}
	}

	return mask;
}

mcp::MVoxelData MHoneyCombFilter::createAndConvoluteFilter(
	const mcp::MVoxelData& aData, const std::vector<std::vector<bool>>& aMask )
{
	LOG_INFO( "createAndConvoluteFilter" );

	mcp::MVoxelData filter( m_DimX, m_DimY, 1 );

	transformVoxelData( [&]( mint i, mint j ) {
		filter.at( i, j ) = m_Avg + aData.at( i, j ) * (float)aMask[j][i];
	} );

	mcp::MVoxelData convedFilter;
	conv2d( filter, arrayToVoxel( MFspecial::disk<1>() ), convedFilter );

	transformVoxelData( [&convedFilter]( mint i, mint j ) {
		convedFilter.at( i, j ) = std::exp( convedFilter.at( i, j ) );
	} );

	maskBorderEdgesAndFinalizeFilter( convedFilter );

	return convedFilter;
}

void MHoneyCombFilter::maskBorderEdgesAndFinalizeFilter( mcp::MVoxelData& aFilter )
{
	LOG_INFO( "maskBorderEdgesAndFinalizeFilter" );

	auto maskYcords = [&aFilter]( mint i, mint j, mint end ) {
		const float val = aFilter.at( 1, 1 );

		for ( ; j < end; ++j )
		{
			aFilter.at( i, j ) = val;
		}
	};

	maskYcords( 0, 0, m_DimX );
	maskYcords( m_DimX - 1, 0, m_DimX );

	auto maskXcords = [&aFilter]( mint i, mint end, mint j ) {
		const float val = aFilter.at( 1, 1 );

		for ( ; i < end; ++i )
		{
			aFilter.at( i, j ) = val;
		}
	};

	maskXcords( 0, m_DimX, 0 );
	maskXcords( 0, m_DimX, m_DimX - 1 );

	float min = std::numeric_limits<float>::max();

	transformVoxelData( [&aFilter, &min]( mint i, mint j ) {
		if ( aFilter.at( i, j ) < min )
		{
			min = aFilter.at( i, j );
		}
	} );

	transformVoxelData( [&aFilter, &min]( mint i, mint j ) {
		aFilter.at( i, j ) = 1.f / ( aFilter.at( i, j ) - min + 1 );
	} );
}

void MHoneyCombFilter::filterImage( mcp::MVoxelData& aImage, const mcp::MVoxelData& aFilter )
{
	LOG_INFO( "filterImage" );

	std::vector<std::complex<double>> filtImg;

	filtImg.reserve( m_FrequencyDomainData.size() );

	auto fftDataIt = m_FrequencyDomainData.begin();
	transformVoxelData( [&aFilter, &filtImg, &fftDataIt]( mint i, mint j ) {
		filtImg.emplace_back(
			static_cast<double>( aFilter.at( i, j ) ) * ( *( fftDataIt++ ) ) );
	} );

	filtImg = MMathematicalUtils::ifftshift( filtImg, m_DimX, m_DimY );

	MFourierTransform fftObject;
	fftObject.inverse2DInPlace( filtImg, std::log2( m_DimX ), std::log2( m_DimY ) );

	float Si = 0.0;
	transformVoxelData(
		[&aImage, &Si]( mint i, mint j ) { Si += aImage.at( i, j ); } );

	float Sf = 0.0;

	for ( const auto& elem : filtImg )
	{
		Sf += std::abs( elem );
	}

	auto filtIt = filtImg.begin();
	transformVoxelData( [&aImage, &filtIt, &Si, &Sf]( mint i, mint j ) {
		aImage.at( i, j ) = float( std::abs( *filtIt++ ) ) * Si / Sf;
	} ); // normalization
}

} // namespace msc
