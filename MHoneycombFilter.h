/*!
 * \file
 * MHoneCombFilter class definition. This file is part of DataFilters module of MSC SDK.
 *
 * \remarks
 * Copyright (C) 2020 Mediso. All rights reserved.
 *
 * \authors
 * pmark
 */
#pragma once

#include <DataFilters2/include/MFilter.h>
#include <MDR/include/MIndex.h>
#include <complex>
#include <functional>
#include <memory>

namespace msc
{

/*!
 * \brief Class for Honey Comb filtering.
 * \todo optimize for speed; separeate bigger functions into classes
 */
class DataFilters2_API MHoneyCombFilter : public MFilter
{
public:
	/*!
	 * \brief Constructor.
	 */
	MHoneyCombFilter();

	/*!
	 * \brief Filter the shapes of the photomultiplier (PMT) from the 2D projection.
	 * This is done by masking the high-frequency regions with a bandpass filter.
	 * \param aData Source voxel data to apply the filter on.
	 * \return True if the filtering was possible and successful, otherwise false.
	 *
	 * See the detailed description of the algorithm : @ref applyFilter
	 */
	bool applyFilter( mcp::MVoxelData& aData ) override;

private:
	/*!
	 * \brief Prepare the data for the creation of the filter.
	 * \param aData The data to be preprocessed
	 * \return The preprocessed data
	 *
	 * See the detailed description of the algorithm : @ref preprocessData
	 */
	std::unique_ptr<mcp::MVoxelData> preprocessData( const mcp::MVoxelData& aData );

	/*!
	 * \brief Apply 2D Fast Fourier transform on data
	 * \param aData The data to be used in the fft2
	 * \return the fft2 data in a 1D column aligned representation of a 2D matrix
	 *
	 * See the detailed description of the algorithm : @ref transformDataToFrequencyDomain
	 */
	void transformDataToFrequencyDomain( const mcp::MVoxelData& aData );

	/*!
	 * \brief Apply logartihm on the fft data so that the high frequency regions
	 * can be more distiungished from the noise
	 \return The data smoothed out
	 *
	 * See the detailed description of the algorithm : @ref dataSmoothing
	 */
	mcp::MVoxelData dataSmoothing();

	/*!
	 * \brief Find and mask out the low frequency region in the center with a
	 * square. This is an important step because it's eases the process of
	 * finding relevant local maxima values. \param aData [in, out] in:
	 * peprocessed data; out: preprocessed data with low frequency region masked
	 * with average
	 *
	 * See the detailed description of the algorithm : @ref findAndMaskLowFrequencyRegion
	 */
	bool findAndMaskLowFrequencyRegion( mcp::MVoxelData& aData );

	/*!
	 * \brief Finds the range of the low frequency region from the center
	 * \param aData peprocessed data
	 *
	 * See the detailed description of the algorithm : @ref findLowFrequencyRegionIn
	 */
	bool findLowFrequencyRegionIn( const mcp::MVoxelData& aData );

	/*!
	 * \brief Mask the low frequency region in the center with the average pixel value in a square with side
	 * length of `m_Dimc * 2`
	 * \param aData [in, out] in: preprocessed data; out preprocessed data with low
	 * frequency region masked with average
	 *
	 * See the detailed description of the algorithm : @ref maskLowFrequencyRegionIn
	 */
	void maskLowFrequencyRegionIn( mcp::MVoxelData& aData );

	/*!
	 * \brief Calculates the average pixel value of the voxelData
	 * \param aData voxel data
	 * \return average of the voxel data
	 *
	 * See the detailed description of the algorithm : @ref averageOfVoxelData
	 */
	float averageOfVoxelData( const mcp::MVoxelData& aData );

	/*!
	 * \brief Colletcting the high frequency regions center corresponding to the shapes of the PMT
	 * \param aExtremaMaxValues the extrema maximum values
	 * \return the 2D indices of the high frequency regions centers
	 *
	 * See the detailed description of the algorithm : @ref collectHighFrequencyRegionsCenter
	 */
	std::unique_ptr<std::vector<mcp::MIndex2D<std::uint64_t>>>
		collectHighFrequencyRegionsCenter( const std::vector<float>& aExtremaMaxValues );

	/*!
	 * \brief Finding of the high frequency regions center corresponding to the shapes of the PMT
	 * \param aCenterIndices the 2D indices of the high frequency regions center
	 * \param aExtremaMaxValues the extrema maximum values
	 * \param aCenterValueLowerBound region center lower bound
	 * \return how many high frequency regions were found?
	 *
	 * See the detailed description of the algorithm : @ref findHighFrequencyRegionsCenter
	 */
	int findHighFrequencyRegionsCenter(
		std::vector<mcp::MIndex2D<std::uint64_t>>& aCenterIndices,
		const std::vector<float>& aExtremaMaxValues, float aCenterValueLowerBound );

	/*!
	 * \brief Drawing around the high frequency regions
	 * \param aHoneyCombIndices 2D indicies of the high frequency regions center
	 * \param aData preprocessed data
	 * \return 2D bandpass filter
	 *
	 * See the detailed description of the algorithm : @ref drawAroundHighFrequencyRegions
	 */
	std::unique_ptr<std::vector<std::vector<bool>>> drawAroundHighFrequencyRegions(
		const std::vector<mcp::MIndex2D<std::uint64_t>>& aHoneyCombIndices,
		const mcp::MVoxelData&							 aData );

	/*!
	 * \brief Create filter and postprocess it
	 * \param aData preprocessed data
	 * \param aMask 2D mask
	 * \return convolutaed filter
	 *
	 * See the detailed description of the algorithm : @ref createAndConvoluteFilter
	 */
	mcp::MVoxelData createAndConvoluteFilter(
		const mcp::MVoxelData& aData, const std::vector<std::vector<bool>>& aMask );

	/*!
	 * \brief Mask border edges and finalize filter
	 * \param aFilter Source voxel data to apply the filter on.
	 * \return True if the filtering was possible and successful, false otherwise.
	 *
	 * See the detailed description of the algorithm : @ref maskBorderEdgesAndFinalizeFilter
	 */
	void maskBorderEdgesAndFinalizeFilter( mcp::MVoxelData& aFilter );

	/*!
	 * \brief Apply filter on the original image
	 * \param aData Source voxel data to apply the filter on.
	 * \param aFilter create filter
	 * \return True if the filtering was possible and successful, false otherwise.
	 *
	 * See the detailed description of the algorithm : @ref filterImage
	 */
	void filterImage( mcp::MVoxelData& aImage, const mcp::MVoxelData& aFilter );

	/*!
	 * \brief Converts voxleData to vector<vector> matrix representation
	 * \param aData source data
	 * \return vector<vector> representation of the voxelData
	 *
	 * See the detailed description of the algorithm : @ref voxelToMatrix
	 */
	std::vector<std::vector<float>> voxelToMatrix( const mcp::MVoxelData& aData );

	/*!
	 * \brief Converts vector<vector> matrix representation to voxelData
	 * \param aData Source data
	 * \return voxelData representation of the data
	 *
	 * See the detailed description of the algorithm : @ref matrixToVoxel
	 */
	mcp::MVoxelData matrixToVoxel( const std::vector<std::vector<float>>& aData );

	/*!
	 * \brief convert array<array> representation to voxelData
	 * \param aData Source data
	 * \return voxelData representation
	 *
	 * See the detailed description of the algorithm : @ref arrayToVoxel
	 */
	template <std::uint64_t tSize>
	mcp::MVoxelData
		arrayToVoxel( const std::array<std::array<double, tSize>, tSize>& aData );

	/*!
	 * \brief Iterates on the plane generated by m_DimX, m_DimY and apply
	 * aFunction on it's indeces
	 *  \param aFunction the function to be applied
	 *
	 * See the detailed description of the algorithm : @ref transformVoxelData
	 */
	void transformVoxelData( const std::function<void( mint, mint )>& aFunction );

private:
	mint   m_DimX{ 0 };
	mint   m_DimY{ 0 };
	mint   m_Dimc{ 0 };
	float m_Avg{ 0.0 };

	std::vector<std::complex<double>> m_FrequencyDomainData;
	std::vector<std::uint64_t>		  m_ExtremaIndicesX;
	std::vector<std::uint64_t>		  m_ExtremaIndicesY;
};

template <std::uint64_t tSize>
mcp::MVoxelData MHoneyCombFilter::arrayToVoxel(
	const std::array<std::array<double, tSize>, tSize>& aData )
{
	const mint dimX = aData[0].size();
	const mint dimY = aData.size();

	mcp::MVoxelData voxelData( dimX, dimY, 1 );

	for ( mint j = 0; j < dimY; ++j )
	{
		for ( mint i = 0; i < dimX; ++i )
		{
			voxelData.at( i, j ) = (float)aData[i][j];
		}
	}

	return voxelData;
}

} // namespace msc