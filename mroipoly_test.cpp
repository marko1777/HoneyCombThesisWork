#include <ScientificCore/include/MRoipoly.h>
#include <gtest/gtest.h>
#include <iostream>

TEST( MRoipoly, roipoly_emptyPoints )
{
	const std::uint64_t dimX = 3;
	const std::uint64_t dimY = 3;

	msc::MRoipoly roi( dimX, dimY );

	std::vector<std::uint64_t> x;
	std::vector<std::uint64_t> y;

	auto actual = roi.mask( x, y );

	EXPECT_TRUE( actual.empty() );
}

TEST( MRoipoly, roipoly_zeroSpace )
{
	const std::uint64_t dimX = 0;
	const std::uint64_t dimY = 0;
	msc::MRoipoly		roi( dimX, dimY );

	std::vector<std::uint64_t> x{ 1, 2, 3 };
	std::vector<std::uint64_t> y{ 2, 3, 4 };

	auto actual = roi.mask( x, y );

	EXPECT_TRUE( actual.empty() );
}

TEST( MRoipoly, roipoly_cube )
{
	const std::uint64_t dimX = 6;
	const std::uint64_t dimY = 6;
	msc::MRoipoly		roi( dimX, dimY );

	std::vector<std::uint64_t> x{ 1, 4, 4, 1 };
	std::vector<std::uint64_t> y{ 1, 1, 4, 4 };

	auto actual = roi.mask( x, y );

	std::vector<mcp::MIndex2D<std::uint64_t>> expected;
	expected.reserve( 16 );

	for ( int i = 1; i <= 4; ++i )
	{
		for ( int j = 1; j <= 4; ++j )
		{
			expected.emplace_back( i, j );
		}
	}

	EXPECT_EQ( actual, expected );
}
