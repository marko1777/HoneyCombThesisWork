#include <ScientificCore/include/MMathematicalUtils.h>
#include <array>
#include <gtest/gtest.h>
#include <iostream>
#include <vector>

TEST( MMathematicalUtils, sub2ind_empty )
{
	std::vector<std::uint64_t> row;
	std::vector<std::uint64_t> col;

	auto result = msc::MMathematicalUtils::sub2ind( { 0, 0 }, row, col );

	EXPECT_TRUE( result.empty() );
}

TEST( MMathematicalUtils, sub2ind_diffSize )
{
	std::vector<std::uint64_t> row( 3 );
	std::vector<std::uint64_t> col( 4 );

	auto result = msc::MMathematicalUtils::sub2ind( { 3, 4 }, row, col );

	EXPECT_TRUE( result.empty() );
}

TEST( MMathematicalUtils, sub2ind_basic )
{
	std::vector<std::uint64_t> row{ 0, 1, 2, 0 };
	std::vector<std::uint64_t> col{ 1, 1, 1, 2 };

	auto actual = msc::MMathematicalUtils::sub2ind( { 3, 3 }, row, col );

	std::vector<std::uint64_t> expected{ 3, 4, 5, 6 };

	EXPECT_EQ( actual, expected );
}

////////

TEST( MMathematicalUtils, ind2sub_empty )
{
	std::vector<std::uint64_t> ind;

	auto result = msc::MMathematicalUtils::ind2sub( { 3, 3 }, ind );

	EXPECT_TRUE( result.first.empty() && result.second.empty() );
}

TEST( MMathematicalUtils, ind2sub_basic )
{
	std::vector<std::uint64_t> ind{ 2, 3, 4, 5 };

	std::vector<std::uint64_t> actualRow;
	std::vector<std::uint64_t> actualCol;

	std::tie( actualRow, actualCol ) = msc::MMathematicalUtils::ind2sub( { 3, 3 }, ind );

	std::vector<std::uint64_t> expectedRow{ 2, 0, 1, 2 };
	std::vector<std::uint64_t> expectedCol{ 0, 1, 1, 1 };

	EXPECT_EQ( actualRow, expectedRow );
	EXPECT_EQ( actualCol, expectedCol );
}

////////

TEST( MMathematicalUtils, fftshift_empty )
{
	std::vector<std::uint64_t> ind;

	auto actual = msc::MMathematicalUtils::fftshift( ind, 3, 3 );

	EXPECT_TRUE( actual.empty() );
}

TEST( MMathematicalUtils, fftshift_basic )
{
	std::vector<std::uint64_t> ind{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };

	auto actual = msc::MMathematicalUtils::fftshift( ind, 3, 3 );

	std::vector<std::uint64_t> expected{ 9, 7, 8, 3, 1, 2, 6, 4, 5 };

	EXPECT_EQ( actual, expected );
}

////////

TEST( MMathematicalUtils, ifftshift_empty )
{
	std::vector<std::uint64_t> ind;

	auto actual = msc::MMathematicalUtils::ifftshift( ind, 3, 3 );

	EXPECT_TRUE( actual.empty() );
}

TEST( MMathematicalUtils, ifftshift_basic )
{
	std::vector<std::uint64_t> ind{ 9, 7, 8, 3, 1, 2, 6, 4, 5 };

	auto actual = msc::MMathematicalUtils::ifftshift( ind, 3, 3 );

	std::vector<std::uint64_t> expected{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };

	EXPECT_EQ( actual, expected );
}

////////

TEST( MMathematicalUtils, linspace_zero )
{
	auto actual = msc::MMathematicalUtils::linspace<double>( -5, 5, 0 );

	std::vector<double> expected{};

	EXPECT_EQ( actual, expected );
}

TEST( MMathematicalUtils, linspace_one )
{
	auto actual = msc::MMathematicalUtils::linspace<double>( -5, 5, 1 );

	std::vector<double> expected{ 5 };

	EXPECT_EQ( actual, expected );
}

TEST( MMathematicalUtils, linspace_basic )
{
	auto actual = msc::MMathematicalUtils::linspace<double>( -5, 5, 7 );

	std::vector<double> expected{ -5.0, -3.33333, -1.66667, 0, 1.66667, 3.33333, 5 };

	EXPECT_TRUE( actual.size() == expected.size() );

	for ( int i = 0; i < actual.size(); ++i )
	{
		EXPECT_NEAR( actual[i], expected[i], 0.00001 );
	}
}

TEST( MMathematicalUtils, linspace_diff_types )
{
	auto actual = msc::MMathematicalUtils::linspace<double>( (float)-5, (std::uint64_t)5, 7 );

	std::vector<double> expected{ -5.0, -3.33333, -1.66667, 0, 1.66667, 3.33333, 5 };

	EXPECT_TRUE( actual.size() == expected.size() );

	for ( int i = 0; i < actual.size(); ++i )
	{
		EXPECT_NEAR( actual[i], expected[i], 0.00001 );
	}
}

////////

TEST( MMathematicalUtils, linerIndexToSubsrcipt )
{
	std::uint64_t x, y;
	std::tie( x, y ) = msc::MMathematicalUtils::linerIndexToSubsrcipt( 6, 3 );

	EXPECT_EQ( x, 0 );
	EXPECT_EQ( y, 2 );
}

////////

TEST( MMathematicalUtils, meshgrid_lhs_empty )
{
	std::array<double, 0> x;
	std::array<double, 2> y{ 1, 2 };

	auto result = msc::MMathematicalUtils::meshgrid( x, y );

	bool empty = true;

	for ( const auto& ar : result.first )
		empty &= ar.empty();

	for ( const auto& ar : result.second )
		empty &= ar.empty();

	EXPECT_TRUE( empty );
}

TEST( MMathematicalUtils, meshgrid_rhs_empty )
{
	std::array<double, 2> x{ 1, 2 };
	std::array<double, 0> y;

	auto result = msc::MMathematicalUtils::meshgrid( x, y );

	bool empty = true;

	for ( const auto& ar : result.first )
		empty &= ar.empty();

	for ( const auto& ar : result.second )
		empty &= ar.empty();
}

TEST( MMathematicalUtils, meshgrid_both_empty )
{
	std::array<double, 0> x;
	std::array<double, 0> y;

	auto result = msc::MMathematicalUtils::meshgrid( x, y );

	bool empty = true;

	for ( const auto& ar : result.first )
		empty &= ar.empty();

	for ( const auto& ar : result.second )
		empty &= ar.empty();

	EXPECT_TRUE( empty );
}

TEST( MMathematicalUtils, meshgrid_same_size )
{
	std::array<double, 3> x{ 0, 1, 2 };
	std::array<double, 3> y{ 1, 2, 3 };

	std::array<std::array<double, 3>, 3> expectedLhs = {
		std::array<double, 3>{ 0, 1, 2 },
		std::array<double, 3>{ 0, 1, 2 },
		std::array<double, 3>{ 0, 1, 2 } };

	std::array<std::array<double, 3>, 3> expectedRhs = {
		std::array<double, 3>{ 1, 1, 1 },
		std::array<double, 3>{ 2, 2, 2 },
		std::array<double, 3>{ 3, 3, 3 } };

	auto actual = msc::MMathematicalUtils::meshgrid( x, y );

	EXPECT_EQ( expectedLhs, actual.first );
	EXPECT_EQ( expectedRhs, actual.second );
}

TEST( MMathematicalUtils, meshgrid_lhs_diff_size )
{
	std::array<double, 4> x{ 0, 1, 2, 3 };
	std::array<double, 3> y{ 1, 2, 3 };

	std::array<std::array<double, 4>, 3> expectedLhs = {
		std::array<double, 4>{ 0, 1, 2, 3 },
		std::array<double, 4>{ 0, 1, 2, 3 },
		std::array<double, 4>{ 0, 1, 2, 3 } };

	std::array<std::array<double, 4>, 3> expectedRhs = {
		std::array<double, 4>{ 1, 1, 1, 1 },
		std::array<double, 4>{ 2, 2, 2, 2 },
		std::array<double, 4>{ 3, 3, 3, 3 } };

	auto actual = msc::MMathematicalUtils::meshgrid( x, y );

	EXPECT_EQ( expectedLhs, actual.first );
	EXPECT_EQ( expectedRhs, actual.second );
}

TEST( MMathematicalUtils, meshgrid_rhs_diff_size )
{
	std::array<double, 3> x{ 0, 1, 2 };
	std::array<double, 4> y{ 1, 2, 3, 4 };

	std::array<std::array<double, 3>, 4> expectedLhs = {
		std::array<double, 3>{ 0, 1, 2 },
		std::array<double, 3>{ 0, 1, 2 },
		std::array<double, 3>{ 0, 1, 2 },
		std::array<double, 3>{ 0, 1, 2 } };

	std::array<std::array<double, 3>, 4> expectedRhs = {
		std::array<double, 3>{ 1, 1, 1 },
		std::array<double, 3>{ 2, 2, 2 },
		std::array<double, 3>{ 3, 3, 3 },
		std::array<double, 3>{ 4, 4, 4 } };

	auto actual = msc::MMathematicalUtils::meshgrid( x, y );

	EXPECT_EQ( expectedLhs, actual.first );
	EXPECT_EQ( expectedRhs, actual.second );
}