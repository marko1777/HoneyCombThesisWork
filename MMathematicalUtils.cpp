#include <ScientificCore/include/MMathematicalUtils.h>

namespace msc
{

std::vector<std::uint64_t> MMathematicalUtils::sub2ind(
	const std::pair<int, int>& aDim, const std::vector<std::uint64_t>& aRow,
	const std::vector<std::uint64_t>& aCol )
{
	if ( aRow.size() != aCol.size() )
	{
		return {};
	}

	std::vector<std::uint64_t> rv;
	rv.reserve( aRow.size() );

	for ( int i = 0; i < aRow.size(); ++i )
	{
		rv.emplace_back( aRow[i] + ( aCol[i] * aDim.first ) );
	}

	return rv;
}

} // namespace msc