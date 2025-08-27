#include <Hadrons/Modules/MAction/ImprovedStaggered.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MAction;

template class HADRONS_NAMESPACE::MAction::TImprovedStaggered<STAGIMPL>;
#ifdef GRID_DEFAULT_PRECISION_DOUBLE
template class HADRONS_NAMESPACE::MAction::TImprovedStaggered<STAGIMPLF>;
#endif
