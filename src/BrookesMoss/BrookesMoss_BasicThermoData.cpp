/*-----------------------------------------------------------------------*\
|   MOM Library — BrookesMoss<BasicThermoData> explicit instantiation     |
|   Compiled unconditionally. Requires no external dependencies.          |
\*-----------------------------------------------------------------------*/

#if defined(MOM_COMPILED_LIBRARY)
#error "Do not define MOM_COMPILED_LIBRARY when compiling library sources"
#endif

#include "BrookesMoss/BrookesMoss.hpp"

namespace MOM
{
template class BrookesMoss<BasicThermoData>;
}

#if defined(MOM_USE_DICTIONARY)

#include "Dictionary.h"
#include "DictionaryManager.h"
#include "DictionaryGrammar.h"
#include "DictionaryKeyWord.h"

namespace MOM
{
template std::expected<void, std::string> BrookesMoss<BasicThermoData>::SetupFromDictionary<
    OpenSMOKEpp::Dictionary>(OpenSMOKEpp::Dictionary& dict);
}

#endif // MOM_USE_DICTIONARY expected
