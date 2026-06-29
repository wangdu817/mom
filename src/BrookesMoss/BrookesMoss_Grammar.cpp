#include "BrookesMoss/BrookesMoss_Grammar.h"

namespace MOM
{
void BrookesMoss_Grammar::DefineRules()
{

    AddKeyWord(OpenSMOKEpp::DictionaryKeyWord("@BrookesMoss",
                                              OpenSMOKEpp::SINGLE_BOOL,
                                              "Three equations model: on/off (default: true)",
                                              true));

    // ----------------------------------------------------------------------------------------------------------- //
    // Soot models
    // ----------------------------------------------------------------------------------------------------------- //

    AddKeyWord(OpenSMOKEpp::DictionaryKeyWord(
        "@NucleationModel",
        OpenSMOKEpp::SINGLE_STRING,
        "Nucleation model: 0=none, 1=BrookesMoss, 2=BrookesMossHall (default: 1)",
        false));

    AddKeyWord(OpenSMOKEpp::DictionaryKeyWord("@SurfaceGrowthModel",
                                              OpenSMOKEpp::SINGLE_INT,
                                              "Surface growth model: 0=off, 1=on (default: 1)",
                                              false));

    AddKeyWord(OpenSMOKEpp::DictionaryKeyWord(
        "@OxidationModel",
        OpenSMOKEpp::SINGLE_STRING,
        "Oxidation model: 0=none, 1=BrookesMoss, 2=BrookesMossHall (default: 1)",
        false));

    AddKeyWord(OpenSMOKEpp::DictionaryKeyWord("@CoagulationModel",
                                              OpenSMOKEpp::SINGLE_INT,
                                              "Coagulation model: 0=off, 1=on (default: 1)",
                                              false));

    AddKeyWord(OpenSMOKEpp::DictionaryKeyWord("@ThermophoreticModel",
                                              OpenSMOKEpp::SINGLE_INT,
                                              "Thermophoretic model: 0=off, 1=on (default: 1)",
                                              false));

    // ----------------------------------------------------------------------------------------------------------- //
    // PAHs
    // ----------------------------------------------------------------------------------------------------------- //

    AddKeyWord(OpenSMOKEpp::DictionaryKeyWord(
        "@Precursors",
        OpenSMOKEpp::SINGLE_STRING,
        "Species to be assumed as soot precursors (example: @Precursors C2H2;)",
        true));

    AddKeyWord(OpenSMOKEpp::DictionaryKeyWord(
        "@SurfaceGrowthSpecies",
        OpenSMOKEpp::SINGLE_STRING,
        "Species to be assumed as partecipating surface growth species (example: @SurfaceGrowthSpecies C2H2;)",
        true));

    AddKeyWord(OpenSMOKEpp::DictionaryKeyWord(
        "@SootParticleDiameter", OpenSMOKEpp::SINGLE_MEASURE, "Soot particle diameter", true));

    AddKeyWord(OpenSMOKEpp::DictionaryKeyWord("@SootParticleMolecularWeight",
                                              OpenSMOKEpp::SINGLE_MEASURE,
                                              "Soot particle molecular weight",
                                              true));

    AddKeyWord(OpenSMOKEpp::DictionaryKeyWord(
        "@GasClosureDummySpecies",
        OpenSMOKEpp::SINGLE_STRING,
        "Species to be assumed as gaseous dummy species (example: @GasClosureDummySpecies TIO2RU;)",
        true));

    AddKeyWord(OpenSMOKEpp::DictionaryKeyWord(
        "@GasConsumption",
        OpenSMOKEpp::SINGLE_BOOL,
        "Consumption of gaseous species is accounted for (default: false)",
        true));

    // ----------------------------------------------------------------------------------------------------------- //
    // Additional sub-models
    // ----------------------------------------------------------------------------------------------------------- //

    AddKeyWord(OpenSMOKEpp::DictionaryKeyWord("@RadiativeHeatTransfer",
                                              OpenSMOKEpp::SINGLE_BOOL,
                                              "Radiative heat transfer (default: true)",
                                              false));

    AddKeyWord(OpenSMOKEpp::DictionaryKeyWord(
        "@PlanckCoefficient",
        OpenSMOKEpp::SINGLE_STRING,
        "Calculation of Planck Coefficient: Smooke (default) | Kent | Sazhin | none",
        false));

    AddKeyWord(OpenSMOKEpp::DictionaryKeyWord(
        "@SchmidtNumber", OpenSMOKEpp::SINGLE_DOUBLE, "Schmidt number (default: 50)", false));

    // ----------------------------------------------------------------------------------------------------------- //
    // Soot properties
    // ----------------------------------------------------------------------------------------------------------- //

    AddKeyWord(OpenSMOKEpp::DictionaryKeyWord("@SootDensity",
                                              OpenSMOKEpp::SINGLE_MEASURE,
                                              "Density of soot particles (default: 1800 kg/m3)",
                                              false));

    // ----------------------------------------------------------------------------------------------------------- //
    // Model constants
    // ----------------------------------------------------------------------------------------------------------- //

    AddKeyWord(OpenSMOKEpp::DictionaryKeyWord("@Calpha",
                                              OpenSMOKEpp::SINGLE_MEASURE,
                                              "Model constant for soot incipent rate (default: 54 1/s)",
                                              false));

    AddKeyWord(
        OpenSMOKEpp::DictionaryKeyWord("@Talpha",
                                       OpenSMOKEpp::SINGLE_MEASURE,
                                       "Activation temperature of soot inception (default: 21000 K)",
                                       false));

    AddKeyWord(OpenSMOKEpp::DictionaryKeyWord("@Cbeta",
                                              OpenSMOKEpp::SINGLE_DOUBLE,
                                              "Model constant for coagulation rate (default: 1)",
                                              false));

    AddKeyWord(
        OpenSMOKEpp::DictionaryKeyWord("@Cgamma",
                                       OpenSMOKEpp::SINGLE_MEASURE,
                                       "Surface growth rate scaling factor (default: 11700 kgm/kmol/s)",
                                       false));

    AddKeyWord(
        OpenSMOKEpp::DictionaryKeyWord("@Tgamma",
                                       OpenSMOKEpp::SINGLE_MEASURE,
                                       "Activation temperature of surface growth (default: 12100 K)",
                                       false));

    AddKeyWord(
        OpenSMOKEpp::DictionaryKeyWord("@Comega",
                                       OpenSMOKEpp::SINGLE_MEASURE,
                                       "Oxidation model constant (default: 105.8125 kgm/kmol/sqrt(K)/s)",
                                       false));

    AddKeyWord(OpenSMOKEpp::DictionaryKeyWord("@EtaColl",
                                              OpenSMOKEpp::SINGLE_DOUBLE,
                                              "Collisional efficiency parameter (default: 0.04)",
                                              false));

    AddKeyWord(OpenSMOKEpp::DictionaryKeyWord("@Coxid",
                                              OpenSMOKEpp::SINGLE_DOUBLE,
                                              "Oxidation rate scaling parameter (default: 0.015)",
                                              false));

    AddKeyWord(OpenSMOKEpp::DictionaryKeyWord("@NucleationExponent",
                                              OpenSMOKEpp::SINGLE_DOUBLE,
                                              "Exponent of concentration in nucleation (default: 1)",
                                              false));

    AddKeyWord(
        OpenSMOKEpp::DictionaryKeyWord("@SurfaceGrowthExponent1",
                                       OpenSMOKEpp::SINGLE_DOUBLE,
                                       "Exponent of concentration in surface growth (default: 1)",
                                       false));

    AddKeyWord(
        OpenSMOKEpp::DictionaryKeyWord("@SurfaceGrowthExponent2",
                                       OpenSMOKEpp::SINGLE_DOUBLE,
                                       "Exponent of mass concentration in surface growth (default: 1)",
                                       false));

    // ----------------------------------------------------------------------------------------------------------- //
    // Minimum values for properties calculation
    // ----------------------------------------------------------------------------------------------------------- //

    AddKeyWord(OpenSMOKEpp::DictionaryKeyWord(
        "@NsNorm", OpenSMOKEpp::SINGLE_MEASURE, "Normalization factor (default: 1e15 #/m3)", false));
}
} // namespace MOM
