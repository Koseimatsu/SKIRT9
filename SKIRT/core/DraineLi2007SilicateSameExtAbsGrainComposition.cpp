/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DraineLi2007SilicateSameExtAbsGrainComposition.hpp"

//////////////////////////////////////////////////////////////////////

DraineLi2007SilicateSameExtAbsGrainComposition::DraineLi2007SilicateSameExtAbsGrainComposition(SimulationItem* parent)
{
    parent->addChild(this);
    setup();
}

//////////////////////////////////////////////////////////////////////

string DraineLi2007SilicateSameExtAbsGrainComposition::name() const
{
    return "Draine_Silicate";
}

//////////////////////////////////////////////////////////////////////

double DraineLi2007SilicateSameExtAbsGrainComposition::bulkDensity() const
{
    return 3e3;
}

//////////////////////////////////////////////////////////////////////

string DraineLi2007SilicateSameExtAbsGrainComposition::resourceNameForOpticalProps() const
{
    return "DustEM_aSil_OpticalProps_SameExtAbs";
}

//////////////////////////////////////////////////////////////////////

string DraineLi2007SilicateSameExtAbsGrainComposition::resourceNameForEnthalpies() const
{
    return "DustEM_aSil_Enthalpies";
}

//////////////////////////////////////////////////////////////////////
