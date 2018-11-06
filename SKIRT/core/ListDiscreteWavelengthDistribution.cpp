/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ListDiscreteWavelengthDistribution.hpp"
#include "FatalError.hpp"
#include "NR.hpp"

//////////////////////////////////////////////////////////////////////

Array ListDiscreteWavelengthDistribution::getWavelengths() const
{
    return NR::array(_wavelengths);
}

//////////////////////////////////////////////////////////////////////
