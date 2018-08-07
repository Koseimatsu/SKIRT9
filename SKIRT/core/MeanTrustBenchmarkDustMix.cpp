/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MeanTrustBenchmarkDustMix.hpp"

//////////////////////////////////////////////////////////////////////

string MeanTrustBenchmarkDustMix::resourceNameForOpticalProps() const
{
    return "MeanTrustBenchmarkOpticalProps";
}

//////////////////////////////////////////////////////////////////////

string MeanTrustBenchmarkDustMix::resourceNameForMuellerMatrix() const
{
    return includePolarization() ? "MeanTrustBenchmarkMuellerMatrix" : "";
}

//////////////////////////////////////////////////////////////////////
