/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "AllCellsLibrary.hpp"
#include "MediumSystem.hpp"

////////////////////////////////////////////////////////////////////

AllCellsLibrary::AllCellsLibrary(SimulationItem* parent)
{
    parent->addChild(this);
    setup();
}

////////////////////////////////////////////////////////////////////

int AllCellsLibrary::numEntries() const
{
    return find<MediumSystem>()->numCells();
}

////////////////////////////////////////////////////////////////////

vector<int> AllCellsLibrary::mapping() const
{
    int numCells = numEntries();
    vector<int> nv(numCells);
    for (int m=0; m!=numCells; ++m) nv[m] = m;
    return nv;
}
