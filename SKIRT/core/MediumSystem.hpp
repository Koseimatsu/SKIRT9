/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MEDIUMSYSTEM_HPP
#define MEDIUMSYSTEM_HPP

#include "SimulationItem.hpp"
#include "Array.hpp"
#include "Medium.hpp"
#include "MaterialMix.hpp"
#include "SpatialGrid.hpp"
class PhotonPacket;
class Random;

//////////////////////////////////////////////////////////////////////

/** An instance of the MediumSystem class represents a complete medium system, which is the
    superposition of one or more transfer media. Each individual medium represents a spatial
    density distribution and defines the material properties of the medium at each location. While
    the specific material properties may vary with location, the fundamental material type must be
    the same throughout the spatial domain for each medium.

    In addition to the media input model, the MediumSystem class includes the spatial grid that
    tessellates the spatial domain of the simulation into cells, and manages the medium state for
    each spatial cell in this grid.

    TODO: add more info on managing the medium state. */
class MediumSystem : public SimulationItem
{
    ITEM_CONCRETE(MediumSystem, SimulationItem, "a medium system")
        ATTRIBUTE_ALLOWED_IF(MediumSystem, "ExtinctionOnlyMode")

    PROPERTY_ITEM_LIST(media, Medium, "the transfer media")
        ATTRIBUTE_DEFAULT_VALUE(media, "GeometricMedium")
        ATTRIBUTE_OPTIONAL(media)

    PROPERTY_ITEM(grid, SpatialGrid, "the spatial grid")
        ATTRIBUTE_DEFAULT_VALUE(grid, "CartesianSpatialGrid")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function calculates and stores initial state information for each cell, including the
        cell volume and the number density for each medium as defined by the input model. */
    void setupSelfAfter() override;

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the medium system, which depends on the (lack of)
        symmetry in the geometries of the media it contains (\em not including the spatial grid). A
        value of 1 means spherical symmetry, 2 means axial symmetry and 3 means none of these
        symmetries. The medium with the least symmetry (i.e. the highest dimension) determines the
        result for the whole system. */
    int dimension() const;

    /** This function returns the dimension of the spatial grid held by the medium system. A value
        of 1 means spherical symmetry, 2 means axial symmetry and 3 means none of these symmetries.
        */
    int gridDimension() const;

    /** This function returns the number of media in the medium system. */
    int numMedia() const;

    /** This function returns the number of cells in the spatial grid held by the medium system. */
    int numCells() const;

    /** This function returns the volume of the spatial cell with index \f$m\f$. */
    double volume(int m) const;

    /** This function returns true if at least one of the media in the medium system has the
        specified fundamental material type (i.e. dust, electrons, or gas). */
    bool hasMaterialType(MaterialMix::MaterialType type) const;

    /** This function returns true if at least one of the media in the medium system contains dust.
        */
    bool hasDust() const { return hasMaterialType(MaterialMix::MaterialType::Dust); }

    /** This function returns true if at least one of the media in the medium system contains
        electrons. */
    bool hasElectrons() const { return hasMaterialType(MaterialMix::MaterialType::Electrons); }

    /** This function returns true if at least one of the media in the medium system contains gas.
        */
    bool hasGas() const { return hasMaterialType(MaterialMix::MaterialType::Gas); }

    /** This function returns true if the medium component with index \f$h\f$ has the specified
        fundamental material type (i.e. dust, electrons, or gas). */
    bool isMaterialType(MaterialMix::MaterialType type, int h) const;

    /** This function returns true if the medium component with index \f$h\f$ contains dust. */
    bool isDust(int h) const { return isMaterialType(MaterialMix::MaterialType::Dust, h); }

    /** This function returns true if the medium component with index \f$h\f$ contains electrons.
        */
    bool isElectrons(int h) const { return isMaterialType(MaterialMix::MaterialType::Electrons, h); }

    /** This function returns true if the medium component with index \f$h\f$ contains gas. */
    bool isGas(int h) const { return isMaterialType(MaterialMix::MaterialType::Gas, h); }

    /** This function returns the number density of the medium component with index \f$h\f$ in
        spatial cell with index \f$m\f$. */
    double numberDensity(int m, int h) const;

    /** This function returns the mass density of the medium component with index \f$h\f$ in
        spatial cell with index \f$m\f$. */
    double massDensity(int m, int h) const;

    /** This function returns the material mix corresponding to the medium component with index
        \f$h\f$ in spatial cell with index \f$m\f$. */
    const MaterialMix* mix(int m, int h) const;

    /** This function randomly returns a material mix corresponding to one of the medium components
        in spatial cell with index \f$m\f$. The sampling is weighted by the scattering opacity
        \f$k=n_h\sigma_h^\text{sca}\f$ at wavelength \f$\lambda\f$ of each medium component with
        index \f$h\f$ in the spatial cell with index \f$m\f$. */
    const MaterialMix* randomMixForScattering(Random* random, double lambda, int m) const;

    /** This function returns the scattering opacity \f$k=n_h\sigma_h^\text{sca}\f$ at wavelength
        \f$\lambda\f$ of the medium component with index \f$h\f$ in spatial cell with index
        \f$m\f$. */
    double opacitySca(double lambda, int m, int h) const;

    /** This function returns the scattering opacity \f$k=\sum_h n_h\sigma_h^\text{sca}\f$ summed
        over all medium components at wavelength \f$\lambda\f$ in spatial cell with index \f$m\f$.
        */
    double opacitySca(double lambda, int m) const;

    /** This function returns the extinction opacity \f$k=n_h\sigma_h^\text{ext}\f$ at wavelength
        \f$\lambda\f$ of the medium component with index \f$h\f$ in spatial cell with index
        \f$m\f$. */
    double opacityExt(double lambda, int m, int h) const;

    /** This function returns the extinction opacity \f$k=\sum_h n_h\sigma_h^\text{ext}\f$ summed
        over all medium components at wavelength \f$\lambda\f$ in spatial cell with index \f$m\f$.
        */
    double opacityExt(double lambda, int m) const;

    /** This function returns the extinction opacity \f$k=\sum_h n_h\sigma_h^\text{ext}\f$ summed
        over all medium components with the specified material type at wavelength \f$\lambda\f$ in
        spatial cell with index \f$m\f$. */
    double opacityExt(double lambda, int m, MaterialMix::MaterialType type) const;

    /** This function returns the scattering albedo \f$\sigma_h^\text{sca}/\sigma_h^\text{ext}\f$
        at wavelength \f$\lambda\f$ of the medium component with index \f$h\f$ in spatial cell with
        index \f$m\f$. */
    double albedo(double lambda, int m, int h) const;

    /** This function returns the weighted scattering albedo \f[\frac{\sum_h
        n_h\sigma_h^\text{sca}} {\sum_h n_h\sigma_h^\text{ext}}\f] over all medium components at
        wavelength \f$\lambda\f$ in spatial cell with index \f$m\f$. */
    double albedo(double lambda, int m) const;

    /** This function calculates the optical depth
        \f$\tau_{\lambda,{\text{path}}}({\boldsymbol{r}},{\boldsymbol{k}})\f$ at wavelength
        \f$\lambda\f$ along a path through the media system starting at the position
        \f${\boldsymbol{r}}\f$ into the direction \f${\boldsymbol{k}}\f$, where \f$\lambda\f$,
        \f${\boldsymbol{r}}\f$ and \f${\boldsymbol{k}}\f$ are obtained from the specified
        PhotonPacket, and it stores the resulting details back into the photon packet.

        The hard work is done by calling the SpatialGrid::path() function which stores the
        geometrical information on the path through the spatial grid into the photon packet: the
        cell indices \f$m\f$ of the cells that are crossed by the path, the path length \f$(\Delta
        s)_m\f$ covered in that particular cell and a total path length counter \f$s_m\f$ that
        gives the total path length covered between the starting point \f${\boldsymbol{r}}\f$ and
        the boundary of the cell.

        With this information given, the calculation of the optical depth is rather
        straightforward. It is calculated as \f[
        \tau_{\lambda,{\text{path}}}({\boldsymbol{r}},{\boldsymbol{k}}) = \sum_m (\Delta s)_m
        \sum_h \varsigma_{\lambda,h}^{\text{ext}}\, n_m, \f] where
        \f$\varsigma_{\lambda,h}^{\text{abs}}\f$ is the extinction cross section corresponding to
        the \f$h\f$'th medium component at wavelength \f$\lambda\f$ and \f$n_{m,h}\f$ the number
        density in the cell with index \f$m\f$ corresponding to the \f$h\f$'th medium component.
        The function also stores the details on the calculation of the optical depth in the photon
        packet, specifically it stores the optical depth covered within the \f$m\f$'th spatial
        cell, \f[ (\Delta\tau_\lambda)_m = (\Delta s)_m \sum_h \varsigma_{\lambda,h}^{\text{ext}}\,
        n_m, \f] and the total optical depth \f$\tau_{\lambda,m}\f$ covered between the starting
        point \f${\boldsymbol{r}}\f$ and the boundary of the cell. */
    void fillOpticalDepth(PhotonPacket* pp);

    //======================== Private Types ========================

private:
    /** This private data structure holds the information maintained per cell and per medium. */
    struct State
    {
        double n;                   // the number density
        const MaterialMix* mix;     // pointer to the material mix
    };

    /** This function returns a writable reference to the state data structure for the given cell
        and medium indices. */
    State& state(int m, int h) { return _Svv[m*_numMedia+h]; }

    /** This function returns a read-only reference to the state data structure for the given cell
        and medium indices. */
    const State& state(int m, int h) const { return _Svv[m*_numMedia+h]; }

    //======================== Data Members ========================

private:
    // initialized during setup
    int _numCells{0};       // index m
    int _numMedia{0};       // index h
    Array _Vv;              // volume of each cell (indexed on m)
    vector<State> _Svv;     // state info for each cell and each medium (indexed on m,h)
};

////////////////////////////////////////////////////////////////

#endif