/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef GKLINEGASMIX_HPP
#define GKLINEGASMIX_HPP

#include "EmittingGasMix.hpp"
#include "GasLineEmission.hpp"

////////////////////////////////////////////////////////////////////

/** The GKLineGasMix class describes rotational line transitions of selected molecules including
    the Goldreich-Kylafis (GK) effect: the linear polarization of molecular lines caused by the
    alignment of the rotational levels in an anisotropic radiation field and a magnetic field.

    <b>Molecular data</b>

    The class supports the molecular species listed in the \em species property. It reads exactly
    the same resource files as the NonLTELineGasMix class (energy levels and statistical weights,
    Einstein A coefficients, collisional coefficients for each interaction partner, and molecular
    mass), using the shared GasLineEmission loader. The rotational quantum number \f$J\f$ of each
    level is derived from its statistical weight, \f$g=s\,(2J+1)\f$, where \f$s\f$ is the nuclear
    spin degeneracy included in the resource data for the species (\f$s=3\f$ for HCN, \f$s=1\f$
    otherwise).

    <b>Polarized statistical equilibrium</b>

    In each spatial cell, the class solves the polarized statistical equilibrium equations
    (pol-SEE) in the strong magnetic field approximation, retaining for each level the irreducible
    density matrix elements \f$\rho^0_0(J)\f$ (population) and \f$\rho^2_0(J)\f$ (alignment, for
    \f$J\ge1\f$) with the quantization axis along the local magnetic field, following Landi
    Degl'Innocenti & Landolfi (2004) and Lankhaar & Vlemmings (2020). The radiative rates depend
    on the line-profile averaged isotropic radiation tensor element \f$\bar{J}^0_0\f$ (the mean
    intensity) and on the anisotropic element \f$[\bar{J}^2_0]_\mathrm{MF}\f$ in the magnetic field
    frame. The latter is obtained by rotating the five global-frame components
    \f$[\bar{J}^2_m]_\mathrm{GF}\f$ that are accumulated during the photon packet life cycle (only
    for simulations that include a GKLineGasMix medium). Collisions are assumed to be isotropic.
    The equations and their implementation are documented in the separate note "Goldreich-Kylafis
    effect in SKIRT". The level populations, level alignments, and the radiation field tensor
    elements for each line are stored as custom state variables, so that they can be output with
    a CustomStateProbe.

    <b>Polarized emission and extinction</b>

    Line photon packets are emitted with the angular distribution and linear polarization of the
    aligned upper level (the actual photon packets are launched isotropically; the anisotropy is
    applied to the emission peel-off photon packets). During emission peel-off towards a distant
    instrument, the simulation applies the dichroic extinction of the aligned molecules along the
    path, propagating the two linear polarization modes parallel and perpendicular to the
    projected magnetic field in each cell. For the calculation of the radiation field, the
    isotropic opacity is used.

    <b>Requirements</b>

    The simulation must have a magnetic field (configured in one of the media) and polarization
    support in all media. As for NonLTELineGasMix, the radiation field wavelength grid must resolve
    the line profiles. The number density of the species is imported or configured through the
    medium; the kinetic temperature, collision partner densities and turbulence velocity are
    imported or taken from the default values configured for this material mix. */
class GKLineGasMix : public EmittingGasMix
{
    /** The enumeration type indicating the molecular species represented by a given GKLineGasMix
        instance. */
    ENUM_DEF(Species, Test, Formyl, HydrogenCyanide, CarbonMonoxide)
        ENUM_VAL(Species, Test, "Fictive two-level test molecule (TT)")
        ENUM_VAL(Species, Formyl, "Formyl cation (HCO+)")
        ENUM_VAL(Species, HydrogenCyanide, "Hydrogen cyanide (HCN)")
        ENUM_VAL(Species, CarbonMonoxide, "Carbon monoxide (CO)")
    ENUM_END()

    ITEM_CONCRETE(GKLineGasMix, EmittingGasMix,
                  "A gas mix supporting polarized rotational transitions (Goldreich-Kylafis effect)")
        ATTRIBUTE_TYPE_INSERT(GKLineGasMix, "CustomMediumState,DynamicState")

        PROPERTY_ENUM(species, Species, "the molecular species being represented")
        ATTRIBUTE_DEFAULT_VALUE(species, "CarbonMonoxide")

        PROPERTY_INT(numEnergyLevels, "the number of energy levels used (or 999 for all supported)")
        ATTRIBUTE_MIN_VALUE(numEnergyLevels, "2")
        ATTRIBUTE_MAX_VALUE(numEnergyLevels, "999")
        ATTRIBUTE_DEFAULT_VALUE(numEnergyLevels, "999")
        ATTRIBUTE_DISPLAYED_IF(numEnergyLevels, "Level2")

        PROPERTY_DOUBLE(defaultTemperature, "the default temperature of the gas")
        ATTRIBUTE_QUANTITY(defaultTemperature, "temperature")
        ATTRIBUTE_MIN_VALUE(defaultTemperature, "[0")
        ATTRIBUTE_MAX_VALUE(defaultTemperature, "1e9]")
        ATTRIBUTE_DEFAULT_VALUE(defaultTemperature, "1000")
        ATTRIBUTE_DISPLAYED_IF(defaultTemperature, "Level2")

        PROPERTY_DOUBLE_LIST(defaultCollisionPartnerRatios,
                             "the default relative abundances of the collisional partners")
        ATTRIBUTE_MIN_VALUE(defaultCollisionPartnerRatios, "[0")
        ATTRIBUTE_MAX_VALUE(defaultCollisionPartnerRatios, "1e50]")
        ATTRIBUTE_DEFAULT_VALUE(defaultCollisionPartnerRatios, "1e4")
        ATTRIBUTE_DISPLAYED_IF(defaultCollisionPartnerRatios, "Level2")

        PROPERTY_DOUBLE(defaultTurbulenceVelocity, "the default (non-thermal) turbulence velocity")
        ATTRIBUTE_QUANTITY(defaultTurbulenceVelocity, "velocity")
        ATTRIBUTE_MIN_VALUE(defaultTurbulenceVelocity, "[0 km/s")
        ATTRIBUTE_MAX_VALUE(defaultTurbulenceVelocity, "100000 km/s]")
        ATTRIBUTE_DEFAULT_VALUE(defaultTurbulenceVelocity, "0 km/s")
        ATTRIBUTE_DISPLAYED_IF(defaultTurbulenceVelocity, "Level2")

        PROPERTY_DOUBLE(maxChangeInLevelPopulations,
                        "the maximum relative change for the level populations in a cell to be considered converged")
        ATTRIBUTE_MIN_VALUE(maxChangeInLevelPopulations, "[0")
        ATTRIBUTE_MAX_VALUE(maxChangeInLevelPopulations, "1]")
        ATTRIBUTE_DEFAULT_VALUE(maxChangeInLevelPopulations, "0.05")
        ATTRIBUTE_DISPLAYED_IF(maxChangeInLevelPopulations, "Level2")

        PROPERTY_DOUBLE(maxFractionNotConvergedCells,
                        "the maximum fraction of not-converged cells for all cells to be considered converged")
        ATTRIBUTE_MIN_VALUE(maxFractionNotConvergedCells, "[0")
        ATTRIBUTE_MAX_VALUE(maxFractionNotConvergedCells, "1]")
        ATTRIBUTE_DEFAULT_VALUE(maxFractionNotConvergedCells, "0.01")
        ATTRIBUTE_DISPLAYED_IF(maxFractionNotConvergedCells, "Level2")

        PROPERTY_DOUBLE(maxChangeInGlobalLevelPopulations,
                        "the maximum relative change for the global level populations to be considered converged")
        ATTRIBUTE_MIN_VALUE(maxChangeInGlobalLevelPopulations, "[0")
        ATTRIBUTE_MAX_VALUE(maxChangeInGlobalLevelPopulations, "1]")
        ATTRIBUTE_DEFAULT_VALUE(maxChangeInGlobalLevelPopulations, "0.05")
        ATTRIBUTE_DISPLAYED_IF(maxChangeInGlobalLevelPopulations, "Level2")

        PROPERTY_DOUBLE(lowestOpticalDepth, "Lower limit of (negative) optical depth along a cell diagonal")
        ATTRIBUTE_MIN_VALUE(lowestOpticalDepth, "[-10")
        ATTRIBUTE_MAX_VALUE(lowestOpticalDepth, "0]")
        ATTRIBUTE_DEFAULT_VALUE(lowestOpticalDepth, "-2")
        ATTRIBUTE_DISPLAYED_IF(lowestOpticalDepth, "Level3")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function loads the molecular data for the configured species from the resource files
        (see NonLTELineGasMix), derives the rotational quantum numbers, precomputes the angular
        momentum coupling coefficients of the pol-SEE for each line, and verifies that the
        simulation has a magnetic field and a radiation field wavelength grid covering the lines. */
    void setupSelfBefore() override;

    //======== Capabilities =======

public:
    /** This function returns true, indicating that this material may have a negative absorption
        cross section (net stimulated emission). */
    bool hasNegativeExtinction() const override;

    /** This function returns true, indicating that the cross sections depend on specific state
        variables other than the number density. */
    bool hasExtraSpecificState() const override;

    /** This function returns DynamicStateType::PrimaryIfMergedIterations (as NonLTELineGasMix). */
    DynamicStateType hasDynamicMediumState() const override;

    /** This function returns true, indicating that this material supports secondary line emission
        from gas. */
    bool hasLineEmission() const override;

    /** This function returns true. The lines do not scatter, but the class supports polarization:
        the emitted line radiation is polarized and the extinction is dichroic. The latter two
        effects are handled through the dedicated functions of this class rather than through the
        hasPolarizedEmission() and hasPolarizedAbsorption() mechanisms for spheroidal grains. */
    bool hasPolarizedScattering() const override;

    /** This function returns true, indicating that this material mix tracks the alignment of its
        energy levels with respect to the local magnetic field. */
    bool hasLevelAlignment() const override;

    //======== Medium state setup =======

public:
    /** This function returns descriptors for the import parameters: the number densities of the
        collisional partners and the turbulence velocity (as for NonLTELineGasMix). */
    vector<SnapshotParameter> parameterInfo() const override;

    /** This function returns the descriptors of the specific state variables. In addition to the
        number density and effective temperature, the custom variables are (in this order): the
        kinetic temperature, the number density of each collision partner, the population of each
        level, the relative alignment \f$\sigma^2_0=\rho^2_0/\rho^0_0\f$ of each level, and for
        each line \f$k\f$ the line-profile averaged mean intensity \f$\bar{J}^0_0\f$, the five
        global-frame components \f$[\bar{J}^2_0]_\mathrm{GF}\f$, \f$\mathrm{Re},\mathrm{Im}\,
        [\bar{J}^2_1]_\mathrm{GF}\f$, \f$\mathrm{Re},\mathrm{Im}\,[\bar{J}^2_2]_\mathrm{GF}\f$,
        and the magnetic-frame element \f$[\bar{J}^2_0]_\mathrm{MF}\f$. */
    vector<StateVariable> specificStateVariableInfo() const override;

    /** This function initializes the specific state variables: the imported or default kinetic
        temperature, collision partner densities and turbulence, LTE level populations, and zero
        alignments and radiation tensor elements. */
    void initializeSpecificState(MaterialState* state, double metallicity, double temperature,
                                 const Array& params) const override;

    //======== Medium state updates =======

    /** This function updates the level populations assuming an isotropic radiation field (zero
        anisotropy). Because hasLevelAlignment() returns true, the medium system calls
        updateSpecificStateWithAnisotropy() instead. */
    UpdateStatus updateSpecificState(MaterialState* state, const Array& Jv) const override;

    /** This function solves the pol-SEE for the given material state, using the mean intensity
        \f$J_\lambda\f$ on the radiation field wavelength grid (\em Jv) and the five global-frame
        rank-2 components on the same grid (\em J2v, laid out as J2v[5*ell+c] with
        c=0..4 for \f$J^2_0\f$, Re \f$J^2_1\f$, Im \f$J^2_1\f$, Re \f$J^2_2\f$, Im \f$J^2_2\f$).
        The magnetic field is taken from the material state. */
    UpdateStatus updateSpecificStateWithAnisotropy(MaterialState* state, const Array& Jv,
                                                   const Array& J2v) const override;

    /** This function returns true if the medium state can be considered converged (same criteria
        as NonLTELineGasMix). */
    bool isSpecificStateConverged(int numCells, int numUpdated, int numNotConverged, MaterialState* currentAggregate,
                                  MaterialState* previousAggregate) const override;

    //======== Low-level material properties =======

public:
    /** This function returns the mass of a molecule of the species. */
    double mass() const override;

    /** This function returns zero (cross sections depend on the level populations). */
    double sectionAbs(double lambda) const override;

    /** This function returns zero. */
    double sectionSca(double lambda) const override;

    /** This function returns zero. */
    double sectionExt(double lambda) const override;

    //======== High-level photon life cycle =======

    /** This function returns the isotropic absorption opacity (as NonLTELineGasMix, i.e. for the
        total populations; the small direction-dependent contribution of the alignment is only
        taken into account during emission peel-off, see polarizedOpacitiesExt()). */
    double opacityAbs(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function returns zero. */
    double opacitySca(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function returns the isotropic extinction opacity, equal to the absorption opacity. */
    double opacityExt(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function does nothing because the lines do not scatter. */
    bool peeloffScattering(double& I, double& Q, double& U, double& V, double& lambda, Direction bfkobs, Direction bfky,
                           const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function does nothing because the lines do not scatter. */
    void performScattering(double lambda, const MaterialState* state, PhotonPacket* pp) const override;

    /** This function returns the opacities for radiation linearly polarized parallel
        (\em kpar) and perpendicular (\em kper) to the projection of the magnetic field on the
        plane perpendicular to the propagation direction, for a propagation direction making an
        angle \f$\vartheta\f$ with the magnetic field (\em cosTheta \f$=\cos\vartheta\f$):
        \f$k_{\parallel,\perp} = k_0 + k_2 f_{\parallel,\perp}\f$ with
        \f$f_\parallel=(3\cos^2\vartheta-2)/\sqrt{2}\f$ and \f$f_\perp=1/\sqrt{2}\f$. If the cell
        has no magnetic field, both equal the isotropic opacity. */
    void polarizedOpacitiesExt(double lambda, const MaterialState* state, const PhotonPacket* pp, double cosTheta,
                               double& kpar, double& kper) const override;

    //======== Secondary emission =======

    /** This function returns the line centers of the supported transitions. */
    Array lineEmissionCenters() const override;

    /** This function returns the particle mass for each line. */
    Array lineEmissionMasses() const override;

    /** This function returns the line luminosities (integrated over all directions) in the cell. */
    Array lineEmissionSpectrum(const MaterialState* state, const Array& Jv) const override;

    /** This function returns the product \f$a=w^{(2)}_{J_u,J_l}\,\sigma^2_0(J_u)\f$ for line
        \em k in the cell, which determines the angular distribution and polarization of the
        spontaneous emission: relative to an isotropic distribution, the emitted intensity is
        \f$1+a(3\cos^2\vartheta-1)/(2\sqrt{2})\f$ and the linear polarization (w.r.t. the projected
        magnetic field) is \f$Q/I = -3a\sin^2\vartheta/\sqrt{2}\,/\,(2+a(3\cos^2\vartheta-1)/\sqrt{2})\f$.
        The function returns zero if the cell has no magnetic field. */
    double lineEmissionAlignmentFactor(const MaterialState* state, int k) const override;

    //======== Temperature =======

    /** This function returns the (effective) temperature stored in the specific state. */
    double indicativeTemperature(const MaterialState* state, const Array& Jv) const override;

    //======== Stand-alone pol-SEE (also used for testing) =======

public:
    /** This function solves the pol-SEE for a molecule described by \em model (as loaded by
        GasLineEmission; energies in J, Einstein B coefficients in the per-wavelength convention),
        with twice the rotational quantum number of each level in \em twoJ, kinetic temperature
        \em Tkin, collision partner densities \em nPartner (m-3), total number density \em nTotal,
        and for each line the line-profile averaged isotropic and magnetic-frame anisotropic
        radiation tensor elements \em J0 and \em J2 (per-wavelength units, W/m3/sr). It returns the
        level populations in \em n (summing to nTotal) and the relative alignments in \em sigma
        (zero for levels with J < 1). The function throws a FatalError if the system is singular. */
    static void solvePolarizedSEE(const GasLineEmission::AtomicModel& model, const vector<int>& twoJ, double Tkin,
                                  const vector<double>& nPartner, double nTotal, const vector<double>& J0,
                                  const vector<double>& J2, vector<double>& n, vector<double>& sigma);

    /** This function returns the coupling factor \f$w^{(2)}_{J,J'} = (-1)^{1+J+J'}\sqrt{3(2J+1)}
        \{1\,1\,2;\,J\,J\,J'\}\f$ for twice the angular momenta \em twoJ and \em twoJp; it vanishes
        for \f$J<1\f$. */
    static double couplingFactorW2(int twoJ, int twoJp);

private:
    /** This function performs the actual state update for the given radiation tensor elements per
        line and returns the average relative change in the level populations. */
    double solveAndStore(MaterialState* state, const vector<double>& J0, const vector<double>& J2GF) const;

    //======================== Data Members ========================

private:
    string _name;  // species resource name
    GasLineEmission _gasLineEmission;
    GasLineEmission::AtomicModel _model;
    vector<int> _twoJ;     // twice the rotational quantum number of each level
    vector<double> _wUL;   // w^(2)_{Ju,Jl} per line (emission, stimulated emission)
    vector<double> _wLU;   // w^(2)_{Jl,Ju} per line (absorption)
    vector<double> _coef;  // pol-SEE coupling coefficients per line (see the implementation file)

    // the radiation field wavelength grid for this simulation
    int _numWavelengths{0};
    Array _lambdav;
    Array _dlambdav;

    // custom variable indices; initialized in specificStateVariableInfo()
    int _indexKineticTemperature{0};
    int _indexFirstColPartnerDensity{0};
    int _indexFirstLevelPopulation{0};
    int _indexFirstLevelAlignment{0};
    int _indexFirstLineTensor{0};  // 7 variables per line: J00, J20GF, J21re, J21im, J22re, J22im, J20MF
};

////////////////////////////////////////////////////////////////////

#endif
