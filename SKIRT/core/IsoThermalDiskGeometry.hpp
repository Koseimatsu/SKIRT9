/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ISOTHERMALDISKGEOMETRY_HPP
#define ISOTHERMALDISKGEOMETRY_HPP

#include "SepAxGeometry.hpp"

////////////////////////////////////////////////////////////////////

/** The IsoThermalDiskGeometry class is a subclass of the SepAxGeometry class, and describes axisymmetric
    geometries characterized by a double-exponential profile, in which the density decreases
    exponentially in the radial and the vertical directions. In formula form \f[ \rho(R,z) = \rho_0\,
    {\text{e}}^{-\frac{R}{h_R}-\frac{|z|}{h_z}}, \f]. The model contains two free parameters:
    the scale length \f$h_R\f$ and the vertical scale height \f$h_z\f$. */
class IsoThermalDiskGeometry : public SepAxGeometry
{
    ITEM_CONCRETE(IsoThermalDiskGeometry, SepAxGeometry, "an iso-thermal disk geometry")

        PROPERTY_DOUBLE(scaleLength, "the scale length")
        ATTRIBUTE_QUANTITY(scaleLength, "length")
        ATTRIBUTE_MIN_VALUE(scaleLength, "]0")

        PROPERTY_DOUBLE(scaleHeight, "the scale height")
        ATTRIBUTE_QUANTITY(scaleHeight, "length")
        ATTRIBUTE_MIN_VALUE(scaleHeight, "]0")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function verifies the validity of the parameters. The central density \f$\rho_0\f$ is
        set by the normalization condition that the total mass equals one. The central density
        is expressed as \f[ \rho_0 = \frac{1}{ 4\pi\, h_R^2\, h_z }. \f] */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

    /** This function returns the density \f$\rho(R,z)\f$ at the cylindrical radius \f$R\f$ and
        height \f$z\f$. It just implements the analytical formula. */
    double density(double R, double z) const override;

    /** This function returns the cylindrical radius \f$R\f$ of a random position drawn from the
        geometry, by picking a uniform deviate \f${\cal{X}}\f$ and solving the equation \f[
        {\cal{X}} = 2\pi \int_0^R \rho_R(R')\, R'\, {\text{d}}R' \f] for \f$R\f$. Substituting the
        exponential radial profile (without truncation) into this equation, we obtain \f[ {\cal{X}}
        = 1 - \left( 1+\frac{R}{h_R} \right) \exp \left( -\frac{R}{h_R} \right). \f] This equation
        can be solved by means of the Lambert function of order \f$-1\f$, yielding \f[ R = h_R
        \left[ -1-W_{-1} \left( \frac{ {\cal{X}}-1}{\text{e}} \right) \right]. \f] The Lambert
        function \f$W_{-1}(z)\f$ is implemented in the function SpecialFunctions::LambertW1. */
    double randomCylRadius() const override;

    /** This function returns the height \f$z\f$ of a random position drawn from the geometry, by
        picking a uniform deviate \f${\cal{X}}\f$ and solving the equation \f[ {\cal{X}} =
        \int_{-\infty}^z \rho_z(z')\, {\text{d}}z' \f] for \f$z\f$. For the exponential disk
        geometry, this integration is simple, and the inversion results in \f[ z = \begin{cases} \;
        h_z\,\ln(2{\cal{X}}) & \text{if $0<{\cal{X}}<\tfrac{1}{2}$,} \\ \;-h_z\,\ln[2(1-{\cal{X}})]
        & \text{if $\tfrac{1}{2}<{\cal{X}}<1$.} \end{cases} \f]  */
    double randomZ() const override;

    /** This function returns the surface density along a line in the equatorial plane starting at
        the centre of the coordinate system, i.e. \f[ \Sigma_R = \int_0\infty \rho(R,0)\,
        {\text{d}}R. \f] For the exponential disc geometry we find \f$ \Sigma_R = \rho_0 h_R \f$. */
    double SigmaR() const override;

    /** This function returns the surface density along the Z-axis, i.e. the integration of the
        density along the entire Z-axis, \f[ \Sigma_Z = \int_{-\infty}^\infty \rho(0,0,z)\,
        {\text{d}}z.\f] For the exponential disc geometry we find \f$ \Sigma_Z =
        2\,\rho_0 h_z \f$. */
    double SigmaZ() const override;

    //======================== Data Members ========================

private:
    // aliases to discoverable data members for ease of notation and backwards compatibility
    const double& _hR{_scaleLength};
    const double& _hz{_scaleHeight};

    // data members initialized during setup
    double _rho0{0.};
};

////////////////////////////////////////////////////////////////////

#endif
