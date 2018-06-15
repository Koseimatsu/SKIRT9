/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "PhotonPacket.hpp"
#include "AngularDistributionInterface.hpp"
#include "PolarizationStateInterface.hpp"
#include "RedshiftInterface.hpp"

////////////////////////////////////////////////////////////////////

PhotonPacket::PhotonPacket()
{
}

////////////////////////////////////////////////////////////////////

void PhotonPacket::launch(size_t historyIndex, double lambda, double L, Position bfr, Direction bfk,
                          RedshiftInterface* rsi,
                          AngularDistributionInterface* adi,
                          PolarizationStateInterface* psi)
{
    _lambda = lambda;
    _W = L * lambda;
    _lambda0 = lambda;
    _rsi = rsi;
    _adi = adi;
    _psi = psi;
    _compIndex = 0;
    _historyIndex = historyIndex;
    _nscatt = 0;
    _bfr = bfr;
    _bfk = bfk;
    if (rsi) applyRedshift(rsi->redshiftForDirection(bfk));
    if (psi) setPolarized(psi->polarizationForDirection(bfk));
    else setUnpolarized();
}

////////////////////////////////////////////////////////////////////

void PhotonPacket::setPrimaryOrigin(int sourceCompIndex)
{
    _compIndex = sourceCompIndex + 1;
}

////////////////////////////////////////////////////////////////////

void PhotonPacket::setSecondaryOrigin(int mediumCompIndex)
{
    _compIndex = - (mediumCompIndex + 1);
}

////////////////////////////////////////////////////////////////////

void PhotonPacket::launchEmissionPeelOff(const PhotonPacket* pp, Direction bfk)
{
    _lambda = pp->_lambda;
    _W = pp->_W;
    _compIndex = pp->_compIndex;
    _historyIndex = pp->_historyIndex;
    _nscatt = 0;
    _bfr = pp->_bfr;
    _bfk = bfk;

    if (pp->_rsi)
    {
        _lambda = pp->_lambda0;  // recover the source's rest-frame wavelength
        applyRedshift(pp->_rsi->redshiftForDirection(bfk));
    }
    if (pp->_adi) applyBias(pp->_adi->probabilityForDirection(bfk));
    if (pp->_psi) setPolarized(pp->_psi->polarizationForDirection(bfk));
    else setUnpolarized();
}

////////////////////////////////////////////////////////////////////

void PhotonPacket::launchScatteringPeelOff(const PhotonPacket* pp, Direction bfk, double w)
{
    _lambda = pp->_lambda;
    _W = pp->_W * w;
    _compIndex = pp->_compIndex;
    _historyIndex = pp->_historyIndex;
    _nscatt = pp->_nscatt + 1;
    _bfr = pp->_bfr;
    _bfk = bfk;
    setUnpolarized();
}

////////////////////////////////////////////////////////////////////

void PhotonPacket::propagate(double s)
{
    _bfr += s*_bfk;
}

////////////////////////////////////////////////////////////////////

void PhotonPacket::scatter(Direction bfk)
{
    _nscatt++;
    _bfk = bfk;
}

////////////////////////////////////////////////////////////////////

void PhotonPacket::applyBias(double w)
{
    _W *= w;
}

////////////////////////////////////////////////////////////////////

void PhotonPacket::applyRedshift(double z)
{
    _lambda *= 1+z;
}

////////////////////////////////////////////////////////////////////
