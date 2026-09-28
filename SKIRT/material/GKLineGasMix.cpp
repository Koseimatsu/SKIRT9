/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "GKLineGasMix.hpp"
#include "Configuration.hpp"
#include "Constants.hpp"
#include "DisjointWavelengthGrid.hpp"
#include "FatalError.hpp"
#include "Log.hpp"
#include "MaterialState.hpp"
#include "NR.hpp"
#include "StringUtils.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // ============== Wigner symbols ==============

    // All angular momentum arguments are given as TWICE their value, so that half-integer values
    // can be represented. The symbols are evaluated with the Racah formulas using logarithms of
    // factorials; for the symbols needed here (with at least one small argument) the sums contain
    // only a few terms, so that the relative accuracy is close to machine precision.

    double logFactorial(int n)
    {
        return std::lgamma(n + 1.);
    }

    int phase(int n)
    {
        return (n % 2 == 0) ? 1 : -1;
    }  // (-1)^n for integer n (also negative)

    bool triangle(int a, int b, int c)
    {
        return (a + b + c) % 2 == 0 && c <= a + b && c >= std::abs(a - b);
    }

    double logDelta(int a, int b, int c)
    {
        return 0.5
               * (logFactorial((a + b - c) / 2) + logFactorial((a - b + c) / 2) + logFactorial((-a + b + c) / 2)
                  - logFactorial((a + b + c) / 2 + 1));
    }

    double wigner3j(int j1, int j2, int j3, int m1, int m2, int m3)
    {
        if (m1 + m2 + m3 != 0 || !triangle(j1, j2, j3)) return 0.;
        if (std::abs(m1) > j1 || std::abs(m2) > j2 || std::abs(m3) > j3) return 0.;
        if ((j1 + m1) % 2 || (j2 + m2) % 2 || (j3 + m3) % 2) return 0.;
        int kmin = std::max(0, std::max((j2 - j3 - m1) / 2, (j1 - j3 + m2) / 2));
        int kmax = std::min((j1 + j2 - j3) / 2, std::min((j1 - m1) / 2, (j2 + m2) / 2));
        double pre =
            logDelta(j1, j2, j3)
            + 0.5
                  * (logFactorial((j1 + m1) / 2) + logFactorial((j1 - m1) / 2) + logFactorial((j2 + m2) / 2)
                     + logFactorial((j2 - m2) / 2) + logFactorial((j3 + m3) / 2) + logFactorial((j3 - m3) / 2));
        double sum = 0.;
        for (int k = kmin; k <= kmax; ++k)
        {
            double den = logFactorial(k) + logFactorial((j3 - j2 + m1) / 2 + k) + logFactorial((j3 - j1 - m2) / 2 + k)
                         + logFactorial((j1 + j2 - j3) / 2 - k) + logFactorial((j1 - m1) / 2 - k)
                         + logFactorial((j2 + m2) / 2 - k);
            sum += phase(k) * std::exp(pre - den);
        }
        return phase((j1 - j2 - m3) / 2) * sum;
    }

    double wigner6j(int a, int b, int c, int d, int e, int f)
    {
        if (!triangle(a, b, c) || !triangle(a, e, f) || !triangle(d, b, f) || !triangle(d, e, c)) return 0.;
        int t1 = (a + b + c) / 2, t2 = (a + e + f) / 2, t3 = (d + b + f) / 2, t4 = (d + e + c) / 2;
        int p1 = (a + b + d + e) / 2, p2 = (b + c + e + f) / 2, p3 = (c + a + f + d) / 2;
        int tmin = std::max(std::max(t1, t2), std::max(t3, t4));
        int tmax = std::min(p1, std::min(p2, p3));
        double pre = logDelta(a, b, c) + logDelta(a, e, f) + logDelta(d, b, f) + logDelta(d, e, c);
        double sum = 0.;
        for (int t = tmin; t <= tmax; ++t)
        {
            double den = logFactorial(t - t1) + logFactorial(t - t2) + logFactorial(t - t3) + logFactorial(t - t4)
                         + logFactorial(p1 - t) + logFactorial(p2 - t) + logFactorial(p3 - t);
            sum += phase(t) * std::exp(pre + logFactorial(t + 1) - den);
        }
        return sum;
    }

    // {a b c; d e f; g h i} = sum_x (-1)^{2x} (2x+1) {a b c; f i x} {d e f; b x h} {g h i; x a d}
    double wigner9j(int a, int b, int c, int d, int e, int f, int g, int h, int i)
    {
        int xmin = std::max(std::abs(a - i), std::max(std::abs(d - h), std::abs(b - f)));
        int xmax = std::min(a + i, std::min(d + h, b + f));
        double sum = 0.;
        for (int x = xmin; x <= xmax; x += 2)
            sum += phase(x) * (x + 1) * wigner6j(a, b, c, f, i, x) * wigner6j(d, e, f, b, x, h)
                   * wigner6j(g, h, i, x, a, d);
        return sum;
    }

    // ============== pol-SEE coupling coefficients ==============

    // For each line, the coupling coefficients are stored in a flat array with NUMCOEF values:
    // four blocks of 2x2x2 values indexed by (rank index of the row level, rank index of the
    // column level, rank index of the radiation tensor), with rank index 0 for K=0 and 1 for K=2,
    //   block 0: rp for absorption into the upper level from the lower level
    //   block 1: rm for stimulated emission out of the upper level
    //   block 2: rp for stimulated emission into the lower level from the upper level
    //   block 3: rm for absorption out of the lower level
    // followed by two values tp (rank index 0, 1) for spontaneous emission into the lower level.
    // The rp and rm values still need to be multiplied by [J_u] B_ul J^K, and tp by A_ul.
    // These are the coefficients of SKIRT_PORTAL/SEE_POL (Landi Degl'Innocenti & Landolfi 2004, Eq. 7.20).
    constexpr int NUMCOEF = 34;

    int coefIndex(int block, int ik1, int ik2, int iKr)
    {
        return block * 8 + ik1 * 4 + ik2 * 2 + iKr;
    }

    double fk(int k1, int k2, int K)
    {
        return std::sqrt(3. * (2 * k1 + 1) * (2 * k2 + 1) * (2 * K + 1));
    }

    double couplingRp(int j1, int j2, int k1, int k2, int K)
    {
        return fk(k1, k2, K) * wigner9j(2, j1, j2, 2, j1, j2, 2 * K, 2 * k1, 2 * k2)
               * wigner3j(2 * k1, 2 * k2, 2 * K, 0, 0, 0);
    }

    double couplingRm(int j1, int j2, int k1, int k2, int K)
    {
        return fk(k1, k2, K) * phase(1 + (j1 - j2) / 2) * wigner6j(2, 2, 2 * K, j1, j1, j2)
               * wigner6j(2 * k1, 2 * k2, 2 * K, j1, j1, j1) * wigner3j(2 * k1, 2 * k2, 2 * K, 0, 0, 0);
    }

    double couplingTp(int j1, int j2, int k1)
    {
        return (j2 + 1.) * phase(1 + (j1 + j2) / 2) * wigner6j(j1, j1, 2 * k1, j2, j2, 2);
    }

    vector<double> buildCoefficients(const GasLineEmission::AtomicModel& model, const vector<int>& twoJ)
    {
        int numLines = model.numLines();
        vector<double> coef(numLines * NUMCOEF, 0.);
        for (int t = 0; t != numLines; ++t)
        {
            int ju = twoJ[model.indexUpRad[t]];
            int jl = twoJ[model.indexLowRad[t]];
            double* c = &coef[t * NUMCOEF];
            for (int ik1 = 0; ik1 != 2; ++ik1)
                for (int ik2 = 0; ik2 != 2; ++ik2)
                    for (int iK = 0; iK != 2; ++iK)
                    {
                        int k1 = 2 * ik1, k2 = 2 * ik2, K = 2 * iK;
                        c[coefIndex(0, ik1, ik2, iK)] = couplingRp(ju, jl, k1, k2, K);
                        c[coefIndex(1, ik1, ik2, iK)] = couplingRm(ju, jl, k1, k2, K);
                        c[coefIndex(2, ik1, ik2, iK)] = couplingRp(jl, ju, k1, k2, K);
                        c[coefIndex(3, ik1, ik2, iK)] = couplingRm(jl, ju, k1, k2, K);
                    }
            c[32] = couplingTp(jl, ju, 0);
            c[33] = couplingTp(jl, ju, 2);
        }
        return coef;
    }

    // clamped log-log interpolation (as used by GasLineEmission for the collisional coefficients)
    double clampedLogLog(double x, const vector<double>& xv, const vector<double>& yv)
    {
        int n = static_cast<int>(xv.size());
        if (x < xv[0]) return yv[0];
        if (x >= xv[n - 1]) return yv[n - 1];
        int i = static_cast<int>(std::upper_bound(xv.begin(), xv.end(), x) - xv.begin()) - 1;
        double x1 = xv[i], x2 = xv[i + 1], f1 = yv[i], f2 = yv[i + 1];
        if (f1 <= 0 || f2 <= 0)
        {
            if (x == x1) return f1;
            if (x == x2) return f2;
            return 0.;
        }
        return f1 * exp(log(x / x1) / log(x2 / x1) * (log(f2 / f1)));
    }

    // solves the dense linear system A x = b in place with LU decomposition and partial pivoting;
    // returns false if the matrix is singular
    bool solveLinearSystem(vector<double>& A, vector<double>& b, int n)
    {
        for (int col = 0; col != n; ++col)
        {
            int pivot = col;
            double maxval = std::abs(A[col * n + col]);
            for (int row = col + 1; row < n; ++row)
            {
                double val = std::abs(A[row * n + col]);
                if (val > maxval)
                {
                    maxval = val;
                    pivot = row;
                }
            }
            if (maxval == 0.) return false;
            if (pivot != col)
            {
                for (int k = 0; k != n; ++k) std::swap(A[col * n + k], A[pivot * n + k]);
                std::swap(b[col], b[pivot]);
            }
            double diag = A[col * n + col];
            for (int row = col + 1; row < n; ++row)
            {
                double factor = A[row * n + col] / diag;
                if (factor == 0.) continue;
                for (int k = col; k != n; ++k) A[row * n + k] -= factor * A[col * n + k];
                b[row] -= factor * b[col];
            }
        }
        for (int row = n - 1; row >= 0; --row)
        {
            double sum = b[row];
            for (int k = row + 1; k != n; ++k) sum -= A[row * n + k] * b[k];
            b[row] = sum / A[row * n + row];
        }
        return true;
    }

    // solves the pol-SEE; see GKLineGasMix::solvePolarizedSEE()
    void solvePolSEE(const GasLineEmission::AtomicModel& model, const vector<int>& twoJ, const vector<double>& coef,
                     double Tkin, const vector<double>& nPartner, double nTotal, const vector<double>& J0,
                     const vector<double>& J2, vector<double>& n, vector<double>& sigma)
    {
        // bookkeeping of the unknowns: rho_J0 for each level, and rho_J2 for levels with J >= 1
        int numLevels = model.numLevels();
        vector<int> a0(numLevels), a2(numLevels, -1);
        int numAlign = 0;
        for (int p = 0; p != numLevels; ++p)
        {
            a0[p] = numAlign++;
            if (twoJ[p] >= 2) a2[p] = numAlign++;
        }
        auto index = [&a0, &a2](int p, int ik) { return ik == 0 ? a0[p] : a2[p]; };
        auto numRanks = [&twoJ](int p) { return twoJ[p] >= 2 ? 2 : 1; };

        vector<double> M(numAlign * numAlign, 0.);
        auto add = [&M, numAlign](int row, int col, double value) { M[row * numAlign + col] += value; };

        // radiative transitions
        for (int t = 0; t != model.numLines(); ++t)
        {
            int up = model.indexUpRad[t];
            int low = model.indexLowRad[t];
            double A = model.einsteinA[t];
            double Bt = model.einsteinBul[t] * (twoJ[up] + 1);
            double BJ[2] = {Bt * J0[t], Bt * J2[t]};
            const double* c = &coef[t * NUMCOEF];
            int nu = numRanks(up), nl = numRanks(low);

            // equations for the upper level
            for (int ik1 = 0; ik1 != nu; ++ik1)
            {
                int row = index(up, ik1);
                add(row, row, -A);
                for (int iK = 0; iK != 2; ++iK)
                {
                    for (int ik2 = 0; ik2 != nu; ++ik2)
                        add(row, index(up, ik2), -c[coefIndex(1, ik1, ik2, iK)] * BJ[iK]);
                    for (int ik2 = 0; ik2 != nl; ++ik2)
                        add(row, index(low, ik2), c[coefIndex(0, ik1, ik2, iK)] * BJ[iK]);
                }
            }

            // equations for the lower level
            for (int ik1 = 0; ik1 != nl; ++ik1)
            {
                int row = index(low, ik1);
                for (int iK = 0; iK != 2; ++iK)
                {
                    for (int ik2 = 0; ik2 != nl; ++ik2)
                        add(row, index(low, ik2), -c[coefIndex(3, ik1, ik2, iK)] * BJ[iK]);
                    for (int ik2 = 0; ik2 != nu; ++ik2)
                        add(row, index(up, ik2), c[coefIndex(2, ik1, ik2, iK)] * BJ[iK]);
                }
                // spontaneous emission conserves the rank
                if (ik1 < nu) add(row, index(up, ik1), c[32 + ik1] * A);
            }
        }

        // collisional transitions (isotropic collisions: only rank 0 is transferred, all ranks decay);
        // the upward rates follow from detailed balance without flooring
        for (size_t cp = 0; cp != model.colPartner.size(); ++cp)
        {
            const auto& partner = model.colPartner[cp];
            double Trep = std::max(partner.T.front(), std::min(Tkin, partner.T.back()));
            double np = nPartner[cp];
            for (size_t t = 0; t != partner.indexUpCol.size(); ++t)
            {
                int up = partner.indexUpCol[t];
                int low = partner.indexLowCol[t];
                double Kul = clampedLogLog(Trep, partner.T, partner.Kul[t]);
                double Klu = Kul * model.weight[up] / model.weight[low]
                             * exp(-(model.energy[up] - model.energy[low]) / Constants::k() / Tkin);
                double down = Kul * np;
                double upr = Klu * np;
                double g = std::sqrt((twoJ[low] + 1.) / (twoJ[up] + 1.));
                add(a0[up], a0[low], g * upr);
                add(a0[low], a0[up], down / g);
                for (int ik = 0; ik != numRanks(up); ++ik) add(index(up, ik), index(up, ik), -down);
                for (int ik = 0; ik != numRanks(low); ++ik) add(index(low, ik), index(low, ik), -upr);
            }
        }

        // replace the first equation by the normalization sum_J sqrt(2J+1) rho_J0 = nTotal
        vector<double> rhs(numAlign, 0.);
        for (int col = 0; col != numAlign; ++col) M[col] = 0.;
        for (int p = 0; p != numLevels; ++p) M[a0[p]] = std::sqrt(twoJ[p] + 1.);
        rhs[0] = nTotal;

        if (!solveLinearSystem(M, rhs, numAlign)) throw FATALERROR("The polarized statistical equilibrium is singular");

        n.resize(numLevels);
        sigma.resize(numLevels);
        for (int p = 0; p != numLevels; ++p)
        {
            double rho0 = rhs[a0[p]];
            n[p] = std::sqrt(twoJ[p] + 1.) * rho0;
            sigma[p] = (a2[p] >= 0 && rho0 != 0.) ? rhs[a2[p]] / rho0 : 0.;
            if (!std::isfinite(n[p]) || !std::isfinite(sigma[p]))
                throw FATALERROR("The polarized statistical equilibrium has a non-finite solution");
        }
    }

    // ============== Line profile ==============

    // dispersion of a line profile in wavelength space
    double sigmaForLine(double center, double temperature, double mass)
    {
        return center / Constants::c() * sqrt(Constants::k() * temperature / mass);
    }

    // normalized Gaussian with center mu and dispersion sigma
    double gaussian(double x, double mu, double sigma)
    {
        double u = (x - mu) / sigma;
        constexpr double front = 0.25 * M_SQRT2 * M_2_SQRTPI;
        return front / sigma * exp(-0.5 * u * u);
    }

    // line profile range considered in the calculations, in units of the Gaussian dispersion
    constexpr double PROFILE_RANGE = 4.;

    // allowed fractional errors on the integration of the line profile over the radiation field grid
    constexpr double MAX_GAUSS_ERROR_WARN = 0.01;
    constexpr double MAX_GAUSS_ERROR_FAIL = 0.10;

    // number of radiation tensor variables stored per line
    constexpr int NUMLINETENSOR = 7;
}

////////////////////////////////////////////////////////////////////

void GKLineGasMix::solvePolarizedSEE(const GasLineEmission::AtomicModel& model, const vector<int>& twoJ, double Tkin,
                                     const vector<double>& nPartner, double nTotal, const vector<double>& J0,
                                     const vector<double>& J2, vector<double>& n, vector<double>& sigma)
{
    solvePolSEE(model, twoJ, buildCoefficients(model, twoJ), Tkin, nPartner, nTotal, J0, J2, n, sigma);
}

////////////////////////////////////////////////////////////////////

double GKLineGasMix::couplingFactorW2(int twoJ, int twoJp)
{
    if (twoJ < 2) return 0.;
    return phase(1 + (twoJ + twoJp) / 2) * std::sqrt(3. * (twoJ + 1)) * wigner6j(2, 2, 4, twoJ, twoJ, twoJp);
}

////////////////////////////////////////////////////////////////////

void GKLineGasMix::setupSelfBefore()
{
    EmittingGasMix::setupSelfBefore();

    // get the resource name of the configured species, the names of its collision partners,
    // and the nuclear spin degeneracy included in the statistical weights
    vector<string> colNames{"H2"};
    int spin = 1;
    switch (species())
    {
        case Species::Test: _name = "TT"; break;
        case Species::Formyl: _name = "HCO+"; break;
        case Species::HydrogenCyanide:
            _name = "HCN";
            colNames = {"H2", "e-"};
            spin = 3;
            break;
        case Species::CarbonMonoxide: _name = "CO"; break;
    }

    // load the molecular model with the loader shared with NonLTELineGasMix
    _gasLineEmission.initialize(this);
    _gasLineEmission.loadAtomicModel(_name, colNames, numEnergyLevels(), _model);
    auto log = find<Log>();
    log->info("Loaded molecular model for " + _name + " from " + std::to_string(3 + 2 * colNames.size())
              + " resource files");
    int numLines = _model.numLines();
    if (numLines == 0) throw FATALERROR("There are no radiative transitions; increase the number of energy levels");

    // derive the rotational quantum number of each level from its statistical weight g = spin*(2J+1)
    int numLevels = _model.numLevels();
    _twoJ.resize(numLevels);
    for (int p = 0; p != numLevels; ++p)
    {
        double twoJplus1 = _model.weight[p] / spin;
        int rounded = static_cast<int>(std::lround(twoJplus1));
        if (rounded < 1 || std::abs(twoJplus1 - rounded) > 1e-6)
            throw FATALERROR("Statistical weight of level " + std::to_string(p) + " of " + _name
                             + " does not correspond to a rotational level");
        _twoJ[p] = rounded - 1;
    }

    // precompute the coupling coefficients
    _wUL.resize(numLines);
    _wLU.resize(numLines);
    for (int k = 0; k != numLines; ++k)
    {
        int ju = _twoJ[_model.indexUpRad[k]];
        int jl = _twoJ[_model.indexLowRad[k]];
        _wUL[k] = couplingFactorW2(ju, jl);
        _wLU[k] = couplingFactorW2(jl, ju);
    }
    _coef = buildCoefficients(_model, _twoJ);

    // log the lines
    auto units = find<Units>();
    log->info("Radiative lines for " + _name + " (Goldreich-Kylafis):");
    for (int k = 0; k != numLines; ++k)
        log->info("  (" + std::to_string(_model.indexUpRad[k]) + "-" + std::to_string(_model.indexLowRad[k]) + ") "
                  + StringUtils::toString(units->owavelength(_model.center[k])) + " " + units->uwavelength()
                  + ", J=" + StringUtils::toString(0.5 * _twoJ[_model.indexUpRad[k]]) + "-"
                  + StringUtils::toString(0.5 * _twoJ[_model.indexLowRad[k]]));
    log->info("Collisional partner(s) for " + _name + ": " + StringUtils::join(colNames, ", "));

    // the alignment is defined with respect to the local magnetic field
    auto config = find<Configuration>();
    if (!config->hasMagneticField())
        throw FATALERROR("GKLineGasMix requires a magnetic field (configure one in one of the media)");

    // verify that the radiation field wavelength grid covers the line centers and cache it
    auto rfwlg = config->radiationFieldWLG();
    if (rfwlg)
    {
        rfwlg->setup();
        for (int k = 0; k != numLines; ++k)
            if (rfwlg->bin(_model.center[k]) < 0)
                throw FATALERROR("Radiation field wavelength grid does not cover the central line for transition ("
                                 + std::to_string(_model.indexUpRad[k]) + "-" + std::to_string(_model.indexLowRad[k])
                                 + ")");
        _numWavelengths = rfwlg->numBins();
        _lambdav = rfwlg->lambdav();
        _dlambdav = rfwlg->dlambdav();
    }
}

////////////////////////////////////////////////////////////////////

bool GKLineGasMix::hasNegativeExtinction() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

bool GKLineGasMix::hasExtraSpecificState() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

MaterialMix::DynamicStateType GKLineGasMix::hasDynamicMediumState() const
{
    return DynamicStateType::PrimaryIfMergedIterations;
}

////////////////////////////////////////////////////////////////////

bool GKLineGasMix::hasLineEmission() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

bool GKLineGasMix::hasPolarizedScattering() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

bool GKLineGasMix::hasLevelAlignment() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

vector<SnapshotParameter> GKLineGasMix::parameterInfo() const
{
    vector<SnapshotParameter> result;
    for (const auto& partner : _model.colPartner)
        result.push_back(SnapshotParameter::custom(partner.name + " number density", "numbervolumedensity", "1/cm3"));
    result.push_back(SnapshotParameter::custom("turbulence velocity", "velocity", "km/s"));
    return result;
}

////////////////////////////////////////////////////////////////////

vector<StateVariable> GKLineGasMix::specificStateVariableInfo() const
{
    vector<StateVariable> result{StateVariable::numberDensity(), StateVariable::temperature()};
    auto self = const_cast<GKLineGasMix*>(this);
    int index = 0;

    self->_indexKineticTemperature = index;
    result.push_back(StateVariable::custom(index++, "kinetic gas temperature", "temperature"));

    self->_indexFirstColPartnerDensity = index;
    for (const auto& partner : _model.colPartner)
        result.push_back(StateVariable::custom(index++, partner.name + " number density", "numbervolumedensity"));

    self->_indexFirstLevelPopulation = index;
    for (int p = 0; p != _model.numLevels(); ++p)
        result.push_back(
            StateVariable::custom(index++, "population of level " + std::to_string(p), "numbervolumedensity"));

    self->_indexFirstLevelAlignment = index;
    for (int p = 0; p != _model.numLevels(); ++p)
        result.push_back(StateVariable::custom(index++, "alignment of level " + std::to_string(p), ""));

    self->_indexFirstLineTensor = index;
    for (int k = 0; k != _model.numLines(); ++k)
    {
        string line = " at line " + std::to_string(k);
        result.push_back(StateVariable::custom(index++, "mean intensity" + line, "wavelengthmeanintensity"));
        result.push_back(
            StateVariable::custom(index++, "anisotropic mean intensity J20" + line, "wavelengthmeanintensity"));
        result.push_back(StateVariable::custom(index++, "anisotropic mean intensity J21 real part" + line,
                                               "wavelengthmeanintensity"));
        result.push_back(StateVariable::custom(index++, "anisotropic mean intensity J21 imaginary part" + line,
                                               "wavelengthmeanintensity"));
        result.push_back(StateVariable::custom(index++, "anisotropic mean intensity J22 real part" + line,
                                               "wavelengthmeanintensity"));
        result.push_back(StateVariable::custom(index++, "anisotropic mean intensity J22 imaginary part" + line,
                                               "wavelengthmeanintensity"));
        result.push_back(StateVariable::custom(index++, "anisotropic mean intensity J20 in magnetic frame" + line,
                                               "wavelengthmeanintensity"));
    }
    return result;
}

////////////////////////////////////////////////////////////////////

// Macro's for accessing custom variables in the material state (as in NonLTELineGasMix)
#define setKineticTemperature(value) setCustom(_indexKineticTemperature, (value))
#define kineticTemperature() custom(_indexKineticTemperature)
#define setColPartnerDensity(index, value) setCustom(_indexFirstColPartnerDensity + (index), (value))
#define colPartnerDensity(index) custom(_indexFirstColPartnerDensity + (index))
#define setLevelPopulation(index, value) setCustom(_indexFirstLevelPopulation + (index), (value))
#define levelPopulation(index) custom(_indexFirstLevelPopulation + (index))
#define setLevelAlignment(index, value) setCustom(_indexFirstLevelAlignment + (index), (value))
#define levelAlignment(index) custom(_indexFirstLevelAlignment + (index))
#define setLineTensor(line, c, value) setCustom(_indexFirstLineTensor + NUMLINETENSOR * (line) + (c), (value))

////////////////////////////////////////////////////////////////////

void GKLineGasMix::initializeSpecificState(MaterialState* state, double /*metallicity*/, double temperature,
                                           const Array& params) const
{
    if (state->numberDensity() > 0.)
    {
        int numColPartners = _model.numColPartners();
        int numLevels = _model.numLevels();

        // kinetic temperature from import or default, effective temperature including turbulence
        double Tkin = temperature >= 0. ? temperature : defaultTemperature();
        state->setKineticTemperature(Tkin);
        double vturb = params.size() ? params[numColPartners] : defaultTurbulenceVelocity();
        state->setTemperature(Tkin + 0.5 * vturb * vturb * _model.mass / Constants::k());

        // collision partner densities from import or default
        if (params.size())
        {
            for (int c = 0; c != numColPartners; ++c) state->setColPartnerDensity(c, params[c]);
        }
        else
        {
            const auto& ratios = defaultCollisionPartnerRatios();
            if (static_cast<int>(ratios.size()) < numColPartners)
                throw FATALERROR("The number of collision partners exceeds the number of default ratios");
            for (int c = 0; c != numColPartners; ++c)
                state->setColPartnerDensity(c, state->numberDensity() * ratios[c]);
        }

        // LTE level populations (ground state for a non-positive temperature), no alignment
        Array pops(numLevels);
        for (int p = 0; p != numLevels; ++p)
            pops[p] = _model.weight[p] * exp(-_model.energy[p] / Constants::k() / Tkin);
        double sum = pops.sum();
        if (sum > 0.)
            pops *= state->numberDensity() / sum;
        else
        {
            pops = 0.;
            pops[0] = state->numberDensity();
        }
        for (int p = 0; p != numLevels; ++p)
        {
            state->setLevelPopulation(p, pops[p]);
            state->setLevelAlignment(p, 0.);
        }
    }
}

////////////////////////////////////////////////////////////////////

UpdateStatus GKLineGasMix::updateSpecificState(MaterialState* state, const Array& Jv) const
{
    Array J2v(5 * Jv.size());
    return updateSpecificStateWithAnisotropy(state, Jv, J2v);
}

////////////////////////////////////////////////////////////////////

UpdateStatus GKLineGasMix::updateSpecificStateWithAnisotropy(MaterialState* state, const Array& Jv,
                                                             const Array& J2v) const
{
    UpdateStatus status;
    if (state->numberDensity() > 0)
    {
        int numLines = _model.numLines();
        vector<double> J0(numLines, 0.);
        vector<double> J2GF(5 * numLines, 0.);

        // average the radiation tensor elements over the normalized line profile of each transition
        // (as NonLTELineGasMix does for the mean intensity)
        for (int k = 0; k != numLines; ++k)
        {
            double center = _model.center[k];
            double sigma = sigmaForLine(center, state->temperature(), _model.mass);
            double lambdamin = center - PROFILE_RANGE * sigma;
            double lambdamax = center + PROFILE_RANGE * sigma;
            int ellmin = std::lower_bound(begin(_lambdav), end(_lambdav), lambdamin) - begin(_lambdav);
            int ellmax = std::upper_bound(begin(_lambdav), end(_lambdav), lambdamax) - begin(_lambdav);
            double gsum = 0.;
            double sums[6] = {0., 0., 0., 0., 0., 0.};
            for (int ell = ellmin; ell < ellmax; ++ell)
            {
                double gdlambda = gaussian(_lambdav[ell], center, sigma) * _dlambdav[ell];
                gsum += gdlambda;
                sums[0] += Jv[ell] * gdlambda;
                for (int c = 0; c != 5; ++c) sums[c + 1] += J2v[5 * ell + c] * gdlambda;
            }
            if (abs(gsum - 1.) > MAX_GAUSS_ERROR_WARN)
            {
                string message = "Integral of Gaussian line profile over radiation field grid equals "
                                 + StringUtils::toString(gsum) + " rather than unity for " + _name + " transition ("
                                 + std::to_string(_model.indexUpRad[k]) + "-" + std::to_string(_model.indexLowRad[k])
                                 + "); refine the radiation field wavelength grid around the line";
                if (abs(gsum - 1.) > MAX_GAUSS_ERROR_FAIL) throw FATALERROR(message);
                find<Log>()->warning(message);
            }
            double norm = gsum > 0. ? 1. / gsum : 1.;
            J0[k] = sums[0] * norm;
            for (int c = 0; c != 5; ++c) J2GF[5 * k + c] = sums[c + 1] * norm;
        }

        double change = solveAndStore(state, J0, J2GF);
        if (change > maxChangeInLevelPopulations())
            status.updateNotConverged();
        else
            status.updateConverged();
    }
    return status;
}

////////////////////////////////////////////////////////////////////

double GKLineGasMix::solveAndStore(MaterialState* state, const vector<double>& J0, const vector<double>& J2GF) const
{
    int numLines = _model.numLines();
    int numLevels = _model.numLevels();

    // rotate the global-frame rank-2 components to the magnetic field frame:
    //   [J20]_MF = C20(b) J20 + 2 sum_{m=1,2} Re(C2m*(b) J2m)
    // with the Racah-normalized spherical harmonics C2m of the magnetic field direction b;
    // without magnetic field there is no preferred axis and the anisotropy is not used
    Vec B = state->magneticField();
    double Bnorm = B.norm();
    vector<double> J2MF(numLines, 0.);
    if (Bnorm > 0.)
    {
        double bx = B.x() / Bnorm, by = B.y() / Bnorm, bz = B.z() / Bnorm;
        double s32 = std::sqrt(1.5);
        double c20 = 0.5 * (3. * bz * bz - 1.);
        double c21r = -s32 * bx * bz, c21i = -s32 * by * bz;
        double c22r = 0.5 * s32 * (bx * bx - by * by), c22i = s32 * bx * by;
        for (int k = 0; k != numLines; ++k)
        {
            const double* j = &J2GF[5 * k];
            J2MF[k] = c20 * j[0] + 2. * (c21r * j[1] + c21i * j[2]) + 2. * (c22r * j[3] + c22i * j[4]);
        }
    }

    // store the radiation tensor elements
    for (int k = 0; k != numLines; ++k)
    {
        state->setLineTensor(k, 0, J0[k]);
        for (int c = 0; c != 5; ++c) state->setLineTensor(k, c + 1, J2GF[5 * k + c]);
        state->setLineTensor(k, 6, J2MF[k]);
    }

    // solve the pol-SEE
    int numColPartners = _model.numColPartners();
    vector<double> nPartner(numColPartners);
    for (int c = 0; c != numColPartners; ++c) nPartner[c] = state->colPartnerDensity(c);
    vector<double> n, sigma;
    solvePolSEE(_model, _twoJ, _coef, state->kineticTemperature(), nPartner, state->numberDensity(), J0, J2MF, n,
                sigma);

    // store the results, keeping track of the change in the populations
    double change = 0.;
    for (int p = 0; p != numLevels; ++p)
    {
        double oldPop = state->levelPopulation(p);
        state->setLevelPopulation(p, n[p]);
        state->setLevelAlignment(p, sigma[p]);
        if (n[p] > 0.)
            change += abs(oldPop / n[p] - 1.);
        else if (oldPop > 0.)
            change += 1.;
    }
    return change / numLevels;
}

////////////////////////////////////////////////////////////////////

bool GKLineGasMix::isSpecificStateConverged(int numCells, int /*numUpdated*/, int numNotConverged,
                                            MaterialState* currentAggregate, MaterialState* previousAggregate) const
{
    double fractionNotConverged = static_cast<double>(numNotConverged) / static_cast<double>(numCells);
    double changeInGlobalLevelPops = 0.;
    for (int p = 0; p != _model.numLevels(); ++p)
    {
        double currentPop = currentAggregate->levelPopulation(p);
        double previousPop = previousAggregate->levelPopulation(p);
        double diff = previousPop > 0. ? abs((currentPop - previousPop) / previousPop) : (currentPop > 0. ? 1. : 0.);
        if (diff > changeInGlobalLevelPops) changeInGlobalLevelPops = diff;
    }

    auto log = find<Log>();
    log->info("GKLineGasMix convergence info:");
    log->info("  Fraction of not converged cells is " + StringUtils::toString(fractionNotConverged * 100., 'f', 2)
              + "% (convergence criterion is " + StringUtils::toString(maxFractionNotConvergedCells() * 100., 'f', 2)
              + "%)");
    log->info("  Global level populations changed by " + StringUtils::toString(changeInGlobalLevelPops * 100., 'f', 2)
              + "% compared to previous iteration (convergence criterion is "
              + StringUtils::toString(maxChangeInGlobalLevelPopulations() * 100., 'f', 2) + "%)");

    return fractionNotConverged <= maxFractionNotConvergedCells()
           && changeInGlobalLevelPops <= maxChangeInGlobalLevelPopulations();
}

////////////////////////////////////////////////////////////////////

double GKLineGasMix::mass() const
{
    return _model.mass;
}

////////////////////////////////////////////////////////////////////

double GKLineGasMix::sectionAbs(double /*lambda*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////

double GKLineGasMix::sectionSca(double /*lambda*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////

double GKLineGasMix::sectionExt(double /*lambda*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////

namespace
{
    // applies the lower limit to a (negative) opacity, as in NonLTELineGasMix
    double limitNegativeOpacity(double opacity, const MaterialState* state, double lowestOpticalDepth)
    {
        if (opacity < 0.)
        {
            double diagonal = 1.7320508 * cbrt(state->volume());  // correct only for cubical cell
            if (opacity * diagonal < lowestOpticalDepth) opacity = lowestOpticalDepth / diagonal;
        }
        return opacity;
    }
}

////////////////////////////////////////////////////////////////////

double GKLineGasMix::opacityAbs(double lambda, const MaterialState* state, const PhotonPacket* /*pp*/) const
{
    double opacity = 0.;
    if (state->numberDensity() > 0.)
    {
        constexpr double front = Constants::h() * Constants::c() / 4. / M_PI;
        for (int k = 0; k != _model.numLines(); ++k)
        {
            double center = _model.center[k];
            double sigma = sigmaForLine(center, state->temperature(), _model.mass);
            if (std::abs(lambda - center) <= PROFILE_RANGE * sigma)
            {
                double transrate = state->levelPopulation(_model.indexLowRad[k]) * _model.einsteinBlu[k]
                                   - state->levelPopulation(_model.indexUpRad[k]) * _model.einsteinBul[k];
                opacity += front / center * transrate * gaussian(lambda, center, sigma);
            }
        }
        opacity = limitNegativeOpacity(opacity, state, lowestOpticalDepth());
    }
    return opacity;
}

////////////////////////////////////////////////////////////////////

double GKLineGasMix::opacitySca(double /*lambda*/, const MaterialState* /*state*/, const PhotonPacket* /*pp*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////

double GKLineGasMix::opacityExt(double lambda, const MaterialState* state, const PhotonPacket* pp) const
{
    return opacityAbs(lambda, state, pp);
}

////////////////////////////////////////////////////////////////////

bool GKLineGasMix::peeloffScattering(double& /*I*/, double& /*Q*/, double& /*U*/, double& /*V*/, double& /*lambda*/,
                                     Direction /*bfkobs*/, Direction /*bfky*/, const MaterialState* /*state*/,
                                     const PhotonPacket* /*pp*/) const
{
    return false;
}

////////////////////////////////////////////////////////////////////

void GKLineGasMix::performScattering(double /*lambda*/, const MaterialState* /*state*/, PhotonPacket* /*pp*/) const {}

////////////////////////////////////////////////////////////////////

void GKLineGasMix::polarizedOpacitiesExt(double lambda, const MaterialState* state, const PhotonPacket* /*pp*/,
                                         double cosTheta, double& kpar, double& kper) const
{
    kpar = 0.;
    kper = 0.;
    if (state->numberDensity() > 0.)
    {
        constexpr double front = Constants::h() * Constants::c() / 4. / M_PI;
        double fpar = (3. * cosTheta * cosTheta - 2.) / M_SQRT2;
        double fper = 1. / M_SQRT2;
        bool aligned = !state->magneticField().isNull();
        for (int k = 0; k != _model.numLines(); ++k)
        {
            double center = _model.center[k];
            double sigma = sigmaForLine(center, state->temperature(), _model.mass);
            if (std::abs(lambda - center) <= PROFILE_RANGE * sigma)
            {
                int up = _model.indexUpRad[k];
                int low = _model.indexLowRad[k];
                double kl = state->levelPopulation(low) * _model.einsteinBlu[k];
                double ku = state->levelPopulation(up) * _model.einsteinBul[k];
                double al = aligned ? _wLU[k] * state->levelAlignment(low) : 0.;
                double au = aligned ? _wUL[k] * state->levelAlignment(up) : 0.;
                double profile = front / center * gaussian(lambda, center, sigma);
                kpar += profile * (kl * (1. + al * fpar) - ku * (1. + au * fpar));
                kper += profile * (kl * (1. + al * fper) - ku * (1. + au * fper));
            }
        }
        kpar = limitNegativeOpacity(kpar, state, lowestOpticalDepth());
        kper = limitNegativeOpacity(kper, state, lowestOpticalDepth());
    }
}

////////////////////////////////////////////////////////////////////

Array GKLineGasMix::lineEmissionCenters() const
{
    return NR::array(_model.center);
}

////////////////////////////////////////////////////////////////////

Array GKLineGasMix::lineEmissionMasses() const
{
    Array masses(_model.numLines());
    masses = _model.mass;
    return masses;
}

////////////////////////////////////////////////////////////////////

Array GKLineGasMix::lineEmissionSpectrum(const MaterialState* state, const Array& /*Jv*/) const
{
    Array luminosities(_model.numLines());
    if (state->numberDensity() > 0.)
    {
        double front = Constants::h() * Constants::c() * state->volume();
        for (int k = 0; k != _model.numLines(); ++k)
            luminosities[k] =
                front / _model.center[k] * _model.einsteinA[k] * state->levelPopulation(_model.indexUpRad[k]);
    }
    return luminosities;
}

////////////////////////////////////////////////////////////////////

double GKLineGasMix::lineEmissionAlignmentFactor(const MaterialState* state, int k) const
{
    if (state->numberDensity() <= 0. || state->magneticField().isNull()) return 0.;
    return _wUL[k] * state->levelAlignment(_model.indexUpRad[k]);
}

////////////////////////////////////////////////////////////////////

double GKLineGasMix::indicativeTemperature(const MaterialState* state, const Array& /*Jv*/) const
{
    return state->temperature();
}

////////////////////////////////////////////////////////////////////
