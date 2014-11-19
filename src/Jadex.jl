module Jadex
### Sections
# * Constants
# * Data file parser
# * Grid calculator
### TODO
# *
###

import Base: show, showcompact


##############################################################################
# Constants
##############################################################################
# Physical constants in CGS
const clight  = 2.99792458e10   # speed of light     (cm/s)
const hplanck = 6.6260963e-27   # Planck constant    (erg/Hz)
const kboltz  = 1.3806505e-16   # Boltzmann constant (erg/K)
const amu     = 1.67262171e-24  # atomic mass unit   (g)
const fk      = hplanck * clight / kboltz
const thc     = 2 * hplanck * clight

# Mathematical constants
const fgauss  = 1.0645 * 8.0 * pi

# Computational constants
const maxpart = 9      # max number of collision partners
const maxtemp = 99     # max number of collision temperatures
const maxlev  = 2999   # max number of energy levels
const mxline  = 99999  # max number of energy levels
const maxcoll = 99999  # max number of collisional transitions

const miniter = 10     # minimum number of iterations
const maxiter = 9999   # maximum number of iterations

const ccrit   = 1e-6   # relative tolerance of solution
const eps     = 1e-30  # round-off error
const minpop  = 1e-20  # minimum level population

# Directories
const datadir = "../../data/"


##############################################################################
# Data file parser
##############################################################################
# Parse molecular data files and create input types.

immutable Molecule
    # Header
    specref::String  # molecule
    amass::Integer  # molecular weight
    # Energy levels
    nlev::Integer  # number of energy levels
    eterm::Array{FloatingPoint,1}  # energy levels, in cm^-1
    gstat::Array{FloatingPoint,1}  # statistical weights
    qnum::Array{ASCIIString,1}  # quantum numbers
    # Transitions
    nline::Integer  # number of radiative transitions
    iupp::Array{Integer,1}  # upper state
    ilow::Array{Integer,1}  # lower state
    aeinst::Array{FloatingPoint,1}  # Einstein A
    spfreq::Array{FloatingPoint,1}  # spectral line frequencies
    eup::Array{FloatingPoint,1}  # upper state energy, E_u / k
    xnu::Array{FloatingPoint,1}  # difference in energy levels between up/low
    # Collision rates
    npart::Integer  # number of collision partners
    colliders::Array{CollisionPartner,1}  # list of colliders
end


immutable CollisionPartner
    collref::String  # collision partner reference name
    ncoll::Integer  # number of collisional transitions
    ntemp::Integer  # number of collisional temperatures
    temp::Array{FloatingPoint,1}  # temperatures
    lcu::Array{FloatingPoint,1}  # upper state of collision
    lcl::Array{FloatingPoint,1}  # lower state of collision
    coll::Array{FloatingPoint,2}  # collision rates, cm^3 s^-1
end
const valid_partners = ["h2", "p-h2", "o-h2", "e", "h", "he", "h+"]


function readdata(specref::String)
    # Read in data file
    f = datadir + specref + ".dat" |> open |> readlines
    f = [strip(l) for l in f]
    # Parse header
    amass = f[4] |> float
    # Parse energies
    nlev = f[6] |> int
    eterm = Array(FloatingPoint, nlev)
    gstat = Array(FloatingPoint, nlev)
    qnum = Array(String, nlev)
    for ii=8:7+nlev
        l = f[ii] |> split
        eterm[ii] = l[2] |> float
        gstat[ii] = l[3] |> float
        qnum[ii] = l[4]
    end
    # Parse transitions
    nline = f[9+nlev] |> int
    iupp = Array(Integer, nline)
    ilow = Array(Integer, nline)
    aeinst = Array(FloatingPoint, nline)
    spfreq = Array(FloatingPoint, nline)
    eup = Array(FloatingPoint, nline)
    xnu = Array(FloatingPoint, nline)
    for ii=11+nlev:10+nlev+nline
        l = f[ii] |> split
        iupp[ii] = l[2] |> int
        ilow[ii] = l[3] |> int
        aeinst[ii] = l[4] |> float
        spfreq[ii] = l[5] |> float
        eup[ii] = l[6] |> float
        xnu[ii] = eterm[iupp[ii]] - eterm[ilow[ii]]
    end
    # Parse collision partners
    npart = f[13+nlev+nline] |> int
    colliders = Array(CollisionPartner, npart)
    for ii, line in enumerate(f)
        if line == "!COLLISIONS BETWEEN"
            collref = valid_partners[int(f[ii+1])]
            # TODO
            colp = CollisionPartner(collref, ncoll, ntemp, temp, lcu, lcl, coll)
            colliders[ii] = colp
        end
    end
    Molecule(specref, amass, nlev, eterm, gstat, qnum, nline, iupp, ilow,
        aeinst, spfreq, eup, xnu, npart, colliders)
end


##############################################################################
# Run descriptor
##############################################################################
# Describe a transfer run

immutable RunDef
    mol::Molecule
    outn::String
    freq::(Real, Real)
    tkin::Real
    colps::Array{CollisionPartner, 1}
    tbg::Real
    coldens::Real
    linewidth::Real
    geometry::String
end


##############################################################################
# Matrix
##############################################################################
# Compute matrix

function matrix(rundef::RunDef, niter::Integer, conv::Bool)
    nlev = rundef.mol.nlev
    totdens

    reducem = false
    # Initialize rate matrix
    rhs = zeros(nlev+1)
    yrate = zeros(nlev+1, nlev+1)
    for ilev=1:nlev
        for jlev=1:nlev
            yrate[ilev, jlev] = -eps * totdens
        end
        yrate[nlev+1, ilev] = 1.0
        rhs[ilev] = eps * totdens
        yrate[ilev, nlev+1] = eps * totdens
    end
    rhs[nlev+1] = eps * totdens

    if niter == 0
        for iline=1:nline
            m = iupp[iline]
            n = ilow[iline]
            etr = fk * xnu[iline] / trj[iline]
            if etr >= 160.0
                exr = 0.0
            else
                exr = 1.0 / (exp(etr) - 1.0)
            end
        end
    else
        # Subsequent iterations: use escape probability
        cddv = cdmol / deltav
        # Count optically thick lines
        nthick = 0
        nfat = 0

        for iline=1:nline
            xt = xnu[iline]^3.0
            m  = iupp[iline]
            n  = ilow[iline]
            # Calculate source fn
            hnu = fk * xnu[iline] / tex[iline]
            if hnu >= 160.0
                bnutex = 0.0
            else
                bnutex = thc * xt / (exp(fk * xnu[iline] / tex[iline]) - 1.0)
            end
            # Calculate line optical depth
            taul[iline] = cddv * (xpop[n] * gstat[m] / gstat[n] - xpop[m]) / (fgaus * xt / aeinst[iline])
            if taul[iline] > 1e-2 nthick += 1 end
            if taul[iline] > 1e5  nfat   += 1 end
            # Use escape probability approx for internal intensity
            beta = escprob[taul[iline]]
            bnu  = totalb[iline] * beta
            exr  = bnu / (thc * xt)
            # Radiative contribution to the rate matrix
            yrate[m,m] = yrate[m,m] + aeinst[iline] * (beta + exr)
            yrate[n,n] = yrate[n,n] + aeinst[iline] * (gstat[m] * exr / gstat[n])
            yrate[m,n] = yrate[m,n] - aeinst[iline] * (gstat[m] / gstat[n] * exr)
            yrate[n,m] = yrate[n,m] - aeinst[iline] * (beta + exr)
        end
    end

    # Warn user if convergence problems expected
    if niter == 1 && nfat > 0
        print("*** Warning: Some lines have very high optical depth")
    end

    # Contribution for collisional processes to the rate matrix
    for ilev=1:nlev
        yrate[ilev,ilev] = yrate[ilev,ilev] + ctot[ilev]
        for jlev=1:nlev
            if ilev != jlev
                yrate[ilev,jlev] -= crate[jlev,ilev]
            end
        end
    end

    if reducem
        # An auxillary array is passed to the linear equation solver after
        # renormalization. The array `yrate` retains the original matrix
        # elements.
        uarray = copy(yrate)
        # Test whether the matrix should be reduced to exclude the radiatively
        # coupled levels.
        redcrit = 10.0 * tkin / fk
        nreduce = 0
        for ilev=1:nlev
            if eterm[ilev] <= redcrit
                nreduce += 1
            end
        end
        # We now separate the collisionally coupled levels from thsoe that are
        # now coupled mainly by radiatve processes, compute an effective
        # cascade matrix for transfer rates from one low-lying level to another
        # and then solve this reduced system of equations explicitly for the
        # low-lying levels only.
        for jlev=1:nreduce
            for ilev=1:nreduce
                for klev=nreduce+1:nlev
                    uarray[ilev,jlev] = abs(yrate[klev,jlev] * yrate[ilev,klev] / yrate[klev,klev]) + uarray[ilev,jlev]
                end
            end
        end
        # Invert the reduced matrix `uarray`
        # TODO ccall uarray
    else
        # TODO ccall yrate
    end

    # Level populations are the normalized RHS components
    total = sum(rhs)

    # Limit population to minpop
    for ilev=1:nlev
        xpopold[ilev] = xpop[ilev]
        xpop[ilev] = max(minpop, rhs[ilev] / total)
        # if first iteration, no old population
        if niter == 0
            xpopold[ilev] = xpop[ilev]
        end
    end

    # Compute excitation temperatures of the lines
    tsum = 0.0
    for iline=1:nline
        m  = iupp[iline]
        n  = ilow[iline]
        xt = xnu[iline]^3
        if niter == 0
            if xpop[n] <= minpop || xpop[m] <= minpop
                tex[iline] = totalb[iline]
            else
                tex[iline] = fk * xnu[iline] / (log(xpop[n] * gstat[m] / (xpop[m] * gstat[n])))
            end
        else
            if xpop[n] <= minpop || xpop[m] <= minpop
                thistex = tex[iline]
            else
                thistex = fk * xnu[iline]
            end
            if taul[iline] > 0.01
                tsum += abs((thistex - tex[iline]) / thistex)
            end
            tex[iline] = 0.5 * (thistex + tex[iline])
            taul[iline] = cddv * (xpop[n] * gstat[m] / gstat[n] - xpop[m]) / (fgauss * xt / aeinst[iline])
        end
    end

    # Introduce a minimum number of iterations
    if niter >= miniter
        if nthick == 0; conv = true end
        if tsum / nthick < ccrit; conv = true end
    end

    # Now do the underrelaxation
    for ilev=1:nlev
        xpop[ilev] = 0.3 * xpop[ilev] + 0.7 * xpopold[ilev]
    end
end


function escprob(τ::Real, geometry::String="sphere")
    τr = τ / 2.0
    # Uniform sphere formula from Osterbrock (Astrophysics of Gaseous Nebulae
    # and Active Galactic Nuclei) Appendix 2 with power law approximations for
    # large and small tau
    if geometry == "sphere"
        if abs(τr) < 0.1
            β = 1.0 - 0.75 * τr + τr^2 / 2.5 - τr^3 / 6.0 + τr^4 / 17.5
        elseif abs(τr) > 50
            β = 0.75 / τr
        else
            β = 0.75 / τr * (1.0 - 1.0 / 2τr^2) + (1.0 / τr + 1.0 / 2τr^2) * exp(-τ)
        end
    # Expanding sphere, Large Velocity Gradient, or Sobolev case. Formular from
    # De Jong, Boland, and Dalgarno (1980, A&A 91, 68).
    # Corrected by factor of 2 in order to return 1 for tau=1.
    elseif geometry == "lvg"
        if abs(τr) < 0.01
            β = 1.0
        elseif abs(τr) < 7.0
            β = 2.0 * (1.0 - exp(-2.34τr)) / 4.65τr
        else
            β = 4τr * (sqrt(log(τr / sqrt(π)))) \ 2.0
        end
    # Slab geometry (e.g. shocks): de Jong, Dalgarno, and Chu (1975), ApJ 199,
    # 69 (again with power law approximations)
    elseif geometry == "slab"
        if abs(3τ) < 0.1
            β = 1.0 - 1.5 * (τ + τ^2)
        elseif abs(3τ) > 50.0
            β = 1.0 / 3τ
        else
            β = (1.0 - exp(-3τ)) / 3τ
        end
    else
        throw(ArgumentError())
    end
    β
end

##############################################################################
# Grid calculator
##############################################################################
# Create grids of `RunDef`s and add definitions to relevant functions to
# accept grid input.
