# ------------------------------------------------------------------------------
# Bogoliubov - de Gennes solver
#
# At the moment, only 2D systems are implemented (though the Shape type stores
# 3 floats, the first two are ignored).
# ------------------------------------------------------------------------------
module BdgSolver
using Constants
#= using Roots =#
export Material, Shape, Parameters, System, Hamiltonian, print_parameters

const NKSI = 50

immutable Material
    name::AbstractString
    ρ::Float64
    ħω::Float64
    λ::Float64
end


type Shape
    Lx::Float64
    Ly::Float64
    Lz::Float64
end


# ------------------------------------------------------------------------------
# Parameters type definition and constructor
# ------------------------------------------------------------------------------
type Parameters
    μ::Float64
    ν::Int32
    kmax::Float64
end


function Parameters(material, shape)
    μ = chemical_potential(material.ρ, material.ħω, shape.Lz)
    ν = get_ν(μ, material.ħω, shape.Lz)
    kmax = get_kmax(material.ħω, μ, shape.Lz)

    Parameters(μ, ν, kmax)
end


""" Calculate the number of bands below the Fermi energy.
    Input:
        - μ: Fermi energy (chemical potential)
        - ħω: Debye energy
        - Lz: Thickness of the system
    Output:
        - ν: Number of bands
"""
get_ν(μ, ħω, Lz) = floor(Integer, Lz/π * sqrt((μ + ħω)/h22m))


""" Calculate the maximal wavevector below the Fermi energy.
    Input:
        - μ: Fermi energy (chemical potentia)
        - ħω: Debye energy
        - Lz: Film thickness
    Output:
        - Maximal wavevector
"""
get_kmax(μ, ħω, Lz) = sqrt((μ + ħω - π^2 / Lz^2) / h22m)


""" Calculates the chemical potential self-consistently.

The chemical potential determines the number of occupied bands ν. ν, however,
also determines the position of the chemical potential μ. Hence, we must solve
for μ self-consistently. [See reference]

    Input:
        - ρ: Charge carrier density, a material dependent parameter
        - ħω: Debye energy
        - Lz: The thickness of the film.
    Output:
        - μ: The chemical potential.
"""
function chemical_potential(ρ, ħω, Lz)
    ϵ = 0.000001  # Tolerance for the self-consistent loop
    μ_old = 0.0
    μ = 2*h22m*π*ρ  # Naive definition of μ

    while abs(μ - μ_old) > ϵ
        μ_old = μ
        ν = get_ν(μ_old, ħω, Lz)
        μ = calculate_μ(ρ, ν, Lz)
    end
    μ
end


""" Calculates the chemical potential from ν, for a 2D system, according to
    [ref].

    Input:
        - ρ: Charge carrier density
        - ν: Number of occupied bands
        - Lz: Thickness of the system
    Output:
        - μ: The chemical potential for given ρ, ν, Lz.
"""
function calculate_μ(ρ, ν, Lz)
    2 * h22m * π * Lz * (ρ + π/(6*Lz^3) * ν * (ν + 0.5) * (ν + 1)) / ν
end



# -----------------------------------------------------------------------------
# System and Hamiltonian type 
# -----------------------------------------------------------------------------
type DensityOfStates
    ξ
    N  # The DOS as an array.
    #corrections::Array<T>  # An array of functions.
end


type Hamiltonian
    DOS::DensityOfStates
end


function Hamiltonian(material::Material, shape::Shape, parameters::Parameters)
    """ Set the DOS to the standard 2D DOS (step functions). 
    
    Because we want to be able to apply corrections to this, we store the DOS as
    a 2D array.

    Input:
        - material: The material, needed to get ħω
        - shape: Needed for Lz
        - parameters: Needed for μ and ν

    """
    DOS(i, ξ) = θ( ξ .- h22m * π^2 * (i.' + 1).^2 / shape.Lz^2)/ 2 / π / h22m

    is = collect(range(0, parameters.ν))
    ξs = collect(linspace(2*(parameters.μ - material.ħω),
                          2*(parameters.μ + material.ħω), NKSI))
    N = DOS(is, ξs)
    Hamiltonian(DensityOfStates(ξs, N))
end


type System
    material::Material
    shape::Shape
    parameters::Parameters
    H::Hamiltonian
end


function System(material::Material, shape::Shape, parameters::Parameters)
    H = Hamiltonian(material, shape, parameters)
    System(material, shape, parameters, H)
end


""" Get Tc by solving the determinantal equation described in [Ref].

    We calculate the thermal integral, multiply by the overlaps, and calculate the
    determinant of this matrix for several values of T. The determinant should 
    vanish for T = Tc.
        
        Input:
            - system
        Output:
            - Tc: obviously.

    This is also where we apply possible corrections to the DOS.
"""


function get_Tc(system::System)
    χ_DW(ξ) = θ(ξ - system.parameters.μ + system.material.ħω) -
              θ(ξ - system.parameters.μ - system.material.ħω)

    ξlims = (system.parameters.μ - 2*system.material.ħω,
             system.parameters.μ + 2*system.material.ħω)
    ξs = collect(linspace(ξlims[1], ξlims[2], NKSI))
    Φ = calculate_overlaps(system.parameters, system.shape) 
    N_mu = evaluate_at(system.H.DOS, system.parameters.μ)
    g = system.material.λ/N_mu
    
    println("N(μ) = $N_mu")
    println("g = $g")

    N = system.H.DOS.N .* χ_DW(ξs)
    D(β) = thermal_det(β, N, ξs, Φ, g)
    
    βs = collect(linspace(0.0, 100.0, 100))
    Ds = map(D, βs)
    (βs, Ds)
    #= Roots.fzero(D, [0.1, 10]) =#
end


function thermal_det(β, N, ξs, Φ, g)
    """ Calculate the thermal determinant that discriminates between
    superconducting and non-superconducting regimes. Tc is determined by the
    condition det(M - I) = 0.

        Input:
            - β: Inverse temperature
            - N: A (discretely sampled) Density of States
            - ξs: Energy samples along which we integrate.
    """
    dξ = (ξs[end] - ξs[1])/NKSI  # Or, simply ξs[2]-ξs[1] ...
    weight = F(ξs, β)
    I = dξ * weight.' * N  # Thermal integral
    det( g/4 * Φ .* diag(I) - eye(length(I))) 
end

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

""" Heaviside function """
θ(x) = 0.5 * (1 + sign(x))


""" Integrate function (actually functional :) )

    Take in a function of a single variable (probably a closure), sample it on
    a range, (that is passed as well), and integrate!

    Input:
        - f: Function (of a single variable!) to be integrated
        - I: Integration interval as a tuple
    Output:
        - The integral of f over I.

    The plan is to make this as dimension-agnostic as possible.

    CAVEAT: If the function returns a multidimensional array, we always assume
        the index to be integrated over to be the FIRST index.
"""
function integrate(f, I)
    dx = 0.001
    xmin, xmax = I
    xs = collect(linspace(xmin, xmax, floor(Integer, ((xmax - xmin)/dx))))
    sum(f(xs), 1)*dx
end


""" Get the Schrödinger wavefunction overlaps """
function calculate_overlaps(parameters::Parameters, shape::Shape)
    (1 + eye(parameters.ν)/2)/shape.Lz
end


""" Define the standard thermal weight function.

    In absence of corrections, this is simply
    F(ξ) = tanh(βξ / 2) / ξ

    Input:
        - ξ: energy (could be a range, hopefully)
        - β: Inverse temperature (1/(kB T) )
    Output:
        - F(ξ), clearly...
"""
F(ξ, β) = tanh(β .* ξ / 2) ./ ξ


function evaluate_at(DOS::DensityOfStates, E)
    """ Return the DOS (summed over all bands i) at a given energy E

    Since we store the DOS as a discrete array, the evaluation is only
    approximate. We look for the nearest match to E in the ξ-array stored in 
    the DOS, and return the corresponding value of N (again, summed over the 
    second axis).

    Input:
        - DOS: A density of states 
        - E: the energy at which to evaluate the DOS.
    Output:
        - The density of states at E.
    """
    index = indmin(abs(DOS.ξ - E))
    sum(DOS.N, 2)[index]
end


# -----------------------------------------------------------------------------
# Printing routines
# -----------------------------------------------------------------------------

function print_parameters(material::Material)
    str="""
--------------------------------------------------------------------------------
--- Material properties --------------------------------------------------------
--------------------------------------------------------------------------------
--- Carrier density:            ρ   = $(material.ρ)         Bohr^{-1},
--- Debye energy:               ħω  = $(material.ħω)        Ha,
--- Electron phonon coupling:   λ   = $(material.λ)
--------------------------------------------------------------------------------
"""
    print(str)
end


function print_parameters(shape::Shape)
    str="""
--------------------------------------------------------------------------------
--- Shape dimensions -----------------------------------------------------------
--------------------------------------------------------------------------------
--- Lx  =   $(shape.Lx) Bohr,   Ly  =   $(shape.Ly) Bohr,   Lz = $(shape.Lz) Bohr
--------------------------------------------------------------------------------
"""
    print(str)
end


function print_parameters(parameters::Parameters)
    str="""
--------------------------------------------------------------------------------
--- Derived parameters  --------------------------------------------------------
--------------------------------------------------------------------------------
--- Chemical potential          μ       = $(parameters.μ)       Ha,
--- Maximum band index          ν       = $(parameters.ν),
--- Maximum wavevector          kmax    = $(parameters.kmax)    Bohr^{-1}
--------------------------------------------------------------------------------
"""
    print(str)
end


function print_parameters(system::System)
    print_parameters(system.shape)
    print_parameters(system.material)
    print_parameters(system.parameters)
end


end
