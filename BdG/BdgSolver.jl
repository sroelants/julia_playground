# ------------------------------------------------------------------------------
# Bogoliubov - de Gennes solver
#
# At the moment, only 2D systems are implemented (though the Shape type stores
# 3 floats, the first two are ignored).
# ------------------------------------------------------------------------------
module BdgSolver
using Constants
export Material, Shape, Parameters, System, Hamiltonian, integrate

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
    μ = 2*h22m*pi*ρ  # Naive definition of μ

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
    2 * h22m * pi * Lz * (ρ + pi/(6*Lz^3) * ν * (ν + 0.5) * (ν + 1)) / ν
end


# -----------------------------------------------------------------------------
# System and Hamiltonian type 
# -----------------------------------------------------------------------------
type Hamiltonian
    DOS
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
    DOS(i, ξ) = 1/h22m * θ( ξ .- h22m * π^2 * (i.' + 1).^2 / shape.Lz^2)

    is = collect(range(0, parameters.ν))
    ξs = collect(linspace(2*(parameters.μ - material.ħω),
                          2*(parameters.μ + material.ħω), NKSI))
    N = DOS(is, ξs)
    Hamiltonian(N)
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

    N = system.H.DOS .* χ_DW(ξs)
    println(N)
    D(β) = thermal_det(β, N, ξs)
    D(1.1)
end

function thermal_det(β::Float64, N, ξs)
    dξ = (ξs[end] - ξs[1])/NKSI
    weight = F(ξs, β)

    I = dξ * weight.' * N
    prod(I)
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
F(ξ, β) = tanh(β * ξ / 2) / ξ

end
