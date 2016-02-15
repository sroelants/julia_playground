#------------------------------------------------------------------------------
# List of physical constants, IN ATOMIC UNITS, that make it easier to express 
# certain relationships and convert quantities.
#
# Basically, this means:
#   - m_e = 1.0         (electron mass)
#   - e = 1             (electron charge)
#   - ħ = 1.0           (Planck's constant)
#   - 1/(4*pi*ϵ) = 1.0  (Coulomb's constant)
#
# Distances are expressed in Bohr,
# Energies are expressed in Hartrees.
#-------------------------------------------------------------------------------

module Constants
export ħ, m_e, h22m, eV, eV_Ha, nm, kB_eV, kB


ħ = 1.0
m_e = 1.0
h22m = ħ^2/(2*m_e)
eV = 1/27.21138505
eV_Ha = eV
nm = 18.89726124565
kB_eV = 8.6173324e-5
kB = kB_eV * eV_Ha
end
