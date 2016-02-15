
FFTW.set_num_threads(4)
# define a potential function
# v(x) = x.^4 - x.^2 # Mexican hat 
# v(x) = x.^2 # Harmonic oscillator
v(x) = 100.0*(sign(-(x-0.5)) + sign(x-0.5)) # Square well
#v(x) = 0.0*x

# K-space definition of the Hamiltonian
hamiltonian(k, v) = 0.5*diagm(k.^2) + v

# Define the problem
L = 2.0 # Size of system
N = 1000 # Number of samples
xs = collect(linspace(-L/2.0, L/2.0, N))
ns = 1:N
ks = Array(1:N)

# Find the potential in k-space
vxs = v(xs)
vkn = fft(vxs)
vmn = [ m > n ? vkn[m-n+1] : conj(vkn[n-m+1]) for m=ns, n=ns ].'
H = convert(Array{Complex, 2}, hamiltonian(ks, vmn))

# Diagonalize H
E = eigfact(H)
ψ1x = ifft(E[:vectors][:,1])
#= ψ2x = ifft(E[:vectors][:,2]) =#
#= ψ3x = ifft(E[:vectors][:,3]) =#
#= ψ4x = ifft(E[:vectors][:,4]) =#

θs = angle(ψ1x)

using PyPlot
#= plot(xs, θs, color="purple", linewidth=1.0, linestyle="--") =#
#= plot(xs, vxs/1000, color="black", linewidth=1.0, linestyle="--") =#
#= plot(xs, E[:values], color="red", linewidth=2.0, linestyle="--") =#
plot(xs[1:30], abs(E[:vectors][:,1])[1:30], color="red", linewidth=2.0, linestyle="--")
#= plot(xs, imag(ψ1x), color="blue", linewidth=2.0, linestyle="--") =#
#= plot(xs, imag(ψ2x), color="blue", linewidth=2.0, linestyle="--") =#
#= plot(xs, abs2(ψ3x), color="pink", linewidth=2.0, linestyle="--") =#
#= plot(xs, abs2(ψ4x), color="orange", linewidth=2.0, linestyle="--") =#
savefig("plot.svg")
