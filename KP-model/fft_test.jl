L = 2.0
θ(x) = 0.5*(1+sign(x))
hat(x) = 0.5*(θ(x+0.5) - θ(x-0.5))

xs = collect(linspace(-L/2, L/2, 10))
components = fft(hat(xs))
orig = ifft(components)

using PyPlot
plot(xs, abs2(orig))
savefig("hat.svg")
