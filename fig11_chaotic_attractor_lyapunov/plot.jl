# Figure 11: chaotic attractor of system (10) at (mu, eta) = (-2.128e-3, -3.99), projected onto the (y_E, y_N)-plane,
# and time series of K_N and Psi, with convective shutdown events
using CairoMakie, DelimitedFiles

d  = readdlm("attractor.dat")
ev = readdlm("events.dat")
eq = readdlm("equilibrium.dat")
t, yN, yE, KN, Psi = d[:, 1] .* 1e-5, d[:, 2], d[:, 3], d[:, 6], d[:, 7] .* 1e3
cr = extrema(Psi)

fig = Figure(size = (500, 800))
ax1 = Axis(fig[1, 1], xlabel = "y_E", ylabel = "y_N", limits = (1.8, 3.1, 0.6, 2.6))
lines!(ax1, yE, yN, color = Psi, colormap = Reverse(:roma), colorrange = cr, linewidth = 1)
scatter!(ax1, ev[:, 3], ev[:, 2], color = "#f27234", strokecolor = :black, strokewidth = 0.75, markersize = 7)
scatter!(ax1, [eq[2, 1]], [eq[1, 1]], color = :black, strokecolor = :black, strokewidth = 0.75, markersize = 7)   # equilibrium
scatter!(ax1, [eq[2, 2]], [eq[1, 2]], color = :white, markersize = 3)                                              # shooting point

ax2 = Axis(fig[2, 1], ylabel = "K_N", yticks = ([1e-2/2.1e4, 8/2.1e4], ["κ_d", "κ_N"]), limits = (0, 4, -2e-5, 4e-4))
lines!(ax2, t, KN, color = Psi, colormap = Reverse(:roma), colorrange = cr, linewidth = 1.5)
scatter!(ax2, ev[:, 1] .* 1e-5, ev[:, 6], color = "#f27234", strokecolor = :black, strokewidth = 0.75, markersize = 7)

ax3 = Axis(fig[3, 1], xlabel = "t (×10⁵)", ylabel = "Ψ (×10⁻³)", limits = (0, 4, 1.25, 2.5))
lines!(ax3, t, Psi, color = Psi, colormap = Reverse(:roma), colorrange = cr, linewidth = 1.5)
scatter!(ax3, ev[:, 1] .* 1e-5, ev[:, 7] .* 1e3, color = "#f27234", strokecolor = :black, strokewidth = 0.75, markersize = 7)

rowsize!(fig.layout, 1, Aspect(1, 1))
save("fig11_chaotic_attractor_lyapunov.png", fig, px_per_unit = 3)
