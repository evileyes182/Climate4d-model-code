# Figure 12: boundary-crisis curve BC of system (10), a fold of homoclinic orbits of a saddle periodic orbit,
# in the (mu, eta)-plane
# Run: julia plot.jl (after auto run.py); writes fig12_boundary_crisis_lin.png
using CairoMakie, DelimitedFiles

fig = Figure(size = (500, 500))
ax = Axis(fig[1, 1], xlabel = "μ (×10⁻³)", ylabel = "η", limits = (-2.069, -2.054, -4.17, -4.085))

d = readdlm("bc.dat")
lines!(ax, d[:, 1] .* 1e3, d[:, 2], color = "#b40005", linewidth = 2)

save("fig12_boundary_crisis_lin.png", fig, px_per_unit = 3)
