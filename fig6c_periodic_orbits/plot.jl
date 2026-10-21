# Figure 6(c): equilibria and periodic orbits of system (10) in mu for eta = -4.08,
# plotted as min(Psi) and max(Psi)
# Run: julia plot.jl (after python3 run.py); writes fig6c_periodic_orbits.png
using CairoMakie, DelimitedFiles

fig = Figure(size = (500, 350))
ax = Axis(fig[1, 1], xlabel = "μ (×10⁻³)", ylabel = "min(Ψ) & max(Ψ) (×10⁻³)", limits = (-2.165, -2.025, 1.0, 2.45))

for name in ["eq_1", "eq_2"]
    d = readdlm(name * ".dat")
    mu, psi, stable = d[:, 1] .* 1e3, d[:, 2] .* 1e3, d[:, 4] .== 1
    lines!(ax, mu, psi, color = "#231f20", linewidth = 2)                                         # unstable
    lines!(ax, mu, ifelse.(stable .& (psi .> 0), psi, NaN), color = "#4471b8", linewidth = 2)  # stable northward
    lines!(ax, mu, ifelse.(stable .& (psi .< 0), psi, NaN), color = "#68cce1", linewidth = 2)  # stable southward
end

d = readdlm("po_1.dat")
mu, psimin, psimax, stable = d[:, 1] .* 1e3, d[:, 2] .* 1e3, d[:, 3] .* 1e3, d[:, 4] .== 1
for psi in (psimin, psimax)
    lines!(ax, mu, psi, color = "#a298cb", linewidth = 2)                                # unstable
    lines!(ax, mu, ifelse.(stable, psi, NaN), color = "#f6861f", linewidth = 2)          # stable
end

for (name, ty, col) in [("eq_1", "HB", "#69aa44"), ("eq_1", "LP", "#231f20"), ("po_1", "LP", "#2682c5"), ("po_1", "PD", "#a64c9d"), ("po_1", "TR", "#4dc0ac")]
    s = readdlm(name * "_special.dat")
    k = s[:, 4] .== ty
    for psi in (s[k, 2], s[k, 3])
        scatter!(ax, s[k, 1] .* 1e3, psi .* 1e3, color = col, strokecolor = :black, strokewidth = 0.75, markersize = 9)
    end
end

save("fig6c_periodic_orbits.png", fig, px_per_unit = 3)
