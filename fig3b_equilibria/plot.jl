# Figure 3(b): equilibria of system (10) in mu for eta = -1.5, plotted as Psi
# Run: julia plot.jl (after auto run.py); writes fig3b_equilibria.png
using CairoMakie, DelimitedFiles

fig = Figure(size = (500, 300))
ax = Axis(fig[1, 1], xlabel = "μ (×10⁻³)", ylabel = "Ψ (×10⁻³)", limits = (-5, 0, -1.5, 2.5))

for name in ["eq_1", "eq_2"]
    d = readdlm(name * ".dat")
    mu, psi, stable = d[:, 1] .* 1e3, d[:, 2] .* 1e3, d[:, 3] .== 1
    lines!(ax, mu, psi, color = "#939598", linewidth = 2)                                         # unstable
    lines!(ax, mu, ifelse.(stable .& (psi .> 0), psi, NaN), color = "#4471b8", linewidth = 2)  # stable northward
    lines!(ax, mu, ifelse.(stable .& (psi .< 0), psi, NaN), color = "#68cce1", linewidth = 2)  # stable southward

    s = readdlm(name * "_special.dat")
    hb, lp = s[:, 3] .== "HB", s[:, 3] .== "LP"
    scatter!(ax, s[hb, 1] .* 1e3, s[hb, 2] .* 1e3, color = "#69aa44", strokecolor = :black, strokewidth = 0.75, markersize = 9)
    scatter!(ax, s[lp, 1] .* 1e3, s[lp, 2] .* 1e3, color = :black, markersize = 9)
end

save("fig3b_equilibria.png", fig, px_per_unit = 3)
