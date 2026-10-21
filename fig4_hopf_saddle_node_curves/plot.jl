# Figure 4: one Hopf curve H (solid supercritical, dashed subcritical) and one saddle-node curve SN
# of system (10) in the (mu, eta)-plane
# Run: julia plot.jl (after python3 run.py); writes fig4_hopf_saddle_node_curves.png
using CairoMakie, DelimitedFiles

fig = Figure(size = (500, 400))
ax = Axis(fig[1, 1], xlabel = "μ (×10⁻³)", ylabel = "η", limits = (-5, 0, -6.2, 0))

for name in ["lp_1", "lp_2"]
    d = readdlm(name * ".dat")
    lines!(ax, d[:, 1] .* 1e3, d[:, 2], color = "#231f20", linewidth = 2)
end

for name in ["hb_1", "hb_2"]
    d = readdlm(name * ".dat")
    mu, eta, l1 = d[:, 1] .* 1e3, d[:, 2], d[:, 3]
    lines!(ax, mu, ifelse.(l1 .< 0, eta, NaN), color = "#69aa44", linewidth = 2)                      # supercritical
    lines!(ax, mu, ifelse.(l1 .> 0, eta, NaN), color = "#69aa44", linewidth = 2, linestyle = :dash)  # subcritical
end

save("fig4_hopf_saddle_node_curves.png", fig, px_per_unit = 3)
