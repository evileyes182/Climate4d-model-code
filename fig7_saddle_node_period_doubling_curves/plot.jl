# Figure 7: saddle-node curve S and period-doubling curve PD of periodic orbits of system (10)
# in the (mu, eta)-plane
using CairoMakie, DelimitedFiles

fig = Figure(size = (500, 500))
ax = Axis(fig[1, 1], xlabel = "μ (×10⁻³)", ylabel = "η", limits = (-2.17, -2.045, -4.29, -3.68))

for (name, col) in [("s_1", "#1384d4"), ("s_2", "#1384d4"), ("pd_1", "#bb5bbd"), ("pd_2", "#bb5bbd")]
    d = readdlm(name * ".dat")
    lines!(ax, d[:, 1] .* 1e3, d[:, 2], color = col, linewidth = 2)
end

save("fig7_saddle_node_period_doubling_curves.png", fig, px_per_unit = 3)
