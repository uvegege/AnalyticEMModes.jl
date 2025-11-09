Pkg.add("SciPy")
using SciPy


begin
    xse = range(0.0, 2π, 300)
    m = 1
    q = 3.4
    y_py = SciPy.special.mathieu_sem(m, q, xse*180/pi)
    y_ju = se_m(m, q, xse)
    y_py, yp_py = y_py[1], y_py[2]
    y_ju, yp_ju = getindex.(y_ju, 1), getindex.(y_ju, 2)

    vpy = maximum(abs.(y_py))
    vju = maximum(abs.(y_ju))
    y_py, yp_py = y_py/vpy, yp_py/vpy 
    y_ju, yp_ju = y_ju/vju, yp_ju/vju
    fig = Figure()
    ax1 = Axis(fig[1,1], subtitle = "se_m")
    lines!(ax1, xse, y_py)
    lines!(ax1, xse, y_ju)
    ax2 = Axis(fig[1,2], subtitle = "se'_m")
    lines!(ax2, xse, yp_py)
    lines!(ax2, xse, yp_ju)
    @show norm(y_py - y_ju)
    @show norm(yp_py - yp_ju)
    fig
end


begin
    xse = range(0.0, 2π, 300)
    m = 7
    q = 17.4
    y_py = SciPy.special.mathieu_cem(m, q, xse*180/pi)
    y_ju = ce_m(m, q, xse)
    y_py, yp_py = y_py[1], y_py[2]
    y_ju, yp_ju = getindex.(y_ju, 1), getindex.(y_ju, 2)

    vpy = maximum(abs.(y_py))
    vju = maximum(abs.(y_ju))
    y_py, yp_py = y_py/vpy, yp_py/vpy 
    y_ju, yp_ju = y_ju/vju, yp_ju/vju
    fig = Figure()
    ax1 = Axis(fig[1,1], subtitle = "ce_m")
    lines!(ax1, xse, y_py)
    lines!(ax1, xse, y_ju)
    ax2 = Axis(fig[1,2], subtitle = "ce'_m")
    lines!(ax2, xse, yp_py)
    lines!(ax2, xse, yp_ju)
    @show norm(y_py - y_ju)
    @show norm(yp_py - yp_ju)
    fig
end


begin
    ξs = range(0.0, 2.0, 300)
    m = 1
    q = 32.4
    y_py = SciPy.special.mathieu_modcem1(m, q, ξs)
    y_ju = Ce_m(m, q, ξs)
    y_py, yp_py = y_py[1], y_py[2]
    y_ju, yp_ju = getindex.(y_ju, 1), getindex.(y_ju, 2)

    vpy = maximum(abs.(y_py))
    vju = maximum(abs.(y_ju))
    y_py, yp_py = y_py/vpy, yp_py/vpy 
    y_ju, yp_ju = y_ju/vju, yp_ju/vju
    fig = Figure()
    ax1 = Axis(fig[1,1], subtitle = "Ce_m")
    lines!(ax1, ξs, y_py)
    lines!(ax1, ξs, y_ju)
    ax2 = Axis(fig[1,2], subtitle = "Ce'_m")
    lines!(ax2, ξs, yp_py)
    lines!(ax2, ξs, yp_ju)
    @show norm(y_py - y_ju)
    @show norm(yp_py - yp_ju)
    fig
end


begin
    ξs = range(0.0, 2.0, 300)
    m = 2
    q = 121.4
    y_py = SciPy.special.mathieu_modsem1(m, q, ξs)
    y_ju = Se_m(m, q, ξs)
    y_py, yp_py = y_py[1], y_py[2]
    y_ju, yp_ju = getindex.(y_ju, 1), getindex.(y_ju, 2)

    vpy = maximum(abs.(y_py))
    vju = maximum(abs.(y_ju))
    y_py, yp_py = y_py/vpy, yp_py/vpy 
    y_ju, yp_ju = y_ju/vju, yp_ju/vju
    fig = Figure()
    ax1 = Axis(fig[1,1], subtitle = "Se_m")
    lines!(ax1, ξs, y_py)
    lines!(ax1, ξs, y_ju)
    ax2 = Axis(fig[1,2], subtitle = "Se'_m")
    lines!(ax2, ξs, yp_py)
    lines!(ax2, ξs, yp_ju)
    @show norm(y_py - y_ju)
    @show norm(yp_py - yp_ju)
    fig
end


begin
    qs = range(0.1, 450.0, 300)
    ξi = 0.76
    m = 3
    y_py = SciPy.special.mathieu_modcem1(m, qs, ξi)
    y_ju = Ce_m.(m, qs, ξi)
    y_py, yp_py = y_py[1], y_py[2]
    y_ju, yp_ju = getindex.(y_ju, 1), getindex.(y_ju, 2)

    vpy = maximum(abs.(y_py))
    vju = maximum(abs.(y_ju))
    y_py, yp_py = y_py/vpy, yp_py/vpy 
    y_ju, yp_ju = y_ju/vju, yp_ju/vju
    fig = Figure()
    ax1 = Axis(fig[1,1], subtitle = "Ce_m")
    lines!(ax1, qs, y_py)
    lines!(ax1, qs, y_ju)
    ax2 = Axis(fig[1,2], subtitle = "Ce'_m")
    lines!(ax2, qs, yp_py)
    lines!(ax2, qs, yp_ju)
    @show norm(y_py - y_ju)
    @show norm(yp_py - yp_ju)
    fig
end


begin
    qs = range(0.1, 450.0, 300)
    ξi = 0.6
    m = 2
    y_py = SciPy.special.mathieu_modsem1(m, qs, ξi)
    y_ju = Se_m.(m, qs, ξi)
    y_py, yp_py = y_py[1], y_py[2]
    y_ju, yp_ju = getindex.(y_ju, 1), getindex.(y_ju, 2)

    vpy = maximum(abs.(y_py))
    vju = maximum(abs.(y_ju))
    y_py, yp_py = y_py/vpy, yp_py/vpy 
    y_ju, yp_ju = y_ju/vju, yp_ju/vju
    fig = Figure()
    ax1 = Axis(fig[1,1], subtitle = "Se_m")
    lines!(ax1, qs, y_py)
    lines!(ax1, qs, y_ju)
    ax2 = Axis(fig[1,2], subtitle = "Se'_m")
    lines!(ax2, qs, yp_py)
    lines!(ax2, qs, yp_ju)
    @show norm(y_py - y_ju)
    @show norm(yp_py - yp_ju)
    fig
end


fig = Figure()
ax1 = Axis(fig[1,1], subtitle = "Ce_m")
qs = range(0.1, 450.0, 3000)
ξi = 0.6931471805599453
y_ju = getindex.(Ce_m.(1, qs, ξi),1)
lines!(ax1, qs, y_ju)
