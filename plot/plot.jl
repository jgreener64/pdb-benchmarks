# Plot benchmark results

using DataFrames
using Gadfly

df = readtable("benchmarks.csv")

p = plot(df,
    x=:package,
    y=:time,
    color=:benchmark,
    Scale.y_log10,
    Guide.xlabel(nothing),
    Guide.ylabel("time / s"),
    Geom.point,
    Geom.line,
)

draw(PNG("plot/plot.png", 9inch, 6inch), p)
