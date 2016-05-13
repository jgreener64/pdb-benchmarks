# Plot benchmark results

using DataFrames
using Gadfly

df = readtable("benchmarks.csv",
    names=[:package, :benchmark, :time],
    header=false,
)

p = plot(df,
    x=:package,
    y=:time,
    color=:benchmark,
    Scale.y_log10,
    Guide.xlabel(nothing),
    Guide.ylabel("time / s"),
)

draw(SVG("plot/plot.svg", 5inch, 3inch), p)
