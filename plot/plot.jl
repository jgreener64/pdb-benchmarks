# Plot benchmark results

using CSV
using Gadfly
using Cairo
using Fontconfig

df = CSV.read("C:\\Users\\Joe\\software\\pdb-benchmarks\\benchmarks.csv")

theme = Theme(
    background_color="white",
    panel_stroke="black",
    major_label_color="black",
    minor_label_color="black",
)

p = Gadfly.with_theme(theme) do
    plot(df,
        x=:Package,
        y=:Runtime,
        color=:Benchmark,
        Scale.y_log10,
        Guide.xlabel(nothing),
        Guide.ylabel("Runtime / s"),
        Geom.point,
        Geom.line,
    )
end

draw(PNG("C:\\Users\\Joe\\software\\pdb-benchmarks\\plot\\plot.png", dpi=300), p)
