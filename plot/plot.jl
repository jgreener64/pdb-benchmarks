# Plot benchmark results

using CSV
using DataFrames
using Gadfly
using Cairo
using Fontconfig

df = CSV.read("benchmarks.csv")

benchmarks = [
    "Parse PDB 1CRN",
    "Parse PDB 1HTQ",
    "Parse mmCIF 1CRN",
    "Parse mmCIF 1HTQ",
    "Parse MMTF 1CRN",
    "Parse MMTF 1HTQ",
    "Count",
    "Distance",
    "Ramachandran",
]
benchind(b) = findfirst(x -> x == b, benchmarks)
df_sorted = sort(df, (order(:Benchmark, by=benchind), :Runtime))

theme = Theme(
    background_color="white",
    panel_stroke="black",
    major_label_color="black",
    minor_label_color="black",
    highlight_width=0mm,
)

p = Gadfly.with_theme(theme) do
    plot(df_sorted,
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

draw(PNG("plot/plot.png", dpi=300), p)
