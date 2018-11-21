using Documenter, ForneyLab

makedocs(modules = [ForneyLab],
    clean = true,
    sitename = "ForneyLab.jl",
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting-started/introduction.md",
        "API" => "API/algorithms.md",
        "Internals" => "internals/algorithms1.md"
            # ,
            # "Compositing" => "man/compositing.md",
            # "Backends" => "man/backends.md",
            # "Themes" => "man/themes.md",
        # "Gallery" => Any[
        #     "Geometries" => "gallery/geometries.md",
        #     "Guides" => "gallery/guides.md",
        #     "Statistics" => "gallery/statistics.md",
        #     "Coordinates" => "gallery/coordinates.md",
        #     "Scales" => "gallery/scales.md",
        #     "Shapes" => "gallery/shapes.md",
        # ],
        # "Library" => Any[
        #     "Gadfly" => "lib/gadfly.md",
        #     "Geometries" => "lib/geometries.md",
        #     "Guides" => "lib/guides.md",
        #     "Statistics" => "lib/statistics.md",
        #     "Coordinates" => "lib/coordinates.md",
        #     "Scales" => "lib/scales.md",
        #     "Shapes" => "lib/shapes.md",
        # ],
        # "Development" => Any[
        #     "Rendering Pipeline" => "dev/pipeline.md",
        #     "Regression Testing" => "dev/regression.md",
        #     "Relationship with Compose.jl" => "dev/compose.md",
        # ]
    ]
)

deploydocs(
    repo   = "github.com/GiovineItalia/Gadfly.jl.git",
)
