using Documenter, GrowthMaps, Weave, IJulia

basedir = dirname(@__FILE__)

example = joinpath(basedir, "src/example.jmd")

mdpath = joinpath(basedir, "src/example.md")
notebookdir = joinpath(basedir, "build/notebook")
pdfdir = joinpath(basedir, "build/pdf")

mkpath(joinpath(basedir, "build/assets"))
mkpath.((notebookdir, pdfdir))

# Generate examples pdf
# weave(example, out_path=pdfdir, doctype="pandoc2pdf")

# Generate examples markdown and images
weave(example, out_path=mdpath, doctype="github")

# Generate examples notebook
# convert_doc(example, joinpath(notebookdir, "example.ipynb"))

# Generate HTML docs
makedocs(
    modules = [GrowthMaps],
    sitename = "GrowthMaps.jl",
    pages = [
        "Home" => "index.md",
        "Examples" => "example.md",
    ],
    clean = false,
)

deploydocs(
    repo = "github.com/cesaraustralia/GrowthMaps.jl.git",
)
