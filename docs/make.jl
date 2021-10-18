using Documenter, GrowthMaps, Literate

# Is a file markdown
ismd(f) = splitext(f)[2] == ".md" 

# Generate a title from a file name
# "establishment_models.md" -> "Establishment Models" 
function title(filename) 
    stem, ext = splitext(filename)
    join(map(titlecase, split(stem, "_")), " ")
end

# Generate title => filename pairs for markdown files in `path`
pages(path) = [title(fn) => fn for fn in readdir(path) if ismd(fn)]

# Define directories
docsdir = realpath(joinpath(dirname(pathof(GrowthMaps)), "../docs"))
builddir = joinpath(docsdir, "build")
lit = joinpath(docsdir, "lit")
output = joinpath(docsdir, "src")
jupyterdir = joinpath(builddir, "jupyter")
plutodir = joinpath(builddir, "pluto")

# Generate examples notebookfe/pluto
mkpath(jupyterdir)
mkpath(plutodir)

# Build markdown and notebook docs
for (root, _, files) in walkdir(lit), file in files
    splitext(file)[2] == ".jl" || continue
    in_path = joinpath(root, file)
    out_path = splitdir(replace(in_path, lit => output))[1]
    Literate.markdown(in_path, out_path; documenter=true)
    Literate.notebook(in_path, plutodir; execute=false, flavor=:pluto)
    Literate.notebook(in_path, jupyterdir; execute=false, flavor=:jupyter)
end

# Copy Project and Manifest toml files so that
# e.g. binder can use them to build the notebook project.
# cp(joinpath(basedir, "Project.toml"), joinpath(builddir, "Project.toml"); force=true)
# cp(joinpath(basedir, "Manifest.toml"), joinpath(builddir, "Manifest.toml"); force=true)

# cd(joinpath(docsdir, "build"))

# Generate HTML docs
makedocs(
    modules = [GrowthMaps],
    sitename = "GrowthMaps.jl",
    pages = [
        "Home" => "index.md",
        pages(output)...
    ],
    clean = false,
    strict = true,
)

deploydocs(
    repo = "github.com/cesaraustralia/GrowthMaps.jl.git",
)

