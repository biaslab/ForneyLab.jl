# This file validates that all demos run without errors

orig_pwd = pwd()
demo_dir = normpath(joinpath(dirname(@__FILE__()), "..", "demo"))
cd(demo_dir)

# Convert all iJulia notebooks in the demo directory to regular .jl files
run(`ipython nbconvert --profile julia --to python *.ipynb`)
map(file -> mv(file, file[1:end-3]*".jl"), filter!(r"\.py$", readdir(demo_dir))) # Rename files from .py to .jl

# Run the demos
try
    for demo_file in filter!(r"\.jl$", readdir(demo_dir))
        println("\n>> Checking demo: $(basename(demo_file)[1:end-3])\n")
        try
            # Run the demo
            readall(`julia $(demo_file)`)
            # Delete julia script if the demo was executed successfully
            rm(demo_file)
        catch exception
            error("The demo generated an error. Have a look at the .jl file that generated the error.")
        end
    end
finally
    # Restore original working directory
    cd(orig_pwd)
end