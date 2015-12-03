# This file validates that all demos run without errors
# Simply run this file to check the demos

orig_pwd = pwd()
demo_dir = normpath(joinpath(dirname(@__FILE__()), "..", "demo"))
forneylab_main = joinpath(normpath(joinpath(dirname(@__FILE__()), "..", "src")), "ForneyLab.jl")
cd(demo_dir)

# Convert all iJulia notebooks in the demo directory to regular .jl files
run(`ipython nbconvert --profile julia --to python *.ipynb`)
for file in  filter!(r"\.py$", readdir(demo_dir))
    # Generate Julia script from .py file
    open(file[1:end-3]*".jl", "w") do f
        # Insert a line to include ForneyLab
        write(f, "include(\"$(forneylab_main)\")\n")
        # Write the actual Julia script
        write(f, readall(file))
    end
    # Remove .py file
    rm(file)
end

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
