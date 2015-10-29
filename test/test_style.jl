#####################
# Source style tests
#####################

function findAllJuliaFiles(path::AbstractString, files::Array{AbstractString,1}=Array(AbstractString,0)) 
    dircontents = readdir(path)
    for content in dircontents
        content_path = joinpath(path, content)
        if isdir(content_path)
            findAllJuliaFiles(content_path, files)
        elseif content_path[end-2:end] == ".jl"
            push!(files, content_path)
        end
    end

    return files
end

facts("Source style tests") do
    context("Source files should not contain tabs as indents, but 4 spaces") do
        forneylab_root = realpath(joinpath(dirname(@__FILE__()),".."))
        paths = [joinpath(forneylab_root, "src"), joinpath(forneylab_root, "test")]
        for path in paths
            julia_files = findAllJuliaFiles(path)
            for julia_file in julia_files
                context("Source file should not contain tabs: $julia_file") do
                    @fact search(open(readall, julia_file), "\t") --> 0:-1
                end
            end
        end
    end
end
