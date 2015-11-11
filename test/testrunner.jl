using ForneyLab
using FactCheck

function collectTestFiles(path::AbstractString, files::Set{AbstractString}, exclude::Set{AbstractString}, regex=r".*")
    path = abspath(path)
    if isdir(path)
        object_names = readdir(path)
        dirs = filter(x -> isdir(joinpath(path, x)), object_names)
        for dir_name in dirs
            files = union(files, collectTestFiles(joinpath(path, dir_name), files, exclude, regex))
        end
        dir_files = Set{AbstractString}([joinpath(path, x) for x in filter(x -> (isfile(joinpath(path, x))
                                                                                 && contains(x, "test_")
                                                                                 && ismatch(regex, x)
                                                                                 && !(x in exclude)),
                                                                           object_names)])
        return union(files, dir_files)
    elseif isfile(path)
        return Set{AbstractString}(path)
    end
end

function collectTestFilesInPaths(paths, files::Set{AbstractString}, exclude::Set{AbstractString}, regex=r".*")
    for path in paths
        test_files_in_path = collectTestFiles(path, files, exclude, regex)
        if test_files_in_path != nothing
            files = union(files, test_files_in_path)
        end
    end
    return files
end

abstract TestRunner

function setUp(runner::TestRunner)
end

function tearDown(runner::TestRunner)
end

type ForneyLabTestRunner<:TestRunner
    files::Set{AbstractString}
end

function setUp(runner::ForneyLabTestRunner)
    current_graph = FactorGraph()
    current_algorithm = nothing
end

function tearDown(runner::ForneyLabTestRunner)
    current_graph = nothing
    current_algorithm = nothing
end

function runTests(runner::ForneyLabTestRunner)
    for test_file in runner.files
        setUp(runner)
        include(test_file)
        tearDown(runner)
    end
    exitstatus()
end

using ArgParse
include("integration_helpers.jl")

function main(args)
    s = ArgParseSettings(description = "Test runner settings:")

    @add_arg_table s begin
        "--paths", "-p"
        nargs = '*'
        default = Any["./"]
        arg_type = AbstractString
        help = "Paths to look for test files in"

        "--regex", "-r"
        default = ".*"
        arg_type = AbstractString
        help = "Regular expression for filenames"

        "--exclude", "-e"
        nargs = '*'
        default = Any["test_forneylab.jl", "test_demos.jl"]
        arg_type = AbstractString
        help = "Files to exclude (for example, test_forneylab.jl)"
    end

    parsed_args = parse_args(s)

    regex = parsed_args["regex"]
    paths = parsed_args["paths"]
    exclude = Set{AbstractString}(parsed_args["exclude"])

    test_files = collectTestFilesInPaths(paths, Set{AbstractString}(), exclude, Regex(regex))

    fl_test_runner = ForneyLabTestRunner(test_files)
    runTests(fl_test_runner)
end

main(ARGS)
