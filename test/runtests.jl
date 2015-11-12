# Run the tests for ForneyLab
# This file is automatically executed if you call Pkg.test("ForneyLab")

include("testrunner.jl")

test_files = collectTestFilesInPaths(["./"], Set{AbstractString}(), Set{AbstractString}(["test_demos.jl"]), Regex(".*"))
fl_test_runner = ForneyLabTestRunner(test_files)
runTests(fl_test_runner)
