module HelpersTest

using Test
import ForneyLab: ensureMatrix, isApproxEqual, isRoundedPosDef, huge, tiny, format, leaftypes, isValid, invalidate!, inv, diageye, *, ^
import LinearAlgebra: Diagonal, isposdef, I, Hermitian
import PDMats: AbstractPDMat, PDMat, PDiagMat

@testset "Helpers" begin
    @testset "ensureMatrix" begin
        # should convert input argument to a Matrix or Nothing
        @test ensureMatrix([1.0, 2.0]) == Diagonal([1.0, 2.0])
        @test ensureMatrix(Diagonal([1.0, 2.0])) == Diagonal([1.0, 2.0])
        @test ensureMatrix(Matrix(1.0I,2,2)) == Matrix(1.0I,2,2)
        @test ensureMatrix(1.0) == Matrix(1.0I,1,1)
        @test ensureMatrix(nothing) == nothing
    end

    @testset "isApproxEqual" begin
        # should work for scalars, vectors and matrices
        @test isApproxEqual(1.0, 1.0+1e-15) == true
        @test isApproxEqual(1.0, 1.0+1e-9) == false
        @test isApproxEqual([1.0, 1.0], [1.0, 1.0].+1e-15) == true
        @test isApproxEqual([1.0, 1.0], [1.0, 1.0].+1e-9) == false
        @test isApproxEqual(Matrix(1.0I,3,3), Matrix(1.0I,3,3).+1e-15) == true
        @test isApproxEqual(Matrix(1.0I,3,3), Matrix(1.0I,3,3).+1e-9) == false
        @test isApproxEqual(PDMat(Matrix(1.0I,3,3)), PDMat(Matrix(1.0I,3,3).+1e-15)) == true
        @test isApproxEqual(PDMat(Matrix(1.0I,3,3)), PDMat(Matrix(1.0I,3,3).+1e-9)) == false
        @test isApproxEqual(PDiagMat([3.0,3.0]), PDiagMat([3.0,3.0].+1e-15)) == true
        @test isApproxEqual(PDiagMat([3.0,3.0]), PDiagMat([3.0,3.0].+1e-9)) == false
    end

    @testset "isRoundedPosDef" begin
        # should check positive definiteness with robustness for numerical errors
        K = inv([4.0 3.0 2.0;
                 3.0 4.0 3.0;
                 2.0 3.0 4.0])
        @test isposdef(K) == false
        @test isRoundedPosDef(K) == true
        @test isRoundedPosDef(Diagonal([1.0, 2.0, 3.0])) == true
    end

    @testset "inv" begin
        # Matrix inversion on a positive (semi)definite matrix (Array)
        A = Diagonal([2.0, 3.0])
        @test inv(A) == Diagonal([1/2.0, 1/3.0])
        B = [2.0 1.0; 1.0 1.0]
        @test isApproxEqual(inv(B), [1.0 -1.0; -1.0 2.0])
        # Matrix inversion on a positive (semi)definite matrix (PDMats objects)
        D = PDiagMat([1.0, 2.0, 4.0])
        @test isApproxEqual(inv(D), Diagonal([1.0, 0.5, 0.25]))
        E = PDMat([3.0 1.0; 1.0 2.0])
        @test isApproxEqual(inv(E), Array{Float64,2}([0.4 -0.2; -0.2 0.6]))
    end

    @testset "diageye" begin
        # should be shorthand for Diagonal(eye(M))
        M = diageye(3)
        @test typeof(M) == Diagonal{Float64, Array{Float64,1}}
        @test M == Diagonal(ones(3))
    end

    @testset "operators for diagonal matrices" begin
        # should yield Diagonal
        D = Diagonal([1.0, 2.0])
        M = [2.0 1.0; 1.0 2.0]
        @test D.*D == Diagonal([1.0, 4.0])
        @test D.*M == Diagonal([2.0, 4.0])
        @test M.*D == Diagonal([2.0, 4.0])
        @test D^0.5 == Diagonal([1.0, sqrt(2.0)])
        @test sqrt(D) == Diagonal([1.0, sqrt(2.0)])
    end

    @testset "format" begin
        # should return formatted strings for various input types
        @test format(0.0000001) == "1.00e-07"
        @test format(0.0) == "0.00"
        @test format(true) == "true"
        @test format(:a) == "a"
        @test format([7.345345456456456464564645645645, 0.00005345, -0.000145, -108.0]) == "[7.35, 5.34e-05, -1.45e-04, -1.08e+02]"
        # @test format([7.345345456456456464564645645645 0.00005345; -0.000145 -108.0]) == "[[7.35, 5.34e-05], [-1.45e-04, -1.08e+02]]"
    end

    @testset "leaftypes" begin
        # should return all subtypes that are leafs on the type tree
        @test Set(leaftypes(Real)) == Set([BigInt, Bool, UInt128, UInt16, UInt32, UInt64, UInt8, Int128, Int16, Int32, Int64, Int8, BigFloat, Float16, Float32, Float64])
    end
end

end # module
