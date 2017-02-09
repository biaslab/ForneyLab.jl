module HelpersTest

using Base.Test
import ForneyLab: ensureMatrix, isApproxEqual, isRoundedPosDef, truncate, pad, huge, tiny, format, isValid, invalidate!, cholinv, diageye, *, .*, ^
 
@testset "Helpers" begin
    @testset "ensureMatrix" begin
        # should convert input argument to a Matrix or Void
        @test ensureMatrix([1.0, 2.0]) == Diagonal([1.0, 2.0])
        @test ensureMatrix(Diagonal([1.0, 2.0])) == Diagonal([1.0, 2.0])
        @test ensureMatrix(eye(2)) == eye(2)
        @test ensureMatrix(1.0) == eye(1)
        @test ensureMatrix(nothing) == nothing
    end

    @testset "isApproxEqual" begin
        # should work for scalars, vectors and matrices
        @test isApproxEqual(1.0, 1.0+1e-15) == true
        @test isApproxEqual(1.0, 1.0+1e-9) == false
        @test isApproxEqual([1.0, 1.0], [1.0, 1.0]+1e-15) == true
        @test isApproxEqual([1.0, 1.0], [1.0, 1.0]+1e-9) == false
        @test isApproxEqual(eye(3,3), eye(3,3)+1e-15) == true
        @test isApproxEqual(eye(3,3), eye(3,3)+1e-9) == false
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

    @testset "isValid" begin
        # should check validity of scalars, vectors and matrices
        @test isValid(1.0) == true
        @test isValid(NaN) == false
        @test isValid([1.0, 2.0]) == true
        @test isValid([1.0, NaN]) == true
        @test isValid([NaN, 2.0]) == false
        @test isValid(eye(2)) == true
        @test isValid([1.0 NaN; 2.0 3.0]) == true
        @test isValid([NaN 1.0; 2.0 3.0]) == false
        @test isValid(Diagonal([1.0, 2.0])) == true
        @test isValid(Diagonal([NaN, 2.0])) == false
    end

    @testset "invalidate!" begin
        # should invalidate vectors and matrices
        A = [1.0, 2.0]
        @test isValid(invalidate!(A)) == false
        @test isValid(A) == false
        A = [1.0 2.0; 3.0 4.0]
        @test isValid(invalidate!(A)) == false
        @test isValid(A) == false
        A = Diagonal([1.0, 2.0])
        @test isValid(invalidate!(A)) == false
        @test isValid(A) == false
    end

    @testset "cholinv" begin
        # should perform a matrix inversion on a positive (semi)definite matrix
        A = [2.0 1.0; 1.0 2.0]
        @test_approx_eq cholinv(A) inv(A)
        A = Diagonal([2.0, 3.0])
        @test cholinv(A) == inv(A)
    end

    @testset "diageye" begin
        # should be shorthand for Diagonal(eye(M))
        M = diageye(3)
        @test typeof(M) == Diagonal{Float64}
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

    @testset "truncate" begin
        # should truncate a string to a specified length
        @test truncate("spsbrats", 9) == "spsbrats"
        @test truncate("spsbrats", 7) == "spsb..."
    end

    @testset "pad" begin
        # should pad a string with spaces to a specified length
        @test pad("spsbrats", 9) == "spsbrats "
        @test pad("spsbrats", 7) == "spsb..."
    end

    @testset "format" begin
        " should return formatted strings for various input types"
        @test format(0.0000001) == "1.00e-07"
        @test format(0.0) == "0.00"
        @test format(true) == "true"
        @test format(:a) == "a"
        @test format([7.345345456456456464564645645645, 0.00005345, -0.000145, -108.0]) == "[7.35, 5.34e-05, -1.45e-04, -1.08e+02]"
        @test format([7.345345456456456464564645645645 0.00005345; -0.000145 -108.0]) == "[[7.35, 5.34e-05][-1.45e-04, -1.08e+02]]"
    end
end

end # module