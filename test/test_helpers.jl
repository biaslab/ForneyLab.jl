#####################
# Unit tests
#####################

facts("Helper function unit tests") do
    context("ensureMatrix should convert an array with one element to a matrix type") do
        @fact ForneyLab.ensureMatrix([1.0, 2.0]) --> Diagonal([1.0, 2.0])
        @fact ForneyLab.ensureMatrix(Diagonal([1.0, 2.0])) --> Diagonal([1.0, 2.0])
        @fact ForneyLab.ensureMatrix(eye(2)) --> eye(2)
        @fact ForneyLab.ensureMatrix(1.0) --> eye(1)
        @fact ForneyLab.ensureMatrix(nothing) --> nothing
    end

    context("isApproxEqual should work for scalars, vectors and matrices") do
        @fact isApproxEqual(1.0, 1.0+1e-15) --> true
        @fact isApproxEqual(1.0, 1.0+1e-9) --> false
        @fact isApproxEqual([1.0, 1.0], [1.0, 1.0]+1e-15) --> true
        @fact isApproxEqual([1.0, 1.0], [1.0, 1.0]+1e-9) --> false
        @fact isApproxEqual(eye(3,3), eye(3,3)+1e-15) --> true
        @fact isApproxEqual(eye(3,3), eye(3,3)+1e-9) --> false
    end

    context("isRoundedPosDef should check positive difiniteness with robustness for numerical errors") do
        K = inv([4.0 3.0 2.0;
                 3.0 4.0 3.0;
                 2.0 3.0 4.0])
        @fact isposdef(K) --> false
        @fact ForneyLab.isRoundedPosDef(K) --> true
        @fact ForneyLab.isRoundedPosDef(Diagonal([1.0, 2.0, 3.0])) --> true
    end

    context("isValid should check validity of scalars, vectors and matrices") do
        @fact isValid(1.0) --> true
        @fact isValid(NaN) --> false
        @fact isValid([1.0, 2.0]) --> true
        @fact isValid([1.0, NaN]) --> true
        @fact isValid([NaN, 2.0]) --> false
        @fact isValid(eye(2)) --> true
        @fact isValid([1.0 NaN; 2.0 3.0]) --> true
        @fact isValid([NaN 1.0; 2.0 3.0]) --> false
        @fact isValid(Diagonal([1.0, 2.0])) --> true
        @fact isValid(Diagonal([NaN, 2.0])) --> false
    end

    context("invalidate! should invalidate vectors and matrices") do
        A = [1.0, 2.0]
        @fact isValid(invalidate!(A)) --> false
        @fact isValid(A) --> false
        A = [1.0 2.0; 3.0 4.0]
        @fact isValid(invalidate!(A)) --> false
        @fact isValid(A) --> false
        A = Diagonal([1.0, 2.0])
        @fact isValid(invalidate!(A)) --> false
        @fact isValid(A) --> false
    end

    context("cholinv() should perform a matrix inversion on a positive semidefinite matrix") do
        A = [2.0 1.0; 1.0 2.0]
        @fact cholinv(A) --> roughly(inv(A))
        A = Diagonal([2.0, 3.0])
        @fact cholinv(A) --> inv(A)
    end

    context("diageye() should return a diagonal eye matrix") do
        M = diageye(3)
        @fact typeof(M) --> Diagonal{Float64}
        @fact M --> Diagonal(ones(3))
    end

    context("The following arithmetic operations on diagonal matrix should yield diagonal matrices") do
        D = Diagonal([1.0, 2.0])
        M = [2.0 1.0; 1.0 2.0]
        @fact D.*D --> Diagonal([1.0, 4.0])
        @fact D.*M --> Diagonal([2.0, 4.0])
        @fact M.*D --> Diagonal([2.0, 4.0])
        @fact D^0.5 --> Diagonal([1.0, sqrt(2.0)])
        @fact sqrt(D) --> Diagonal([1.0, sqrt(2.0)])
    end

    context("ensureMessage! should assign a message to an interface if there is none") do
        node = TerminalNode(Gaussian(m=5.0, V=1.0))
        @fact node.i[:out].message --> nothing
        @fact ForneyLab.ensureMessage!(node.i[:out], Gaussian) --> Message(vague(Gaussian))
        @fact node.i[:out].message --> Message(vague(Gaussian))
        @fact ForneyLab.ensureMessage!(node.i[:out], MvGaussian{2}) --> Message(vague(MvGaussian{2}))
        @fact node.i[:out].message --> Message(vague(MvGaussian{2}))
        @fact typeof(ForneyLab.ensureMessage!(node.i[:out], Delta{Bool}).payload) --> Delta{Bool}
        @fact typeof(ForneyLab.ensureMessage!(node.i[:out], Delta{Float64}).payload) --> Delta{Float64}
        @fact typeof(ForneyLab.ensureMessage!(node.i[:out], MvDelta{Float64, 2}).payload) --> MvDelta{Float64, 2}
        @fact typeof(ForneyLab.ensureMessage!(node.i[:out], MatrixDelta{Float64, 2, 3}).payload) --> MatrixDelta{Float64, 2, 3}
    end

    context("ensureMessage! should return the present message is there is one and the type matches") do
        node = TerminalNode(Gaussian(m=5.0, V=1.0))
        node.i[:out].message = Message(Gaussian(m=5.0, V=1.0))
        @fact ForneyLab.ensureMessage!(node.i[:out], Gaussian) --> Message(Gaussian(m=5.0, V=1.0))
        @fact node.i[:out].message --> Message(Gaussian(m=5.0, V=1.0))
    end

    context("ensureMessage! should assign a message to an interface if the type does not match") do
        node = TerminalNode(Gaussian(m=5.0, V=1.0))
        node.i[:out].message = Message(Gaussian(m=5.0, V=1.0))
        @fact typeof(ForneyLab.ensureMessage!(node.i[:out], Delta{Float64})) --> Message{Delta{Float64}}
    end

    context("truncate() should truncate a string to a specified length") do
        @fact ForneyLab.truncate("spsbrats", 9) --> "spsbrats"
        @fact ForneyLab.truncate("spsbrats", 7) --> "spsb..."
    end

    context("pad() should pad a string with spaces to a specified length") do
        @fact ForneyLab.pad("spsbrats", 9) --> "spsbrats "
        @fact ForneyLab.pad("spsbrats", 7) --> "spsb..."
    end

    context("format() should format distributions") do
        @fact format(0.0000001) --> "1.00e-07"
        @fact format(0.0) --> "0.00"
        @fact format(true) --> "true"
        @fact format(:a) --> "a"
        @fact format([7.345345456456456464564645645645, 0.00005345, -0.000145, -108.0]) --> "[7.35, 5.34e-05, -1.45e-04, -1.08e+02]"
        @fact format([7.345345456456456464564645645645 0.00005345; -0.000145 -108.0]) --> "[[7.35, 5.34e-05][-1.45e-04, -1.08e+02]]"
        @fact format(Delta()) --> "δ(m=1.00)"
        @fact format(Gamma()) --> "Gam(a=1.00, b=1.00)"
        @fact format(Gaussian()) --> "N(m=0.00, V=1.00)"
        @fact format(InverseGamma()) --> "Ig(a=3.00, b=2.00)"
        @fact format(StudentsT()) --> "St(μ=0.00, λ=1.00, ν=1.00e+12)"
        @fact format(Beta()) --> "Bet(a=1.00, b=1.00)"

        @fact format(MvDelta()) --> "δ(m=[1.00])"
        @fact format(MvGaussian()) --> "N(m=[0.00], V=[[1.00]])"
        @fact format(MvGaussian(m=zeros(2), V=diageye(2))) --> "N(m=[0.00, 0.00], V=diag[1.00, 1.00])"
        @fact format(NormalGamma()) --> "Ng(m=0.00, β=1.00, a=1.00, b=1.00)"
    end

    context("expand() should expand a dictionary indexed by arrays into separate entries") do
        d = Dict([1 2; 3 4] => 'a', [5 6] => 'b', [7] => 'c', 8 => 'd')
        @fact ForneyLab.expand(d) --> Dict(1 => 'a', 2 => 'a', 3 => 'a', 4 => 'a', 5 => 'b', 6 => 'b', 7 => 'c', 8 => 'd')
    end
end
