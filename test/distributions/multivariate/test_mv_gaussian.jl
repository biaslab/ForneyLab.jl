#####################
# Unit tests
#####################

facts("MvGaussian unit tests") do
    context("MvGaussian() should initialize a Gaussian distribution") do
        @fact MvGaussian().V --> reshape([1.0],1,1)
        @fact MvGaussian().m --> [0.0]
        @fact typeof(MvGaussian(m=[0.0], W=eye(1))) --> MvGaussian{1}
        @fact MvGaussian(m=[0.0], W=eye(1)).W --> eye(1)
        @fact MvGaussian(m=[0.0], W=Diagonal([1.0])).W --> Diagonal([1.0])
        @fact MvGaussian(xi=[0.0], W=eye(1)).W --> eye(1)
        @fact typeof(MvGaussian(m=[1.0], V=eye(1)).V) --> Array{Float64, 2} # cast single value to matrix
        @fact typeof(MvGaussian(m=[0.0, 0.0], V=eye(2))) --> MvGaussian{2} # multivariate
        @fact MvGaussian(m=[0.0, 0.0], V=eye(2)).V --> eye(2) # multivariate
        @fact MvGaussian(m=[0.0, 0.0], V=diageye(2)).V --> diageye(2) # multivariate
        @fact_throws MvGaussian(V=eye(1), W=eye(1))
        @fact_throws MvGaussian(m=[0.0], xi=[0.0])
        @fact_throws MvGaussian(xi=[0.0])
        @fact_throws MvGaussian(m=[0.0, 0.0], V=[0. 0.;0. 0.])  # V should be positive definite
        @fact_throws MvGaussian(m=[0.0, 0.0], V=[Inf 0.;0. 1.]) # V should be bounded
        @fact_throws MvGaussian(m=[0.0, 0.0], W=[0. 0.;0. 0.])  # W should be positive definite
        @fact_throws MvGaussian(m=[0.0, 0.0], W=[Inf 0.;0. 1.]) # W should be bounded
        @fact pdf(MvGaussian(m=ones(2), W=eye(2)), [0.0;0.0]) --> roughly(0.05854983152431917, atol=1e-8)
        @fact pdf(MvGaussian(xi=ones(2), V=eye(2)), [0.0;0.0]) --> roughly(0.05854983152431917, atol=1e-8)
    end

    context("vague() should initialize a vague (almost uninformative) Gaussian distribution") do
        dist = vague(MvGaussian{1})
        @fact dist.m --> [0.0]
        @fact dist.V --> Diagonal([huge])

        dist = vague(MvGaussian{2})
        @fact dist.m --> zeros(2)
        @fact dist.V --> Diagonal([huge, huge])
    end

    context("isProper() should indicate whether a Gaussian distribution is proper") do
        @fact isProper(MvGaussian()) --> true
        @fact isProper(MvGaussian(m=[0.0, 0.0], V=diageye(2))) --> true
        @fact isProper(MvGaussian(m=[0.0], V= -1*eye(1))) --> false
        @fact isProper(MvGaussian(m=[0.0], W= -1*eye(1))) --> false
    end

    context("MvGaussian can be sampled") do
        @fact typeof(sample(MvGaussian(m=[1.2, 2.7], V=[2.0 -0.5; -0.5 1.5]))) --> Array{Float64, 1}
        @fact typeof(sample(MvGaussian(m=[1.2, 2.7], V=Diagonal([2.0, 1.5])))) --> Array{Float64, 1}
    end

    context("Underdetermined MvGaussian should be detected by isWellDefined()") do
        @fact isWellDefined(MvGaussian()) --> true
        @fact isWellDefined(MvGaussian(m=[0.0], V=eye(1))) --> true
        @fact isWellDefined(MvGaussian(m=[0.0], W=eye(1))) --> true
        @fact isWellDefined(MvGaussian(m=[0.0, 0.0], W=diageye(2))) --> true
        @fact isWellDefined(MvGaussian(xi=[0.0], V=eye(1))) --> true
        @fact isWellDefined(MvGaussian(xi=[0.0], W=eye(1))) --> true
        @fact isWellDefined(MvGaussian(m=[0.0], xi=[0.0], W=eye(1), V=eye(1))) --> true

        dist = MvGaussian(m=[0.0], V=eye(1))
        invalidate!(dist.m)
        @fact isWellDefined(dist) --> false

        dist = MvGaussian(m=[0.0], V=eye(1))
        invalidate!(dist.V)
        @fact isWellDefined(dist) --> false

        dist = MvGaussian(xi=[0.0], W=eye(1))
        invalidate!(dist.xi)
        @fact isWellDefined(dist) --> false

        dist = MvGaussian(xi=[0.0], W=eye(1))
        invalidate!(dist.W)
        @fact isWellDefined(dist) --> false

        dist = MvGaussian(m=[0.0], V=eye(1), W=eye(1))
        invalidate!(dist.m)
        @fact isWellDefined(dist) --> false

        dist = MvGaussian(m=[0.0], xi=[0.0], V=eye(1))
        invalidate!(dist.V)
        @fact isWellDefined(dist) --> false
    end
    context("Conversions between valid parametrizations of a MvGaussian should be consistent") do
        # Defined as (m,V)
        @fact isConsistent(ensureParameters!(MvGaussian(m=[0.0], V=eye(1)), (:m, :W))) --> true
        @fact isConsistent(ensureParameters!(MvGaussian(m=[0.0], V=eye(1)), (:xi, :V))) --> true
        @fact isConsistent(ensureParameters!(MvGaussian(m=[0.0], V=eye(1)), (:xi, :W))) --> true
        # Defined as (m,W)
        @fact isConsistent(ensureParameters!(MvGaussian(m=[0.0], W=eye(1)), (:m, :V))) --> true
        @fact isConsistent(ensureParameters!(MvGaussian(m=[0.0], W=eye(1)), (:xi, :V))) --> true
        @fact isConsistent(ensureParameters!(MvGaussian(m=[0.0], W=eye(1)), (:xi, :W))) --> true
        # Defined as (xi,V)
        @fact isConsistent(ensureParameters!(MvGaussian(xi=[2.0], V=eye(1)), (:m, :V))) --> true
        @fact isConsistent(ensureParameters!(MvGaussian(xi=[2.0], V=eye(1)), (:m, :W))) --> true
        @fact isConsistent(ensureParameters!(MvGaussian(xi=[2.0], V=eye(1)), (:xi, :W))) --> true
        # Defined as (xi,W)
        @fact isConsistent(ensureParameters!(MvGaussian(xi=[2.0], W=eye(1)), (:m, :V))) --> true
        @fact isConsistent(ensureParameters!(MvGaussian(xi=[2.0], W=eye(1)), (:m, :W))) --> true
        @fact isConsistent(ensureParameters!(MvGaussian(xi=[2.0], W=eye(1)), (:xi, :V))) --> true
    end
    context("mean(d), var(d), cov(d) should return correct results") do
        @fact mean(MvGaussian(m=[1.0, 2.0], V=eye(2))) --> [1.0, 2.0]
        @fact mean(MvGaussian(xi=[1.0, 2.0], V=2.0*eye(2))) --> [2.0, 4.0]
        @fact isProper(MvGaussian(xi=[1.0, 2.0], V=-2.0*eye(2))) --> false
        @fact var(MvGaussian(m=[1.0, 2.0], V=diagm([2.0, 4.0]))) --> [2.0, 4.0]
        @fact var(MvGaussian(m=[1.0, 2.0], W=diagm([2.0, 4.0]))) --> roughly([0.5, 0.25])
        @fact isProper(MvGaussian(m=[1.0, 2.0], W=diagm([-2.0, 4.0]))) --> false
        @fact cov(MvGaussian(m=[1.0, 2.0], V=[2.0 1.0;1.0 2.0])) --> [2.0 1.0;1.0 2.0]
        @fact cov(MvGaussian(m=[1.0, 2.0], W=[2.0 1.0;1.0 2.0])) --> roughly(inv([2.0 1.0;1.0 2.0]))
    end
    context("Inconsistent overdetermined MvGaussian should be detected by isConsistent()") do
        @fact isConsistent(MvGaussian(m=[0.0], xi=[1.0], W=eye(1))) --> false
        @fact isConsistent(MvGaussian(m=[0.0], V=eye(1), W=2*eye(1))) --> false
        @fact isConsistent(MvGaussian(m=[0.0], xi=[1.0], V=eye(1), W=2*eye(1))) --> false
    end

    context("differentialEntropy() should evaluate the differential entropy") do
        @fact ForneyLab.differentialEntropy(MvGaussian()) --> roughly(1.4189385332046727)
    end
end

facts("Marginal calculations for the MvGaussian") do
    context("calculateMarginal!(edge) should give correct result and save the marginal to the edge") do
        initializePairOfTerminalNodes(MvGaussian(m=[0.0], V=eye(1)), MvGaussian(m=[0.0], V=eye(1)))
        edge = n(:t1).i[:out].edge
        n(:t1).i[:out].message = Message(MvGaussian(m=[0.0], V=eye(1)))
        n(:t2).i[:out].message = Message(MvGaussian(m=[0.0], V=eye(1)))
        marginal_dist = calculateMarginal!(edge)
        @fact edge.marginal --> marginal_dist
        ensureParameters!(marginal_dist, (:m, :V))
        @fact edge.marginal.m --> [0.0]
        @fact isApproxEqual(edge.marginal.V, reshape([0.5], 1, 1)) --> true
    end

    context("prod! for MvGaussian") do
        marg = ensureParameters!(MvGaussian(m=[0.0], V=eye(1)) * MvGaussian(m=[0.0], V=eye(1)), (:m, :V))
        @fact marg.m --> [0.0]
        @fact isApproxEqual(marg.V, reshape([0.5], 1, 1)) --> true
        @fact MvGaussian(m=zeros(2), V=eye(2)) * MvDelta(ones(2)) --> MvDelta(ones(2))
        @fact MvDelta(ones(2)) * MvGaussian(m=zeros(2), V=eye(2)) --> MvDelta(ones(2))
        @fact_throws MethodError MvDelta(ones(3)) * MvGaussian(m=zeros(2), V=eye(2))
        @fact ForneyLab.prod!(MvGaussian(m=zeros(2), V=eye(2)), MvDelta(ones(2)), MvGaussian(m=zeros(2), V=eye(2))) --> MvGaussian(m=ones(2), V=tiny*eye(2))
        @fact ForneyLab.prod!(MvDelta(ones(2)), MvGaussian(m=zeros(2), V=eye(2)), MvGaussian(m=zeros(2), V=eye(2))) --> MvGaussian(m=ones(2), V=tiny*eye(2))
    end
end

facts("MvGaussian converts") do
    context("MvDelta should be convertible to MvGaussian with tiny variance") do
        @fact convert(MvGaussian, MvDelta([1.0, 2.0])) --> MvGaussian(m=[1.0, 2.0], V=tiny*diageye(2)) # Vector
    end

    context("Message{MvDelta} should be convertible to Message{MvGaussian} with tiny variance") do
        @fact convert(Message{MvGaussian}, Message(MvDelta([1.0, 2.0]))) --> Message(MvGaussian(m=[1.0, 2.0], V=tiny*diageye(2))) # Vector
    end

    context("Gaussian should be convertible to MvGaussian") do
        @fact convert(MvGaussian, Gaussian(m=1.0, V=2.0)) --> MvGaussian(m=[1.0], V=2.0*eye(1))
    end

    context("MvGaussian should be convertible to Gaussian") do
        @fact convert(Gaussian, MvGaussian(m=[1.0], V=2.0*eye(1))) --> Gaussian(m=1.0, V=2.0)
        @fact_throws convert(Gaussian, MvGaussian(m=ones(2), V=eye(2)))
    end
end
