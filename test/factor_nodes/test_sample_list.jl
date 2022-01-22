module SampleListTest

using Test
using ForneyLab
using ForneyLab: outboundType, isApplicable, prod!, unsafeMean, unsafeCov, unsafeVar, dims, bootstrap
using ForneyLab: SPSampleListOutNPP, VBSampleListOut

@testset "dims" begin
    @test dims(Distribution(Univariate, SampleList, s=[0.0, 1.0], w=[0.5, 0.5])) == ()
    @test dims(Distribution(Multivariate, SampleList, s=[[0.0], [1.0]], w=[0.5, 0.5])) == (1,)
    @test dims(Distribution(MatrixVariate, SampleList, s=[mat(0.0), mat(1.0)], w=[0.5, 0.5])) == (1,1)
end

@testset "unsafeMean" begin
    @test unsafeMean(Distribution(Univariate, SampleList, s=[1.0, 1.0], w=[0.5, 0.5])) == 1.0
    @test unsafeMean(Distribution(Multivariate, SampleList, s=[[1.0], [1.0]], w=[0.5, 0.5])) == [1.0]
    @test unsafeMean(Distribution(MatrixVariate, SampleList, s=[mat(1.0), mat(1.0)], w=[0.5, 0.5])) == mat(1.0)
end

@testset "unsafeVar" begin
    @test unsafeVar(Distribution(Univariate, SampleList, s=[1.0, 1.0], w=[0.5, 0.5])) == 0.0
    @test unsafeVar(Distribution(Multivariate, SampleList, s=[[1.0], [1.0]], w=[0.5, 0.5])) == [0.0]
end

@testset "unsafeCov" begin
    @test unsafeCov(Distribution(Univariate, SampleList, s=[1.0, 1.0], w=[0.5, 0.5])) == 0.0
    @test unsafeCov(Distribution(Multivariate, SampleList, s=[[1.0], [1.0]], w=[0.5, 0.5])) == mat(0.0)
    @test unsafeCov(Distribution(MatrixVariate, SampleList, s=[eye(2), eye(2)], w=[0.5, 0.5])) == zeros(4,4)
end

@testset "Univariate SampleList prod!" begin
    dist_gaussian = Distribution(Univariate, GaussianMeanVariance, m=0.0, v=1.0)
    dist_prod = dist_gaussian * Distribution(Univariate, SampleList, s=[0.0, 1.0], w=[0.5, 0.5])
    dist_true_prod = Distribution(Univariate, SampleList, s=[0.0, 1.0], w=[0.6224593312018546, 0.37754066879814546])

    @test dist_prod.params[:s] == dist_true_prod.params[:s]
    @test dist_prod.params[:w] == dist_true_prod.params[:w]
    @test dist_prod.params[:logintegrand]([-0.7, 1.5]) == logPdf.([dist_gaussian],[-0.7, 1.5])
end

@testset "Multivariate SampleList prod!" begin
    dist_gaussian = Distribution(Multivariate, GaussianMeanVariance, m=[0.0], v=mat(1.0))
    dist_prod = dist_gaussian * Distribution(Multivariate, SampleList, s=[[0.0], [1.0]], w=[0.5, 0.5])
    dist_true_prod = Distribution(Multivariate, SampleList, s=[[0.0], [1.0]], w=[0.6224593312018546, 0.37754066879814546])

    @test dist_prod.params[:s] == dist_true_prod.params[:s]
    @test dist_prod.params[:w] == dist_true_prod.params[:w]
    @test dist_prod.params[:logintegrand]([[-0.7], [1.5]]) == logPdf.([dist_gaussian],[[-0.7], [1.5]])
end

@testset "Initialization of logproposal, logintegrand, unnormalizedweights" begin
    dist_beta = Distribution(Univariate, Beta, a=2.0, b=1.0)
    dist_gamma = Distribution(Univariate, Gamma, a=2.0, b=1.0)
    dist_sample = Distribution(Univariate, SampleList, s=[0.0, 1.0], w=[0.5, 0.5])
    dist_betagamma = dist_beta*dist_gamma
    dist_gammasample = dist_gamma*dist_sample
    dist_betagg = dist_betagamma*dist_gamma

    @test haskey(dist_betagamma.params, :logproposal)
    @test haskey(dist_betagamma.params, :logintegrand)
    @test haskey(dist_betagamma.params, :unnormalizedweights)

    @test !haskey(dist_gammasample.params, :logproposal)
    @test haskey(dist_gammasample.params, :logintegrand)
    @test !haskey(dist_gammasample.params, :unnormalizedweights)

    @test haskey(dist_betagg.params, :logproposal)
    @test haskey(dist_betagg.params, :logintegrand)
    @test haskey(dist_betagg.params, :unnormalizedweights)
end

@testset "bootstrap" begin
    p1 = Distribution(Univariate, SampleList, s=[2.0], w=[1.0])
    p2 = Distribution(Univariate, PointMass, m=0.0)
    @test bootstrap(p1, p2) == [2.0]

    p1 = Distribution(Multivariate, SampleList, s=[[2.0]], w=[1.0])
    p2 = Distribution(MatrixVariate, PointMass, m=mat(tiny))
    @test isapprox(bootstrap(p1, p2)[1][1], 2.0, atol=1e-4)

    p1 = Distribution(Univariate, GaussianMeanVariance, m=2.0, v=0.0)
    p2 = Distribution(Univariate, SampleList, s=[0.0], w=[1.0])
    @test bootstrap(p1, p2) == [2.0]

    p1 = Distribution(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(0.0))
    p2 = Distribution(MatrixVariate, SampleList, s=[mat(tiny)], w=[1.0])
    @test isapprox(bootstrap(p1, p2)[1][1], 2.0, atol=1e-4)
end


#-------------
# Update rules
#-------------

@testset "SPSampleListOutNPP" begin
    @test SPSampleListOutNPP <: SumProductRule{SampleList}
    @test outboundType(SPSampleListOutNPP) == Message{SampleList}
    @test isApplicable(SPSampleListOutNPP, [Nothing, Message{PointMass}, Message{PointMass}])

    @test ruleSPSampleListOutNPP(nothing, Message(Multivariate, PointMass, m=[0.0,2.0,5.2]), Message(Multivariate, PointMass, m=[0.4,0.5,0.1])) == Message(Univariate, SampleList, s=[0.0,2.0,5.2], w=[0.4,0.5,0.1])
    @test ruleSPSampleListOutNPP(nothing, Message(Multivariate, PointMass, m=[[0.0,2.0],[5.2,-0.6]]), Message(Multivariate, PointMass, m=[0.4,0.6])) == Message(Multivariate, SampleList, s=[[0.0,2.0],[5.2,-0.6]], w=[0.4,0.6])
    @test ruleSPSampleListOutNPP(nothing, Message(Multivariate, PointMass, m=[eye(2),eye(2)]), Message(Multivariate, PointMass, m=[0.4,0.6])) == Message(MatrixVariate, SampleList, s=[eye(2),eye(2)], w=[0.4,0.6])
end

@testset "VBSampleListOut" begin
    @test VBSampleListOut <: NaiveVariationalRule{SampleList}
    @test outboundType(VBSampleListOut) == Message{SampleList}
    @test isApplicable(VBSampleListOut, [Nothing, Distribution, Distribution])

    @test ruleVBSampleListOut(nothing, Distribution(Multivariate, PointMass, m=[0.0,2.0,5.2]), Distribution(Multivariate, PointMass, m=[0.4,0.5,0.1])) == Message(Univariate, SampleList, s=[0.0,2.0,5.2], w=[0.4,0.5,0.1])
    @test ruleVBSampleListOut(nothing, Distribution(Multivariate, PointMass, m=[[0.0,2.0],[5.2,-0.6]]), Distribution(Multivariate, PointMass, m=[0.4,0.6])) == Message(Multivariate, SampleList, s=[[0.0,2.0],[5.2,-0.6]], w=[0.4,0.6])
    @test ruleVBSampleListOut(nothing, Distribution(Multivariate, PointMass, m=[eye(2),eye(2)]), Distribution(Multivariate, PointMass, m=[0.4,0.6])) == Message(MatrixVariate, SampleList, s=[eye(2),eye(2)], w=[0.4,0.6])
end


#------------
# Integration
#------------

g(x) = x

@testset "Sampling code generation" begin
    # Define a model
    fg = FactorGraph()
    N = 4
    @RV s ~ SampleList(collect(1.0:N), ones(N)/N)
    @RV x ~ GaussianMeanVariance(s, 0.0)
    @RV y ~ Nonlinear{Sampling}(x, g=g)

    # Define an algorithm
    pfz = PosteriorFactorization(fg)
    algo = messagePassingAlgorithm(y)
    code = algorithmSourceCode(algo)

    @test occursin("ruleSPSampleListOutNPP", code)
    @test occursin("ruleSPGaussianMeanVarianceOutNSP", code)
    @test occursin("ruleSPNonlinearSOutNM", code)
end

# Generated algorithm code
function step!(data::Dict, marginals::Dict=Dict(), messages::Vector{Message}=Array{Message}(undef, 3))
    messages[1] = ruleSPSampleListOutNPP(nothing, Message(Multivariate, PointMass, m=[1.0, 2.0, 3.0, 4.0]), Message(Multivariate, PointMass, m=[0.25, 0.25, 0.25, 0.25]))
    messages[2] = ruleSPGaussianMeanVarianceOutNSP(nothing, messages[1], Message(Univariate, PointMass, m=0.0))
    messages[3] = ruleSPNonlinearSOutNM(g, nothing, messages[2])

    marginals[:y] = messages[3].dist

    return marginals
end

@testset "Sampling code execution" begin
    # Execute the algorithm code
    marginals = step!(Dict())

    @test marginals[:y] == Distribution(Univariate, SampleList, s=collect(1.0:4.0), w=ones(4)/4)
end

@testset "Differential Entropy estimate" begin
    fg = FactorGraph()

    @RV x ~ Gamma(3.4, 1.2)
    @RV s1 ~ Nonlinear{Sampling}(x, g=g)
    @RV y1 ~ Poisson(s1)
    @RV y2 ~ Poisson(x)

    @RV z ~ Gamma(3.4, 1.2)
    @RV s2 ~ Nonlinear{Sampling}(z, g=g)
    @RV y3 ~ Poisson(z)
    @RV y4 ~ Poisson(s2)

    @RV w ~ Gamma(3.4, 1.2)
    @RV y5 ~ Poisson(w)
    @RV y6 ~ Poisson(w)

    placeholder(y1,:y1)
    placeholder(y2,:y2)
    placeholder(y3,:y3)
    placeholder(y4,:y4)
    placeholder(y5,:y5)
    placeholder(y6,:y6)

    pfz = PosteriorFactorization(fg)
    algo = messagePassingAlgorithm([x,z,w])
    source_code = algorithmSourceCode(algo)
    eval(Meta.parse(source_code));

    v_1, v_2 = 7, 5
    data = Dict(:y1 => v_1, :y2 => v_2, :y3 => v_1, :y4 => v_2, :y5 => v_1, :y6 => v_2)
    marginals = step!(data)

    @test isa(marginals[:x].params[:entropy], Float64)
    @test isa(marginals[:z].params[:entropy], Float64)
end

end
