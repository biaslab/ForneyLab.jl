module NonconjugateTest

using Test
using ForneyLab
import ForneyLab: outboundType, isApplicable, prod!, logPdf, unsafeMean, unsafeVar
import ForneyLab: SPNonconjugateInFN, SPNonconjugateOutNG

@testset "prod!" begin
    #@test ProbabilityDistribution(Bernoulli, p=0.2) * ProbabilityDistribution(Bernoulli, p=0.8) == ProbabilityDistribution(Bernoulli, p=0.5000000000000001)
    #@test_throws Exception ProbabilityDistribution(Bernoulli, p=0.0) * ProbabilityDistribution(Bernoulli, p=1.0)
    println(convert(ProbabilityDistribution{Univariate, GaussianMeanVariance}, ProbabilityDistribution(Univariate, GaussianMeanVariance, m=1.0, v=2.0)*ProbabilityDistribution(Univariate, GaussianMeanVariance, m=2.0, v=3.0)))
    f_dummy(x) = x
    println(ruleSPNonconjugateInFN(Message(Univariate, GaussianMeanVariance, m=1.0, v=2.0),nothing,f_dummy).dist*ProbabilityDistribution(Univariate, GaussianMeanVariance, m=2.0, v=3.0))
end

#-------------
# Update rules
#-------------

@testset "SPNonconjugateInFN" begin
    f_dummy(x) = x
    @test SPNonconjugateInFN <: SumProductRule{Nonconjugate}
    @test outboundType(SPNonconjugateInFN) == Message{Function}
    @test isApplicable(SPNonconjugateInFN, [Message{Union{Bernoulli, Beta, Categorical, Dirichlet, Gaussian, Gamma, LogNormal, Poisson, Wishart}}, Nothing])
    f(x) = ruleSPNonconjugateInFN(Message(Univariate, GaussianMeanVariance, m=2.0, v=1.0),nothing,f_dummy).dist.params[:log_pdf](x)
    @test f(1.5) == logPdf(ProbabilityDistribution(Univariate, GaussianMeanVariance, m=2.0, v=1.0),1.5)
end

 @testset "SPNonconjugateOutNG" begin
     f_dummy(x) = x
     @test SPNonconjugateOutNG <: SumProductRule{Nonconjugate}
     @test outboundType(SPNonconjugateOutNG) == Message{Abstract_dist}
     @test isApplicable(SPNonconjugateOutNG, [Nothing, Message{Gaussian}])
     @test abs(ruleSPNonconjugateOutNG(nothing,Message(Univariate, GaussianMeanVariance, m=2.0, v=1.0),f_dummy).dist.params[:m] - unsafeMean(ProbabilityDistribution(Univariate, Abstract_dist, m=2.0, v=1.0, f=f_dummy))) < 0.05
     @test abs(ruleSPNonconjugateOutNG(nothing,Message(Univariate, GaussianMeanVariance, m=2.0, v=1.0),f_dummy).dist.params[:v] - unsafeVar(ProbabilityDistribution(Univariate, Abstract_dist, m=2.0, v=1.0, f=f_dummy))) < 0.05
 end

end # module
