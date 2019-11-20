module NonconjugateTest

using Test
using ForneyLab
import ForneyLab: outboundType, isApplicable, prod!, logPdf, unsafeMean, unsafeVar
import ForneyLab: SPNonlinearPTInMN, SPNonlinearPTOutNG

@testset "prod!" begin
    f_dummy(x) = x
    @test abs((convert(ProbabilityDistribution{Univariate, GaussianMeanVariance}, ProbabilityDistribution(Univariate, GaussianMeanVariance, m=4.0, v=2.0)
            *ProbabilityDistribution(Univariate, GaussianMeanVariance, m=2.0, v=3.0))).params[:m]
            -  (ruleSPNonlinearPTInMN(Message(Univariate, GaussianMeanVariance, m=4.0, v=2.0),nothing,f_dummy).dist*ProbabilityDistribution(Univariate, GaussianMeanVariance, m=2.0, v=3.0)).params[:m]) < 0.1
    @test abs((convert(ProbabilityDistribution{Univariate, GaussianMeanVariance}, ProbabilityDistribution(Univariate, GaussianMeanVariance, m=4.0, v=2.0)
            *ProbabilityDistribution(Univariate, GaussianMeanVariance, m=2.0, v=3.0))).params[:v]
            -  (ruleSPNonlinearPTInMN(Message(Univariate, GaussianMeanVariance, m=4.0, v=2.0),nothing,f_dummy).dist*ProbabilityDistribution(Univariate, GaussianMeanVariance, m=2.0, v=3.0)).params[:v]) < 0.1
    @test abs((convert(ProbabilityDistribution{Univariate, GaussianMeanVariance}, ProbabilityDistribution(Univariate, GaussianMeanVariance, m=1.2, v=1.0)
            *ProbabilityDistribution(Univariate, GaussianMeanVariance, m=2.6, v=0.5))).params[:m]
            -  (ruleSPNonlinearPTInMN(Message(Univariate, GaussianMeanVariance, m=1.2, v=1.0),nothing,f_dummy).dist*ProbabilityDistribution(Univariate, GaussianMeanVariance, m=2.6, v=0.5)).params[:m]) < 0.1
    @test abs((convert(ProbabilityDistribution{Univariate, GaussianMeanVariance}, ProbabilityDistribution(Univariate, GaussianMeanVariance, m=6.5, v=4.1)
            *ProbabilityDistribution(Univariate, GaussianMeanVariance, m=12.0, v=3.0))).params[:v]
            -  (ruleSPNonlinearPTInMN(Message(Univariate, GaussianMeanVariance, m=6.5, v=4.1),nothing,f_dummy).dist*ProbabilityDistribution(Univariate, GaussianMeanVariance, m=12.0, v=3.0)).params[:v]) < 0.1
end

#-------------
# Update rules
#-------------

@testset "SPNonlinearPTInMN" begin
    f_dummy(x) = x
    @test SPNonlinearPTInMN <: SumProductRule{Nonconjugate}
    @test outboundType(SPNonlinearPTInMN) == Message{Function}
    @test isApplicable(SPNonlinearPTInMN, [Message{Union{Bernoulli, Beta, Categorical, Dirichlet, Gaussian, Gamma, LogNormal, Poisson, Wishart}}, Nothing])
    f(x) = ruleSPNonlinearPTInMN(Message(Univariate, GaussianMeanVariance, m=2.0, v=1.0),nothing,f_dummy).dist.params[:log_pdf](x)
    @test f(1.5) == logPdf(ProbabilityDistribution(Univariate, GaussianMeanVariance, m=2.0, v=1.0),1.5)
end

 @testset "SPNonlinearPTOutNG" begin
     f_dummy(x) = x
     @test SPNonlinearPTOutNG <: SumProductRule{Nonconjugate}
     @test outboundType(SPNonlinearPTOutNG) == Message{Abstract_dist}
     @test isApplicable(SPNonlinearPTOutNG, [Nothing, Message{Gaussian}])
     @test abs(ruleSPNonlinearPTOutNG(nothing,Message(Univariate, GaussianMeanVariance, m=2.0, v=1.0),f_dummy).dist.params[:m] - unsafeMean(ProbabilityDistribution(Univariate, Abstract_dist, m=2.0, v=1.0, f=f_dummy))) < 0.05
     @test abs(ruleSPNonlinearPTOutNG(nothing,Message(Univariate, GaussianMeanVariance, m=2.0, v=1.0),f_dummy).dist.params[:v] - unsafeVar(ProbabilityDistribution(Univariate, Abstract_dist, m=2.0, v=1.0, f=f_dummy))) < 0.05
 end

end # module
