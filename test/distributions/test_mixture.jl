facts("Mixture unit tests") do
    context("Construction") do
        dd = Mixture([Gaussian(); Gaussian()], [0.4; 0.6])
        @fact typeof(dd) --> Mixture{Gaussian}
        @fact dd.components --> [Gaussian(), Gaussian()]
        @fact dd.weights --> [0.4; 0.6]
        @fact_throws Mixture([Gaussian(), Gaussian()])
        @fact_throws Mixture([Gaussian(), Gamma()], [0.4; 0.6])
        @fact_throws Mixture([MvGaussian(m=zeros(2),V=eye(2)), MvGaussian(m=zeros(3),V=eye(3))], [0.4; 0.6])
        dd = Mixture([Gaussian(m=1.0, V=2.0); Gaussian(m=2.0, V=3.0)], [0.4; 0.6])
        @fact pdf(dd, 1.5) --> roughly(0.4*pdf(Gaussian(m=1.0, V=2.0), 1.5) + 0.6*pdf(Gaussian(m=2.0, V=3.0), 1.5), atol=1e-6)
    end

    context("vague! and vague should be implemented") do
        dd = Mixture([Gaussian(); Gaussian()], [0.4; 0.6])
        ForneyLab.vague!(dd)
        @fact dd.components[1] --> vague(Gaussian)
        @fact dd.components[2] --> vague(Gaussian)
        vague_d = vague(Mixture{MvGaussian{2}})
        @fact vague_d.components[1] --> vague(MvGaussian{2})
    end

    context("mean should be implemented") do
        dd = Mixture([Gaussian(m=3.0; V=2.0); Gaussian(m=1.0; V=3.0)], [0.4; 0.6])
        @fact mean(dd) --> 0.4*3.0 + 0.6*1.0
    end

    context("== operator") do
        d1 = Mixture([Gaussian(m=1.0,V=2.0); Gaussian(m=3.0,V=2.0)], [0.3;0.7])
        d2 = Mixture([Gaussian(m=1.0,V=2.0); Gaussian(m=3.0,V=2.0)], [0.3;0.7])
        d3 = Mixture([Gaussian(m=1.0,V=2.0); Gaussian(m=3.1,V=2.0)], [0.3;0.7])
        d4 =
        @fact d1 --> d2
        @fact (d1==d3) --> false
        @fact (d1 == Mixture([Gaussian(m=1.0,V=2.0); Gaussian(m=3.0,V=2.0)], [0.7;0.3])) --> false
        @fact d1 --> Mixture([Gaussian(m=3.0,V=2.0); Gaussian(m=1.0,V=2.0)], [0.7;0.3])
    end

    context("dimensions(::Mixture)") do
        @fact dimensions(Mixture([Gaussian(m=1.0,V=2.0); Gaussian(m=3.0,V=2.0)], [0.3;0.7])) --> 1
        @fact dimensions(Mixture{MvGaussian{4}}) --> 4
    end

    context("prod() should yield correct result") do
        # Test product of mixture and single component
        d = Gaussian(m=0.0, V=1.0)
        dm = Mixture([Gaussian(m=1.0, V=0.1); Gaussian(m=2.0, V=0.1)], [0.4; 0.6])
        dp = dm * d
        x_test = [1.0; 1.5; 0.0] # evaluate the PDF of the product at these points, and check that it is proportional to pdf(d,x)*pdf(dm,x)
        scaling_factors = [pdf(dp,x) / (pdf(d, x) * pdf(dm, x)) for x in x_test]
        @fact scaling_factors[2] --> roughly(scaling_factors[1], atol=1e-6)
        @fact scaling_factors[3] --> roughly(scaling_factors[1], atol=1e-6)

        # Test product of two mixtures
        d1 = Mixture([Gaussian(m=1.0, V=0.1); Gaussian(m=2.0, V=1.0)], [0.4; 0.6])
        d2 = Mixture([Gaussian(m=3.0, V=0.2); Gaussian(m=-4.0, V=3.0)], [0.1; 0.9])
        d3 = d1 * d2
        @fact length(d3.components) --> 4
        x_test = [1.0; 1.5; 0.0] # evaluate the PDF of the product at these points, and check that it is proportional to pdf(d,x)*pdf(dm,x)
        scaling_factors = [pdf(d3,x) / (pdf(d1, x) * pdf(d2, x)) for x in x_test]
        @fact scaling_factors[2] --> roughly(scaling_factors[1], atol=1e-6)
        @fact scaling_factors[3] --> roughly(scaling_factors[1], atol=1e-6)
    end
end