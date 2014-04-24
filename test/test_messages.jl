facts("GaussionMessage") do
    context("GaussianMessage() should initialize a Gaussian message") do
        @fact GaussianMessage().V => ones(1, 1)
        @fact GaussianMessage().m => [0.0]
        @fact GaussianMessage(m=[0.0], W=[1.0]).W => ones(1, 1)
        @fact GaussianMessage(xi=[0.0], W=[1.0]).W => ones(1, 1)
        @fact typeof(GaussianMessage(m=[1.0], V=[1.0]).V) => Array{Float64, 2} # cast single value to matrix
        @fact_throws GaussianMessage(V=[1.0], W=[1.0])
        @fact_throws GaussianMessage(m=[0.0], xi=[0.0])
        @fact_throws GaussianMessage(xi=[0.0], V=[1.0])
    end
end

facts("GeneralMessage") do
    context("GeneralMessage() should initiatize a multiplication parameter as message") do
        @fact GeneralMessage(1.0).value => 1.0
        @fact GeneralMessage([1.0, 2.0]).value => [1.0, 2.0]
    end
end