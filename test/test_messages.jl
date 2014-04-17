context("Specific message types") do
	facts("GaussianMessage should initialize a Gaussian message") do
		@fact GaussianMessage().V => [1.0]
		@fact GaussianMessage().m => [0.0]
		@fact GaussianMessage(m=[0.0],W=[1.0]).W => [1.0]
		@fact GaussianMessage(xi=[0.0],W=[1.0]).W => [1.0]
		@fact_throws GaussianMessage(V=[1.0],W=[1.0])
		@fact_throws GaussianMessage(m=[0.0],xi=[0.0])
		@fact_throws GaussianMessage(xi=[0.0],V=[1.0]) 
	end

	facts("GeneralMessage should initiatize a multiplication parameter as message") do
		@fact GeneralMessage(1.0).value => 1.0
		@fact GeneralMessage([1.0, 2.0]).value => [1.0, 2.0]
	end
end