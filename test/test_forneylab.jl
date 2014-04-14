module TestForneyLab

using FactCheck
using ForneyLab

context("Messages") do
	facts("It should initialize gaussian messages") do
		@fact ForneyLab.GaussianMessage().V => [1.0]
		@fact ForneyLab.GaussianMessage().m => [0.0]
	end

	facts("It should initiatize a multiplication parameter as message") do
		@fact ForneyLab.ScalarParameterMessage(1.0).value => 1.0
		@fact ForneyLab.ScalarParameterMessage([1.0, 2.0]).value => [1.0, 2.0]
	end
end

context("Nodes") do
	facts("It should initialize multiplication nodes") do
		@fact typeof(ForneyLab.MultiplicationNode()) => ForneyLab.MultiplicationNode
		@fact ForneyLab.MultiplicationNode().multiplier.message => nothing
	end
end

context("Edges") do
	facts("It should initialize interfaces as half-edges on a node") do
	end

	facts("It should initiatize edges as a combination of two interfaces") do
	end
end

context("calculatemessage") do
end

end # module TestForneyLab