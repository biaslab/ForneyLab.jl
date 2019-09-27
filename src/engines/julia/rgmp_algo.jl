mutable struct rgmpAlgorithm
    g::Function # Pdf function of output variable
    num_epochs1::Int64 # Epochs to update the variational parameters
    num_epochs2::Int64 # Number of epochs for for monte carlo gradients which can be set to 1 for fully stochastic optimization
    lr::Float64 # Learning rate
    alg::String

    function rgmpAlgorithm()

    function RGMP(out, in1, g::Function, num_epochs1::Int64, num_epochs2::Int64, lr::Float64, id=ForneyLab.generateId(RGMP))
        @ensureVariables(out, in1)
        self = new(id, Array{Interface}(undef, 2), Dict{Symbol,Interface}(), g, num_epochs1, num_epochs2, lr)
        ForneyLab.addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:in1] = self.interfaces[2] = associate!(Interface(self), in1)

        return self
    end
end
