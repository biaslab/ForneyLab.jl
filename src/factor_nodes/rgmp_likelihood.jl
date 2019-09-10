export RGMP_likelihood

"""
Description:

    Non-conjugate likelihood functions can be incorporated into
    the existing factor graph with rgmp_likelihood node. in1 interface
    expects Gaussian message and sends an approximated message back.
    Primitive implementation requires input-output relation as log_pdf function.
    Optimizer by default is stochastic gradient descent. Initial points of posterior
    are the statistics of the incoming Gaussian message.

Interfaces:

    1. out
    2. in1

Construction:

    RGMP_likelihood(out, in1, g::Function, num_epochs1::Int64, num_epochs2::Int64, lr::Float64)
"""
mutable struct RGMP_likelihood <: SoftFactor
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol, Interface}

    g::Function # Pdf function of output variable
    num_epochs1::Int64 # Epochs to update the variational parameters
    num_epochs2::Int64 # Number of epochs for for monte carlo gradients which can be set to 1 for fully stochastic optimization
    lr::Float64 # Learning rate

    function RGMP_likelihood(out, in1, g::Function, num_epochs1::Int64, num_epochs2::Int64, lr::Float64, id=ForneyLab.generateId(RGMP_likelihood))
        @ensureVariables(out, in1)
        self = new(id, Array{Interface}(undef, 2), Dict{Symbol,Interface}(), g, num_epochs1, num_epochs2, lr)
        ForneyLab.addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:in1] = self.interfaces[2] = associate!(Interface(self), in1)

        return self
    end
end

slug(::Type{RGMP_likelihood}) = "RGMP"
