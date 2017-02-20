type GaussianMeanVariance
    id::Symbol
    interfaces::Vector{Interface}

    function GaussianMeanVariance(out::Variable, mean::Variable, variance::Variable; id=generateNodeId(GaussianMeanVariance))
        self = new(id, Array(Interface, 3))
        self.interfaces[1] = connect!(Interface(self), mean)
        self.interfaces[2] = connect!(Interface(self), variance)
        self.interfaces[3] = connect!(Interface(self), out)

        addNode!(currentGraph(), self)

        return self
    end
end

# TODO: rules id dictionary