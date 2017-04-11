function rule(::Type{ForneyLab.SPGaussianMeanVariancePPV}, inputs)
    m = unsafeMean(inputs[1].payload)
    v = unsafeMean(inputs[2].payload)

    return Message{GaussianMeanVariance}(m, v)
end