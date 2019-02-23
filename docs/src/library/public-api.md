# User API

Documentation for `ForneyLab.jl`'s public API.

See [Internal API](@ref) for internal package docs.

## Contents
```@contents
Pages = ["public-api.md"]
Depth = 5
```

## Index
```@index
Modules = [ForneyLab]
Pages = ["public-api.md"]
Order = [:macro, :module, :constant, :type, :function]
```

## Model specification
```@docs
ForneyLab.@RV
ForneyLab.FactorGraph
ForneyLab.Variable
ForneyLab.currentGraph
```

### Factor nodes
```@autodocs
Modules = [ForneyLab]
Private = false
Pages = collect(Iterators.flatten([[joinpath(root[4:end], file) for file in files] for (root, dirs, files) in walkdir("../src/factor_nodes/")]))
Order = [:macro, :module, :constant, :type, :function]
```

## Scheduling
```@docs
ForneyLab.MarginalSchedule
ForneyLab.RecognitionFactorization
ForneyLab.Schedule
ForneyLab.currentRecognitionFactorization
ForneyLab.expectationPropagationSchedule
ForneyLab.sumProductSchedule
ForneyLab.variationalExpectationPropagationSchedule
ForneyLab.variationalSchedule
```

## Algorithm generation
```@docs
ForneyLab.expectationPropagationAlgorithm
ForneyLab.freeEnergyAlgorithm
ForneyLab.sumProductAlgorithm
ForneyLab.variationalAlgorithm
ForneyLab.variationalExpectationPropagationAlgorithm
```

## Algorithm execution
```@docs
ForneyLab.Message
ForneyLab.PointMass
ForneyLab.ProbabilityDistribution
```

## Helper

```@autodocs
Modules = [ForneyLab]
Private = true
Pages = ["src/helpers.jl"]
Order = [:macro, :module, :constant, :type, :function]
```
