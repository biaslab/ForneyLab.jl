# User API

Documentation for `ForneyLab.jl`'s user API.

If you want to know how you can extend `ForneyLab.jl` (e.g. register new update rules), see [Developer API](@ref).

## Contents
```@contents
Pages = ["user-api.md"]
Depth = 5
```

## Index
```@index
Modules = [ForneyLab]
Pages = ["user-api.md"]
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
ForneyLab.MarginalTable
ForneyLab.PosteriorFactorization
ForneyLab.InferenceAlgorithm
ForneyLab.Schedule
ForneyLab.currentInferenceAlgorithm
```

## Algorithm assembly
```@docs
ForneyLab.messagePassingAlgorithm
```

## Algorithm code generation
```@docs
ForneyLab.algorithmSourceCode
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
