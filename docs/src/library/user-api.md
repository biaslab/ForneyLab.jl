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
ForneyLab.Algorithm
ForneyLab.Schedule
ForneyLab.currentAlgorithm
ForneyLab.expectationPropagationSchedule
ForneyLab.sumProductSchedule
ForneyLab.variationalExpectationPropagationSchedule
ForneyLab.variationalSchedule
```

## Algorithm assembly
```@docs
ForneyLab.expectationPropagationAlgorithm
ForneyLab.freeEnergyAlgorithm
ForneyLab.sumProductAlgorithm
ForneyLab.variationalAlgorithm
ForneyLab.variationalExpectationPropagationAlgorithm
```

## Algorithm code generation
```@docs
ForneyLab.algorithmString
ForneyLab.freeEnergyString
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
