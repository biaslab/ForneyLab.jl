# Developer API

Documentation for `ForneyLab.jl`'s developer API.

## Contents
```@contents
Pages = ["developer-api.md"]
Depth = 5
```

## Index
```@index
Modules = [ForneyLab]
Pages = ["developer-api.md"]
Order = [:macro, :module, :constant, :type, :function]
```

## Extended rules registration
```@docs
ForneyLab.@expectationPropagationRule
ForneyLab.@marginalRule
ForneyLab.@naiveVariationalRule
ForneyLab.@structuredVariationalRule
ForneyLab.@sumProductRule
```

## Graph (low-level)
```@docs
ForneyLab.Edge
ForneyLab.Interface
ForneyLab.Terminal
```

## Scheduler (low-level)
```@docs
ForneyLab.MarginalRule
ForneyLab.MarginalEntry
ForneyLab.MarginalUpdateRule
ForneyLab.MessageUpdateRule
ForneyLab.PosteriorFactor
ForneyLab.ScheduleEntry
```
