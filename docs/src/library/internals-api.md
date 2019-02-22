# Internals API

## Index
```@index
Modules = [ForneyLab]
Pages = ["user-api.md"]
Order = [:macro, :module, :constant, :type, :function]
```

## Extended factor nodes
```@docs
ForneyLab.@composite
```

## Extended rules
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
ForneyLab.MarginalScheduleEntry
ForneyLab.MarginalUpdateRule
ForneyLab.MessageUpdateRule
ForneyLab.RecognitionFactor
ForneyLab.ScheduleEntry
```
