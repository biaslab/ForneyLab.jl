# Construct a model

### Index
```@index
Modules = [ForneyLab]
Pages = ["construct-a-model.md"]
Order = [:macro, :module, :constant, :type, :function]
```


### Random variable
```@docs
ForneyLab.@RV
```

### Factor graph
```@docs
ForneyLab.FactorGraph
```

### Variable
```@docs
ForneyLab.Variable
```

### Current graph
```@docs
ForneyLab.currentGraph
```

### Node types

```@autodocs
Modules = [ForneyLab]
Private = false
Pages = collect(Iterators.flatten([[joinpath(root[4:end], file) for file in files] for (root, dirs, files) in walkdir("../src/factor_nodes/")]))
Order = [:macro, :module, :constant, :type, :function]
```
