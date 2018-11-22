# Algorithms

### Index
```@index
Modules = [ForneyLab]
Pages = ["algorithms.md"]
Order = [:macro, :module, :constant, :type, :function]
```

### Description
```@autodocs
Modules = [ForneyLab]
Private = false
Pages = collect(Iterators.flatten([[joinpath(root[4:end], file) for file in files] for (root, dirs, files) in walkdir("../src/algorithms/")]))
Order = [:macro, :module, :constant, :type, :function]
```
