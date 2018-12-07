# General

### Index
```@index
Modules = [ForneyLab]
Pages = ["general.md"]
Order = [:macro, :module, :constant, :type, :function]
```

### Description
```@autodocs
Modules = [ForneyLab]
Private = false
Pages = [[joinpath(root[4:end], file) for file in files] for (root,dirs,files) in walkdir("../src/")][1]
Order = [:macro, :module, :constant, :type, :function]
```
