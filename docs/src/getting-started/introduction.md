## Factor Graphs and Variables

A central concept in ForneyLab is the (random) `Variable` type. After including ForneyLab and indicating that we start a new `FactorGraph`, we can declare a `Variable` by calling its constructor function:


```julia
using ForneyLab

# Declare a new graph
g = FactorGraph()

# Declare a variable
x = Variable(id=:x)

g.variables
```




    Dict{Symbol,Variable} with 1 entry:
      :x => Variable(:x, Edges:…



The call to `FactorGraph()` creates a factor graph type and registers the graph as the currently active graph. Note that the variable has been associated with an edge in the currently active graph.

ForneyLab comes equipped with the `@RV` macro to define random variables. For instance, defining a new variable `y` with identifier `:y` and associating the variable to the current graph can also be accomplished by executing `@RV y`:


```julia
@RV y

g.variables
```




    Dict{Symbol,Variable} with 2 entries:
      :y => Variable(:y, Edges:…
      :x => Variable(:x, Edges:…



We can assign a probability distribution to a random variable by the `~` operator:


```julia
@RV z ~ GaussianMeanVariance(0.0, 1.0)
g.variables
```




    Dict{Symbol,Variable} with 5 entries:
      :clamp_2 => Variable(:clamp_2, Edges:…
      :y       => Variable(:y, Edges:…
      :clamp_1 => Variable(:clamp_1, Edges:…
      :z       => Variable(:z, Edges:…
      :x       => Variable(:x, Edges:…



Note that the graph now also includes two variables with id `:clamp_1` and `:clamp_2`. These two variables correspond to the mean and variance parameters for the Gaussian and are clamped to values `0.0` and `1.0` respectively.

If you have [graphviz](https://www.graphviz.org/) installed, then you can draw the factor graph. (Edges (variables) that are not constrained by any factor are not drawn):


```julia
ForneyLab.draw(g) # draw the graph
```



In case you don't like the automatically generated id's `:clamp_1` and `:clamp_2`, you could have declared the parameters of the Gaussian distribution through the `@RV` macro and associated a `Clamp` distribution with these variables:


```julia
g2 = FactorGraph()
@RV m ~ Clamp(0.0)
@RV v ~ Clamp(1.0)
@RV z ~ GaussianMeanVariance(m, v)
g2.variables
```




    Dict{Symbol,Variable} with 3 entries:
      :m => Variable(:m, Edges:…
      :v => Variable(:v, Edges:…
      :z => Variable(:z, Edges:…




```julia
ForneyLab.draw(g2)
```




The graph stores the identifier of each variable. This is useful because now we can retrieve a variable from a graph by its identifier, .e.g.,


```julia
g2.variables[:m]
```




    Variable(:m, Edges:
    Edge belonging to variable m: ( clamp_1.i[out] )----( gaussianmeanvariance_1.i[m] ).
    )



Let's build another simple factor graph for
$$\begin{align*}
p(x,y,z) &= p(z|x,y)\,p(x)\,p(y) \\
  &= \delta(z-x-y)\,\mathcal{N}(x\,|\,0.0,1.0)\,\mathcal{N}(y\,|\,2.0,3.0)
\end{align*}$$


```julia
g3 = FactorGraph()
@RV x ~ GaussianMeanVariance(0.0, 1.0)
@RV y ~ GaussianMeanVariance(2.0, 3.0)
@RV z = x + y
ForneyLab.draw(g3)
```



Next, we could be interested in inferring a property of `z`, e.g., the mean parameter for the distribution over `z`. This process can be automated by message passing in the graph. The next set of demo's will expand on the various methods for describing graphs and message passing inference methods with ForneyLab.
