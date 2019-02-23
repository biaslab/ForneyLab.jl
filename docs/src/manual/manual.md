# Manual

## Creating a new factor graph

## Adding `Variable`s to a graph

### Using the `Variable` constructor

### Using the `@RV` macro

#### Custom ids
`@RV[...]`

## Adding factor nodes
Factor nodes are used to define the relationship between different random variables. They assign probability distributions to a random variable as a function of other random variables. See [Factor nodes](@ref) for a complete list of the available factor nodes.

## Drawing a factor graph

## Overloaded operators

## Online vs. batch processing

### Factor Graphs and Variables

A central concept in ForneyLab is the (random) `Variable` type. After including ForneyLab and indicating that we start a new `FactorGraph`, we can declare a `Variable` by calling its constructor function:

```@example 1
using ForneyLab

# Declare a new graph
g = FactorGraph()

# Declare a variable
x = Variable(id=:x)

g.variables
```

The call to `FactorGraph()` creates a factor graph type and registers the graph as the currently active graph. Note that the variable has been associated with an edge in the currently active graph.

ForneyLab comes equipped with the `@RV` macro to define random variables. For instance, defining a new variable `y` with identifier `:y` and associating the variable to the current graph can also be accomplished by executing `@RV y`:

```@example 1
@RV y
g.variables
```

We can assign a probability distribution to a random variable using the `~` operator:

```@example 1
@RV z ~ GaussianMeanVariance(0.0, 1.0)
g.variables
```

Note that the graph now has two new variables with ids `:clamp_1` and `:clamp_2`. These two variables correspond to the mean and variance parameters for the Gaussian and are clamped to values `0.0` and `1.0` respectively.

If you have [graphviz](https://www.graphviz.org/) installed, then you can draw the factor graph. (Edges (variables) that are not constrained by any factor are not drawn):

```@example 1
ForneyLab.draw(g) # draw the graph
```

In case you want to define your own ids in place of the automatically generated `:clamp_1` and `:clamp_2`, you can declare the parameters of the Gaussian distribution through the `@RV` macro and associate a `Clamp` distribution with these variables:

```@example 1
g2 = FactorGraph()
@RV m ~ Clamp(0.0)
@RV v ~ Clamp(1.0)
@RV z ~ GaussianMeanVariance(m, v)
g2.variables
```

```@example 1
ForneyLab.draw(g2)
```

The graph stores the identifier of each variable. This is useful because now we can retrieve a variable from a graph by its identifier, .e.g.,

```@example 1
g2.variables[:m]
```

Let's build another simple factor graph for

$$\begin{align*}
p(x,y,z) &= p(z|x,y)\,p(x)\,p(y) \\
  &= \delta(z-x-y)\,\mathcal{N}(x\,|\,0.0,1.0)\,\mathcal{N}(y\,|\,2.0,3.0)
\end{align*}$$

```@example 1
g3 = FactorGraph()
@RV x ~ GaussianMeanVariance(0.0, 1.0)
@RV y ~ GaussianMeanVariance(2.0, 3.0)
@RV z = x + y
ForneyLab.draw(g3)
```

Now suppose we are interested in inferring a property of `z`, e.g., the mean parameter for the distribution over `z`. This process can be automated by message passing in the graph. The next set of demo's will expand on the various methods for describing graphs and message passing inference methods with ForneyLab.
