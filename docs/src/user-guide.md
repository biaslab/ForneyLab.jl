# User guide

This user guide describes the basic usage of ForneyLab for solving inference problems. The main content is divided in three parts:
- [Specifying a model](@ref)
- [Generating an algorithm](@ref)
- [Executing an algorithm](@ref)

## [Installation](@id user-guide-installation)
ForneyLab can be installed using the Julia package manager. From the Julia REPL, type ] to enter the Pkg REPL mode and run

```julia
pkg> add ForneyLab
```

!!! note
    If you want to use the graph visualization functions, you need to have [GraphViz](http://www.graphviz.org/) installed.

To import ForneyLab into the active Julia session run  

```@example 1
using ForneyLab
```

## Specifying a model

Probabilistic models incorporate the element of randomness to describe an event or phenomenon. They do this by using random variables and probability distributions. These models can be represented diagrammatically using probabilistic graphical models (PGMs). A factor graph is a type of PGM that is well suited to cast inference tasks in terms of graphical manipulations.

### Creating a new factor graph

Factor graphs are represented by the `FactorGraph` composite type (struct). To instantiate a new (empty) factor graph we use its constructor function that takes no arguments
```@example 1
g = FactorGraph()
nothing # hide
```

ForneyLab keeps track of an *active* factor graph at all times. Future operations on the graph, such as adding variables or nodes, affect the active graph. The call to `FactorGraph()` creates a new instance of a factor graph and registers it as the active graph.

To get the active graph run
```@example 1
fg = currentGraph()
nothing # hide
```

To set the active graph run
```@example 1
setCurrentGraph(g)
nothing # hide
```

### Adding random variables (edges)

Random variables are represented as edges on Forney-style factor graphs. You can add a random variable to the active factor graph by instantiating the `Variable` composite type. The constructor function takes an `id` of type `Symbol` as argument. For example, running
```@example 1
x = Variable(id=:x)
nothing # hide
```
associates the variable `x` to an edge of the active factor graph.

Alternatively, the `@RV` macro can be used for the same purpose in a more compact form. Executing the following line is equivalent to the previous approach.

```@example 1
@RV x
nothing # hide
```
By default, the `@RV` macro uses the variable's name to create a `Symbol` that is assigned to the `id` field of the `Variable` object (`:x` in this example). However, if a variable with that `id` already exists, then ForneyLab will create a name of the form `variable_x`, where `x` is a number that increments. In case you want to provide a custom `id`, the `@RV` macro accepts an optional argument between square brackets for this purpose. For example,
```@example 1
@RV [id=:my_id] x
nothing # hide
```
adds the variable `x` to the active graph and assigns the `:my_id` symbol to its `id` field. Later we will see that this is useful once we start visualizing the factor graph.

### Adding factor nodes
Factor nodes are used to define the relationship between different random variables. They assign probability distributions to a random variable as a function of other variables. See [Factor nodes](@ref) for a complete list of the available factor nodes in ForneyLab.

We can assign a probability distribution to a random variable using the `~` operator together with the `@RV` macro. For example, to create a Gaussian distributed random variable `y`, where its mean and variance are controlled by the random variables `m` and `v` respectively, we could run
```@example 1
@RV m
@RV v
@RV y ~ GaussianMeanVariance(m, v)
nothing # hide
```

### Visualizing a factor graph

ForneyLab provides a `draw` function to visualize a factor graph. It takes a `FactorGraph` object as argument. Let's visualize the factor graph that we defined in the previous section.
```@example 1
ForneyLab.draw(g)
```
Edges that are not connected to any factor node are not drawn.

### Clamping

Suppose that we know that the variance of the random variable `y`, of the previous model, is fixed to a certain value. ForneyLab provides a special factor node to impose this kind of constraint called a `Clamp`. Clamp factor nodes can be implicitly defined by using literals like in the following example
```@example 1
g = FactorGraph() # create a new factor graph
@RV m
@RV y ~ GaussianMeanVariance(m, 1.0)
ForneyLab.draw(g)
```
Here, the literal `1.0` passed as the second argument to `GaussianMeanVariance` creates a clamp node implicitly. Clamp factor nodes are visualized with a gray background.

On the other hand, if you want to assign a custom `id` to the `Clamp` factor node, then you have to instantiate them explicitly, i.e.
```@example 1
g = FactorGraph() # create a new factor graph
@RV m
@RV v ~ Clamp(1.0)
@RV y ~ GaussianMeanVariance(m, v)
ForneyLab.draw(g)
```

### Placeholders

Placeholders are a kind of `Clamp` factor nodes that act as entry points for data. They associate a given random variable with a buffer through which data is fed at a later point. This buffer has an `id`, a dimensionality and a data type. Placeholders are created with the `placeholder` function. Suppose that we observed a series of one-dimensional floating-point data points that we plan to feed to the `y` random variable of the model of the previous section. We would then need to define `y` as a placeholder as follows.
```@example 1
g = FactorGraph() # create a new factor graph
@RV m
@RV v ~ Clamp(1.0)
@RV y ~ GaussianMeanVariance(m, v)
placeholder(y, :y)
ForneyLab.draw(g)
```
Placeholders default to one-dimensional floating-point data. In case we want to override this with, for example, 3-dimensional integer data, we would need to specify the `dims` and `datatpye` parameters of the `placeholder` function as follows
```julia
placeholder(y, :y, dims=(3,), datatype=Int)
```
In section [Executing an algorithm](@ref) we will see how the data is fed to the placeholders.


### Overloaded operators

ForneyLab supports the use of the `+`, `-` and `*` operators between random variables that have certain types of probability distributions. This is known as *operator overloading*. These operators are internally modelled as a factor nodes in ForneyLab. As an example, a two-component Gaussian mixture can be defined as follows  

```@example 1
g = FactorGraph() # create a new factor graph
@RV x ~ GaussianMeanVariance(0.0, 1.0)
@RV y ~ GaussianMeanVariance(2.0, 3.0)
@RV z = x + y
placeholder(z, :z)
ForneyLab.draw(g)
```

### Online vs. offline learning

Online learning involves processing observations one at a time. In a Bayesian setting, this reduces to applying Bayes rule in a recursive fashion, i.e. the posterior distribution for a given random variable, becomes the prior for the next processing step. Since we are feeding one observation at each time step, the factor graph will have *one* placeholder for every observed variable. All of the factor graphs that we have seen so far, were specified to process data in this fashion. Let's look at a simple example in order to contrast it with its offline counterpart.
```@example 1
g = FactorGraph() # create a new factor graph
@RV x ~ GaussianMeanVariance(0.0, 1.0)
@RV y ~ GaussianMeanVariance(x, 1.0)
placeholder(y, :y)
ForneyLab.draw(g)
```
Just like we have seen before, there is one placeholder that accepts one observation at a time.

On the other hand, offline learning involves processing a batch of samples at every time step. This translates in a factor graph that has one placeholder linked to a random variable for *each* sample in the batch. We can specify this type of models using a `for` loop like in the following example.
```@example 1
g = FactorGraph()   # create a new factor graph
N = 3               # number of observations
y = Vector{Variable}(undef, N)
@RV x ~ GaussianMeanVariance(0.0, 1.0)
for i = 1:N
    @RV y[i] ~ GaussianMeanVariance(x, 1.0)
    placeholder(y[i], :y, index=i)
end
ForneyLab.draw(g)
```
The important thing to note here is that we need an array of `N` observed random variables in which every element is linked to a dedicated index of the placeholder's buffer. This buffer can be thought of as a `N` dimensional array of `Clamp` factor nodes. We achieve this link by means of the `index` parameter of the `placeholder` function.

In section [Executing an algorithm](@ref) we will see examples of how the data is fed to the placeholders in each of this scenarios.

!!! note
    Batch processing does not perform well with large datasets at the moment. We are working on this issue.


## Generating an algorithm

ForneyLab supports code generation for three different types of message-passing algorithms:
- [Belief propagation](https://en.wikipedia.org/wiki/Belief_propagation)
- [Variational message passing](https://en.wikipedia.org/wiki/Variational_message_passing)
- [Expectation propagation](https://en.wikipedia.org/wiki/Expectation_propagation)

Whereas belief propagation computes exact inference for the random variables of interest, variational message passing (VMP) and expectation propagation (EP) algorithms are approximation methods that can be applied to a larger range of models.

### Exact inference

The way to instruct ForneyLab to generate a belief propagation algorithm (also known as a sum-product algorithm) is by using the `sumProductAlgorithm` function. This function takes as argument(s) the random variable(s) for which we want to infer the posterior distribution. As an example, consider the following hierarchical model where the mean of a Gaussian distribution is represented by another Gaussian distribution.  
```@example 1
g = FactorGraph() # create a new factor graph
@RV m2 ~ GaussianMeanVariance(0.0, 1.0)
@RV m1 ~ GaussianMeanVariance(m2, 1.0)
@RV y ~ GaussianMeanVariance(m1, 1.0)
placeholder(y, :y)
ForneyLab.draw(g)
```
If we were only interested in inferring the posterior distribution of `m1` then we would run
```@example 1
algorithm_string = sumProductAlgorithm(m1)
nothing # hide
```
On the other hand, if we were interested in the posterior distributions of both `m1` and `m2` we would then need to pass them as elements of an array, i.e.
```@example 1
algorithm_string = sumProductAlgorithm([m1, m2])
nothing # hide
```

Note that the message-passing algorithm returned by the `sumProductAlgorithm` function is a `String` that contains the definition of a Julia function. In order to be able to execute this function, we first need to parse this string as Julia expression to then evaluate it in the current scope that the program is running on. This can be done as follows
```@example 1
algorithm_expr = Meta.parse(algorithm_string)
nothing # hide
```
```julia
:(function step!(data::Dict, marginals::Dict=Dict(), messages::Vector{Message}=Array{Message}(undef, 4))
      #= none:3 =#
      messages[1] = ruleSPGaussianMeanVarianceOutVPP(nothing, Message(Univariate, PointMass, m=0.0), Message(Univariate, PointMass, m=1.0))
      #= none:4 =#
      messages[2] = ruleSPGaussianMeanVarianceOutVGP(nothing, messages[1], Message(Univariate, PointMass, m=1.0))
      #= none:5 =#
      messages[3] = ruleSPGaussianMeanVarianceMPVP(Message(Univariate, PointMass, m=data[:y]), nothing, Message(Univariate, PointMass, m=1.0))
      #= none:6 =#
      messages[4] = ruleSPGaussianMeanVarianceMGVP(messages[3], nothing, Message(Univariate, PointMass, m=1.0))
      #= none:8 =#
      marginals[:m1] = (messages[2]).dist * (messages[3]).dist
      #= none:9 =#
      marginals[:m2] = (messages[1]).dist * (messages[4]).dist
      #= none:11 =#
      return marginals
  end)
```

```@example 1
eval(algorithm_expr)
```

At this point a new function named `step!` becomes available in the current scope. This function contains a message-passing algorithm that infers `m1` and `m2` given one or more `y` observations. In the section [Executing an algorithm](@ref) we will see how this function is used.

### Approximate inference

Variational message passing (VMP) and expectation propagation (EP) algorithms are generated much in the same way as the belief propagation algorithm we saw earlier. There are several differences though. Let's look at an example of VMP. In this example we will generate a VMP algorithm that infers both the mean and the precision of a Gaussian distributed variable. The main difference with respect to belief propagation is that for VMP we need to define a recognition factorization function that dictates our assumption of how the joint probability of our model factorizes. This assumption is what makes VMP an approximation method. A common approach is to assume that all random variables of interest (also known as hidden variables) factorize with respect each other. This is known as the *mean field* assumption.

```@example 1
g = FactorGraph() # create a new factor graph
@RV m ~ GaussianMeanVariance(0, 10)
@RV w ~ Gamma(0.1, 0.1)
@RV y ~ GaussianMeanPrecision(m, w)
placeholder(y, :y)
draw(g)
```


#### Computing free energy


## Executing an algorithm



Mention Online vs. batch processing

### Online processing

### Batch processing


!!! note
    Batch processing does not perform well with large datasets at the moment. We are working on this issue.
