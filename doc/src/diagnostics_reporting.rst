***************************
 Diagnostics and reporting
***************************

Several helper functions for inspecting constructed objects are available. This can be useful for debugging or inpecting (large) graphs.

Inspect objects
---------------

.. function:: show()

    Can be called on any object to display a summary of that instance and sometimes provides more detailed show suggestions.


Draw graphs and subgraphs
-------------------------

The `draw` function visualizes the argument object. Messages that are present at the moment of rendering are denoted with filled black circles. Setting the optional argument `external_viewer=true` renders the drawing in an external viewer. A related helper function `drawPdf(*, filename)` behaves similar to `draw` but writes the visualization to a pdf file.

If you want to use graph visualizations, GraphViz needs to be installed. 

.. function:: draw(::FactorGraph)

	Visualizes the argument graph. A dotted green edge indicates a `wrap`.

.. function:: draw(::Subraph)

	This method is only relevant when working with variational message passing and visualizes a subgraph. A dashed red edge indicates an external edge from the perspective of the argument subgraph.

.. function:: draw(::CompositeNode)

	Draws the internal graph of a composite node.

.. function:: draw(::Set{Node})

	Visualizes a set of nodes and all edges inbetween.


Display implemented update rules
--------------------------------

The function `rules` is used to inspect the message computation rules that are implemented in ForneyLab.

.. function:: rules(::Node)

	Tabulates all message update rules for an argument node, together with references. Passing an optional argument `format=:list` shows the actual update functions.


Verify algorithm execution
--------------------------

Sometimes you want an overview of all the computation rules that were called. This can be accomplished by calling `setVerbosity` before algorithm execution. ForneyLab now prints a table of all update calls.

.. function:: setVerbosity(::Bool)

	Sets the verbosity.