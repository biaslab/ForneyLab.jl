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

The ``draw(object)`` function visualizes the argument object. Messages that are present at the moment of rendering are denoted with filled black circles. Setting the optional keyword argument ``external_viewer`` to one of the possible values (``:Default``, ``:iTerm``) will render the drawing in an external viewer of your choice. The related function ``drawPdf(object, filename)`` is similar to ``draw`` but writes the visualization to a pdf file.

If you want to use graph visualizations, the GraphViz package needs to be installed.

.. function:: draw(::FactorGraph)

	Visualizes the argument graph. A dotted green edge indicates a ``wrap``.

.. function:: draw(::Subgraph)

	This method is only relevant when working with variational message passing and visualizes a subgraph. A dashed red edge indicates an external edge from the perspective of the argument subgraph.

.. function:: draw(::CompositeNode)

	Draws the internal graph of a composite node.

.. function:: draw(::Set{Node})

	Visualizes a set of nodes and all edges inbetween.


Display implemented update rules
--------------------------------

The function ``rules`` is used to inspect the message computation rules that are implemented in ForneyLab.

.. function:: rules(::Node)

	Tabulates all message update rules for an argument node, together with references. Pass the optional keyword argument ``format=:list`` to show the actual update functions.


Verify algorithm execution
--------------------------

Sometimes you want an overview of all the computation rules that were executed. This can be accomplished by enabling the vebose mode using ``setVerbosity()`` before algorithm execution. ForneyLab will print a table of all update calls in verbose mode.

.. function:: setVerbosity(::Bool)

	Enable/disable the verbose mode.
