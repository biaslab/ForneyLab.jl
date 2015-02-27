ForneyLab.jl Documentation README
=================================

The documentation is written (just like Julia's documentation) in reStructuredText, a good reference for which
is the [Documenting Python](http://docs.python.org/devguide/documenting.html)
chapter of the Python Developer's Guide.


Prerequisites for building the documentation
--------------------------------------------

The documentation is built using [Sphinx](http://sphinx.pocoo.org/) and LaTeX.
On Ubuntu, you'll need the following packages installed:

    latex-cjk-all
    texlive
    texlive-lang-cjk
    texlive-latex-extra

On OS X, you can install install MacTex using the GUI installer


Building the documentation
--------------------------

Build the documentation by running

    $ sudo make deps
    $ make helpdb.jl
    $ make html
    $ make latexpdf

Sphinx extensions and theme
---------------------------
The extensions to Sphinx and the theme are in the
https://github.com/JuliaLang/JuliaDoc repository.
