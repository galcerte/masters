# The Scotogenic Model of Neutrino Masses and Dark Matter in Warped Extra Dimensions
This is my master's thesis in physics. This repository includes the LaTeX source file
that compiles to the finished version of the thesis itself, as well as the software
and other files required for numerical calculations and visualizations.

Given that the Python scripts were made and then extensively modified on short
notice, they are not in the best shape they could be. However, I will try to
make them a bit more friendly to anyone else that's not me when I find the
time.

`thesis.py` deals with what is shown in the thesis itself, and `slides.py` generates
various graphs exclusive to its presentation.

## Compiling everything
In order to compile the thesis' PDF, it's necessary to come up with the various
figures shown in it. So, first execute `thesis.py`,

```
$ chmod +x thesis.py
$ ./thesis.py
```

then compile the LaTeX to a PDF. I used LuaLaTeX given that I make heavy use
of TikZ-Feynman, as this package is also written in Lua. To include the bibliography,
I used Biber, so whenever I'd compile the thesis, I'd execute this,

```
$ lualatex Memoria.tex
$ biber
$ lualatex Memoria.tex
```

If references are not changed in the thesis itself, I believe it wasn't necessary
to call Biber again like this.
