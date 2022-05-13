# The Scotogenic Model of Neutrino Masses and Dark Matter in Warped Extra Dimensions
This is my master's thesis in physics. This repository includes the LaTeX source file
that compiles to the finished version of the thesis itself, as well as the software
and other files required for numerical calculations and visualizations contained
therein.

## Summary
This Master's Thesis has, as a subject of study, the embedding of the scotogenic model, which
explains the generation of neutrino masses and dark matter, in the context of the Randall-Sundrum
scenario, where an extra compactified spacelike dimension is added to the usual Minkowski space. We
shall start by briefly reviewing both of these known models as seen in the literature, before
combining them and start to explore some of the implications of such an embedding.

## File contents
`Memoria.tex` is the main course, the thesis itself. It's probably going to be
the most polished file of this repository, as people actually had to read this one.

`Notas.tex` is just that, the file I used to note explanations, demonstrations
and clarifications on concepts used in this work.

`TFM.bib` contains the references that were used or could have been used in both
of the previous documents.

`thesis.py` and `slides.py` make various numerical calculations and generate graphs
shown in the thesis and the defense's slides, respectively. Given that the
Python scripts were made and then extensively modified on short notice, they
are not in the best shape they could be. However, I will try to make them a bit
more friendly to anyone else that's not me when I find the time.

## Dependencies
### Memoria.tex
You will need LuaLaTeX, as well as the following packages:

- `geometry`
- `mathtools`
- `amssymb`
- `graphicx`
- `babel`
- `float`
- `subfigure`
- `pdfpages`
- `relsize`
- `biblatex`
- `csquotes`
- `dirtytalk`
- `multicol`
- `tikz-feynman`
- `slashed`
- `physics`
- `hyperref`

### thesis.py
You will need Python 3+, as well as the following libraries:

- NumPy
- SciPy
- Matplotlib

## Compiling
In order to compile the thesis' PDF, it's necessary to come up with the various
figures shown in it. So, first execute `thesis.py`,

```
$ chmod +x thesis.py
$ ./thesis.py
```

then compile the LaTeX to a PDF. I used LuaLaTeX since I make heavy use of
TikZ-Feynman, as this package is mostly written in Lua. To include the
bibliography, I used Biber, so whenever I'd compile the thesis, I'd execute
this,

```
$ lualatex Memoria.tex
$ biber
$ lualatex Memoria.tex
```

If references are not changed in the thesis itself, I believe it isn't necessary
to call Biber again like this.
