# The Scotogenic Model of Neutrino Masses and Dark Matter in Warped Extra Dimensions
## Abstract
This master's thesis has, as a subject of study, the embedding of the scotogenic model, which
explains the generation of neutrino masses and dark matter, in the context of the Randall-Sundrum
scenario, where an extra compactified spacelike dimension is added to the usual Minkowski space. We
shall start by briefly reviewing both of these known models as seen in the literature, before
combining them and start to explore some of the implications of such an embedding.

ELI12: What if we added a few new particles as well as an extra (spatial) dimension?

## File contents
`Memoria.tex` is the main course, the thesis itself. It's probably going to be
the most polished file of this repository, as people actually had to read this one.

`Notas.tex` is just that, the file I used to note explanations, demonstrations
and clarifications on concepts used in this work.

`TFM.bib` contains the references that were used or could have been used in both
of the previous documents.

`thesis.py` and `slides.py` make various numerical calculations and generate graphs
shown in the thesis and the defense's slides, respectively.

## Dependencies
### `Memoria.tex`
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

### `thesis.py`
You will need Python 3+, and the following libraries:

- NumPy
- SciPy
- Matplotlib

## Compiling
In order to compile the thesis' PDF, it's necessary to come up with the various
figures shown in it. So, first execute `thesis.py`
(these commands are for Linux/BSD/Mac systems, but it's trivial to do the
same on Windows)

```
$ cd masters
$ chmod +x thesis.py
$ ./thesis.py
```

then compile the LaTeX to a PDF. You need to use LuaLaTeX in order to do this,
since I made heavy use of [TikZ-Feynman](https://jpellis.me/projects/tikz-feynman/),
and this package is mostly written in Lua. To include the bibliography, I used
Biber, so whenever I'd compile the thesis, I'd execute this,

```
$ lualatex Memoria.tex
$ biber
$ lualatex Memoria.tex
```

If references are not changed in `Memoria.tex`, I believe it isn't necessary
to call Biber and LuaLaTeX again like this.

## Future plans
Now, the thesis itself is finished, so I will not add any new actual content
(i.e. science). Further research has its place in additional papers down the
line, however I decided that, in the near future, I will not pursue a career in
academia. If you wish to follow work on this or similar topics, I suggest you
keep an eye on what my supervisors will publish.

I am content with its formatting and the grade I got, however I still feel like
I could do some tinkering with the delivery. By this, I mean that I'm very
intrigued by [Org-mode](https://orgmode.org/index.html), and I might just end
up reformatting this thesis into an org file with all the code I used in
chunks. I also plan to clean up the scripts I used, given that they were made
and then extensively modified on short notice, so they are not in the best
shape they could be.
