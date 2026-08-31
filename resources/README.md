# Shared style for Grad-IO slides

All decks use the shared package `resources/teaching_slides.sty` (Metropolis
theme, course colors/fonts, bibliography setup, common macros). The old
`\input{../resources/preamble.tex}` pattern is retired.

## How to start a new deck

```latex
\documentclass[aspectratio=169,11pt]{beamer}
\usepackage{teaching_slides}

\title{...}
\author{Chris Conlon}
\institute{Grad IO}
\date{Fall 2026}

\begin{document}
...
```

Class options (`handout`, `notes=show`, font size, aspect ratio) go on the
`\documentclass` line as usual.

## How the package is found

`teaching_slides.sty` (and `fixpauseincludegraphics.sty`) live in this
directory; the copy in the personal texmf tree
(`~/Library/texmf/tex/latex/...`) is a **symlink** to the repo file, so the
repo is the single source of truth — edit it here only. On a new machine,
recreate the symlinks (or add this directory to `TEXINPUTS`):

```sh
mkdir -p ~/Library/texmf/tex/latex/teaching_slides
ln -sf "$(pwd)/resources/teaching_slides.sty" ~/Library/texmf/tex/latex/teaching_slides/
mkdir -p ~/Library/texmf/tex/latex/fixpauseincludegraphics
ln -sf "$(pwd)/resources/fixpauseincludegraphics.sty" ~/Library/texmf/tex/latex/fixpauseincludegraphics/
```

## What the package provides

- Metropolis theme configuration, course palette, Fira/TeX Gyre font fallbacks.
- `\graphicspath` covering `./resources/` at several nesting depths.
- apacite/natbib bibliography setup with gray citations and `\citepos`.
- Math macros: `\E`, `\Var`, `\Cov`, `\Pr` (renders $\mathbb{P}$), `\prob{...}`,
  `\argmax`, `\argmin`, `\abs`, `\norm`, `\overbar`; script letters `\calS`,
  `\calD`, `\calJ`, `\calF`, `\calH`.
- Slide helpers: `\cmark`/`\xmark`, `\goto`/`\goback`/`\buttons`,
  `wideitemize`/`wideenumerate`; `listings` configured for R.

Notation conventions for the course live in `NOTATION.md` at the repo root.

## Building

`./build_all.sh` builds every deck that uses `teaching_slides` with
latexmk + XeLaTeX. To build one deck:

```sh
cd "Week X- .../" && latexmk -pdf -pdflatex="xelatex -interaction=nonstopmode" deck.tex
```
