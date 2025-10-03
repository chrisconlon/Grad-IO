Shared preamble for Grad-IO slides

Purpose
- Central canonical preamble: `resources/preamble.tex`.
- Use a single guarded, portable preamble so all lecture slides share look-and-feel.

How to use
- Instead of `\documentclass{beamer}` put at the top of your .tex file:

  - If you need custom beamer class options (handout, notes=show, xcolor, aspectratio, etc.) set them first:

    \def\beamerclassoptions{[notes=show,aspectratio=169]}

  - Then input the canonical preamble (relative path from the file):

    \input{../resources/preamble.tex}

  - After that, write `\begin{document}` and the rest of the slide content as normal.

Notes and conventions
- The canonical preamble uses guarded definitions (\providecommand, \@ifundefined,
  \IfFileExists) to avoid redefinition errors when multiple files include it.
- The preamble adds a flexible `\graphicspath` that includes `./` and several `../resources/` locations so `\includegraphics` works regardless of nesting depth.
- Prefer relative paths for portability inside the repo. Adjust the `..` count based on file depth.

If you need to keep a local package or macro that must be loaded before the preamble, set `\beamerclassoptions` (if needed) and then `\input` the preamble; local macros can be declared after the input but consider using `\providecommand` to avoid redefinition warnings.

Troubleshooting
- "Two \documentclass" errors: ensure you removed any `\documentclass{beamer}` lines; use the wrapper pattern above instead.
- "Command already defined" warnings: use `\providecommand` or `\@ifundefined` in local files when adding macros that may also exist in `resources/preamble.tex`.

Contact
- If uncertain about the correct relative path for `\input{.../resources/preamble.tex}`, ask the repository maintainer or open an issue.
