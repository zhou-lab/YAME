# docs/

The YAME documentation site: <https://zhou-lab.github.io/YAME/>

`index.html` is the site. It is a single self-contained page — styles, scripts
and the format diagrams (inline SVG) all live in the file, with `logo.png` as
the only asset. There is no build step, no theme, and no dependencies: edit the
file, open it in a browser, commit. The same arrangement as the other tools in
the suite (`kycg`, `sesame-cli`, `methscope-cli`, `cinderplot`, `tabl`).

`.github/workflows/pages.yml` uploads this directory to GitHub Pages verbatim
on any push to `main` that touches `docs/`.

## Conventions worth keeping

- **The page answers two questions, in two tables.** *Which format do I need*
  (ordered by what data the reader has, not by format number) and *what can I
  run on it* (18 commands in task groups, each with a seven-format strip).
  Everything else is one keystroke away inside a `<details>` drawer. Resist
  adding a third top-level section: the previous version had seventeen, one
  per converted Jekyll page, and nobody could find anything.
- **A row is a drawer.** The `<summary>` *is* the table row — a CSS grid, not
  a `<table>` — so a closed drawer costs one line and an open one keeps its
  place. No JS is needed to open them.
- **The filter must open what it matches.** A closed `<details>` still holds
  its text, so a filter that only hid sections would claim a hit the reader
  cannot see. There is exactly one `input` listener; keep it that way.
- **Flags are not documented here.** The page points at `yame <command> -h`,
  which is printed by the code that parses them and so cannot drift. A flag
  table would be a second source of truth that silently goes stale — the old
  Jekyll site documented several options that no longer existed.
- **Facts come from `src/`, not from the previous page.** Every format
  claim (packing width, NA representation, which `-f` letter) and every
  command's accepted formats were read out of the C. When behaviour changes,
  re-read the source rather than editing around the old sentence.
- **Terms are defined once**, in the glossary section. A dotted term elsewhere
  on the page is `<span class="tip" data-tip="SLUG">`; the script reads the
  definition out of `#gl-SLUG` at hover time, so a tooltip cannot disagree with
  the glossary entry it points at.
- **The page is theme-aware.** Colours come from CSS variables defined in four
  places — two `prefers-color-scheme` blocks and two `:root[data-theme=...]`
  overrides used by the toggle. A new variable must be added to all four, or
  toggling the theme against the OS preference leaves it at the wrong value.
- **Examples are runnable.** Every `<pre>` gets a copy button, so a block
  should be something a reader can paste, not a fragment.
