# Conda recipe for yame

Builds the `yame` binary and publishes it to the **zhou-lab** channel.

## Install (end users)

```sh
conda install -c zhou-lab -c conda-forge yame
```

Reference data is **not** shipped in the package — it is fetched into the
user's own store on demand:

```sh
yame fetch          # browse the catalogue and choose
```

That split is deliberate. The catalogue runs to hundreds of megabytes and is
versioned independently of the software, so packaging it would make every data
update a software release. What the binary carries is the *registry*: the tag
and digest of every directory it can verify, so a download is checked against a
pin this build holds rather than against whatever the server serves.

The store is shared with the other tools that link `libyame.a`, so a file
fetched here is found by them without a second download.

## Two recipes, on purpose

| | this recipe | [bioconda](https://github.com/bioconda/bioconda-recipes/tree/master/recipes/yame) |
|---|---|---|
| source | the checked-out tree (`path: ..`) | the release tarball for a tag |
| publishes | a commit, as soon as we want it | a tag, after review and the bot |
| channel | `zhou-lab` | `bioconda` |

Keep the two in step where it matters — same dependencies, same tests. If they
differ, the difference should be deliberate and written down.

## Build & upload (maintainers)

Requires `conda-build` and `anaconda-client`:

```sh
conda build -c conda-forge conda-recipe/
anaconda upload -u zhou-lab "$(conda build -c conda-forge conda-recipe/ --output)"
```

In practice you do not run these by hand: pushing a `v*` tag builds every
platform and publishes them together (`.github/workflows/conda-build.yml`),
using the `ANACONDA_TOKEN` repository secret. Publishing is a separate job that
requires all platforms to have built, so a release cannot land for linux-64
alone.

Before tagging, `src/yame_version.h` and this recipe's `version` must agree
with the tag — CI checks all three, because a package labelled 1.29 shipping a
binary that says 1.28 is worse than a failed build.

## libcurl is required

`yame fetch`, the catalogue browser, `summary -b`, and resolving `-R` / `-m`
against the store are all built on libcurl. Since v1.29 the Makefile *fails*
without it rather than dropping the feature silently, so this recipe cannot
accidentally publish a yame that cannot download. The recipe test asserts the
working state positively (`fetch available`) rather than grepping for the
absence of a warning, which would pass on any output at all.
