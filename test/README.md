# The test suite

**Status, 2026-09-11.** The plan below is mostly implemented; this file is its
permanent home, moved out of `tmp/` once the work landed. What exists:

- `test/run.sh` runs every `t_*.sh`; `make test`, and the conda build on both
  platforms, run it. 14 scripts; every one of the 18 subcommands is exercised.
- The three-configuration matrix (release / `-DNDEBUG` / ASan+UBSan) runs in
  CI as the `test` job. `t_valgrind.sh` covers the leak and overrun class.
- `probe.c` + `t_probe.sh` are the library-contract layer.
- `t_corrupt.sh` is the adversarial layer.
- `t_store_info.sh` is the first layer-5 member: real fixtures through the
  registry, skipped cleanly without `YAME_DATA_HOME`.
- CI job `api-surface` diffs the transitive public header set (everything
  reachable from `cdata.h`/`cfile.h`) against the last tag and fails unless
  a commit names the changed header. Its first dry run caught a real break:
  a macro added to `wzmisc.h` collides with a same-named function in
  methscope, so the next methscope YAME bump would not compile.
- Line coverage is measured by `scripts/coverage.sh` and gated in CI against
  the committed badge in either direction. 57.2% as of this date; `format5.c`
  at 0% is correct (obsolete, decodable but not packable).
- The methscope docs harness now scores every command in a block (`set -e`),
  not only the last -- the gap a zero-byte `upscale` shipped through.
- `t_multimask.sh` + `probe_multi.c` / `probe_runs.c` / `probe_index.c` cover
  the summary kernels (the walk and the inverted index, `summary -I`), and
  `t_mrmp_fixture.sh` checks runs against records on a `.cm` methscope
  exported from a `.mrmp`, the one input this repo cannot generate.

Still open: the remaining per-subcommand depth (`hprint` 22%, `format6` 28%,
`rowop` 35%, `format4` 41%, `format7` 40%), and any further layer-5 members.
The seed-regression table at the end is the reason each layer exists.

---

# plan — a test suite for YAME

**Type:** plan. **Date:** 2026-09-11. **Status:** proposed; layer 2 has one
member (`t_pairwise.sh`) and the runner exists.

## goal, and the constraints that shape it

Catch the class of bug that reached us through consumers this month, without
making the repo heavier. Every regression since v1.40 was found by methscope,
the docs harness, or a hand comparison against a second binary -- never by
YAME itself. The suite has to close that gap under four constraints:

- **No framework, no dependencies.** Bash, POSIX awk, and the binary.
- **No committed data.** Fixtures are packed from inline text at run time.
  Anything large comes through the registry and is gitignored.
- **Runs on both CI legs** (linux-64, osx-arm64) inside the conda build, so
  a wrong answer fails the package and the publish job holds the version.
- **Small tree.** One runner, one script per subcommand, one C probe.

## what exists today

| piece                          | what it checks                                      |
| ------------------------------ | --------------------------------------------------- |
| `test/run.sh`                  | runs every `test/t_*.sh`, tallies, exits non-zero   |
| `test/t_pairwise.sh`           | set output, -S summary, names, broadcast, refusals  |
| `make test`                    | builds if needed, runs the runner                   |
| `conda-recipe/build.sh`        | runs the suite after the build, before install      |
| `conda-recipe/meta.yaml` test  | banner version, registry loads, fetch reports ready |
| labjournal docs harness        | runs the methscope docs page blocks (see layer 5)   |

`zhou-lab/YAME_test` (retired July 2026) had no assertions; nothing from it
is worth recovering.

## layout

```
test/
  run.sh              runner (exists)
  t_<subcommand>.sh   one per subcommand, self-contained, layer 1+2
  t_corrupt.sh        layer 3, adversarial input
  probe.c             layer 4, built by run.sh against libyame.a
  t_probe.sh          runs the probe
  t_store_*.sh        layer 5, skipped unless YAME_DATA_HOME is set
  data/               gitignored; registry fetch target for layer 5
  fixtures/           committed inputs a test cannot make itself (small):
                      one real mask; mrmp_export/, a .cm exported from a
                      .mrmp by methscope with its expected means (7 KB),
                      driven by t_mrmp_fixture.sh + probe_fixture.c
```

Each `t_*.sh`: `set -euo pipefail`, `mktemp -d` with a trap, fixtures from
heredocs, expected values recomputed in awk from the same text (an
independent oracle), goldens only as `unpack` text and never as binary. Exit
codes and empty-stdout-on-refusal are assertions, not afterthoughts.

## the five layers

### 1. format round-trips (in each `t_<cmd>.sh`, shared idiom)

For every format a subcommand dispatches on: `pack` then `unpack` returns
the input text; a `subset`/`split`/`cat`/`index -s` cycle returns a store
where `subset` by name gives the original record; `chunk` then `cat` is the
original. Most of the fifteen code-review bugs sat in format branches no
caller reached; a format grid is what reaches them.

Formats to grid: 0, 1, 2, 3, 4, 6, 7. Format 4 (beta, NA negative) and 7
(coordinates) are the ones with the least coverage anywhere today.

### 2. per-subcommand behaviour

| script            | must assert (beyond the round-trip)                                   |
| ----------------- | --------------------------------------------------------------------- |
| `t_pack.sh`       | every `-f`; `-f c` refuses a multi-char line (was: silent truncation) |
| `t_unpack.sh`     | `-f 0/1/2` on fmt3; `-c` at n == k*chunk exits 0 (was: fatal after)   |
| `t_index.sh`      | `-1` per concatenated sample gives distinct offsets; `-s`; dup name   |
| `t_subset.sh`     | by name, raw path and `-z` path identical; missing name leaves no file|
| `t_split.sh`      | one file per record, names from index                                 |
| `t_info.sh`       | NSample/Nrow/Format per record; NA without index                      |
| `t_rowsub.sh`     | `-B`, `-I`, `-m`; out-of-range beg fatal; `-B 4294967306_…` refused   |
| `t_chunk.sh`      | exact multiple of chunk size; `cat` of chunks == input                |
| `t_summary.sh`    | against a mask; fmt6 set/meth flavours                                |
| `t_pairwise.sh`   | exists                                                                |
| `t_binarize.sh`   | threshold edge; output universe == covered sites                      |
| `t_mask.sh`       | default vs `-v`; fmt3 mask means coverage > 0                         |
| `t_dsample.sh`    | fixed seed reproducible; n == k*8 rows (was: heap overrun)            |
| `t_perturb.sh`    | fixed seed reproducible; flip count within tolerance                  |
| `t_rowop.sh`      | binasum, musum, stat columns vs awk; cometh refuses non-fmt3          |
| `t_hprint.sh`     | region and full-dataset on fmt0/3/4/6; alphabet unchanged             |
| `t_fetch.sh`      | `-l` lists; a named single file needs no `-y`; a dir refuses without  |
| `t_chunkchar.sh`  | text chunking, line boundaries                                        |

Each row that names a past bug is a regression test first.

### 3. adversarial input (`t_corrupt.sh`)

The contract is one sentence: **every input either yields the right answer
or exits non-zero with a message; exit 0 with wrong data is the only
forbidden outcome.** Cases, each run through `info`, `unpack` and `subset`:

- truncate a 4-record store at every interior byte (the map from
  2026-09-10: today 369 fatal, 7 clean at member boundaries -- assert that
  map holds, so a change in the reader's leniency is a visible diff)
- 200 bytes of noise; an empty file; a text file
- a store followed by raw bytes (the bundle shape) -- `info` must fail,
  a bounded read through the probe must succeed
- two end markers in a row (the concatenation trap)
- a record whose header promises more bytes than follow
- an `.idx` pointing past end of file, and one whose offsets are all equal

`cx_read_record` is one entry point, so this layer is also where a fuzzer
plugs in later (afl or libFuzzer over the reader) without new structure.

### 4. library contract (`probe.c` + `t_probe.sh`)

One C file linked with `$(yame-config --libs)`, exercising what consumers
call: `open_cfile`, `read_cdata1`, `cx_read_record` with and without a
limit, `free_cdata`, `decompress_in_situ`, `f3_get_mu`. Assertions are the
documented semantics: END at the limit, NOT_CX past it, `n == 0` after
`free_cdata`. The `mask->n`-after-free breakage was a semantics change nobody
tested for; this is the test.

Plus a CI step that diffs `src/cdata.h` and `src/cfile.h` against the last
tag and fails on a changed signature unless the commit message names it.
methscope's `RELEASE.md` does that diff by hand today.

### 5. integration on real data (`t_store_*.sh`)

Fetch through the registry into `test/data/`, skip cleanly when
`YAME_DATA_HOME` is unset, run on release tags rather than every push.
First members: `info` and `summary` on `hg38/data/human_hg38_test.cg`; a
bounded read of `hg38/models/hg38_10k1.updecx`. Also fix the docs
harness (`20260824_docs_runthrough.py`) to score every command in a block,
not the last one -- that scoring is how a zero-byte `upscale` shipped.

## build matrix

Run layers 1–4 under three configurations. Both legs build natively, so this
is three `make` invocations in CI, not new tooling:

```sh
make clean && make -j && make test                          # release flags
make clean && make -j CC="cc -DNDEBUG" && make test         # what conda ships
make clean && make -j CC="cc -fsanitize=address -g -O1" \
  && ASAN_OPTIONS=detect_leaks=1 make test                  # memory errors
```

The assert-side-effect bug existed only in the second; the `dsample` overrun
and the fmt1 overread would have failed the third on first run.

## rules

- Every bug fix ships with the test that would have caught it.
- Independent oracle (awk/python recompute) over stored goldens; goldens as
  `unpack` text when unavoidable; never binary.
- POSIX only: no `md5sum`, no GNU awk, no `readlink -f`.
- A test that cannot fail is not a test: point the suite at the previous
  release once and confirm the new test goes red.
- Nothing under `test/` but scripts and one C file. Data goes through the
  registry.

## order of work

1. **`t_corrupt.sh`** (layer 3) and the **NDEBUG matrix leg**. Together
   these would have caught every regression since v1.40. Half a day.
2. **`probe.c`** (layer 4), seeded from the 2026-09-10 bounded-read probe.
   An hour.
3. **`t_pack.sh`, `t_unpack.sh`, `t_index.sh`, `t_subset.sh`** -- the four
   subcommands every pipeline passes through, with the format grid. A day.
4. The remaining `t_<cmd>.sh`, one per sitting, regression rows first.
5. ASan leg, header-diff gate, layer 5 and the docs-harness fix.

## seed regressions, and the layer that owns each

| bug (2026-08/09)                                  | layer | script          |
| ------------------------------------------------- | ----- | --------------- |
| `unpack -c` fatal at exact chunk multiple          | 2     | t_unpack        |
| `unpack -c` leaked a record per chunk              | 3 mtx | ASan leg        |
| `pack -f c` kept first character silently          | 2     | t_pack          |
| corrupt/truncated read as empty store, rc 0        | 3     | t_corrupt       |
| `float_t` vs 4-byte contract                       | 1     | fmt4 round-trip |
| `dsample` heap overrun at n == k*8                 | 3 mtx | ASan + t_dsample|
| `cdata_nbytes` ignored unit                        | 4     | probe           |
| `free_cdata` double free on aux                    | 4     | probe           |
| FMT6 macro index unparenthesised                   | 1     | fmt6 round-trip |
| `rowsub` parsed uint64 with atoi                   | 2     | t_rowsub        |
| stale unit/aux across records                      | 1     | mixed-fmt store |
| unaligned uint16 stores                            | 3 mtx | ASan (UBSan)    |
| fmt1 overread on short stream                      | 3     | t_corrupt       |
| `subset` left a zero-record file on bad name       | 2     | t_subset        |
| assert side effects under -DNDEBUG (3 sites)       | mtx   | NDEBUG leg      |
| bundle prefix read fatal (updecx)                  | 3+4   | t_corrupt+probe |
| `mask->n` read after `free_cdata` (methscope)      | 4     | probe           |
| `index -1` constant offsets (conda only)           | 2+mtx | t_index + NDEBUG|

## what not to do

No C unit-test framework and no line-coverage target. YAME's behaviour is
defined at the format and command boundary; the property grid buys more
correctness per line than mocking internals. Unit tests are worth it only
for pure helpers (`pack_value`, the `FMT6_*` macros, `MU2beta`), and those
fit in `probe.c`.
