#!/bin/bash
## Layer 5: every mask in the real hg38 knowledgebase, through every path of
## `summary`, against a dense methylome, a single cell and a binary-call store.
## Correctness is byte identity across the three paths; speed and memory are
## checked against a committed baseline.
##
## WHY. `summary` now has three ways to accumulate one table: the per-mask
## path summarize1() (the original arithmetic), the one-pass kernel (the query
## read once, masks 64 rows at a time), and the inverted index (-I). All three
## must print the same bytes, and the fast one must actually be the one that
## runs. Both failures have happened silently: a guard sent every mask down
## the slow path and all tests passed while 6 s became 107 s. Synthetic masks
## did not catch it either -- the kernel's first version was SLOWER than the
## loop on TFBS while 1.9x on synthetic masks. So this runs the real files:
## 32 masks, binary and state, one record and many (the 1,359-record files
## as their first 16), sparse (Blacklist) and claiming every row (state
## masks), against fmt3 dense, fmt3 sparse and fmt6.
##
## SPEED. Each run is timed and its peak RSS taken; both are compared with
## test/fixtures/kb_baseline.tsv. The tolerance is loose -- 2x the seconds
## plus half a second, 1.5x the memory plus 50 MB -- because release machines
## are shared and noisy, and the regressions this exists to catch are 5x to
## 18x, not 20%. YAME_KB_REBASE=1 rewrites the baseline from this run; do that
## on an idle machine, and say so in the commit.
##
## GROUND TRUTH. test/fixtures/kb_expected.tsv holds every table the walk
## printed, keyed by mask and query, without the two file-path columns (a
## table over 1,000 rows as its sha256 and every 200th row). The
## three paths agreeing with each other is not enough: a change they share
## (an accessor, a beta rounding) moves all three together. Every table must
## match its expected rows byte for byte; a miss prints the diff, which names
## the row and the column. The expected file is written with the baseline
## (YAME_KB_REBASE=1) and should change only with a deliberate change to what
## summary prints, said so in the commit. The single-cell query is a seeded
## dsample of the dense one, so a change to dsample's draw moves its 32
## tables too. The store files are dated and carry SHA256SUMS.
##
## ANOTHER BINARY. YAME_KB_REF=/path/to/an/older/yame runs that binary too, on
## every mask and query, and requires its table to be byte-identical to this
## one's; its seconds go into the table as <mask>/<query>/ref, so a release
## can be set beside the one before it. It is not compared with the baseline.
##
## Opt-in like the other layer-5 test: release checks, not every push. About
## two minutes with the store mounted; big masks run as 16-record subsets.
set -euo pipefail
export LC_ALL=C                     # sort and join must agree on the order
YAME=${YAME:?export YAME=/path/to/yame}
here=$(cd "$(dirname "$0")" && pwd)
if [ -z "${YAME_TEST_LAYER5:-}" ]; then
  echo "skip: layer 5 runs with YAME_TEST_LAYER5=1 (release checks)"; exit 0
fi
if [ -z "${YAME_DATA_HOME:-}" ]; then
  echo "skip: YAME_DATA_HOME unset (layer 5 needs the shared store)"; exit 0
fi
KB=$YAME_DATA_HOME/hg38/KYCG
DATA=$YAME_DATA_HOME/hg38/data
[ -d "$KB" ] && [ -f "$DATA/human_hg38_immune_mixture.cg" ] && [ -f "$DATA/human_hg38_test.cg" ] ||
  { echo "skip: hg38 knowledgebase or data not in $YAME_DATA_HOME"; exit 0; }
TIME=/usr/bin/time
"$TIME" -f '%e %M' true 2>/dev/null || { echo "skip: GNU time is needed for the timings"; exit 0; }
base=$here/fixtures/kb_baseline.tsv
expected=$here/fixtures/kb_expected.tsv

sha256() { if command -v sha256sum >/dev/null; then sha256sum "$1"; else shasum -a 256 "$1"; fi | cut -d' ' -f1; }

d=$(mktemp -d); trap 'rm -rf "$d"' EXIT
cd "$d"

## ---- the queries -------------------------------------------------------------
## dense:  a whole-genome fmt3 methylome, 4.2M covered rows
## cell:   the same record downsampled to 200K covered rows, seeded, so it is
##         the same rows every run -- the shape of a single cell
## calls:  a fmt6 binary-call store (set + universe), the other query kernel
cp "$DATA/human_hg38_immune_mixture.cg" dense.cg
"$YAME" dsample -s 7 -N 200000 dense.cg > cell.cg 2>/dev/null
cp "$DATA/human_hg38_test.cg" calls.cg
for q in dense cell calls; do
  [ -s $q.cg ] || { echo "could not prepare $q.cg"; exit 1; }
done

## ---- one run: time it, keep the table, keep the path counters -------------
## run <tag> <out> [env=val ...] -- <summary args>
## Where `run` finds the binary: YAME, or the reference when tag ends in /ref.
run() {
  local tag=$1 out=$2; shift 2
  local envs=() bin=$YAME
  while [ "$1" != "--" ]; do envs+=("$1"); shift; done; shift
  case "$tag" in */ref) bin=$YAME_KB_REF ;; esac
  env "${envs[@]}" YAME_SUMMARY_PATH=1 "$TIME" -f '%e %M' -o "$out.time" \
    "$bin" summary "$@" > "$out" 2> "$out.err" ||
    { echo "$tag: yame summary $* failed"; cat "$out.err"; exit 1; }
  ## (an older reference binary prints no counters: an empty .path is fine)
  { grep 'one-pass' "$out.err" || true; } | tail -1 | sed 's/.*one-pass //' > "$out.path"
  printf '%s\t%s\t%s\n' "$tag" "$(cut -d' ' -f1 "$out.time")" "$(cut -d' ' -f2 "$out.time")" >> measured.tsv
  ## remembered, so a run the baseline flags can be timed once more
  { printf '%s\t%s' "$tag" "$out"; for e in "${envs[@]}"; do printf '\t%s' "$e"; done
    printf '\t--'; for a in "$@"; do printf '\t%s' "$a"; done; printf '\n'; } >> cmds.tsv
}
paths() { cat "$1.path"; }        # "K, per-mask F, inverted-index I"

: > measured.tsv
: > got_expected.tsv
fails=0
bad() { echo "  FAIL: $*"; fails=$((fails + 1)); }

for m in "$KB"/*.cm; do
  name=$(basename "$m" .cm); name=${name%.*}      # TFBS.20220921 -> TFBS
  fmt=$("$YAME" info "$m" 2>/dev/null | sed -n 2p | cut -f5)
  recs=$("$YAME" info "$m" 2>/dev/null | tail -n +2 | grep -c .)
  [ "$recs" -gt 0 ] || { bad "$name: info lists no record"; continue; }

  ## A big file runs as its first 16 records, through subset (the .idx the big
  ## files have). That keeps every variety and every path while the whole
  ## thing stays near two minutes and a quarter of a gigabyte: the full TFBS
  ## alone is 7 s a walk, 22 s and 1 GB an index, and 3 minutes a per-mask
  ## pass, and a regression that shows on 1,359 records shows on 16 -- the
  ## 6 s -> 107 s fallback would read 0.3 s -> 2.4 s here. The full files are
  ## for a hand comparison (YAME_KB_REF) on a quiet machine, not the release.
  if [ "$recs" -gt 16 ]; then
    [ -f "$m.idx" ] || { echo "  note: $name has $recs records and no .idx; skipped"; continue; }
    "$YAME" subset -H 16 "$m" > sub_$name.cm 2>/dev/null ||
      { bad "$name: subset -H 16 failed"; continue; }
    m=$PWD/sub_$name.cm; recs=16
  fi

  for q in dense cell calls; do
    ## the three paths, byte-identical
    run "$name/$q/walk"    w_$name.$q -- -m "$m" $q.cg
    ## and the walk against the committed tables (columns 2 and 4 on: no
    ## paths). A table over 1,000 rows -- Win100k is 29,511 windows -- goes in
    ## as a hash of the whole plus every 200th row, so the file stays small
    ## and a miss still shows real rows.
    tail -n +2 w_$name.$q | cut -f2,4- > t_$name.$q
    if [ "$(grep -c . t_$name.$q)" -gt 1000 ]; then
      printf '%s/%s\t#sha256\t%s\t%s rows\n' "$name" "$q" "$(sha256 t_$name.$q)" "$(grep -c . t_$name.$q)" >> got_expected.tsv
      /usr/bin/awk -v k="$name/$q" 'NR % 200 == 1 { print k "\t" $0 }' t_$name.$q >> got_expected.tsv
    else
      sed "s|^|$name/$q	|" t_$name.$q >> got_expected.tsv
    fi
    run "$name/$q/permask" p_$name.$q YAME_SUMMARY_KERNEL=0 -- -m "$m" $q.cg
    cmp -s w_$name.$q p_$name.$q || bad "$name/$q: the kernel and the per-mask path differ"
    ## -I is for many records against many masks, so it runs on the
    ## multi-record masks and the format 3 queries, which is its scenario
    if [ "$q" != calls ] && [ "$recs" -gt 1 ]; then
      run "$name/$q/index" i_$name.$q -- -I -m "$m" $q.cg
      cmp -s w_$name.$q i_$name.$q || bad "$name/$q: -I and the walk differ"
    fi
    ## the reference binary, when given: the same bytes, whatever it costs
    if [ -n "${YAME_KB_REF:-}" ]; then
      run "$name/$q/ref" r_$name.$q -- -m "$m" $q.cg
      cmp -s w_$name.$q r_$name.$q ||
        { bad "$name/$q: $YAME_KB_REF prints a different table"
          diff w_$name.$q r_$name.$q | head -6 | sed 's/^/    /'; }
    fi
    ## the table's shape: one row per (record, state); universe = the genome
    rows=$(tail -n +2 w_$name.$q | grep -c .)
    if [ "$fmt" = 2 ]; then
      states=$("$YAME" info "$m" 2>/dev/null | sed -n 2p | cut -f7 | sed 's/^N=\([0-9]*\).*/\1/')
      [ "$rows" -eq "$states" ] || bad "$name/$q: $rows rows for a $states-state mask"
      ## a state mask claims every row: its N_mask sum is the universe
      tail -n +2 w_$name.$q | /usr/bin/awk -F'\t' '{ s += $7; u = $5 } END { exit (s == u) ? 0 : 1 }' ||
        bad "$name/$q: the states' N_mask do not sum to N_univ"
    else
      [ "$rows" -eq "$recs" ] || bad "$name/$q: $rows rows for $recs records"
    fi
    [ "$q" = calls ] || tail -n +2 w_$name.$q | /usr/bin/awk -F'\t' '$5 != 29401795 { exit 1 }' ||
      bad "$name/$q: N_univ is not the genome's row count"

    ## which path ran -- the check identical output cannot make
    if [ "$fmt" = 2 ]; then
      paths w_$name.$q | grep -q "^0, per-mask $recs," || bad "$name/$q: a state mask should walk per mask: $(paths w_$name.$q)"
    else
      paths w_$name.$q | grep -q "^$recs, per-mask 0," || bad "$name/$q: the kernel did not take every mask: $(paths w_$name.$q)"
      if [ "$q" != calls ] && [ "$recs" -gt 1 ]; then
        paths i_$name.$q | grep -q "inverted-index $recs\$" || bad "$name/$q: -I did not take every mask: $(paths i_$name.$q)"
      fi
    fi
  done
done

## ---- the tables against the ground truth --------------------------------------
if [ -n "${YAME_KB_REBASE:-}" ] || [ ! -f "$expected" ]; then
  { printf 'mask/query\tQuery\tMask\tN_univ\tN_query\tN_mask\tN_overlap\tLog2OddsRatio\tBeta\tDepth\n'
    cat got_expected.tsv; } > "$expected"
  echo "expected tables written to $expected ($(grep -c . got_expected.tsv) rows); commit it"
elif ! cmp -s <(tail -n +2 "$expected") got_expected.tsv; then
  bad "the tables differ from $expected:"
  diff <(tail -n +2 "$expected") got_expected.tsv | head -20 | sed 's/^/    /'
fi

## ---- speed and memory against the baseline -----------------------------------
sort -o measured.tsv measured.tsv
cp measured.tsv all.tsv
if [ -n "${YAME_KB_REBASE:-}" ] || [ ! -f "$base" ]; then
  { printf 'run\tseconds\trss_kb\n'; grep -v '/ref	' measured.tsv; } > "$base"
  echo "baseline written to $base ($(grep -c . measured.tsv) runs); commit it"
else
  tail -n +2 "$base" | sort > base.tsv
  grep -v '/ref	' measured.tsv > cur.tsv || true
  join -t "$(printf '\t')" -j 1 base.tsv cur.tsv > joined.tsv
  n=$(grep -c . joined.tsv)
  [ "$n" -gt 0 ] || { bad "nothing in common with the baseline"; }
  ## joined: run, baseline seconds, baseline rss, this run's seconds, its rss.
  ## A run over the line is timed ONCE MORE before it counts: on a shared box
  ## a 0.9 s job has taken 3.1 s at load 7, and noise like that does not
  ## repeat, where a real regression does.
  over() {          # <joined.tsv> -> the runs over either line, one per row
    /usr/bin/awk -F'\t' '{ if ($4 > $2 * 2 + 0.5 || $5 > $3 * 1.5 + 51200)
      printf "%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $4, $5 }' "$1"
  }
  over joined.tsv > over1.tsv
  if [ -s over1.tsv ]; then
    echo "  over the baseline once, timing again: $(cut -f1 over1.tsv | paste -sd' ' -)"
    : > again.tsv
    while IFS=$(printf '\t') read -r tag out b_s b_m c_s c_m; do
      line=$(grep "^$tag	" cmds.tsv | head -1)
      [ -n "$line" ] || continue
      ## replay: fields after the out name are envs, then --, then the args
      set -- $(printf '%s' "$line" | cut -f3- | tr '\t' ' ')
      : > measured.tsv
      run "$tag" "$out" "$@"
      cat measured.tsv >> again.tsv
    done < over1.tsv
    join -t "$(printf '\t')" -j 1 base.tsv <(sort again.tsv) > joined2.tsv
    over joined2.tsv | /usr/bin/awk -F'\t' '
      { if ($4 > $2 * 2 + 0.5) printf "  SLOWER %-40s %6.2f s twice, baseline %6.2f s\n", $1, $4, $2
        if ($5 > $3 * 1.5 + 51200) printf "  MEMORY %-40s %6.0f MB twice, baseline %6.0f MB\n", $1, $5/1024, $3/1024
        bad = 1 } END { exit bad }' || fails=$((fails + 1))
  fi
  ## runs the baseline does not know are reported, not failed
  cut -f1 cur.tsv | sort > m.txt; cut -f1 base.tsv | sort > b.txt
  comm -23 m.txt b.txt | sed 's/^/  new run, not in baseline: /'
fi
## the reference beside this binary's walk, when it ran
if [ -n "${YAME_KB_REF:-}" ]; then
  echo
  printf '%-24s %8s %8s %6s\n' "mask/query" ref walk ratio
  /usr/bin/awk -F'\t' '
    { n = split($1, a, "/"); k = a[1] "/" a[2]
      if (a[3] == "ref") r[k] = $2; else if (a[3] == "walk") w[k] = $2 }
    END { for (k in r) if (k in w)
            printf "%-24s %8.2f %8.2f %6.1fx\n", k, r[k], w[k], (w[k] > 0) ? r[k] / w[k] : 0 }' \
    all.tsv | sort
fi

## the table, for the log
printf '%-40s %8s %8s\n' run seconds MB
/usr/bin/awk -F'\t' '{ printf "%-40s %8.2f %8.0f\n", $1, $2, $3/1024 }' all.tsv

[ "$fails" -eq 0 ] || { echo "$fails knowledgebase check(s) failed"; exit 1; }
echo "ok: t_store_kb ($(grep -c . all.tsv) runs over $(ls "$KB"/*.cm | wc -l | tr -d ' ') masks)"
