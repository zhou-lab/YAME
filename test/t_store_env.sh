#!/bin/bash
## Store resolution across the suite. The store is one shared tree, so the
## variable a reader exported for one tool has to move it for all of them:
## with $METHSCOPE_DATA_HOME set and $YAME_DATA_HOME unset, yame must read
## that store, not fall through to ~/.local/share/yame. The docs' Upscale
## block is the case that found this -- `yame hprint` reported the hg38
## coordinate track missing while it sat in the store methscope had resolved.
## Offline: nothing here fetches, it only asks yame where it would look.
set -euo pipefail
YAME=${YAME:?export YAME=/path/to/yame}

d=$(mktemp -d); trap 'rm -rf "$d"' EXIT
ms="$d/msstore"; ys="$d/ystore"
mkdir -p "$ms" "$ys"
cd "$d"

## ---- 1. $METHSCOPE_DATA_HOME alone resolves the store ----------------------
env -u YAME_DATA_HOME -u XDG_DATA_HOME METHSCOPE_DATA_HOME="$ms" \
  "$YAME" fetch -h </dev/null > h.ms 2>&1 || true
grep -q "METHSCOPE_DATA_HOME: $ms" h.ms ||
  { echo "METHSCOPE_DATA_HOME did not resolve the store"; grep DATA_HOME h.ms; exit 1; }

## ---- 2. yame's own variable still wins -------------------------------------
env XDG_DATA_HOME="$d/xdg" METHSCOPE_DATA_HOME="$ms" YAME_DATA_HOME="$ys" \
  "$YAME" fetch -h </dev/null > h.both 2>&1 || true
grep -q "YAME_DATA_HOME: $ys" h.both ||
  { echo "YAME_DATA_HOME did not win over the suite variable"; grep DATA_HOME h.both; exit 1; }

## ---- 3. the order is printed, suite variables last -------------------------
grep -q 'Resolved in order: -d, \$YAME_DATA_HOME, \$METHSCOPE_DATA_HOME,' h.both ||
  { echo "the printed order does not list the suite variable after yame's"; grep Resolved h.both; exit 1; }

## ---- 4. the banner names the variable that actually resolved it ------------
## Saying "YAME_DATA_HOME: <path>" over a path that came from somewhere else
## sends a reader to check a variable they never set.
env -u YAME_DATA_HOME -u XDG_DATA_HOME METHSCOPE_DATA_HOME="$ms" \
  "$YAME" </dev/null > b.ms 2>&1 || true
grep -q "METHSCOPE_DATA_HOME: $ms" b.ms ||
  { echo "the banner does not name the variable that resolved the store"; grep DATA_HOME b.ms; exit 1; }

## ---- 5. neither set: the default, and it says so ---------------------------
env -u YAME_DATA_HOME -u METHSCOPE_DATA_HOME XDG_DATA_HOME="$d/xdg" \
  "$YAME" </dev/null > b.none 2>&1 || true
grep -q "YAME_DATA_HOME: $d/xdg/yame (unset, default)" b.none ||
  { echo "the unset default is wrong or unmarked"; grep DATA_HOME b.none; exit 1; }

echo "t_store_env: ok"
