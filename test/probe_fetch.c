/* A downstream tool, in a dozen lines: `fetch` over ITS OWN registry, under
 * its own name and its own store variable, by linking libyame.a and handing
 * yame_fetch_main() a yame_fetch_cfg_t. This is what methscope does with the
 * registry a downstream tool's own generator produces. Nothing here is
 * yame's table; if fetch.o were still bound to it, this could not link with
 * a different one. */
#include <stdio.h>
#include <string.h>
#include "assets.h"
static const yame_asset_file_t F[] = {
  { "zhou-lab/probe@v2:one.cm", "hg38/probe/one.cm",
    "http://127.0.0.1:1/zhou-lab/probe/v2/one.cm",
    "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855", 0, 1,
    "One", "A probe set.", "probe", "-" },
  { NULL, NULL, NULL, NULL, 0, 0, NULL, NULL, NULL, NULL } };
int main(int argc, char **argv) {
  yame_fetch_cfg_t cfg = { F, YAME_NFILES(F), "methprobe", "METHPROBE_DATA_HOME" };
  return yame_fetch_main(&cfg, argc, argv);
}
