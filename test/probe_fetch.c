/* A downstream tool, in eighteen lines: `fetch` over ITS OWN registry, under
 * its own name and its own store variable, by linking libyame.a and handing
 * yame_fetch_main() a yame_fetch_cfg_t. This is what methscope does with the
 * registry make_registry.sh --tool=methscope generates. Nothing here is
 * yame's catalogue; if fetch.o were still bound to it, this could not link
 * with a different one. */
#include <stdio.h>
#include <string.h>
#include "assets.h"
static const yame_asset_file_t F[] = {
  { "one.cm", "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855", 0, NULL },
  { NULL, NULL, 0, NULL } };
static const yame_asset_reg_t R[] = {
  { "probe", "hg38/probe", "http://127.0.0.1:1/zhou-lab/probe", "v2", "", "hg38/probe",
    "0000000000000000000000000000000000000000000000000000000000000000", F, YAME_NFILES(F), NULL, 0 },
  { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0 } };
int main(int argc, char **argv) {
  yame_fetch_cfg_t cfg = { R, 1, "methprobe", "METHPROBE_DATA_HOME" };
  return yame_fetch_main(&cfg, argc, argv);
}
