#!/usr/bin/env bash
set -euo pipefail

mkdir -p input
mkdir -p output
mkdir -p truth

echo "Generating toy pack inputs in: input/"
echo "Packing outputs into:          output/"
echo

############################################
# Helper: echo command then run it
############################################
run() {
  echo "+ $*"
  "$@"
}

############################################
# FORMAT 0 – Binary (fb)
############################################
cat > input/f0_binary.txt << 'EOF'
1
0
1
1
0
1
EOF

echo "[FMT0] pack input/f0_binary.txt -> output/f0.cg"
run ../yame pack -fb input/f0_binary.txt output/f0.cg
# run diff -u truth/f0.cg output/f0.cg

############################################
# FORMAT 1 – Integer (f1)
############################################
cat > input/f1_int.txt << 'EOF'
0
0
0
2
2
1
EOF

echo "[FMT1] pack input/f1_int.txt -> output/f1.cg"
run ../yame pack -f1 input/f1_int.txt output/f1.cg
# run diff -u truth/f1.cg output/f1.cg

############################################
# FORMAT 2 – States (f2)
############################################
cat > input/f2_states.txt << 'EOF'
Promoter
Promoter
Exon
Intron
Intergenic
Intergenic
EOF

echo "[FMT2] pack input/f2_states.txt -> output/f2.cm"
run ../yame pack -f2 input/f2_states.txt output/f2.cm
# run diff -u truth/f2.cm output/f2.cm

############################################
# FORMAT 3 – MU counts (f3)
############################################
{
  printf '10\t0\n'
  printf '0\t5\n'
  printf '3\t7\n'
  printf '0\t0\n'
  printf '4\t4\n'
  printf '20\t5\n'
} > input/f3_mu.tsv

echo "[FMT3] pack input/f3_mu.tsv -> output/f3.cg"
run ../yame pack -f3 input/f3_mu.tsv output/f3.cg
# run diff -u truth/f3.cg output/f3.cg

############################################
# FORMAT 4 – Beta values (f4)
############################################
cat > input/f4_beta.txt << 'EOF'
0.10
0.50
NA
0.00
1.00
0.75
EOF

echo "[FMT4] pack input/f4_beta.txt -> output/f4.cg"
run ../yame pack -f4 input/f4_beta.txt output/f4.cg
# run diff -u truth/f4.cg output/f4.cg

############################################
# FORMAT 6 – Query/Universe (f6)
############################################
{
  printf '1\t1\n'
  printf '0\t1\n'
  printf '1\t1\n'
  printf '0\t1\n'
  printf '0\t0\n'
  printf '0\t0\n'
} > input/f6_sparse.tsv

echo "[FMT6] pack input/f6_sparse.tsv -> output/f6.cx"
run ../yame pack -f6 input/f6_sparse.tsv output/f6.cx
# run diff -u truth/f6.cx output/f6.cx

############################################
# FORMAT 7 – Coordinates (f7)
############################################
{
  printf "chr1\t100\t101\tCpG_1\n"
  printf "chr1\t200\t201\tCpG_2\n"
  printf "chr1\t300\t301\tCpG_3\n"
  printf "chr1\t400\t401\tCpG_4\n"
  printf "chr1\t500\t501\tCpG_5\n"
  printf "chr1\t600\t601\tCpG_6\n"
} > input/f7_coords.bed

echo "[FMT7] pack input/f7_coords.bed -> output/f7.cr"
run ../yame pack -f7 input/f7_coords.bed output/f7.cr
# run diff -u truth/f7.cr output/f7.cr

echo
echo "Done. Input files in input/, packed files in output/. To enable validation, add truth/*.cg and uncomment diff lines."
