#!/usr/bin/env bash
set -euo pipefail

mkdir -p input
mkdir -p output
mkdir -p truth

run() {
  echo "+ $*"
  "$@"
}

cat > input/f0_binary.txt << 'EOF'
1
0
1
1
0
1
0
1
EOF

cat > input/f1_int.txt << 'EOF'
0
9
a
2
2
1
A
C
EOF

echo "[FMT1] pack input/f1_int.txt -> output/f1.cg"
# run diff -u truth/f1.cg output/f1.cg

cat > input/f2_states.txt << 'EOF'
Promoter
Promoter
Exon
Intron
Intergenic
Intergenic
Promoter
Exon
EOF

echo "[FMT2] pack input/f2_states.txt -> output/f2.cm"
{
  printf '10\t0\n'
  printf '0\t5\n'
  printf '3\t7\n'
  printf '0\t0\n'
  printf '4\t4\n'
  printf '20\t5\n'
  printf '4\t4\n'
  printf '21\t5\n'
} > input/f3_mu.tsv

cat > input/f4_beta.txt << 'EOF'
0.10
0.50
NA
0.00
1.00
0.75
0.62
0.41
EOF

{
  printf '1\t1\n'
  printf '0\t1\n'
  printf '1\t1\n'
  printf '0\t1\n'
  printf '0\t0\n'
  printf '0\t0\n'
  printf '0\t0\n'
  printf '0\t1\n'  
} > input/f6_sparse.tsv

{
  printf "chr1\t100\t101\tCpG_1\n"
  printf "chr1\t200\t201\tCpG_2\n"
  printf "chr1\t300\t301\tCpG_3\n"
  printf "chr1\t400\t401\tCpG_4\n"
  printf "chr1\t500\t501\tCpG_5\n"
  printf "chr1\t600\t601\tCpG_6\n"
  printf "chr2\t500\t501\tCpG_7\n"
  printf "chr2\t600\t601\tCpG_8\n"
} > input/f7_coords.bed

run ../yame pack -f1 input/f1_int.txt output/f1.cg
run ../yame pack -f2 input/f2_states.txt output/f2.cm
run ../yame pack -f3 input/f3_mu.tsv output/f3.cg
run ../yame pack -fb input/f0_binary.txt output/f0.cg
run ../yame pack -f4 input/f4_beta.txt output/f4.cg
run ../yame pack -f6 input/f6_sparse.tsv output/f6.cx
run ../yame pack -f7 input/f7_coords.bed output/f7.cr

# test split
cat output/f1.cg output/f2.cm output/f3.cg output/f6.cx output/f7.cr >output/combined.cg
run ../yame split output/combined.cg output/combine_pre

run ../yame chunk -s 10000000 input/cpg_nocontig.cr output/cpg_nocontig_chunk
