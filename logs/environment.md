# Environment check

- Host: macOS 26.6.2, arm64; Python 3.13.7
- colima: colima version 0.10.3

| Check | Result | Detail |
|---|---|---|
| python 3.13 | PASS | 3.13.7 |
| pkg certifi | PASS | 2026.10.23 |
| pkg gemmi | PASS | 0.7.5 |
| pkg matplotlib | PASS | 3.11.2 |
| pkg mmcif-pdbx | PASS | 2.1.0 |
| pkg numpy | PASS | 2.5.4 |
| pkg openpyxl | PASS | 3.1.5 |
| pkg pandas | PASS | 3.0.6 |
| pkg pytest | PASS | 9.2.0 |
| pkg PyYAML | PASS | 6.0.3 |
| pkg scipy | PASS | 1.18.1 |
| pkg fr3d | PASS | 0.0.2 |
| fr3d import | PASS |  |
| input rfam_seed | PASS | inputs/rfam/15.1/Rfam.seed.gz |
| input rfam_pdb_mapping | PASS | inputs/rfam/15.1/Rfam.pdb.gz |
| input rfam_readme | PASS | inputs/rfam/15.1/README |
| input star3d | PASS | inputs/software/STAR3D_v1.2.tar.gz |
| docker daemon | PASS | x86_64 29.5.2 |
| image star3d-runtime:1 | PASS | sha256:e1602e600f599264c7ca239202433dd009118214db05459744a7d3e59aa61b7e amd64 |
| container java | PASS | openjdk version "1.8.0_504" |
| STAR3D.jar starts | PASS | usage: java -jar STAR3D.jar <options> PDB1 chain1 PDB2 chain2 /  -e <arg>   Gap extension penalty /  -g <arg>   Gap open |
| bundled MC-Annotate executes | PASS | usage: MC-Annotate [-bhlvV] [-e num] [-f <model number>] [-r <residue ids>] <structure file> ... / rc=1 |
| bundled RemovePseudoknots executes | PASS | Incorrect number of required parameters given. / USAGE: RemovePseudoknots <input ct file> <output ct file> [options] / r |
