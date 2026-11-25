# Verification archive

Evidence collected from the independent exact local search of NCBI genome assemblies and transcripts (`step2_verify/`)
and records of the earlier STAR3D production and QC runs (`step3/`). `COLLECT.sha256` lists every file with its sha256.

For redistribution, the STAR3D v1.2 package (source, compiled classes, jars, bundled MC-Annotate and RemovePseudoknots
binaries, libraries and the release tarball) was removed from the QC folder; the smoke test's own outputs
(`out/`, `STAR3D_struct_info/`) are kept. 286 package files were removed and dropped from `COLLECT.sha256`
(365 -> 79 entries). STAR3D is available from its authors (see `config.yaml` for URL and sha256).
