"""Where the Sequence tab's bases come from: the held GFF, a disk cache, AllTheBacteria, or NCBI Datasets.

Two remote stations, one per key kind — a BioSample id goes to AllTheBacteria, a RefSeq/GenBank assembly
stem to NCBI — and never both for one genome. `contig_namespace` reconciles an NCBI assembly's ``NZ_``
headers with the GenBank-style contig names the coordinates were built from.
"""
