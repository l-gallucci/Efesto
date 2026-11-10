# Ln/W reference checks

These reference sequences and labels support regression checks of the shipped official profiles, not an independent accuracy benchmark. Model training sets may contain the same proteins. UniProt sequence data are attributed to UniProt (https://www.uniprot.org/help/license); the three historical Vup-associated sequences were retrieved from NCBI Protein. Exact accessions, URLs, sequence hashes and label limitations are in `reference_sources.json`.

`profile_sources.json` preserves official download URLs, uncompressed HMM hashes, accession versions, seed counts, inherited GA thresholds and interpretation limits. No new HMM was trained and no source threshold was lowered. `reference_expectations.json` defines the regression panel. The larger discovery panel and all failed reference detections are documented in the local review report and `docs/metal_expansion.md`.

The known TupA coverage failures (Q316W1, Q93KD6 and inferred Q749P2) are not included as expected-positive regression cases: their failure is explicitly documented rather than hidden or corrected by an arbitrary cutoff. Absence of a hit is not absence of tungsten acquisition.
