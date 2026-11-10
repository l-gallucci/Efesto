# Custom HMM libraries

Work in a separate library directory while evaluating new or revised profiles.
Do not overwrite source HMMs or lower official GA thresholds merely to recover a
small reference panel.

## Required files

```text
custom_library/
  Category-Name/profile_stem.hmm
  hmm_registry.tsv
  HMM-bitcutoffs.txt
  FeGenie-map.txt
  operon_rules.json                 # optional JSON context rules
```

The filename stem keys the registry, cutoff table, gene map and rule membership.
The internal HMM `NAME` may differ; official imports already use that convention.
Keep these external keys synchronized.

For each model, record the source/version, sequence provenance, NSEQ, category,
cutoff score type, runtime status, references and review limitations. `active`
and `experimental` are enabled states. An unregistered file is not an adequate
way to mark a model as development-only: keep such candidates outside the
runtime library.

## Verify and use

```bash
python scripts/curate_hmm_library.py --verify custom_library/
efesto --faa_dir proteins/ --gff_dir gff/ \
  --hmm_dir custom_library/ --out custom_library_out
```

Use a new output directory when the library changes. `--normalize_hmms` can
convert older HMM formats in the chosen library; it is a format conversion, not
biological calibration.

## Thresholds and new categories

Efesto's cutoff table uses full-sequence bitscores. A domain threshold cannot be
copied into it as though it were a full-sequence threshold. Zero/nonpositive
cutoffs invoke the generic bitscore floor and `low_confidence` interpretation.

Acquisition, resistance and homeostasis should remain separate. Check canonical
selectors in the source before adding categories. Context rules use model stems;
a complete system must not be inferred by counting duplicate copies of one component.

Review biological specificity, aliases, real sequence controls and model
competition before deployment. Similar HMM consensus sequences identify a
redundancy candidate, not proof that models are functionally interchangeable.
See [methodology](hmm_methodology.md) and [curation](hmm_library_curation.md).
