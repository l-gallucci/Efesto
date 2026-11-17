#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2026 Efesto contributors.

"""Verify release archives and run the installed wheel outside its source tree."""
import argparse
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import subprocess
import tarfile
import tempfile
import venv
import zipfile


def sha(data):
    return hashlib.sha256(data).hexdigest()


def check(wheel, sdist, smoke=False):
    with zipfile.ZipFile(wheel) as whl, tarfile.open(sdist) as tar:
        top = tar.getmembers()[0].name.split("/")[0]
        names = whl.namelist()
        source_names = tar.getnames()
        def source(rel):
            return tar.extractfile(top + "/" + rel).read()
        registry_bytes = source("hmm_library/hmm_registry.tsv")
        assert whl.read("efesto/data/hmm_library/hmm_registry.tsv") == registry_bytes
        rows = list(csv.DictReader(io.StringIO(registry_bytes.decode()), delimiter="\t"))
        for row in rows:
            rel = "hmm_library/" + row["hmm_file"]
            assert sha(whl.read("efesto/data/" + rel)) == sha(source(rel)), rel
        for name in ("HMM-bitcutoffs.txt", "FeGenie-map.txt", "category_aliases.tsv",
                     "category_catalog.tsv", "category_curation.tsv", "model_licenses.tsv",
                     "LICENSES.md"):
            rel = "hmm_library/" + name
            assert whl.read("efesto/data/" + rel) == source(rel), rel
        for name in ("LICENSE", "THIRD_PARTY_NOTICES.md", "CITATION.cff"):
            assert whl.read("efesto/data/" + name) == source(name)
        model_licenses = list(csv.DictReader(io.StringIO(
            source("hmm_library/model_licenses.tsv").decode()), delimiter="\t"))
        assert len(model_licenses) == len(rows)
        assert {r["stem"] for r in model_licenses} == {r["stem"] for r in rows}
        for rel in ("LICENSES/MIT.txt", "LICENSES/GPL-3.0.txt", "LICENSES/README.md"):
            assert whl.read("efesto/data/" + rel) == source(rel), rel
            assert any(n.endswith(".dist-info/licenses/" + rel) for n in names), rel
        for seed_alignment in ("MtrA", "MtoA", "CymA", "ImoA", "DA_402", "CbcL"):
            assert source("data/seeds/curated/" + seed_alignment + "_seeds.aln")
        cbcl_sources = json.loads(source("data/seeds/curated/CbcL_source_records.json"))
        assert len(cbcl_sources["records"]) == 8
        assert json.loads(source("data/seeds/curated/CbcL_reconstruction.json"))["same_entire_hmm_except_date"]
        assert whl.read("efesto/data/examples/operon_rules.json") == source("operon_rules.json")
        assert sum(n.endswith(".hmm") for n in names) == len(rows)
        assert sum(n.endswith(".hmm") for n in source_names) == len(rows)
        for name in names + source_names:
            assert not any(part in {"benchmark", "_calibration", "transport_workflow"}
                           for part in Path(name).parts), name
        panel = json.loads(source("data/validation/metal_expansion/reference_expectations.json"))
        seed_id = panel["lanM_NF053806"]["positive"][0]
        fasta = source("data/validation/metal_expansion/references.faa").decode()
        sequence = ""
        selected = False
        for line in fasta.splitlines():
            if line.startswith(">"):
                selected = line[1:].split()[0] == seed_id
            elif selected:
                sequence += line.strip()
        assert sequence
    result = {"version": wheel.name.split("-")[1], "profiles": len(rows),
              "enabled_profiles": sum(r["status"] in {"active", "experimental"} for r in rows),
              "profile_and_metadata_hashes_match": True}
    if smoke:
        with tempfile.TemporaryDirectory(prefix="efesto-package-check-") as temp:
            root = Path(temp)
            env_dir = root / "venv"
            venv.EnvBuilder(with_pip=True).create(env_dir)
            python = env_dir / "bin" / "python"
            env = os.environ.copy()
            env.pop("PYTHONPATH", None)
            env.pop("EFESTO_DEV_MODE", None)
            env["PYTHONNOUSERSITE"] = "1"
            env["PIP_DISABLE_PIP_VERSION_CHECK"] = "1"
            def run(args):
                process = subprocess.run(args, cwd=root, env=env, capture_output=True, text=True)
                if process.returncode:
                    raise RuntimeError(process.stdout + "\n" + process.stderr)
                return process.stdout
            run([str(python), "-m", "pip", "install", "--no-index", "--no-deps", str(wheel.resolve())])
            installed = json.loads(run([str(python), "-I", "-c",
                "import json,efesto;from efesto.resources import bundled_hmm_library;"
                "from efesto.io import load_registry;"
                "p=bundled_hmm_library();print(json.dumps({'version':efesto.__version__,"
                "'library':str(p),'profiles':len(load_registry(p))}))"]))
            assert installed["version"] == result["version"]
            assert installed["profiles"] == len(rows)
            assert Path(installed["library"]).resolve().is_relative_to(env_dir.resolve()), installed
            command = env_dir / "bin" / "efesto"
            assert run([str(command), "--version"]).strip() == "Efesto " + result["version"]
            run([str(command), "--help"])
            faa = root / "proteins"
            faa.mkdir()
            (faa / "sample.faa").write_text(">reference_1\n" + sequence + "\n")
            output = root / "results"
            run([str(command), "--faa_dir", str(faa), "--annotate", "Ln", "--out", str(output),
                 "--all_results", "--threads", "1"])
            with (output / "Efesto-results-long.tsv").open() as handle:
                hits = list(csv.DictReader(handle, delimiter="\t"))
            assert len(hits) == 1 and hits[0]["gene"].lower() == "lanm", hits
            assert hits[0]["threshold_type"] == "model_specific"
            assert "cluster_confidence" not in hits[0]
            assert (output / "Efesto-context-evidence.tsv").is_file()
            assert (output / "Efesto-UniOP-pair-evidence.tsv").is_file()
            result["clean_wheel_install_and_real_hmmer_run"] = True
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("wheel", type=Path)
    parser.add_argument("sdist", type=Path)
    parser.add_argument("--smoke", action="store_true", help="Install and run with real HMMER")
    args = parser.parse_args()
    print(json.dumps(check(args.wheel, args.sdist, args.smoke), indent=2))
