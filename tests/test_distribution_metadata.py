# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2026 Efesto contributors.
"""Release checks reject inconsistent cutoff metadata in actual archives."""
import importlib.util
import io
from pathlib import Path
import tarfile
import zipfile
import pytest


@pytest.mark.parametrize('table,message',[('model\t10.1\n','model'),
                                         ('model\t10.05\nmodel\t10.05\n','Duplicate cutoff')])
def test_distribution_rejects_inconsistent_or_duplicate_cutoffs(tmp_path,table,message):
    root=Path(__file__).resolve().parents[1]
    spec=importlib.util.spec_from_file_location('distribution_check',root/'scripts/check_distribution.py')
    checker=importlib.util.module_from_spec(spec);spec.loader.exec_module(checker)
    registry='stem\tcutoff\thmm_file\nmodel\t10.05\tCategory/model.hmm\n'
    wheel=tmp_path/'efesto-1.0.0-py3-none-any.whl';sdist=tmp_path/'efesto-1.0.0.tar.gz'
    with zipfile.ZipFile(wheel,'w') as archive:
        archive.writestr('efesto/data/hmm_library/hmm_registry.tsv',registry)
    with tarfile.open(sdist,'w:gz') as archive:
        for name,body in [('hmm_registry.tsv',registry),('HMM-bitcutoffs.txt',table)]:
            data=body.encode();item=tarfile.TarInfo('efesto-1.0.0/hmm_library/'+name);item.size=len(data)
            archive.addfile(item,io.BytesIO(data))
    with pytest.raises(AssertionError,match=message):checker.check(wheel,sdist)
