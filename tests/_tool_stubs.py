# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2026 Efesto contributors.
"""Small external-tool stand-ins with observable invocations."""
import json


def write_uniop_stub(path, probability=0.8, membership=True):
    script = '''import json, sys
from pathlib import Path
trace=Path(__file__).with_suffix('.calls.json')
calls=json.loads(trace.read_text()) if trace.exists() else []
calls.append(sys.argv[1:]);trace.write_text(json.dumps(calls))
control=Path(__file__).with_suffix('.control.json')
settings=json.loads(control.read_text()) if control.exists() else {}
out=Path(sys.argv[sys.argv.index('-t')+1]);out.mkdir(parents=True,exist_ok=True)
if '-a' in sys.argv:
    proteins=Path(sys.argv[sys.argv.index('-a')+1]).read_text()
    assert sys.argv[sys.argv.index('--faa_source')+1]=='prodigal_faa'
else:
    proteins='>internal_A\\nMAA\\n>internal_B\\nMKK\\n'
    (out/'internal.faa').write_text(proteins)
(out/'received.txt').write_text(proteins)
(out/'uniop.operon').write_text('idx_genes,idx_op\\n'+MEMBERSHIP)
(out/'uniop.pred').write_text('Gene A\\tGene B\\tPrediction\\n1\\t2\\tPROBABILITY\\n')
sys.exit(settings.get('exit',0))
'''.replace('MEMBERSHIP', repr('"[1, 2]",7\n' if membership else '')).replace('PROBABILITY',str(probability))
    path.write_text(script)
    path.chmod(0o755)
    return path


def invocation_count(script):
    path=script.with_suffix('.calls.json')
    return len(json.loads(path.read_text())) if path.exists() else 0
