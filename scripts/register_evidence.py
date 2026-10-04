"""Append one finished experiment to experiments/EVIDENCE_MANIFEST.json.

Usage (from the repository root): python scripts/register_evidence.py '<json spec>'
Spec keys: id, dir, grade, layer, boundary_commit, inputs, outputs, and optionally
local_inputs, static_only, replay, replay_note. Hashes are computed now. An id that is
already registered is refused: registered evidence is never re-hashed or rewritten.
"""
import hashlib,json,sys
from pathlib import Path
spec=json.loads(sys.argv[1])
p=Path('experiments/EVIDENCE_MANIFEST.json');m=json.loads(p.read_text())
assert spec['id'] not in {e['id'] for e in m['experiments']},'already registered'
h=lambda f:hashlib.sha256(Path(f).read_bytes()).hexdigest()
d=spec['dir']
entry={'id':spec['id'],'dir':d,'grade':spec['grade'],'layer':spec.get('layer'),'recorded_commit':spec.get('boundary_commit'),
 'documents':{f'{d}/{n}':h(f'{d}/{n}') for n in ('BOUNDARY.md','INPUTS.md','RESULTS.md')},
 'inputs':{f:h(f) for f in spec.get('inputs',[])},'outputs':{f:h(f) for f in spec['outputs']}}
for k in ('local_inputs','static_only','replay','replay_note'):
    if k in spec:entry[k]=spec[k]
m['experiments'].append(entry);p.write_text(json.dumps(m,indent=1,ensure_ascii=False)+'\n');print('registered',spec['id'])
