"""Read-only evidence check; NOT a physics validation or a paper reproduction.

Static mode verifies that every registered experiment keeps its BOUNDARY/INPUTS/RESULTS
documents and that frozen documents, inputs and outputs still match
experiments/EVIDENCE_MANIFEST.json. --replay regenerates replayable outputs in a
temporary copy of the working tree and requires byte-identical files; the working
tree itself is never written. A03 replay needs the local paper rasters (--images).

Usage:
  python scripts/check_evidence.py
  python scripts/check_evidence.py --replay [--images ../tmp/pdfs/a03]
  python scripts/check_evidence.py --hash PATH...   # print SHA-256 for a new manifest entry
"""
from pathlib import Path
import argparse,hashlib,json,os,shutil,subprocess,sys,tempfile

ROOT=Path(__file__).resolve().parents[1]
MANIFEST=ROOT/'experiments/EVIDENCE_MANIFEST.json'
REQUIRED_DOCS=('BOUNDARY.md','INPUTS.md','RESULTS.md')
COPY_IGNORE=shutil.ignore_patterns('.venv','.git','__pycache__','*.egg-info','.ipynb_checkpoints','.DS_Store')


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_manifest():
    return json.loads(MANIFEST.read_text())


def static_check(root=ROOT,manifest=None):
    """Return a list of failure strings; empty means every registered file is intact."""
    manifest=manifest or load_manifest()
    failures=[]
    for exp in manifest['experiments']:
        folder=root/exp['dir']
        for doc in REQUIRED_DOCS:
            if not (folder/doc).is_file():
                failures.append(f"{exp['id']}: missing {doc}")
        for group in ('documents','inputs','outputs'):
            for rel,expected in exp.get(group,{}).items():
                path=root/rel
                if not path.is_file():
                    failures.append(f"{exp['id']}: {group} file missing: {rel}")
                elif sha256(path)!=expected:
                    failures.append(f"{exp['id']}: {group} changed: {rel}")
    return failures


def replay_check(images=None):
    """Regenerate outputs in a scratch copy. Returns (failures, report lines)."""
    manifest=load_manifest();failures=[];report=[]
    with tempfile.TemporaryDirectory(prefix='tfln-replay-') as tmp:
        copy=Path(tmp)/'repo'
        shutil.copytree(ROOT,copy,ignore=COPY_IGNORE)
        env=dict(os.environ,PYTHONPATH=str(copy/'src'),MPLBACKEND='Agg')
        for exp in manifest['experiments']:
            spec=exp.get('replay')
            if spec is None:
                report.append(f"SKIP {exp['id']}: not replayable ({exp.get('replay_note','no note')})");continue
            if spec['kind']=='preflight':
                code=('import json,sys;from tfln_mzm.contracts import preflight;'
                      'print(json.dumps(preflight(json.load(open(sys.argv[1])))))')
                run=subprocess.run([sys.executable,'-c',code,str(copy/spec['config'])],cwd=copy,env=env,capture_output=True,text=True)
                if run.returncode!=0:
                    failures.append(f"{exp['id']}: preflight replay failed: {run.stderr.strip()[-300:]}");continue
                stored=json.loads((copy/spec['output']).read_text())
                ok=json.loads(run.stdout)==stored
                report.append(f"{'SAME' if ok else 'DIFF'} {exp['id']}: preflight result")
                if not ok:failures.append(f"{exp['id']}: preflight result differs from {spec['output']}")
                continue
            args=[a.replace('{images}',str(images)) if images else a for a in spec['command']]
            if any('{images}' in a for a in args):
                report.append(f"SKIP {exp['id']}: needs --images (local paper rasters, not in Git)");continue
            run=subprocess.run([sys.executable,*args],cwd=copy,env=env,capture_output=True,text=True)
            if run.returncode!=0:
                failures.append(f"{exp['id']}: replay exited {run.returncode}: {run.stderr.strip()[-300:]}");continue
            for rel,expected in exp['outputs'].items():
                if rel in exp.get('static_only',[]):continue
                same=(copy/rel).is_file() and sha256(copy/rel)==expected
                report.append(f"{'SAME' if same else 'DIFF'} {exp['id']}: {rel}")
                if not same:failures.append(f"{exp['id']}: replay output differs: {rel}")
    return failures,report


def main():
    parser=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--replay',action='store_true')
    parser.add_argument('--images',type=Path,help='folder with p2_2.jpeg and p3_0.png for A03 replay')
    parser.add_argument('--hash',nargs='+',metavar='PATH')
    a=parser.parse_args()
    if a.hash:
        for p in a.hash:print(sha256(p),p)
        return 0
    failures=static_check()
    count=len(load_manifest()['experiments'])
    print(('PASS' if not failures else 'FAIL')+f' static: {count} registered experiments, documents/inputs/outputs against manifest')
    if a.replay:
        images=a.images.resolve() if a.images else None
        replay_failures,report=replay_check(images)
        print('\n'.join(report))
        print(('PASS' if not replay_failures else 'FAIL')+' replay in a temporary copy; working tree untouched')
        failures+=replay_failures
    for f in failures:print('  -',f)
    print('This checks archived evidence only; it does not validate physics or establish Figure 3 reproduction.')
    return 1 if failures else 0


if __name__=='__main__':sys.exit(main())
