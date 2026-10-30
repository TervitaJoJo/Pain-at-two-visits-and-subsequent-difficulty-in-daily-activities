from pathlib import Path
import json,hashlib
ROOT=Path(__file__).resolve().parent.parent
manifest=ROOT/'documentation/FREEZE_MANIFEST.json'
def main():
 obj=json.loads(manifest.read_text(encoding='utf-8'))
 expected={item['path']:item for item in obj['files']}
 actual={p.relative_to(ROOT).as_posix():p for p in ROOT.rglob('*') if p.is_file() and p!=manifest}
 if set(actual)!=set(expected):raise AssertionError({'missing':sorted(set(expected)-set(actual)),'extra':sorted(set(actual)-set(expected))})
 for name,p in actual.items():
  if hashlib.sha256(p.read_bytes()).hexdigest()!=expected[name]['sha256']:raise AssertionError('Changed frozen file: '+name)
 for folder in ROOT.iterdir():
  if not folder.is_dir():raise AssertionError('Files directly below code are prohibited')
  n=sum(p.is_file() for p in folder.rglob('*'))
  if n>100:raise AssertionError(f'{folder.name}: {n} files exceeds 100')
 print(f'PASS: {len(expected)} frozen files verified; every project folder <=100 files.')
if __name__=='__main__':main()
