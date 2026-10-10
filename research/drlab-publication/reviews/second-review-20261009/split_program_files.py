"""Split the verified program into independently extractable bounded ZIPs."""
from pathlib import Path
import hashlib,json,shutil,zipfile
base=Path('/home/user/zhou1800/lbflow/research/drlab-publication/submission/cpc-20261009')
out=base/'program-parts';out.mkdir(exist_ok=False)
limit=80*1024*1024;history='historical/pre-major-revision-reproducibility.zip'
manifest={'input_program_sha256':hashlib.sha256((base/'LBFlow-program-and-data.zip').read_bytes()).hexdigest(),'history_sha256':'33f2b23888193c31bfacebcc4b965d14018fe4c0a380214a324aa213b52260c8','parts':[],'history_chunks':[]}
readme='''LBFlow program and numerical data: extract all program ZIP files into the same fresh directory. Source, samples and raw evidence keep their documented paths. Then run python3 assemble_historical_archive.py to reconstruct the unchanged historical reproduction archive from checksum-verified chunks. README.md and research/drlab-publication/notes/major-revision-reproduction.md give build and reanalysis instructions. The manuscript and supplement are separate submission files.\n'''
def record(p):
 assert p.stat().st_size<100*1024*1024
 manifest['parts'].append({'file':p.name,'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
with zipfile.ZipFile(base/'LBFlow-program-and-data.zip') as src:
 regular=[i for i in src.infolist() if i.filename!=history]
 groups=[[]];size=0
 for info in regular:
  if info.compress_size>limit:raise AssertionError('Single raw file too large')
  if size+info.compress_size>limit and groups[-1]:groups.append([]);size=0
  groups[-1].append(info);size+=info.compress_size
 for idx,group in enumerate(groups,1):
  p=out/f'LBFlow-source-raw-{idx:02d}.zip'
  with zipfile.ZipFile(p,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
   for info in group:
    with src.open(info.filename) as a,z.open(info.filename,'w') as b:shutil.copyfileobj(a,b)
   z.writestr('CPC-program-parts.txt',readme)
  record(p)
 with src.open(history) as stream:
  idx=0;full=hashlib.sha256()
  while True:
   block=stream.read(limit)
   if not block:break
   idx+=1;full.update(block);name=f'historical/chunks/reproduction-{idx:02d}.bin'
   manifest['history_chunks'].append({'file':name,'bytes':len(block),'sha256':hashlib.sha256(block).hexdigest()})
   p=out/f'LBFlow-historical-{idx:02d}.zip'
   with zipfile.ZipFile(p,'w',compression=zipfile.ZIP_STORED) as z:z.writestr(name,block);z.writestr('CPC-program-parts.txt',readme)
   record(p)
  assert full.hexdigest()==manifest['history_sha256']
 # Verify every original regular entry and reconstruct the historical hash.
 entries={}
 for p in sorted(out.glob('*.zip')):
  with zipfile.ZipFile(p) as z:
   assert z.testzip() is None
   for name in z.namelist():
    if name=='CPC-program-parts.txt':continue
    if name in entries:raise AssertionError('Duplicate program entry')
    entries[name]=p
 for info in regular:
  a=hashlib.sha256(src.read(info.filename)).hexdigest()
  with zipfile.ZipFile(entries[info.filename]) as z:b=hashlib.sha256(z.read(info.filename)).hexdigest()
  assert a==b,info.filename
helper='''from pathlib import Path\nimport hashlib,json\nroot=Path(__file__).resolve().parent\nm=json.loads((root/'program-parts-manifest.json').read_text())\ntarget=root/'historical/pre-major-revision-reproducibility.zip'\nif target.exists():\n if hashlib.sha256(target.read_bytes()).hexdigest()!=m['history_sha256']:raise RuntimeError('Existing archive differs')\n print('Historical archive already verified');raise SystemExit\npartial=target.with_suffix('.partial.zip');full=hashlib.sha256()\nwith partial.open('xb') as out:\n for entry in m['history_chunks']:\n  data=(root/entry['file']).read_bytes()\n  if len(data)!=entry['bytes'] or hashlib.sha256(data).hexdigest()!=entry['sha256']:raise RuntimeError('Chunk differs: '+entry['file'])\n  full.update(data);out.write(data)\nif full.hexdigest()!=m['history_sha256']:raise RuntimeError('Historical archive hash differs')\npartial.replace(target);print('Historical archive reconstructed and verified')\n'''
p=out/'LBFlow-reassembly.zip'
with zipfile.ZipFile(p,'w',compression=zipfile.ZIP_DEFLATED) as z:
 z.writestr('assemble_historical_archive.py',helper)
 z.writestr('program-parts-manifest.json',json.dumps(manifest,indent=2)+'\n')
 z.writestr('CPC-program-parts.txt',readme)
record(p)
manifest['verification']='All part CRCs, every original source/raw entry and the complete historical-stream digest verified.'
(out/'parts-verification.json').write_text(json.dumps(manifest,indent=2)+'\n')
print([(x['file'],x['bytes']) for x in manifest['parts']],flush=True)
