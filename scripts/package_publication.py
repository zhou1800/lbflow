#!/usr/bin/env python3
"""Package a clean scientific source snapshot and checksum-indexed evidence.

This only creates local archives. It does not publish, submit or upload them.
The previous full archive is retained unchanged in the new reproduction ZIP.
Frozen manuscripts and output must be outside the public checkout; historical
contact metadata makes these private preparation artifacts.
"""
import argparse
import hashlib
import io
import json
from pathlib import Path
import subprocess
import tarfile
import zipfile


def digest(file):
    with Path(file).open('rb') as stream:return stream_digest(stream)


def stream_digest(stream):
    value=hashlib.sha256()
    for block in iter(lambda:stream.read(1024*1024),b''):value.update(block)
    return value.hexdigest()


def package(repo,raw,manifest_path,previous,validation_path,output,manuscripts):
    repo,output,manuscripts=repo.resolve(),output.resolve(),manuscripts.resolve()
    if output.is_relative_to(repo) or manuscripts.is_relative_to(repo):
        raise ValueError('Private package output and frozen manuscripts must be outside the checkout')
    if subprocess.check_output(['git','status','--porcelain'],cwd=repo,text=True).strip():
        raise RuntimeError('Commit the reviewed source before packaging; the working tree must be clean')
    commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip()
    root=repo/'research/drlab-publication';paper=root/'paper'
    validation=json.loads(validation_path.read_text());manifest=json.loads(manifest_path.read_text())
    previous_sha=validation['backup']['previous_reproducibility_sha256']
    if digest(previous)!=previous_sha:raise AssertionError('Previous-stage backup differs')
    for name,entry in manifest['files'].items():
        path=raw/name
        if not path.is_file() or path.stat().st_size!=entry['bytes'] or digest(path)!=entry['sha256']:
            raise AssertionError('Raw file differs: '+name)
    if digest(paper/'data/major-revision.json')!=validation['dataset_sha256']:raise AssertionError('Canonical dataset changed')
    for name,entry in validation['manuscripts'].items():
        if digest(manuscripts/name)!=entry['sha256']:raise AssertionError('Manuscript changed: '+name)
    snapshot=dict(source_commit=commit,date='2026-10-09',author='Yu Zhou, Independent Researcher, Canada',
        dataset_sha256=validation['dataset_sha256'],raw_manifest_sha256=digest(manifest_path),
        manuscripts=validation['manuscripts'],figures=validation['figures'],new_executions=validation['new_executions'],
        baseline_commit=validation['backup']['commit'],baseline_tag=validation['backup']['tag'],
        previous_stage_sha256=previous_sha,
        status='Major-revision author-review draft; public program archive and final author review pending. No submission, upload or payment.')
    allowed=[]
    for path in sorted(paper.rglob('*')):
        if not path.is_file() or path.suffix not in ('.tex','.bib','.json','.txt','.png','.pdf'):continue
        if path.name in ('main.pdf','supplement.pdf','package-manifest.json','bibliography-metadata.json'):continue
        allowed.append((path,str(path.relative_to(paper))))
    output.mkdir(parents=True,exist_ok=False,mode=0o700)
    project=output/'latex-project.zip';review=output/'cpc-author-review.zip';reproduce=output/'lbflow-cpc-reproducibility.zip'
    targets=[project,review,reproduce];partials=[p.with_suffix('.partial.zip') for p in targets]
    with zipfile.ZipFile(partials[0],'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for path,name in allowed:z.write(path,name)
        z.writestr('snapshot.json',json.dumps(snapshot,indent=2)+'\n')
    with zipfile.ZipFile(partials[1],'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for path,name in allowed:z.write(path,'paper/'+name)
        for name in ('main.pdf','supplement.pdf'):z.write(manuscripts/name,'paper/'+name)
        for path in sorted((root/'notes').glob('*')):
            if path.is_file():z.write(path,'notes/'+path.name)
        z.write(validation_path,'major-revision-validation.json')
        z.write(manifest_path,'major-revision-raw-manifest.json')
        z.writestr('snapshot.json',json.dumps(snapshot,indent=2)+'\n')
    source_archive=subprocess.check_output(['git','archive',commit],cwd=repo)
    with tarfile.open(fileobj=io.BytesIO(source_archive)) as source:
        files={m.name:source.extractfile(m).read() for m in source.getmembers() if m.isfile()}
    with zipfile.ZipFile(partials[2],'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for name,blob in files.items():z.writestr(name,blob)
        for name in ('main.pdf','supplement.pdf'):z.write(manuscripts/name,'research/drlab-publication/paper/'+name)
        for name in manifest['files']:z.write(raw/name,'raw/major-revision/'+name)
        z.write(manifest_path,'raw/major-revision/raw-manifest.json')
        z.write(previous,'historical/pre-major-revision-reproducibility.zip',compress_type=zipfile.ZIP_STORED)
        z.writestr('snapshot.json',json.dumps(snapshot,indent=2)+'\n')
    print('Local archives created; checking CRCs, committed source bytes and every raw digest',flush=True)
    for path in partials:
        with zipfile.ZipFile(path) as z:
            if z.testzip() is not None:raise AssertionError('ZIP CRC failed: '+str(path))
    with zipfile.ZipFile(partials[2]) as z:
        for name,blob in files.items():
            if z.read(name)!=blob:raise AssertionError('Archived source differs: '+name)
        for name,entry in manifest['files'].items():
            with z.open('raw/major-revision/'+name) as stream:
                if stream_digest(stream)!=entry['sha256']:raise AssertionError('Archived raw differs: '+name)
        with z.open('historical/pre-major-revision-reproducibility.zip') as stream:
            if stream_digest(stream)!=previous_sha:raise AssertionError('Archived baseline differs')
    for partial,target in zip(partials,targets):partial.replace(target)
    result=dict(source_commit=commit,packages=[dict(file=p.name,bytes=p.stat().st_size,sha256=digest(p)) for p in targets],
        raw_files=len(manifest['files']),raw_bytes=sum(e['bytes'] for e in manifest['files'].values()),
        previous_archive_bytes=previous.stat().st_size,previous_archive_sha256=previous_sha,
        checks='Every ZIP CRC, every committed source byte, every archived raw digest, both manuscript hashes and unchanged full previous archive verified.',
        source_scope='This commit is the scientific source snapshot; a later receipt-only commit may record the completed archive verification.',
        external_actions='Local packaging only.')
    (output/'package-manifest.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo',type=Path,default=Path(__file__).resolve().parents[1])
    parser.add_argument('--raw-root',type=Path,required=True);parser.add_argument('--raw-manifest',type=Path,required=True)
    parser.add_argument('--previous-stage',type=Path,required=True);parser.add_argument('--validation',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--manuscripts',type=Path,required=True,help='Directory containing unchanged frozen main.pdf and supplement.pdf')
    args=parser.parse_args();package(args.repo.resolve(),args.raw_root.resolve(),args.raw_manifest.resolve(),args.previous_stage.resolve(),args.validation.resolve(),args.output.resolve(),args.manuscripts.resolve())


if __name__=='__main__':main()
