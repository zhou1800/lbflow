#!/usr/bin/env python3
"""Prepare local CPC files from committed source and verified frozen evidence.

Private author contact is inserted only in the manuscript-source package.
This command does not register accounts, upload files, submit or pay.
"""
import argparse
import io
import json
from pathlib import Path
import shutil
import subprocess
import tarfile
import zipfile

from package_publication import digest, stream_digest


def insert_private_contact(paper, contact):
    """Add contact to private copies, leaving public templates intact."""
    main = paper / 'main.tex'
    public_address = r'\address{Independent Researcher, Canada}'
    full_address = ', '.join([contact['affiliation'], contact['street'], contact['city'],
                              contact['state'], contact['postal_code'], contact['country']])
    email = contact['email']
    text = main.read_text()
    if (text.count(public_address) != 1 or '@' not in email or
            any(c in full_address + email for c in '\\{}%$&#_\r\n')):
        raise ValueError('Unexpected contact template or unescaped TeX contact')
    main.write_text(text.replace(public_address, r'\address{' + full_address + '}\n' +
                                r'\ead{' + email + '}'))
    supplement = paper / 'supplement.tex'
    author = r'\author{Yu Zhou\\Independent Researcher, Canada}'
    text = supplement.read_text()
    if text.count(author) != 1:
        raise ValueError('Unexpected supplement author template')
    supplement.write_text(text.replace(author, author[:-1] +
                                      r'\\\texttt{' + email + '}}'))


def prepare(repo, previous, contact_path, samples, output, engine):
    repo, contact_path, output = repo.resolve(), contact_path.resolve(), output.resolve()
    if output.is_relative_to(repo) or contact_path.is_relative_to(repo):
        raise ValueError('Private contact and submission output must be outside the checkout')
    if subprocess.check_output(['git', 'status', '--porcelain'], cwd=repo, text=True).strip():
        raise RuntimeError('Commit candidate source before preparing submission files')
    commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=repo, text=True).strip()
    root = repo / 'research/drlab-publication'
    verified = json.loads((root / 'major-revision-package-verification.json').read_text())
    expected = next(x for x in verified['packages'] if x['file'] == 'lbflow-cpc-reproducibility.zip')
    if digest(previous) != expected['sha256']:
        raise AssertionError('Previously verified reproduction archive changed')
    manifest = json.loads((root / 'major-revision-raw-manifest.json').read_text())
    contact = json.loads(contact_path.read_text())
    if not contact['originality_confirmed'] or not contact['not_under_consideration_elsewhere']:
        raise RuntimeError('Originality and exclusive-submission facts are required')
    sample_check = json.loads((samples / 'verification.json').read_text())
    if {x['case'] for x in sample_check['cases']} != {'channel', 'drop'}:
        raise AssertionError('Both comprehensive sample cases are required')
    for case in sample_check['cases']:
        folder = samples / case['case']
        for name, sha in case['output_sha256'].items():
            if digest(folder / name) != sha:
                raise AssertionError('Sample output changed: ' + name)
        for name, key in [('input.prm', 'input_sha256'), ('geometry.raw', 'geometry_sha256')]:
            if digest(folder / name) != case[key]:
                raise AssertionError('Sample input changed: ' + name)
        if case['exit'] != 0 or not all(x['finite'] for x in case['history']):
            raise AssertionError('Incomplete sample verification')
    output.mkdir(parents=True, exist_ok=False, mode=0o700)
    paper = output / 'manuscript-source'
    paper.mkdir()
    original = root / 'paper'
    for path in original.rglob('*'):
        rel = path.relative_to(original)
        if path.is_file() and (path.suffix in ('.tex', '.bib') or
                               (rel.parts[0] == 'figures' and path.suffix == '.pdf')):
            target = paper / rel
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, target)
    insert_private_contact(paper, contact)
    for name in ('main.tex', 'supplement.tex'):
        subprocess.run([str(engine), '--keep-logs', name], cwd=paper, check=True)
    for src, dst in [('main.pdf', 'manuscript.pdf'), ('supplement.pdf', 'supplement.pdf')]:
        shutil.copy2(paper / src, output / dst)
    for name in ('cover-letter.txt', 'highlights.txt'):
        shutil.copy2(original / name, output / name)
    cover = output / 'cover-letter.txt'
    cover.write_text(cover.read_text().rstrip() + '\n' + contact['email'] + '\n')
    source_zip = output / 'manuscript-source.zip'
    with zipfile.ZipFile(source_zip, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        for p in paper.rglob('*'):
            if p.is_file() and p.suffix in ('.tex', '.bib', '.pdf'):
                z.write(p, str(p.relative_to(paper)))
    archive = subprocess.check_output(['git', 'archive', commit], cwd=repo)
    with tarfile.open(fileobj=io.BytesIO(archive)) as t:
        source = {m.name: t.extractfile(m).read() for m in t.getmembers() if m.isfile()}
    program = output / 'LBFlow-program-and-data.zip'
    program_readme = f'''LBFlow CPC program and research data

Candidate source commit: {commit}
Original frozen scientific source: {verified['source_commit']}

README.md contains build, test and execution instructions. examples/ and
sample-runs/ contain self-contained inputs, verified outputs and run receipts.
The channel and drop samples ran with two MPI processes. The production
executable hash is recorded in sample-runs/verification.json. The numerical
scope, optional controls and original-source provenance remain documented.

raw/major-revision/ contains the complete frozen major-revision evidence.
Its raw-manifest.json checks all {len(manifest['files'])} files. Reanalysis uses
raw/major-revision/analysis-source/collect_major_revision.py. The original
verified historical reproduction ZIP remains unchanged under historical/.
Additional complete-mode response checks and their portable driver are in
research/drlab-publication/reviews/second-review-20261009/.

The manuscript submitted for review is supplied separately. Committed source
contains country-only affiliation; personal correspondence details are used
only in the separate private manuscript package and submission form.
Historical archives preserve original bytes and may contain former contact
metadata. This submission package is private and is not a public release.
No account credentials are included. Code is MIT-licensed; attributed
Bentheimer derivatives retain CC BY 4.0 and third-party notices are preserved.
'''
    with zipfile.ZipFile(previous) as old, zipfile.ZipFile(
            program, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        for name, blob in source.items():
            z.writestr(name, blob)
        for name in manifest['files']:
            archive_name = 'raw/major-revision/' + name
            with old.open(archive_name) as src, z.open(archive_name, 'w') as dst:
                shutil.copyfileobj(src, dst, length=1024*1024)
        z.write(root / 'major-revision-raw-manifest.json', 'raw/major-revision/raw-manifest.json')
        historical = 'historical/pre-major-revision-reproducibility.zip'
        info = zipfile.ZipInfo(historical)
        info.compress_type = zipfile.ZIP_STORED
        info.file_size = old.getinfo(historical).file_size
        with old.open(historical) as src, z.open(info, 'w', force_zip64=True) as dst:
            shutil.copyfileobj(src, dst, length=1024*1024)
        for p in samples.rglob('*'):
            if p.is_file():
                z.write(p, 'sample-runs/' + str(p.relative_to(samples)))
        z.writestr('CPC-program-readme.txt', program_readme)
        z.writestr('submission-source.json', json.dumps(dict(
            source_commit=commit, previous_scientific_source=verified['source_commit'],
            canonical_data_sha256=verified['canonical_data_sha256'],
            raw_manifest_sha256=verified['raw_manifest_sha256'],
            scope='Current source and unchanged frozen evidence; no new full publication LBM/FV trajectories.',
            software_license='MIT', data_license='CC BY 4.0 for Bentheimer derivatives'), indent=2)+'\n')
    for path in (source_zip, program):
        with zipfile.ZipFile(path) as z:
            if z.testzip() is not None:
                raise AssertionError('ZIP CRC failure: ' + str(path))
    with zipfile.ZipFile(program) as z:
        for name, blob in source.items():
            if z.read(name) != blob:
                raise AssertionError('Committed source differs: ' + name)
        for name, entry in manifest['files'].items():
            with z.open('raw/major-revision/' + name) as stream:
                if stream_digest(stream) != entry['sha256']:
                    raise AssertionError('Frozen raw digest differs: ' + name)
        with z.open(historical) as stream:
            if stream_digest(stream) != verified['previous_archive_sha256']:
                raise AssertionError('Historical archive changed')
    result = dict(source_commit=commit, reviewed_commit='f274c1e0555f06cccac8a3d1c60e0dd0ef8f9334',
                  original_scientific_source=verified['source_commit'],
                  packages=[dict(file=p.name, bytes=p.stat().st_size, sha256=digest(p))
                            for p in (source_zip, program)],
                  manuscripts={p.name:dict(bytes=p.stat().st_size, sha256=digest(p))
                               for p in (output/'manuscript.pdf', output/'supplement.pdf')},
                  raw_files=len(manifest['files']), original_archives_preserved=True,
                  private_address_inserted=True, credentials_in_archives=False,
                  author_originality_confirmed=True, author_exclusive_submission_confirmed=True,
                  checks='Both ZIP CRCs, all committed source bytes, all frozen raw digests, historical archive identity and sample input/output hashes verified.',
                  journal_submission_made=False, external_actions='Local preparation only')
    (output/'verification.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2), flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--program-snapshot', type=Path, required=True)
    parser.add_argument('--contact-file', type=Path, required=True)
    parser.add_argument('--samples', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--tex-engine', type=Path, required=True)
    args = parser.parse_args()
    prepare(args.repo.resolve(), args.program_snapshot.resolve(), args.contact_file.resolve(),
            args.samples.resolve(), args.output.resolve(), args.tex_engine.resolve())


if __name__ == '__main__':
    main()
