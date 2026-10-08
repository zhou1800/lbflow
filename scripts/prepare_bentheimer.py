#!/usr/bin/env python3
"""Prepare attributed, checksum-pinned crops of the public Bentheimer image.

The source data retain CC BY 4.0; the repository's MIT code license does not
replace the data license. The fixed crop is a numerical demonstration and
does not establish a representative elementary volume or laboratory property.
"""
import argparse
import hashlib
import json
from pathlib import Path
import urllib.request
import zipfile

import numpy as np
from scipy import ndimage


SOURCE = dict(
    title='Bentheimer Sandstone',
    authors=['Saurabh M. Shah', 'John P. Crawshaw', 'Edo S. Boek'],
    doi='10.6084/m9.figshare.1200096.v1',
    landing_page='https://figshare.com/articles/dataset/Bentheimer_Sandstone/1200096',
    archive_url='https://ndownloader.figshare.com/files/1712678',
    archive_sha256='723f3734ae95978791148135ef8ac96b6d97d8f6b69cac83bd538ec9481a8a27',
    archive_bytes=10918859,
    raw_file='Benth_4.54.raw', raw_bytes=190109375,
    raw_sha256='7123e09e46119d23bd765ff6f2d5172d3387b3867b4dcd725f9b06ba494a76cd',
    shape_zyx=[575, 575, 575], voxel_size_m=4.54e-6,
    license='CC BY 4.0', license_url='https://creativecommons.org/licenses/by/4.0/',
    attribution='Data courtesy of Qatar Carbonates and Carbon Storage Research Centre, '
                'Imperial College London; dataset by Shah, Crawshaw and Boek.',
    labels={'pore': 0, 'grain': 255},
    storage='Unsigned 8-bit; x fastest, then y, then z; zero header offset.')


def digest(path):
    checksum = hashlib.sha256()
    with Path(path).open('rb') as stream:
        while block := stream.read(1 << 20):
            checksum.update(block)
    return checksum.hexdigest()


def verified_source(source=None, cache=None):
    if source:
        path = Path(source).resolve()
    else:
        cache = Path(cache).resolve()
        cache.mkdir(parents=True, exist_ok=True)
        path = cache/SOURCE['raw_file']
        if not path.exists():
            archive = cache/'Bentheimer_segmented.zip'
            if not archive.exists():
                temporary = archive.with_suffix('.download')
                request = urllib.request.Request(SOURCE['archive_url'], headers={'User-Agent': 'LBFlow-data-preparation/0.3'})
                with urllib.request.urlopen(request, timeout=120) as response, temporary.open('wb') as output:
                    while block := response.read(1 << 20):
                        output.write(block)
                if temporary.stat().st_size != SOURCE['archive_bytes'] or digest(temporary) != SOURCE['archive_sha256']:
                    raise ValueError('Downloaded archive does not match the pinned public dataset')
                temporary.rename(archive)
            if digest(archive) != SOURCE['archive_sha256']:
                raise ValueError('Cached archive checksum mismatch')
            with zipfile.ZipFile(archive) as compressed:
                # Extract only the exact member; do not unpack arbitrary archive paths.
                with compressed.open('Bentheimer_segmented/'+SOURCE['raw_file']) as stream, path.open('wb') as output:
                    while block := stream.read(1 << 20):
                        output.write(block)
    if path.stat().st_size != SOURCE['raw_bytes'] or digest(path) != SOURCE['raw_sha256']:
        raise ValueError('Source image checksum or byte count differs from the cited dataset')
    return path


def component_counts(pore, neighbors=6, periodic=False):
    structure = ndimage.generate_binary_structure(3, 1 if neighbors == 6 else 2)
    labels, number = ndimage.label(pore, structure)
    sizes = np.bincount(labels.ravel())
    sizes[0] = 0
    parent = list(range(number+1))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    if periodic:
        # Merge components connected across the periodic faces/edges. All
        # nonperiodic interior connections are already labelled above.
        for delta in np.argwhere(structure)-1:
            if not np.any(delta):
                continue
            opposite = np.roll(labels, tuple(delta), axis=(0, 1, 2))
            joined = (labels > 0) & (opposite > 0) & (labels != opposite)
            for left, right in np.unique(np.stack((labels[joined], opposite[joined]), axis=1), axis=0):
                left, right = find(int(left)), find(int(right))
                if left != right:
                    parent[right] = left
    grouped = {}
    for i in range(1, number+1):
        key = find(i)
        grouped[key] = grouped.get(key, 0)+int(sizes[i])
    return dict(components=len(grouped), largest_component_fraction=(
        max(grouped.values())/int(pore.sum()) if grouped else 0),
        neighborhood=neighbors, periodic=periodic)


def prepare(source, output, size=48, refinement=1):
    if not 8 <= size <= 575 or refinement < 1:
        raise ValueError('Unsupported crop size or refinement')
    source = verified_source(source)
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    begin = (575-size)//2
    image = np.memmap(source, mode='r', dtype=np.uint8, shape=(575,)*3)
    crop = np.array(image[begin:begin+size, begin:begin+size, begin:begin+size])
    if not np.isin(crop, (0, 255)).all():
        raise ValueError('The segmented source contains unexpected labels')
    pore = crop == 0
    # Periodic Euclidean distance to the voxel-center grain mask. This is an
    # initialization aid, not a calibrated geometrical pore/throat radius.
    distance = ndimage.distance_transform_edt(np.tile(pore, (3, 3, 3)))[
        size:2*size, size:2*size, size:2*size]
    nonwetting = pore & (distance >= 3.)
    pure = np.where(pore, 75, 255).astype(np.uint8)
    phases = np.where(pore, np.where(nonwetting, 75, 76), 255).astype(np.uint8)
    # Subdivision preserves the original binary voxel boxes and the initial
    # phase volumes exactly. It adds no new image information.
    if refinement > 1:
        for axis in range(3):
            pure = np.repeat(pure, refinement, axis=axis)
            phases = np.repeat(phases, refinement, axis=axis)
    pure.tofile(output/'single-phase.raw')
    phases.tofile(output/'two-phase.raw')
    metadata = dict(source=SOURCE, crop_origin_xyz=[begin]*3, crop_size_xyz=[size]*3,
        refinement=refinement, simulation_size_xyz=list(pure.shape[::-1]),
        simulation_voxel_size_m=SOURCE['voxel_size_m']/refinement,
        boundaries='Periodic on all axes; halfway bounce-back on every pore/grain link.',
        processing='Centered crop, no removal of isolated pores, no filtering or segmentation '
                   'changes; optional nearest-neighbor voxel subdivision.',
        initial_phase_rule='Phase 0 at periodic EDT >= 3 original voxels; phase 1 in remaining pores.',
        initial_phase0_pore_fraction=float(nonwetting.sum()/pore.sum()),
        porosity=float(pore.mean()),
        periodic_EDT_voxel_quantiles=dict(zip(('q10', 'q50', 'q90', 'q99'),
            np.quantile(distance[pore], (.1, .5, .9, .99)).tolist())),
        connectivity=[component_counts(pore, neighbors, periodic) for neighbors in (6, 18)
                      for periodic in (False, True)],
        files={name:dict(bytes=(output/name).stat().st_size, sha256=digest(output/name))
               for name in ('single-phase.raw', 'two-phase.raw')},
        scope='A cropped periodically repeated numerical specimen, not a laboratory permeability '
              'reference or representative elementary volume. Refinement subdivides image voxels.')
    (output/'manifest.json').write_text(json.dumps(metadata, indent=2)+'\n')
    (output/'DATA-LICENSE.txt').write_text(
        'Bentheimer image data and the derived crops: Creative Commons Attribution 4.0.\n'
        +SOURCE['license_url']+'\n'+SOURCE['landing_page']+'\nDOI: '+SOURCE['doi']+'\n'
        +SOURCE['attribution']+'\nChanges: '+metadata['processing']+'\n')
    return metadata


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', help='Local Benth_4.54.raw; checked against its published checksum')
    parser.add_argument('--cache', default='results/bentheimer-source')
    parser.add_argument('--output', required=True)
    parser.add_argument('--size', type=int, default=48)
    parser.add_argument('--refinement', type=int, default=1)
    args = parser.parse_args()
    source = verified_source(args.source, args.cache)
    print(json.dumps(prepare(source, args.output, args.size, args.refinement), indent=2))


if __name__ == '__main__':
    main()
