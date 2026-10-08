#!/usr/bin/env python3
"""Executable pore-flow regression against independent analytic references."""
import argparse
from pathlib import Path
import uuid

from benchmarks import Study
from porous_media import bundle_study, entry_study


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solver', required=True)
    parser.add_argument('--root', required=True)
    args = parser.parse_args()
    study = Study(args.solver, Path(args.root)/('pores-'+uuid.uuid4().hex[:10]))
    bundle_study(study, scales=(1,), dimensions=(2,))
    entry_study(study, scales=(1,), walls=('legacy',))
    print('Passed analytic phase flux/viscosity and capillary-entry direction checks.')


if __name__ == '__main__':
    main()
