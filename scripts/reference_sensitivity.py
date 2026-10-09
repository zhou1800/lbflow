#!/usr/bin/env python3
"""Quantify continuum reconstruction sensitivity without rerunning the solver.

The envelope is a comparison of specified reconstructions, not a confidence
interval or an estimate of an isolated lattice truncation error. Requires SciPy.
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator

from benchmarks import sha
from interfacial_transport import diffuse_reference, viscosity


def analyze(data):
    records=[]
    for m in data['measurements']:
        height=m['height']
        y=np.arange(height)+.5
        phi=np.array(m['phase0_profile'])
        observed=np.array(m['profile'])
        mu=[(tau-.5)/3 for tau in m['tau']]
        references={}
        for kind,subdivisions in [('linear-256',256),('linear-1024',1024),('pchip-1024',1024)]:
            grid=np.linspace(0,height,height*subdivisions+1)
            fraction=(np.interp(grid,y,phi) if kind.startswith('linear') else
                      PchipInterpolator(y,phi)(np.clip(grid,y[0],y[-1])))
            ref=diffuse_reference(y,grid,viscosity(fraction,mu,m['closure']),m['force'])
            error=float(np.linalg.norm(observed-ref['profile'])/np.linalg.norm(ref['profile']))
            references[kind]=dict(profile=ref['profile'].tolist(),relative_profile_error=error,
                                  continuous_volume_flux=ref['volume_flux'],
                                  phase_fraction_range=[float(fraction.min()),float(fraction.max())])
        if not np.allclose(references['linear-256']['profile'],m['measured_diffuse_profile'],
                           rtol=1e-11,atol=1e-15):
            raise AssertionError('The distributed baseline reconstruction was not reproduced')
        profiles={key:np.array(value['profile']) for key,value in references.items()}
        envelope=[min(r['relative_profile_error'] for r in references.values()),
                  max(r['relative_profile_error'] for r in references.values())]
        records.append(dict(name=m['name'],height=height,beta=m['beta'],closure=m['closure'],
                            references=references,profile_error_envelope=envelope,
                            quadrature_relative_change=float(np.linalg.norm(
                                profiles['linear-256']-profiles['linear-1024'])/np.linalg.norm(observed)),
                            interpolation_relative_change=float(np.linalg.norm(
                                profiles['linear-1024']-profiles['pchip-1024'])/np.linalg.norm(observed))))
    identities=[]
    position=np.linspace(-12,12,20001)
    for ratio in (.25,1.,3.,9.,27.):
        fraction=(1-np.tanh(position))/2
        inverse_a=1/viscosity(fraction,(1.,ratio),'arithmetic')
        shifted=(1-np.tanh(position+.5*math.log(ratio)))/2
        inverse_h=1/viscosity(shifted,(1.,ratio),'harmonic')
        identities.append(dict(ratio=ratio,interface_shift_over_xi=-.5*math.log(ratio),
                               max_absolute_identity_error=float(np.max(abs(inverse_a-inverse_h)))))
    return dict(protocol='Same measured cell-center phase samples; constant wall continuation; linear and shape-preserving PCHIP interpolation; 256/1024 quadrature intervals per cell.',
                scope='Reference sensitivity only. These reconstructions do not bound every possible continuum profile or isolate lattice error. Identity checks are analytic postprocessing, not simulator executions.',
                records=records,shift_identity=identities)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data',default=str(Path(__file__).resolve().parents[1]/
                        'research/drlab-publication/paper/data/interfacial-transport.json'))
    parser.add_argument('--output',required=True)
    args=parser.parse_args()
    result=analyze(json.loads(Path(args.data).read_text()))
    result['input_sha256']=sha(args.data)
    output=Path(args.output)
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(records=len(result['records']),
                          max_quadrature_relative_change=max(r['quadrature_relative_change'] for r in result['records']),
                          max_interpolation_relative_change=max(r['interpolation_relative_change'] for r in result['records']))))


if __name__=='__main__':main()
