#!/usr/bin/env python3
import argparse,uuid
from pathlib import Path
import subprocess,sys
p=argparse.ArgumentParser()
p.add_argument('--driver',required=True);p.add_argument('--root',required=True)
a=p.parse_args()
root=Path(a.root)/('moments-'+uuid.uuid4().hex[:10])
raise SystemExit(subprocess.call([sys.executable,str(Path(__file__).with_name('moment_audit.py')),
                                 '--driver',a.driver,'--root',str(root)]))
