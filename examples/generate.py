#!/usr/bin/env python3
"""Generate a small self-contained channel or two-phase drop input directory."""
import argparse
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('case', choices=('channel', 'drop'))
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    folder = Path(args.output)
    folder.mkdir(parents=True, exist_ok=False)
    channel = args.case == 'channel'
    nx, ny, nz = (8, 10, 4) if channel else (24, 24, 1)
    labels = bytes(255 if y in (0, ny-1) else 75 for z in range(nz)
                   for y in range(ny) for x in range(nx)) if channel else bytes(
        100 if (x-11.5)**2+(y-11.5)**2 <= 25 else 0 for y in range(ny) for x in range(nx))
    (folder/'geometry.raw').write_bytes(labels)
    dimensions = f'{nz} {ny} {nx}' if channel else f'{ny} {nx}'
    model = 'scsp' if channel else 'color'
    tau = '1' if channel else '1 1'
    density = '1' if channel else '1 1'
    indices = ('75 "GV" "0.000001" "1" "0"\n255 "SD" "90" "N/A" "N/A"' if channel else
               '0 "VD" "N/A" "N/A" "0"\n100 "VD" "N/A" "N/A" "1"')
    (folder/'input.prm').write_text(f'''gridFile "geometry.raw";
dimension {dimensions};
directions {19 if channel else 9};
fileType "uint8";
latticeLen 1;
startTime 0;
endTime 1000;
writeInterval 250;
reportInterval 250;
restartInterval 500;
stopError 1e-8;
outPath "field";
logPath ".";
reportDir 0;
reportSlice 1 {nx-2};
fluidModel "{model}";
equilModel "MRT";
rexlationTime {tau};
latticeSpeed 1;
alpha0 {1/3 if channel else 4/9};
density {density};
contactAngle 90;
surfaceTension 0.005;
interfaceThickness 0.7;
wallGradient "legacy";
viscosityInterpolation "arithmetic";
{indices}
''')
    print(f'Generated {nx} x {ny} x {nz} {args.case} in {folder}. Run LBFlow input.prm from that directory.')


if __name__ == '__main__':
    main()
