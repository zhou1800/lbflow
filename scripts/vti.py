"""Read and validate the solver's appended-raw, little-endian VTI output."""
import struct
import xml.etree.ElementTree as ET

def read_vti(path):
    data = path.read_bytes()
    pos = data.index(b'<AppendedData')
    doc = ET.fromstring(data[:pos] + b'</VTKFile>')
    begin = data.index(b'>_', pos) + 2
    im = doc.find('ImageData')
    extent = list(map(int, im.attrib['WholeExtent'].split()))
    nx, ny, nz = [extent[i+1] - extent[i] + 1 for i in (0, 2, 4)]
    total = nx * ny * nz
    fields = {}
    coverage = [0] * total
    for piece in im.findall('Piece'):
        e = list(map(int, piece.attrib['Extent'].split()))
        ids = [(z * ny + y) * nx + x for z in range(e[4], e[5]+1)
               for y in range(e[2], e[3]+1) for x in range(e[0], e[1]+1)]
        for i in ids:
            coverage[i] += 1
        for a in piece.find('PointData').findall('DataArray'):
            name = a.attrib['Name']
            n = int(a.attrib.get('NumberOfComponents', '1'))
            off = begin + int(a.attrib['offset'])
            length = struct.unpack_from('<Q', data, off)[0]
            expected = len(ids) * n * 8
            if length != expected:
                raise ValueError(f'{path.name}: {name}: block {length} != {expected}')
            values = struct.unpack_from('<' + 'd' * (length // 8), data, off + 8)
            field = fields.setdefault(name, {'components': n, 'values': [None] * (total * n)})
            for j, i in enumerate(ids):
                field['values'][i*n:(i+1)*n] = values[j*n:(j+1)*n]
    if any(c != 1 for c in coverage):
        raise ValueError('Pieces do not cover whole grid exactly once')
    if any(v is None for f in fields.values() for v in f['values']):
        raise ValueError('Missing field values')
    return (nx, ny, nz), fields
