#!/usr/bin/env python3
import numpy as np
from pathlib import Path

nx, ny, nz = 32, 32, 32
xmin, xmax = 0.0, 1.0
ymin, ymax = 0.0, 1.0
zmin, zmax = 0.0, 1.0

Bx, By, Bz = 0.0, 0.0, 1.0
Ex, Ey, Ez = 0.1, 0.0, 0.0

x = np.linspace(xmin, xmax, nx)
y = np.linspace(ymin, ymax, ny)
z = np.linspace(zmin, zmax, nz)

out = Path('fields_uniform.txt')
with out.open('w') as f:
    f.write('# x y z Bx By Bz Ex Ey Ez\n')
    # The Fortran reader assumes i fastest, then j, then k.
    for zk in z:
        for yj in y:
            for xi in x:
                f.write(f'{xi:.16e} {yj:.16e} {zk:.16e} {Bx:.16e} {By:.16e} {Bz:.16e} {Ex:.16e} {Ey:.16e} {Ez:.16e}\n')

print(f'Wrote {out} with {nx*ny*nz} rows')
print('Columns: x y z Bx By Bz Ex Ey Ez')
