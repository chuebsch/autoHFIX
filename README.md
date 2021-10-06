# autoHFIX
Adds **AFIX** instructions for hydrogen atoms to **SHELX** input files based on geometry and residual desity automatically

## Usage:

autoHFIX.exe [-options] filename.res

## Options:

**-HKLMAX** maximum hkl index [200] to be read from *filename.hkl* file

## Output files:

* filename.fcf6 fcf file in LIST 6 style
* filename.ins SHELX instructions file with AFIX instructions for hydrogren atoms added.

## Typical Output:

```
"qqqbtp02_a.res"
R1=0.095831 3343 3519 2247.45 23452.2 384.000000
loadFouAndPerform 266 {qqqbtp02_a.fcf6}
qqqbtp02_a.fcf6
sigma 0.105698 max 0.875535 min -0.415202
Uniq Reflections: 3344
Fourier grid dimensions: 180 X 48 X 90 = 777600 grid points.

autoHfix ...74 35

?autoHfix 82 70
test1

```

