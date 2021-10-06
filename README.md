# autoHFIX
Adds **AFIX** instructions for hydrogen atoms to **SHELX** input files based on geometry and residual desity automatically

## Usage:

autoHFIX.exe [-options] filename.res

## Options:

* **-HKLMAX** maximum hkl index [200] to be read from *filename.hkl* file
* **-recheck* loads the resuling ins file to check the effect of the hydrogen atoms added to R value and residual density.

## Output files:

* filename.fcf6 fcf file in LIST 6 style
* filename.ins SHELX instructions file with AFIX instructions for hydrogren atoms added.

## Typical Output:

```
autoHFIX.exe qqqbtp02_a.res
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
or
```
autoHFIX.exe -recheck -HKLMAX 25 qqqbtp02_a.res
"qqqbtp02_a.res"
R1=0.095831 3343 3519 2247.45 23452.2 384.000000
qqqbtp02_a.fcf6
-24<h<24  0<h<8  0<h<14
sigma 0.104637 max 0.865585 min -0.392632
Uniq Reflections: 3249
Fourier grid dimensions: 144 X 48 X 90 = 622080 grid points.

autoHfix ...74 35

?autoHfix 82 70
"qqqbtp02_a.ins"
R1=0.075967 3343 3519 1781.78 23454.6 428.000000
qqqbtp02_a.fcf6
-24<h<24  0<h<8  0<h<14
sigma 0.077720 max 0.602859 min -0.404140
Uniq Reflections: 3249
Fourier grid dimensions: 144 X 48 X 90 = 622080 grid points.

```

