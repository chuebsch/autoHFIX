# autoHFIX
Adds **AFIX** instructions for hydrogen atoms to **SHELX** input files based on geometry and residual desity automatically

## Usage:

autoHFIX.exe [-options] filename.res

## Options:

* **-HKLMAX** maximum hkl index [200] to be read from *filename.hkl* file
* **-recheck** loads the resuling ins file to check the effect of the hydrogen atoms added to R value and residual density.
* **-aula** trys to auto label the structure first no constraits restraints or AFIX instructions may be present before. 

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
or
```
autoHFIX.exe -aula C:\Users\Christian.Huebschle\Desktop\BrukeUsersMeeting\example-1\vitc_a.res
"C:\Users\Christian.Huebschle\Desktop\BrukeUsersMeeting\example-1\vitc_a.res"
R1=0.052446 1239 1243 689.218 13141.4 336.000000
C:\Users\Christian.Huebschle\Desktop\BrukeUsersMeeting\example-1\vitc_a.fcf6
-7<h<7  0<h<7  0<h<19
sigma 0.079392 max 0.583505 min -0.312173
Uniq Reflections: 1240
Fourier grid dimensions: 45 X 45 X 120 = 243000 grid points.

autoHfix ...83 44
It has AFIX = 0 instructions! It has 0 restraints or constraints!
TEST
========Au=La=====================================================================
        .-"-~-"-.
       /.-"-.-"-.\
       ||((o|o))||
       )\__/V\__/(
      / ~ -...- ~ \
     |\` ~. ~ .~ `/|
     | `~ - ^ - ~` |
     | ;  '  :  .  |
      \ . : '  ; '/
       '.   ; ' .'
        `uu---uu`
atoms24 tripel34
Rings 6
Rings 2 are real and complete.
Molecule#0 has 12 atoms
+C00G=O003=C00N=C00M=C00K+
+========================+
Ring #1 size = 5/5 is complete 1
+C00G=O003=C00N=C00M=C00K+
+========================+
AWAYOFCOM 0.523389 1  test-0.498627
Rings 2
+C00G=O003=C00N=C00M=C00K+
+========================+
most far AWAY from COM 0.523389 1
Molecule#1 has 12 atoms
+C00E=O009=C00J=C00L=C00I+
+========================+
Ring #0 size = 5/5 is complete 1
+C00E=O009=C00J=C00L=C00I+
+========================+
AWAYOFCOM 0.555706 0  test-0.54401
Rings 2
+C00E=O009=C00J=C00L=C00I+
+========================+
most far AWAY from COM 0.555706 0
{C00N--C00M--C00K--C00G--O003--C00F--C00H--O00C--O00A--O004--O007--O001--C00J--C00L--C00I--C00E--O009--C00D--C00O--O00B--O008--O005--O006--O002--}
{C1~C2~C3~C4~O1~C5~C6~O2~O3~O4~O5~O6~C7~C8~C9~C10~O7~C11~C12~O8~O9~O10~O11~O12~}
aula ends 24 24
"C:\Users\Christian.Huebschle\Desktop\BrukeUsersMeeting\example-1\vitc_a.res"
R1=0.052446 1239 1243 689.218 13141.4 336.000000
C:\Users\Christian.Huebschle\Desktop\BrukeUsersMeeting\example-1\vitc_a.fcf6
-7<h<7  0<h<7  0<h<19
sigma 0.079392 max 0.583505 min -0.312173
Uniq Reflections: 1240
Fourier grid dimensions: 45 X 45 X 120 = 243000 grid points.

before autoHfix lines 107 atoms 44

 after autoHfix lines 121 atoms 60
```
