# autoHFIX
Adds **AFIX** instructions for hydrogen atoms to **SHELX** input files based on geometry and residual desity automatically

**Usage:**

autoHFIX.exe [-options] filename.res

**Options:**

**-HKLMAX** maximum reflections to be read from *filename.hkl* file

**Output:**

*filename.fcf6* fcf file in LIST 6 style
*filename.ins* SHELX instructions file with AFIX instructions for hydrogren atoms added.


