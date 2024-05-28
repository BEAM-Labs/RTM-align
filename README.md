# RTM-align
Official code for RTM-align: an improved RNA alignment tool with enhanced short sequence performance via post-standardization and fragment alignment

DISCLAIMER:
Permission to use, copy, modify, and distribute this program for any
purpose, with or without fee, is hereby granted, provided that the
notices on the head, the reference information, and this copyright
notice appear in all copies or substantial portions of the Software.
It is provided "as is" without express or implied warranty.

This program is developed based on RNA-align: https://zhanggroup.org/RNA-align/download.html

## Getting Started

### Dataset

All dataset considered can be found at https://doi.org/10.6084/m9.figshare.25903405.v2

### Install

To install RTM-align, follow these steps:

1. **compiler preparation**

requirement for compiling RTM-align: gcc(12.2.0 in our environment), python(3.10.x in our environment for visualizing fragment alignment only)

2. **clone the repo**

```shell
git clone https://github.com/BEAM-Labs/RTM-align.git
cd RTM-align/code
make RTMalign
```

The "-static" flag in MakeFile should be removed on Mac OS, which does not support building static executables.

### Use RTM-align

Briefly, to alignment two single-chain structures (1.pdb and 2.pdb), enter the following:

```shell
./RTMalign 1.pdb 2.pdb
```

If want to visualize fragment alignment (remember to change filex and filey in run.sh):

```shell
bash ./run.sh
```

You can run the program without arguments to obtain a brief instruction.
Full document for all available options can be obained by:

```shell
./RTMalign -h
```