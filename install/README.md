PAT Install Script (install.py)
===============================

This script compiles and installs the required binary tools used by the PAT pipeline.

Requirements
------------
- Python 3.10  (required)
- A working conda environment (recommended)
- g++ (version 8 or newer)
- htslib and zlib development headers available in the environment

A typical setup uses conda:

    conda create -n patenv python=3.10 -y
    conda activate patenv
    conda install -y htslib=1.21 zlib eigen

Usage
-----
From the root directory of the PAT repository, run:

    python install/install.py --prefix /path/to/scripts/

The `--prefix` path should match the value of `"ScriptFolder"` in your Snakemake config JSON.  
For example, if your config contains:

    "ScriptFolder": "/home/user/PAT/scripts/"

then run:
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
    python install/install.py --prefix /home/user/PAT/scripts/

What the Script Does
--------------------
- Detects your active conda environment
- Compiles each required C++ tool located in `scripts/src/`
- Copies the resulting binaries into the folder specified by `--prefix`

You should see output such as:

    [BUILD] kmercounter8 not found, building...
    [DONE] Copied kmercounter8 → /home/user/PAT/scripts/kmercounter8

Troubleshooting
---------------
If you see errors like:

    fatal error: htslib/sam.h: No such file or directory

then install headers in your conda env:

    conda install -y htslib=1.21 zlib

If compilation fails, verify your compiler:

    g++ --version   (must be ≥ 8)

----------------------------------

After installation, you can proceed to run the Snakemake pipeline as described in the main documentation.

