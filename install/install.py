#!/usr/bin/env python3

import sys
import subprocess
import shutil
import re
from pathlib import Path
import os
import shutil as _shutil  # avoid name collision with std shutil above

# =========================
#  Version helpers
# =========================

def ensure_python_version():
    required = (3, 10)
    if sys.version_info < required:
        sys.stderr.write(
            f"[ERROR] Python >= {required[0]}.{required[1]} required, "
            f"found {sys.version_info.major}.{sys.version_info.minor}\n"
        )
        sys.exit(1)
    else:
        print(f"[CHECK] Python version OK: {sys.version_info.major}.{sys.version_info.minor}")


def ensure_gcc_version():
    try:
        out = subprocess.check_output(["gcc", "--version"], text=True)
    except FileNotFoundError:
        sys.stderr.write("[ERROR] gcc not found in PATH. Please install gcc >= 8.\n")
        sys.exit(1)

    first_line = out.splitlines()[0]
    # Try to find something like 11.4.0
    m = re.search(r"(\d+)\.(\d+)\.(\d+)", first_line)
    if not m:
        sys.stderr.write(f"[WARN] Could not parse gcc version from: {first_line}\n")
        return

    major = int(m.group(1))
    if major < 8:
        sys.stderr.write(f"[ERROR] gcc >= 8 required, found {major}.x\n")
        sys.exit(1)
    else:
        print(f"[CHECK] gcc version OK: {first_line.strip()}")


def parse_simple_version(v: str):
    """
    Parse a version like '1.19.1' -> (1, 19, 1).
    If it's weird, fall back to 0.
    """
    parts = v.strip().split(".")
    nums = []
    for p in parts:
        m = re.match(r"(\d+)", p)
        if m:
            nums.append(int(m.group(1)))
        else:
            nums.append(0)
    return tuple(nums)


def version_cmp(a: str, b: str) -> int:
    """
    Compare versions a and b.
    Returns:
      -1 if a < b
       0 if a == b
      +1 if a > b
    """
    ta = parse_simple_version(a)
    tb = parse_simple_version(b)
    n = max(len(ta), len(tb))
    ta = ta + (0,) * (n - len(ta))
    tb = tb + (0,) * (n - len(tb))
    if ta < tb:
        return -1
    elif ta > tb:
        return 1
    else:
        return 0


def version_ge(found: str, required: str) -> bool:
    return version_cmp(found, required) >= 0


def version_in_range(found: str, low: str, high: str) -> bool:
    """
    Return True if low <= found <= high.
    """
    return version_cmp(found, low) >= 0 and version_cmp(found, high) <= 0


# =========================
#  Conda helpers
# =========================

def ensure_conda_available():
    if shutil.which("conda") is None:
        sys.stderr.write(
            "[ERROR] 'conda' not found in PATH.\n"
            "Please install Miniconda/Anaconda and ensure 'conda' is available, "
            "then re-run this script.\n"
        )
        sys.exit(1)
    else:
        print("[CHECK] conda found.")


def conda_has_package(name: str):
    """
    Return (installed, version_str or None).
    """
    try:
        res = subprocess.run(
            ["conda", "list", name],
            capture_output=True,
            text=True,
            check=False,
        )
    except FileNotFoundError:
        return False, None

    if res.returncode != 0:
        return False, None

    for line in res.stdout.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        cols = line.split()
        # expect: name  version  build  channel
        if cols[0] == name and len(cols) >= 2:
            found_ver = cols[1]
            return True, found_ver

    return False, None


def conda_install_packages(specs):
    """
    specs: list like ["htslib=1.21", "eigen>=3.4", "snakemake=6.15.1", "blast", ...]
    Use conda-forge + bioconda explicitly.
    """
    cmd = [
        "conda", "install", "-y",
        "-c", "conda-forge",
        "-c", "bioconda",
    ] + specs
    print(f"[INSTALL] Running: {' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        sys.stderr.write(f"[ERROR] conda install failed: {e}\n")
        sys.exit(1)


def ensure_conda_deps():
    """
    htslib: accept 1.19–1.21; otherwise install 1.21
    eigen: require >= 3.4
    snakemake: exactly 6.15.1
    """
    ensure_conda_available()

    needed_specs = []

    # ---- htslib: require exactly 1.21 ----
    have_hts, hts_ver = conda_has_package("htslib")
    if have_hts and hts_ver is not None:
        if version_cmp(hts_ver, "1.21") == 0:
            print(f"[CHECK] htslib version OK: {hts_ver} == 1.21")
        else:
            print(f"[WARN] htslib {hts_ver} != 1.21; will install htslib=1.21")
            needed_specs.append("htslib=1.21")
    else:
        print("[WARN] htslib not found; will install htslib=1.21")
        needed_specs.append("htslib=1.21")

    # ---- eigen: require >= 3.4 ----
    have_eigen, eigen_ver = conda_has_package("eigen")
    if have_eigen and eigen_ver is not None:
        if version_ge(eigen_ver, "3.4"):
            print(f"[CHECK] eigen version OK: {eigen_ver} >= 3.4")
        else:
            print(f"[WARN] eigen {eigen_ver} < 3.4; will install eigen>=3.4")
            needed_specs.append("eigen>=3.4")
    else:
        print("[WARN] eigen not found; will install eigen>=3.4")
        needed_specs.append("eigen>=3.4")

    # ---- snakemake: must be EXACT 6.15.1 ----
    have_snake, snake_ver = conda_has_package("snakemake")
    if have_snake and snake_ver is not None:
        if version_cmp(snake_ver, "6.15.1") == 0:
            print(f"[CHECK] snakemake version OK: {snake_ver} == 6.15.1")
        else:
            print(f"[WARN] snakemake {snake_ver} != 6.15.1; will install snakemake=6.15.1")
            needed_specs.append("snakemake=6.15.1")
    else:
        print("[WARN] snakemake not found; will install snakemake=6.15.1")
        needed_specs.append("snakemake=6.15.1")

    if needed_specs:
        print(f"[INFO] Missing or incompatible packages, will install via conda: {', '.join(needed_specs)}")
        conda_install_packages(needed_specs)
    else:
        print("[CHECK] htslib, eigen, and snakemake satisfy version requirements.")

    # ---- numpy: required ----
    have_numpy, numpy_ver = conda_has_package("numpy")
    if have_numpy and numpy_ver is not None:
        print(f"[CHECK] numpy already installed (version {numpy_ver})")
    else:
        print("[WARN] numpy not found; will install numpy")
        needed_specs.append("numpy")

    # ---- pandas: required ----
    have_pandas, pandas_ver = conda_has_package("pandas")
    if have_pandas and pandas_ver is not None:
        print(f"[CHECK] pandas already installed (version {pandas_ver})")
    else:
        print("[WARN] pandas not found; will install pandas")
        needed_specs.append("pandas")

    if needed_specs:
        print(f"[INFO] Missing or incompatible packages, will install via conda: {', '.join(needed_specs)}")
        conda_install_packages(needed_specs)
    else:
        print("[CHECK] htslib, eigen, snakemake, numpy, and pandas satisfy requirements.")



def ensure_conda_bio_tools():
    """
    Ensure these are installed via conda (they provide the binaries):
        blast  -> blastn, makeblastdb
        winnowmap
        bedtools
        samtools
        minimap2
    """
    ensure_conda_available()

    # Map binary -> conda package (we don't actually check binaries here,
    # just make sure the packages are installed).
    pkgs = ["blast", "winnowmap", "bedtools", "samtools", "minimap2"]

    needed_specs = []
    for pkg in pkgs:
        have_pkg, ver = conda_has_package(pkg)
        if have_pkg:
            print(f"[CHECK] {pkg} already installed (version {ver})")
        else:
            print(f"[WARN] {pkg} not found; will install {pkg}")
            needed_specs.append(pkg)

    if needed_specs:
        print(f"[INFO] Installing bio tools via conda: {', '.join(needed_specs)}")
        conda_install_packages(needed_specs)
    else:
        print("[CHECK] All bio tool conda packages are present.")


# =========================
#  Build binaries + external tool check
# =========================
def ensure_binaries_exist(ScriptFolder: Path):
    """
    Ensure required binaries are compiled and present in ScriptFolder.
    If missing, build them from their respective subfolders.
    """
    ScriptFolder = Path(ScriptFolder)
    tools = {
        "kmercounter8": "kmercounter",
        "kmer_selector": "kmerselector",
        "kmerpartition": "kmerpartition",
        "kmernorm": "kmernorm",
        "KmerStrd": "kmerstrd",
        "kmertree": "kmertree",
    }

    for binary_name, subfolder_name in tools.items():
        binary_path = ScriptFolder / binary_name

        if binary_path.exists():
            print(f"[CHECK] Binary found: {binary_path}")
            continue

        # Try to build
        source_folder = ScriptFolder / "src" / subfolder_name
        print(f"[BUILD] {binary_name} not found, attempting to build from {source_folder}...")

        if not source_folder.exists():
            sys.stderr.write(f"[ERROR] Source folder not found: {source_folder}\n")
            sys.exit(1)

        # Inherit env and inject conda include/lib paths if available
        env = os.environ.copy()
        conda_prefix = env.get("CONDA_PREFIX")
        if conda_prefix:
            inc = f"{conda_prefix}/include"
            lib = f"{conda_prefix}/lib"

            # Make the compiler see htslib/sam.h and zlib.h in $CONDA_PREFIX/include
            for var in ("CPLUS_INCLUDE_PATH", "CPATH"):
                env[var] = inc + (":" + env[var] if var in env and env[var] else "")

            # Make linker see libhts, libz in $CONDA_PREFIX/lib
            for var in ("LIBRARY_PATH", "LD_LIBRARY_PATH"):
                env[var] = lib + (":" + env[var] if var in env and env[var] else "")

            print(f"[INFO] Using CONDA_PREFIX={conda_prefix} for includes/libs")

        try:
            subprocess.run(["make"], cwd=source_folder, check=True, env=env)
        except subprocess.CalledProcessError:
            sys.stderr.write(f"[ERROR] Failed to build {binary_name} in {source_folder}\n")
            sys.exit(1)

        # After build, move the binary to ScriptFolder if it exists in the subfolder
        built_binary = source_folder / binary_name
        if built_binary.exists():
            _shutil.copy2(built_binary, ScriptFolder)
            print(f"[DONE] Copied {binary_name} → {binary_path}")
        else:
            sys.stderr.write(f"[ERROR] Binary {binary_name} not found after build in {source_folder}\n")
            sys.exit(1)


def precheck():
    # Tools to check (binaries in PATH)
    required_tools = [
        "blastn", "makeblastdb",   # BLAST
        "winnowmap",               # Winnowmap
        "bedtools",                # Bedtools
        "samtools",                # Samtools
        "minimap2",                # Minimap2
    ]

    missing = []
    for tool in required_tools:
        if _shutil.which(tool) is None:
            missing.append(tool)

    if missing:
        sys.stderr.write(f"\n[ERROR] Missing required tools in PATH: {', '.join(missing)}\n")
        sys.stderr.write("They should have been installed via conda; "
                         "please ensure you activated the correct environment.\n")
        sys.exit(1)
    else:
        print("[CHECK] All required tools found in PATH.")


# =========================
#  Main
# =========================

def main():
    # ScriptFolder is ../scripts/ relative to this install.py
    ScriptFolder = (Path(__file__).resolve().parent / ".." / "scripts").resolve()

    print("=== PATs pipeline installer / prechecker ===")
    print(f"[INFO] Using ScriptFolder: {ScriptFolder}")

    ensure_python_version()
    ensure_gcc_version()
    ensure_conda_deps()      # htslib/eigen/snakemake
    ensure_conda_bio_tools() # blast / winnowmap / bedtools / samtools / minimap2
    precheck()               # verify binaries in PATH
    ensure_binaries_exist(ScriptFolder)

    print("\n[OK] Environment looks good and all binaries are present.")
    print("You should be ready to run the pipeline.")


if __name__ == "__main__":
    main()

