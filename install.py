#!/usr/bin/env python3
"""Install and validate the local target-driven PATs pipeline."""

from __future__ import annotations

import argparse
import os
import platform
import shlex
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence


ROOT = Path(__file__).resolve().parent
RUNTIME = ROOT
SCRIPT_DIR = RUNTIME / "scripts"
SOURCE_DIR = SCRIPT_DIR / "src"
REQUIREMENTS = ROOT / "requirements.txt"

PYTHON_MODULES = ("bitarray", "Cython", "numpy", "scipy", "snakemake", "pulp")
SNAKEMAKE_VERSION = "6.5.0"
PULP_VERSION = "2.7.0"
EXTERNAL_COMMANDS = ("samtools", "minimap2", "blastn", "makeblastdb", "winnowmap")
CONDA_PACKAGES = (
    "make",
    "cxx-compiler",
    "htslib",
    "zlib",
    "libcurl",
    "blast",
    "winnowmap",
    "samtools",
    "minimap2",
)


@dataclass(frozen=True)
class BuildTarget:
    source: str
    built_name: str
    installed_name: str
    make_target: str = "all"
    openmp: bool = False


BUILD_TARGETS = (
    BuildTarget("kmerconverter", "kmer_convertor8", "kmer_convertor8", "kmer_convertor8"),
    BuildTarget("kmersearcher", "kmer_searcher", "kmer_searcher"),
    BuildTarget("kmerpartition", "kmerpartition", "kmerpartition", openmp=True),
    BuildTarget("kmerselector", "kmer_selector", "kmer_selector"),
    BuildTarget("kmernorm", "kmernorm", "kmernorm"),
    BuildTarget("kmerstrd", "kmerstrd", "kmerstrd"),
    BuildTarget("kmertree", "kmertree", "kmertree"),
)


def info(message: str) -> None:
    print(f"[PATs install] {message}")


def fail(message: str) -> None:
    raise RuntimeError(message)


def command_path(name: str, conda_prefix: Path | None = None) -> str | None:
    if conda_prefix:
        candidate = conda_prefix / "bin" / name
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return str(candidate)
    return shutil.which(name)


def active_conda_prefix() -> Path | None:
    value = os.environ.get("CONDA_PREFIX")
    return Path(value).resolve() if value else None


def conda_manager() -> str | None:
    return shutil.which("mamba") or shutil.which("conda")


def validate_layout() -> list[str]:
    missing: list[str] = []
    for path in (
        RUNTIME / "pats.py",
        RUNTIME / "Snakefile",
        SCRIPT_DIR,
        SCRIPT_DIR / "matrixformat.py",
        REQUIREMENTS,
    ):
        if not path.exists():
            missing.append(str(path))
    for target in BUILD_TARGETS:
        folder = SOURCE_DIR / target.source
        if not (folder / "Makefile").is_file() or not list((folder / "src").glob("*.cpp")):
            missing.append(str(folder))
    return missing


def install_conda_dependencies(prefix: Path) -> None:
    manager = conda_manager()
    if manager is None:
        fail("conda or mamba is required to install compiled libraries and bioinformatics tools")
    packages = list(CONDA_PACKAGES)
    if platform.system() == "Darwin":
        packages.append("llvm-openmp")
    command = [
        manager,
        "install",
        "--yes",
        "--prefix",
        str(prefix),
        "--channel",
        "conda-forge",
        "--channel",
        "bioconda",
        *packages,
    ]
    info("installing external tools and build libraries with conda/mamba")
    subprocess.run(command, check=True)


def environment_python(conda_prefix: Path | None) -> str:
    if conda_prefix:
        candidate = conda_prefix / "bin" / "python"
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return str(candidate)
    return sys.executable


def install_python_dependencies(python: str) -> None:
    info(f"installing Python packages from {REQUIREMENTS}")
    subprocess.run(
        [python, "-m", "pip", "install", "-r", str(REQUIREMENTS)],
        check=True,
    )


def compiler_candidates(prefix: Path | None) -> list[str]:
    candidates: list[str] = []
    if os.environ.get("CXX"):
        candidates.append(os.environ["CXX"])
    if prefix:
        patterns = ("bin/*-c++", "bin/*-g++", "bin/clang++", "bin/g++")
        for pattern in patterns:
            candidates.extend(str(path) for path in sorted(prefix.glob(pattern)))
    candidates.extend(name for name in ("g++", "c++", "clang++") if shutil.which(name))
    return list(dict.fromkeys(candidates))


def find_compiler(prefix: Path | None) -> str:
    for candidate in compiler_candidates(prefix):
        path = shutil.which(candidate) if os.sep not in candidate else candidate
        if path and Path(path).is_file() and os.access(path, os.X_OK):
            return str(Path(path).resolve())
    fail("no C++ compiler was found; install the conda cxx-compiler package")


def prefix_build_flags(prefix: Path | None) -> tuple[str, str]:
    if prefix is None:
        return "", ""
    include = shlex.quote(str(prefix / "include"))
    library = shlex.quote(str(prefix / "lib"))
    return f"-I{include}", f"-L{library} -Wl,-rpath,{library}"


def openmp_flags(compiler: str, cppflags: str, ldflags: str) -> tuple[str, str]:
    attempts = [("-fopenmp", "-fopenmp")]
    if platform.system() == "Darwin":
        attempts.append(("-Xpreprocessor -fopenmp", "-lomp"))
    source_text = "#include <omp.h>\nint main(){return omp_get_max_threads() < 1;}\n"
    with tempfile.TemporaryDirectory(prefix="pats.openmp.") as temporary:
        source = Path(temporary) / "test.cpp"
        binary = Path(temporary) / "test"
        source.write_text(source_text)
        for compile_flags, link_flags in attempts:
            command = [
                compiler,
                *shlex.split(cppflags),
                *shlex.split(compile_flags),
                str(source),
                "-o",
                str(binary),
                *shlex.split(ldflags),
                *shlex.split(link_flags),
            ]
            result = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            if result.returncode == 0:
                return compile_flags, link_flags
    fail(
        "the available C++ compiler cannot build OpenMP code required by kmerpartition; "
        "install cxx-compiler (and llvm-openmp on macOS) in the active conda environment"
    )


def build_binaries(prefix: Path, jobs: int) -> None:
    compiler = find_compiler(active_conda_prefix())
    cppflags, ldflags = prefix_build_flags(active_conda_prefix())
    omp_compile, omp_link = openmp_flags(compiler, cppflags, ldflags)
    prefix.mkdir(parents=True, exist_ok=True)
    info(f"using C++ compiler: {compiler}")

    for target in BUILD_TARGETS:
        source = SOURCE_DIR / target.source
        make_vars = [
            f"CXX={compiler}",
            f"CPPFLAGS={cppflags}",
            f"LDFLAGS={ldflags}",
        ]
        if target.openmp:
            make_vars.extend(
                [
                    f"OPENMP_CXXFLAGS={omp_compile}",
                    f"OPENMP_LDFLAGS={omp_link}",
                ]
            )
        info(f"building {target.installed_name} from {source.relative_to(ROOT)}")
        subprocess.run(["make", "clean", *make_vars], cwd=source, check=True)
        subprocess.run(
            ["make", f"-j{jobs}", target.make_target, *make_vars],
            cwd=source,
            check=True,
        )
        built = source / target.built_name
        if not built.is_file():
            fail(f"build completed but did not create {built}")
        destination = prefix / target.installed_name
        shutil.copy2(built, destination)
        destination.chmod(destination.stat().st_mode | 0o111)
        subprocess.run(["make", "clean", *make_vars], cwd=source, check=True)
        info(f"installed {destination}")

    for wrapper in ("runblastn", "runmakeblastdb", "runwinnowmap.sh"):
        path = SCRIPT_DIR / wrapper
        if path.is_file():
            path.chmod(path.stat().st_mode | 0o111)


def check_python(python: str) -> list[str]:
    missing: list[str] = []
    for module in PYTHON_MODULES:
        result = subprocess.run(
            [python, "-c", f"import {module}"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        if result.returncode != 0:
            missing.append(module)
    compatibility = subprocess.run(
        [
            python,
            "-c",
            (
                "import pulp, snakemake; "
                f"assert snakemake.__version__ == {SNAKEMAKE_VERSION!r}; "
                f"assert pulp.__version__ == {PULP_VERSION!r}; "
                "assert hasattr(pulp, 'list_solvers')"
            ),
        ],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    if compatibility.returncode != 0:
        missing.append(
            f"snakemake=={SNAKEMAKE_VERSION}/PuLP=={PULP_VERSION} compatibility"
        )
    return missing


def check_external(prefix: Path | None) -> list[str]:
    return [name for name in EXTERNAL_COMMANDS if command_path(name, prefix) is None]


def check_binaries(prefix: Path) -> list[str]:
    return [target.installed_name for target in BUILD_TARGETS if not (prefix / target.installed_name).is_file()]


def report_status(prefix: Path, *, python: bool, external: bool, binaries: bool) -> bool:
    okay = True
    layout_missing = validate_layout()
    if layout_missing:
        okay = False
        info("missing source/layout paths: " + ", ".join(layout_missing))
    else:
        info("all required Python, shell, and C++ source files are present")

    if python:
        python_executable = environment_python(active_conda_prefix())
        missing = check_python(python_executable)
        okay &= not missing
        info(
            f"Python packages ({python_executable}): "
            + ("OK" if not missing else "missing " + ", ".join(missing))
        )
    if external:
        missing = check_external(active_conda_prefix())
        okay &= not missing
        info("external tools: " + ("OK" if not missing else "missing " + ", ".join(missing)))
    if binaries:
        missing = check_binaries(prefix)
        okay &= not missing
        info("PATs binaries: " + ("OK" if not missing else "missing " + ", ".join(missing)))
    return okay


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Install Python packages, external bioinformatics tools, and PATs C++ binaries."
    )
    parser.add_argument(
        "--prefix",
        type=Path,
        default=SCRIPT_DIR,
        help="binary installation directory (default: scripts)",
    )
    parser.add_argument("--jobs", type=int, default=max(1, os.cpu_count() or 1))
    parser.add_argument("--check", action="store_true", help="report status without changing anything")
    parser.add_argument("--skip-python", action="store_true")
    parser.add_argument("--skip-external", action="store_true")
    parser.add_argument("--skip-build", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    prefix = args.prefix.expanduser().resolve()
    if args.jobs < 1:
        fail("--jobs must be greater than zero")
    if not (sys.version_info[:2] == (3, 10)):
        fail(
            "Python 3.10 is required by the pinned Snakemake 6.5 environment; "
            f"found {platform.python_version()}"
        )

    missing_layout = validate_layout()
    if missing_layout:
        fail("required source files are missing: " + ", ".join(missing_layout))

    if args.check:
        return 0 if report_status(
            prefix,
            python=not args.skip_python,
            external=not args.skip_external,
            binaries=not args.skip_build,
        ) else 1

    conda_prefix = active_conda_prefix()
    if not args.skip_external:
        if conda_prefix is None:
            fail("activate the conda environment in which PATs should be installed, then rerun install.py")
        install_conda_dependencies(conda_prefix)
    if not args.skip_python:
        install_python_dependencies(environment_python(conda_prefix))
    if not args.skip_build:
        build_binaries(prefix, args.jobs)

    okay = report_status(
        prefix,
        python=not args.skip_python,
        external=not args.skip_external,
        binaries=not args.skip_build,
    )
    if okay:
        python = environment_python(conda_prefix)
        info(f"installation complete; run: {python} {RUNTIME / 'pats.py'} --help")
        return 0
    return 1


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, RuntimeError, subprocess.CalledProcessError) as exc:
        print(f"[PATs install] ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
