#!/usr/bin/env python3
"""Run one configuration of the validation matrix and compare its output
with the reference results, or store the output as the new reference.

The configurations live in configurations.txt (labels | prefix | args).
A run of `disorder <args> -prefix <prefix>` produces
  <prefix>xsct_*.dat, <prefix>disorder_*.dat   (results, histograms)
and we store its screen output as <prefix without trailing _>.log.
All of these are compared with the files of the same name in the
reference directory:
  - lines with volatile content (timings, dates, library banners, the
    version and the years in the welcome banner) are
    dropped, and the path of the executable in the echoed command line
    is ignored;
  - any NaN or infinity in the output is a failure (also with --generate);
  - the remaining lines must have the same text, and numbers must agree
    within a relative tolerance rtol (default 1e-5; plus an absolute
    tolerance atol for numbers that are zero up to rounding), so that the
    results do not have to be bitwise identical across compilers and
    platforms. (Changing the floating-point code generation, e.g. to
    -O2 -march=native without LTO, moves the cross sections by <~1e-8 and
    the P2B histograms by <~1e-6 relative.)

Usage:
  run_validation.py --disorder BIN --prefix PREFIX --workdir DIR
                    [--refdir DIR] [--generate] [--rtol R] [--atol A]
  run_validation.py --list [--label L]      # print the prefixes
"""
import argparse
import math
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CONFIG = os.path.join(HERE, "configurations.txt")

# Lines containing any of these (case-insensitive) are not compared
VOLATILE = ("total time", "stamped by", "fastjet", "hoppet", "arxiv", "lhapdf",
            "welcome to disorder",          # the version number
            "written by alexander karlberg")  # the years
# The echoed command line starts with the path to the executable
COMMAND_LINE = re.compile(r"#\S*disorder(?=\s)")
# VEGAS's chi^2 per iteration is a diagnostic that is zero up to rounding
# after the first iteration
VEGAS_CHI2 = re.compile(r"(chi\*\*2/it'n =)\s*\S+")
# A NaN or infinity anywhere in the output fails the run (also when
# generating references), since a reference with NaN would compare equal
# to a broken run
NONFINITE = re.compile(rb"(?i)\b(nan|infinity)\b")
# Numbers anywhere in a line (Fortran may glue them to other characters)
NUMBER = re.compile(r"([-+]?(?:\d+\.?\d*|\.\d+)(?:[eEdD][-+]?\d+)?)")


def read_configurations():
    configs = []
    with open(CONFIG, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            labels, prefix, args = (field.strip() for field in line.split("|"))
            configs.append((labels.split(), prefix, args.split()))
    return configs


def output_files(directory, prefix):
    """The files belonging to a run with this prefix."""
    names = []
    for name in sorted(os.listdir(directory)):
        rest = name[len(prefix):] if name.startswith(prefix) else None
        if rest is not None and (rest.startswith("xsct_") or rest.startswith("disorder_")):
            names.append(name)
    log = prefix.rstrip("_") + ".log"
    if os.path.exists(os.path.join(directory, log)):
        names.append(log)
    return names


def drop_fastjet_banner(lines):
    """FastJet prints its banner through C++ std::cout, so where it ends up
    relative to the Fortran output depends on buffering: drop it."""
    out, i = [], 0
    while i < len(lines):
        if lines[i].startswith("#---") and i + 1 < len(lines) and "FastJet release" in lines[i + 1]:
            j = i + 1
            while j < len(lines) and not lines[j].startswith("#---"):
                j += 1
            i = j + 1
            continue
        out.append(lines[i])
        i += 1
    return out


def comparable_lines(path):
    with open(path, encoding="utf-8", errors="replace") as f:
        lines = [COMMAND_LINE.sub("#disorder", l.rstrip("\n")) for l in f]
    lines = drop_fastjet_banner(lines)
    lines = [VEGAS_CHI2.sub(r"\1 *", l) for l in lines]
    return [l for l in lines if not any(v in l.lower() for v in VOLATILE)]


def split_numbers(line):
    """Split a line into its text (compared exactly, ignoring spacing) and
    its numbers (compared with a tolerance)."""
    parts = NUMBER.split(line)
    text = " ".join("".join(parts[0::2]).split())
    numbers = [float(p.replace("D", "E").replace("d", "e")) for p in parts[1::2]]
    return text, numbers


def compare_files(ref, new, rtol, atol):
    """Return (list of differences, largest relative deviation)."""
    ref_lines, new_lines = comparable_lines(ref), comparable_lines(new)
    problems, worst = [], 0.0
    if len(ref_lines) != len(new_lines):
        problems.append(f"different number of lines: {len(ref_lines)} (ref) vs {len(new_lines)}")
    for i, (a, b) in enumerate(zip(ref_lines, new_lines)):
        if a == b:
            continue
        (ta, na), (tb, nb) = split_numbers(a), split_numbers(b)
        bad = ta != tb or len(na) != len(nb)
        for x, y in zip(na, nb):
            if x == y:
                continue
            if not (math.isfinite(x) and math.isfinite(y)):
                bad = True
                continue
            diff, scale = abs(x - y), max(abs(x), abs(y))
            if diff > atol:
                worst = max(worst, diff / scale)
            if diff > max(rtol * scale, atol):
                bad = True
        if bad:
            problems.append(f"line {i+1}:\n    ref: {a}\n    new: {b}")
    return problems, worst


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--disorder", help="path to the disorder executable")
    ap.add_argument("--prefix", help="configuration to run (its prefix)")
    ap.add_argument("--workdir", help="directory to run in (wiped first)")
    ap.add_argument("--refdir", default=os.path.join(HERE, "ref_runs"))
    ap.add_argument("--generate", action="store_true", help="store output as the new reference")
    ap.add_argument("--rtol", type=float, default=1e-5)
    ap.add_argument("--atol", type=float, default=1e-12)
    ap.add_argument("--list", action="store_true", help="list the prefixes and exit")
    ap.add_argument("--label", help="with --list: only configurations with this label")
    opts = ap.parse_args()

    configs = read_configurations()
    if opts.list:
        for labels, prefix, _ in configs:
            if opts.label is None or opts.label in labels:
                print(prefix)
        return 0

    matches = [c for c in configs if c[1] == opts.prefix]
    if len(matches) != 1 or not opts.disorder or not opts.workdir:
        ap.error("need --disorder, --workdir and a --prefix listed in configurations.txt")
    _, prefix, args = matches[0]

    shutil.rmtree(opts.workdir, ignore_errors=True)
    os.makedirs(opts.workdir)
    # Run through a link in the working directory, so that the program
    # name echoed in the output is short and does not depend on the build
    # location.
    os.symlink(os.path.abspath(opts.disorder), os.path.join(opts.workdir, "disorder"))
    command = ["./disorder"] + args + ["-prefix", prefix]
    print("Running in", opts.workdir + ":", " ".join(command), flush=True)
    # The log holds all of stdout followed by all of stderr (as GNU
    # parallel, used to make the original references, wrote it), which
    # does not depend on how the two streams are buffered.
    run = subprocess.run(command, cwd=opts.workdir, capture_output=True)
    log = prefix.rstrip("_") + ".log"
    with open(os.path.join(opts.workdir, log), "wb") as out:
        out.write(run.stdout + run.stderr)
    if run.returncode != 0:
        tail = (run.stdout + run.stderr).decode(errors="replace").splitlines()[-40:]
        print("\n".join(tail))
        print(f"FAILED: disorder exited with status {run.returncode}")
        return 1

    produced = output_files(opts.workdir, prefix)
    nonfinite = False
    for name in produced:
        with open(os.path.join(opts.workdir, name), "rb") as f:
            for number, line in enumerate(f, 1):
                if NONFINITE.search(line):
                    print(f"FAILED: {name}:{number} contains a NaN or infinity:")
                    print("  " + line.decode(errors="replace").rstrip()[:200])
                    nonfinite = True
                    break
    if nonfinite:
        return 1
    if opts.generate:
        os.makedirs(opts.refdir, exist_ok=True)
        for name in output_files(opts.refdir, prefix):
            os.remove(os.path.join(opts.refdir, name))
        for name in produced:
            shutil.copy(os.path.join(opts.workdir, name), opts.refdir)
        print(f"Stored {len(produced)} reference files in {opts.refdir}")
        return 0

    expected = output_files(opts.refdir, prefix)
    failed = False
    for name in sorted(set(expected) - set(produced)):
        print(f"FAILED: {name} was not produced")
        failed = True
    for name in sorted(set(produced) - set(expected)):
        print(f"FAILED: {name} has no reference")
        failed = True
    for name in expected:
        if name not in produced:
            continue
        problems, worst = compare_files(os.path.join(opts.refdir, name),
                                        os.path.join(opts.workdir, name), opts.rtol, opts.atol)
        status = "FAILED" if problems else "passed"
        print(f"{status}: {name} (largest relative deviation {worst:.1e})")
        for p in problems[:20]:
            print("  " + p)
        if len(problems) > 20:
            print(f"  ... and {len(problems) - 20} more differences")
        failed = failed or bool(problems)
    if not expected:
        print(f"FAILED: no reference files for {prefix} in {opts.refdir}")
        failed = True
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
