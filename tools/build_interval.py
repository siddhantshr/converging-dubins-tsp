"""Compile the C++ interval solver (src/dubins/_interval.cpp) in place.

    python tools/build_interval.py

Without the compiled module, src/dubins/interval.py falls back to the
(much slower) pure-Python solver and prints a warning.
"""

from __future__ import annotations

import shlex
import subprocess
import sys
import sysconfig
from pathlib import Path

SRC = Path(__file__).resolve().parent.parent / "src" / "dubins" / "_interval.cpp"


def main() -> int:
    out = SRC.with_name("_interval" + sysconfig.get_config_var("EXT_SUFFIX"))
    # CXX may carry flags, e.g. "clang++ -pthread". No -ffast-math: the solver needs IEEE NaNs.
    cmd = shlex.split(sysconfig.get_config_var("CXX") or "c++")
    cmd += ["-O3", "-std=c++17", "-shared", "-fPIC", f"-I{sysconfig.get_paths()['include']}"]
    if sys.platform == "darwin":
        cmd += ["-undefined", "dynamic_lookup"]  # symbols come from the running Python
    cmd += [str(SRC), "-o", str(out)]
    print(shlex.join(cmd))
    subprocess.run(cmd, check=True)
    print(f"built {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
