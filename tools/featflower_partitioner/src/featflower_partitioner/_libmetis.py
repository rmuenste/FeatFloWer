"""
Prioritised libmetis.so loader.

Search order:
  1. <package_dir>/libmetis.so  (bundled at install time via setup.py)
  2. $FEATFLOWER_BUILD_DIR/extern/libraries/metis-5.1.0/libmetis/libmetis.so
  3. ./libmetis.so  (working directory — legacy behaviour)
  4. ../lib64/libmetis.so  (relative to cwd — legacy behaviour)
  5. libmetis.so  (system path via ldconfig)

Raises RuntimeError listing all tried paths on failure.
Exposes: metis (CDLL handle), metis_func (tuple of METIS partition functions).
"""

import os
from pathlib import Path
from ctypes import CDLL, c_int, c_float, POINTER

_lib_name_posix = "libmetis.so"


def _load_metis():
    candidates = []

    # 1. Bundled copy installed alongside the package (highest priority)
    bundled = Path(__file__).parent / _lib_name_posix
    candidates.append(str(bundled))

    # 2. Build tree pointed to by environment variable
    build_dir = os.environ.get("FEATFLOWER_BUILD_DIR")
    if build_dir:
        candidates.append(
            os.path.join(build_dir, "extern", "libraries", "metis-5.1.0",
                         "libmetis", _lib_name_posix)
        )

    # 3/4. Legacy cwd-relative paths
    candidates.append(os.path.join(os.curdir, _lib_name_posix))
    candidates.append(os.path.join(os.curdir, "..", "lib64", _lib_name_posix))

    # 5. System library path
    candidates.append(_lib_name_posix)

    errors = []
    for path in candidates:
        # For the system-path entry (plain name) skip the exists check so
        # ldconfig gets a chance even when the file isn't visible directly.
        if path != _lib_name_posix and not os.path.exists(path):
            continue
        try:
            return CDLL(path)
        except OSError as exc:
            errors.append((path, str(exc)))

    tried = "\n  ".join(p for p, _ in errors) or "\n  ".join(candidates)
    raise RuntimeError(
        "Could not load the METIS library.  Searched:\n  " + tried + "\n\n"
        "Hint: set FEATFLOWER_BUILD_DIR to your cmake build directory."
    )


_metis = None
_metis_func = None


def _ensure_loaded():
    global _metis, _metis_func
    if _metis is not None:
        return
    _metis = _load_metis()
    _pidx  = POINTER(c_int)
    _preal = POINTER(c_float)
    _PartArgs = (
        _pidx,   # nvtxs
        _pidx,   # ncon
        _pidx,   # xadj
        _pidx,   # adjncy
        _pidx,   # vwgt
        _pidx,   # vsize
        _pidx,   # adjwgt
        _pidx,   # nparts
        _preal,  # tpwgts
        _preal,  # ubvec
        _pidx,   # options
        _pidx,   # edgecut
        _pidx,   # part
    )
    _metis.METIS_PartGraphRecursive.argtypes = _PartArgs
    _metis.METIS_PartGraphKway.argtypes      = _PartArgs
    _metis_func = (
        _metis.METIS_PartGraphRecursive,  # Method 1
        _metis.METIS_PartGraphKway,       # Method 2
        _metis.METIS_PartGraphKway,       # Method 3
    )


def get_metis():
    _ensure_loaded()
    return _metis


def get_metis_func():
    _ensure_loaded()
    return _metis_func
