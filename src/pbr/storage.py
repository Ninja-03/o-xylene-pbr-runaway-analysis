"""(15) Result storage: JSON with provenance (parameter-file path and sha256, in-memory overrides,
source-code hash, package versions, solver settings, run statuses)."""
import datetime
import enum
import hashlib
import json
import math
import pathlib
import platform

import numpy as np

SRC = pathlib.Path(__file__).resolve().parent


def code_hash():
    h = hashlib.sha256()
    for f in sorted(SRC.glob("*.py")):
        h.update(f.name.encode())
        h.update(f.read_bytes())
    return h.hexdigest()


def environment():
    import matplotlib
    import scipy
    import yaml
    return {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__,
            "pyyaml": yaml.__version__, "matplotlib": matplotlib.__version__, "platform": platform.platform()}


def jsonable(o):
    if isinstance(o, enum.Enum):
        return o.value
    if isinstance(o, dict):
        return {str(k): jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return jsonable(o.tolist())
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, (float, np.floating)):
        v = float(o)
        return v if math.isfinite(v) else str(v)
    return o


def save(path, payload, p=None):
    meta = {"created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
            "code_sha256": code_hash(), "environment": environment()}
    if p is not None:
        meta.update(params_path=p.path, params_sha256=p.sha256)
    path = pathlib.Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"meta": meta, "data": jsonable(payload)}, indent=1, allow_nan=False))
    return path


def load(path):
    return json.loads(pathlib.Path(path).read_text())
