"""Run selected verification tests (default: all) and store one JSON per test in results/verification/."""
import pathlib
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parents[1]
for _p in (str(ROOT), str(ROOT / "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from pbr import params, storage  # noqa: E402
from verification import tests  # noqa: E402


def main(ids):
    p = params.load()
    for tid in ids:
        t = time.perf_counter()
        r = getattr(tests, tid)(p)
        r["runtime_s"] = round(time.perf_counter() - t, 2)
        storage.save(ROOT / "results" / "verification" / f"{tid}.json", r, p)
        print(f"{tid:4s} {r['disposition']} ({r['runtime_s']:7.1f} s)  {r['actual']}", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:] or tests.ALL)
