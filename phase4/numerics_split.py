"""Phase 4 numerical robustness near the sensitivity-based boundary, split into parts that each fit one command.
Usage: python phase4/numerics_split.py [phys | tight | outlet]. Writes results/phase4/numerical_robustness_<part>.json
(with this file's own sha256; the study-module hash of the other records is unaffected)."""
import hashlib
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from phase4 import common  # noqa: E402
from phase4.common import read_json, write_json  # noqa: E402
from pbr import params, sensitivity  # noqa: E402

P = params.load()
PD_BRACKET = (100.0, 750.0)   # same documented pressure-drop bracket as studies.py
SELF = hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()


def main(part):
    out = []
    if part == "phys":
        for T0 in (610.0, 625.0, 640.0):
            q = P.with_(T0=T0)
            b, r = sensitivity.find_Ucrit(q), sensitivity.find_Ucrit(q, solver="Radau")
            out.append({"case": f"physical Pi=0 T0={T0:g}", "BDF": b.get("U_crit"), "Radau": r.get("U_crit"),
                        "statuses": [b["status"].value, r["status"].value],
                        "abs_diff_Radau": abs(r.get("U_crit") - b.get("U_crit"))})
            print(out[-1], flush=True)
    elif part == "tight":
        for T0 in (610.0, 625.0, 640.0):
            q = P.with_(T0=T0)
            b = sensitivity.find_Ucrit(q)
            t = sensitivity.find_Ucrit(q, U_bracket=(b["U_crit"] - 10.0, b["U_crit"] + 13.0), n_scan=3,
                                       rtol=1e-11, atol_scaled=1e-13)
            out.append({"case": f"physical Pi=0 T0={T0:g}", "BDF": b.get("U_crit"), "BDF_rtol_1e-11": t.get("U_crit"),
                        "statuses": [b["status"].value, t["status"].value],
                        "abs_diff_tight": abs(t.get("U_crit") - b.get("U_crit"))})
            print(out[-1], flush=True)
    elif part == "outlet":
        ref = [m for m in read_json("boundary_maps_outlet") if m["T0"] == 625.0][0]
        r = sensitivity.find_Ucrit(P, mode="outlet", Pi=0.2, n_scan=3, U_bracket=(ref["U_crit"] - 25.0, ref["U_crit"] + 35.0),
                                   solver="Radau")
        out.append({"case": "physical outlet Pi=0.2 T0=625", "BDF": ref["U_crit"], "Radau": r.get("U_crit"),
                    "statuses": [ref["status"], r["status"].value], "abs_diff_Radau": abs(r.get("U_crit") - ref["U_crit"]),
                    "note": "BDF value from boundary_maps_outlet (bracket [100, 750]); Radau root bracketed near it"})
        print(out[-1], flush=True)
    for o in out:
        o["meta"] = {**common.short_meta(), "numerics_split_sha256": SELF}
    write_json(f"numerical_robustness_{part}", out)


if __name__ == "__main__":
    main(sys.argv[1])
