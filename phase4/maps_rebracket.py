"""Phase 4 map follow-up (separate data-producing module; records carry its own sha256).
(1) Map roots that failed only because the criterion fails at U = 0 (adiabatic event-location failure) are re-run
    with the documented bracket U in [100, 750]. (2) Failed roots whose probes show g jumping from ~0 to ~-1 at the
    root are CLASSIFIED (not overridden): the sign change of g is located to xtol, but it sits at the documented
    discontinuity of g (interior S_T maximum vanishes at S ~ 1), so it is reported as a flagged sign-change location."""
import hashlib
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import numpy as np  # noqa: E402

from phase4 import common  # noqa: E402
from phase4.common import read_json, root_summary, write_json  # noqa: E402
from pbr import params, sensitivity  # noqa: E402

P = params.load()
BR = (100.0, 750.0)
SELF = hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()


def classify_discontinuity(rs):
    lg = rs.get("local_g")
    if (rs["status"] != "PASS" and rs.get("root") is not None and lg is not None and rs.get("sign_changes") == 1
            and rs.get("subscan_sign_changes") == 1):
        lg = np.asarray(lg, float)
        if np.max(np.abs(np.diff(lg))) > 0.5 and np.min(np.abs(lg)) < 1e-3:
            rs["g_discontinuity_at_root"] = True
            rs["sign_change_location"] = rs["root"]
            rs["outcome"] = ("sign change of g located to xtol at a DISCONTINUITY of g (interior S_T maximum vanishes "
                             "at S ~ 1); flagged, not a smooth validated root")
    return rs


out = []
for m in read_json("boundary_maps"):
    if m["status"] == "PASS":
        continue
    if m["detail"].startswith("criterion evaluation failed") and m["kind"] == "U":
        q = P.with_(T0=m["T0"], Tc=m["T0"]) if m["config"] == "normalized" else P.with_(T0=m["T0"])
        rs = root_summary(sensitivity.find_Ucrit(q, config=m["config"], U_bracket=BR), "map",
                          m["label"] + " [U bracket 100-750]", "U", series=m["series"], T0=m["T0"], config=m["config"],
                          mode="none", Pi=0.0, U_bracket=BR, rebracketed=True)
        out.append(classify_discontinuity(rs))
        print(f"rebracketed {m['label']}: {rs['status']} {rs.get('U_crit')} | {rs['outcome'][:70]}", flush=True)
    elif m.get("root") is not None:
        rs = classify_discontinuity(dict(m, label=m["label"] + " [classified]", classified_from_original=True))
        out.append(rs)
        print(f"classified {m['label']}: {rs['outcome'][:70]} at {rs.get('sign_change_location')}", flush=True)
for o in out:
    o["meta"] = {**common.short_meta(), "maps_rebracket_sha256": SELF}
write_json("boundary_maps_rebracket", out)
