"""Phase 4 completion (data-producing; own sha256 in its output): O2-depletion U_crit at T0 = 610 K with the
documented bracket U in [100, 750] (the [0, 750] attempt failed at U = 0, adiabatic event-location limitation).
The base fixed-O2 model is not altered; the variant is applied in memory (o2_variant.o2_depletion)."""
import hashlib
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from phase4 import common, o2_variant  # noqa: E402
from phase4.common import read_json, root_summary, write_json  # noqa: E402
from pbr import params, sensitivity  # noqa: E402

P = params.load()
BR = (100.0, 750.0)
SELF = hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()
with o2_variant.o2_depletion():
    r = sensitivity.find_Ucrit(P.with_(T0=610.0), U_bracket=BR)
rs = root_summary(r, "o2", "O2 Ucrit(T0=610) [U bracket 100-750]", "U", series="O2-depletion variant", T0=610.0,
                  U_bracket=BR, rebracketed=True)
base = [x for x in read_json("o2_roots")["roots"] if x["series"] == "fixed-O2 (base)" and x["T0"] == 610.0][0]
rs["fixed_O2_U_crit"] = base["U_crit"]
rs["diff_vs_fixed_O2"] = rs["U_crit"] - base["U_crit"] if rs["status"] == "PASS" else None
rs["verification"] = {k: rs.get(k) for k in ("converged", "sign_changes", "subscan_sign_changes",
                                             "local_monotonic_at_root", "g_at_root", "bracket", "iterations")}
rs["meta"] = {**common.short_meta(), "o2_rebracket_sha256": SELF}
write_json("o2_roots_rebracket", [rs])
print(f"O2 U_crit(610): {rs['status']} {rs.get('U_crit')} | fixed-O2 {base['U_crit']} | diff {rs['diff_vs_fixed_O2']} | {rs['verification']}")
