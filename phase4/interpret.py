"""Read-only interpretation rules for Phase 4 roots (used by figures/tables/report; produces no data).
Decision (Phase 4 finalization): the PRIMARY interpreted sensitivity-based parametric runaway boundary map
(physical configuration, Pi = 0) is restricted to the validated smooth-root domain T0 <= 650 K, U <= 450 W/(m2 K).
Flagged sign changes at the discontinuity of g are kept as diagnostic data only."""
T0_VALID, U_VALID = 650.0, 450.0
FLAG_TEXT = ("Flagged sign change associated with discontinuity of the sensitivity criterion; "
             "not interpreted as a validated boundary root")


def interpretation(m):
    out = m.get("outcome", "")
    if m.get("g_discontinuity_at_root"):
        return FLAG_TEXT
    if m["status"] == "PASS":
        s = "validated smooth root"
        if "Pi=0" in m.get("series", "") and "physical" in m.get("series", ""):
            T0 = m.get("T0") if m.get("kind") == "U" else m.get("T0_crit")
            U = m.get("U_crit") if m.get("kind") == "U" else m.get("U")
            s += " (primary map domain)" if (T0 <= T0_VALID and U <= U_VALID) else " (outside primary map domain)"
        return s
    if "above the searched range" in out:
        return "U_crit > 750 (above searched range)" if m.get("kind") == "U" else "T0_crit above searched range"
    if "below the searched range" in out:
        return "T0_crit < 600 (below searched range)" if m.get("kind") == "T0" else "U_crit below searched range"
    if m.get("detail", "").startswith("criterion evaluation failed"):
        return "superseded: U = 0 event-location failure; see rebracketed root [100, 750]"
    return "root failure: " + m.get("detail", "")
