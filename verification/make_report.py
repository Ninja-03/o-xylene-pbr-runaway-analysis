"""Assemble results/verification/*.json into results/verification_report.json and .md."""
import hashlib
import json
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]
ORDER = ["A1"] + [f"V{i}" for i in range(1, 14)]


def cell(x):
    return str(x).replace("|", "/").replace("\n", " ")


def main():
    recs = []
    for tid in ORDER:
        f = ROOT / "results" / "verification" / f"{tid}.json"
        if f.exists():
            doc = json.loads(f.read_text())
            d = doc["data"]
            d["meta"] = doc["meta"]
        else:
            d = {"id": tid, "purpose": "", "expected": "", "actual": "NOT RUN", "error": "", "tolerance": "",
                 "pass": False, "comments": "test not run"}
        recs.append(d)
    disp = lambda r: r.get("disposition", "PASS" if r["pass"] else "FAIL")
    summary = {"tests": len(recs), "pass": [r["id"] for r in recs if disp(r) == "PASS"],
               "pass_with_documented_limitation": [r["id"] for r in recs if disp(r) == "PASS WITH DOCUMENTED LIMITATION"],
               "failed": [r["id"] for r in recs if disp(r) == "FAIL"],
               "tolerances_sha256": hashlib.sha256((ROOT / "verification" / "tolerances.yaml").read_bytes()).hexdigest(),
               "params_sha256": hashlib.sha256((ROOT / "params" / "base_case.yaml").read_bytes()).hexdigest()}
    (ROOT / "results" / "verification_report.json").write_text(json.dumps({"summary": summary, "tests": recs}, indent=1))
    lines = ["# Phase 3 verification report", "",
             f"Parameter file sha256 `{summary['params_sha256']}`; tolerances sha256 `{summary['tolerances_sha256']}`.", "",
             f"**PASS: {len(summary['pass'])}; PASS WITH DOCUMENTED LIMITATION: "
             f"{summary['pass_with_documented_limitation'] or 'none'}; FAIL: {summary['failed'] or 'none'}.**", "",
             "Primary runaway criterion: *sensitivity-based parametric runaway boundary* (not a dynamic-stability "
             "boundary, saddle-node or bifurcation point). Morbidelli-Varma generalized criterion: not evaluated "
             "(primary source not accessed).", "",
             "| TEST | PURPOSE | EXPECTED RESULT | ACTUAL RESULT | ERROR | PASS/FAIL | COMMENTS |",
             "|---|---|---|---|---|---|---|"]
    for r in recs:
        lines.append("| " + " | ".join(cell(x) for x in (r["id"], r["purpose"], r["expected"], r["actual"], r["error"],
                                                         disp(r), r["comments"])) + " |")
    (ROOT / "results" / "verification_report.md").write_text("\n".join(lines) + "\n")
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
