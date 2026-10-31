"""(16) Plotting utilities. Figures are diagnostics only (no scientific results in Phase 3)."""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def profiles(runs, p, path, title, note=""):
    """runs: list of (Solution, label). Temperature, conversion and pressure profiles."""
    fig, ax = plt.subplots(1, 3, figsize=(13.0, 3.9))
    for s, lab in runs:
        ax[0].plot(s.z, s.T, label=lab)
        ax[1].plot(s.z, s.X, label=lab)
        ax[2].plot(s.z, s.P / p.P_ref, label=lab)
    ax[0].axhline(p.Tc, color="0.6", ls="--", lw=0.8, label="T_c")
    ax[0].set(xlabel="z (m)", ylabel="T (K)", title="temperature")
    ax[1].set(xlabel="z (m)", ylabel="X (-)", title="conversion")
    ax[2].set(xlabel="z (m)", ylabel="P / P_ref (-)", title="pressure")
    for a in ax:
        a.grid(alpha=0.35)
        a.legend(fontsize=7)
    fig.suptitle(title, fontsize=10)
    if note:
        fig.text(0.5, 0.01, note, ha="center", fontsize=7.5)
    fig.tight_layout(rect=(0.0, 0.06, 1.0, 0.93))
    fig.savefig(path, dpi=130)
    plt.close(fig)


def scan(x, y, path, xlabel, ylabel, title):
    fig, ax = plt.subplots(figsize=(5.2, 3.8))
    ax.plot(x, y, marker="o", ms=3)
    ax.axhline(0.0, color="0.6", lw=0.8)
    ax.set(xlabel=xlabel, ylabel=ylabel, title=title)
    ax.grid(alpha=0.35)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)
