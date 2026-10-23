"""Offline plots of every family-best curve, retaining inadequate alternatives."""

import matplotlib.pyplot as plt


def plot_family_comparison(x, y, fit, data_kind, output):
    families = fit.get("comparison_fits", {})
    if len(families) < 2:
        return
    fig, axes = plt.subplots(
        2, len(families), figsize=(4.2 * len(families), 6.5), layout="constrained", squeeze=False
    )
    for index, (shape, result) in enumerate(families.items()):
        signal, residual = axes[:, index]
        signal.plot(x, y, color="#a0aab1", lw=0.7, label="Original corrected samples")
        signal.plot(x, result["fitted"], color="#137b80", lw=1.3, label="Family-best model")
        for curve in result["components"]:
            signal.plot(x, curve, "--", color="#b77530", lw=0.8)
        signal.set_title(
            f"{shape.capitalize()} | N={len(result['peaks'])}\n{result['status']} (conditional fit)",
            fontsize=11,
        )
        signal.set_ylabel("Corrected intensity (a.u.)")
        signal.legend(frameon=False, fontsize=7)
        residual.plot(x, result["residuals"], color="#435669", lw=0.7)
        residual.axhline(0, color="#444444", lw=0.5)
        residual.set(xlabel=r"Raman shift (cm$^{-1}$)", ylabel="Residual (a.u.)")
    fig.suptitle(
        f"{data_kind.upper()} | same baseline and samples | {fit['status']}\nFamily disagreement is retained; count accuracy is not established",
        fontsize=12,
    )
    for extension in ("png", "pdf"):
        fig.savefig(output / f"shape-comparison.{extension}", dpi=160)
    plt.close(fig)
