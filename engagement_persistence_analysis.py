#!/usr/bin/env python3
"""
Engagement persistence analysis (Python replication of analyze_engaged_bout_survival.m
and the statistics in analyze_reengage_gaps.m).

Starts from the CSVs exported by the MATLAB scripts, so it checks the statistics
independently of the MATLAB code. It does not rebuild bouts or episodes from frame
data; to check that step too, rebuild them from your own state labels and compare
against engaged_bouts_primary.csv / engagement_episodes_primary.csv.

Inputs (in --dir):
  engaged_bouts_primary.csv        one row per engaged bout (from analyze_engaged_bout_survival.m)
  engagement_episodes_primary.csv  one row per engagement episode
  engaged_bout_trials_primary.csv  one row per trial
  reengage_gaps.csv                (optional) one row per gap after a bout (from analyze_reengage_gaps.m)
  reengage_gap_trials.csv          (optional)

Definitions (match the MATLAB scripts):
  bout exit_type: 1 direct disengage, 2 NOTA -> disengage, 3 NOTA -> re-engage,
                  4 capture / window end (censored)
  episode: chain of bouts linked by exit type 3; event = ends in disengagement
           (last exit 1 or 2); censored at capture / window end
  bouts (negative control): event = any exit (1-3); censored = 4
  gap category 1 = contact (min range < 5 cm); re-engage = exit_type 3

Statistics:
  Cox proportional hazards (Efron ties), stratified by mouse, laser as the only
  covariate, fit separately in light and dark; pooled laser x dark interaction.
  Permutation: trial laser labels shuffled within mouse x lighting.
  Bootstrap: mice resampled with replacement.
  RMST: area under the Kaplan-Meier curve to tau, per mouse and condition;
        Wilcoxon signed-rank across mice, laser on vs off.

Usage:
  python engagement_persistence_analysis.py --dir /path/to/csvs --out /path/to/output
Requires: numpy, pandas, scipy, statsmodels, matplotlib
"""
import argparse, json, os
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
from statsmodels.duration.hazard_regression import PHReg
import matplotlib
import matplotlib.ticker
matplotlib.use("Agg")
import matplotlib.pyplot as plt

COND = [(0, 0), (0, 1), (1, 0), (1, 1)]                     # (dark, laser)
COND_NAMES = ["light off", "light on", "dark off", "dark on"]
LIGHT_NAMES = ["Light", "Dark"]
EXIT_NAMES = ["direct disengage", "NOTA -> disengage", "NOTA -> re-engage", "capture / window end"]
C_OFF, C_ON = "black", (0, 0.6, 0.9)
plt.rcParams.update({"font.size": 18})
pd.set_option("display.width", 200)
pd.set_option("display.max_columns", 20)


# ---------------------------------------------------------------- helpers
def cox(df, t, ev, x="laser"):
    """Stratified Cox: returns (beta, se, p) for column x."""
    X = df[[x]].to_numpy(float) if isinstance(x, str) else df[x].to_numpy(float)
    m = PHReg(df[t].to_numpy(float), X, status=df[ev].to_numpy(int),
              strata=df["animal"].to_numpy(), ties="efron").fit(disp=0)
    return m.params, m.bse, m.pvalues


def km(t, e):
    """Kaplan-Meier: distinct event times and S(t) just after each."""
    t, e = np.asarray(t, float), np.asarray(e, bool)
    ut = np.unique(t[e])
    if ut.size == 0:
        return ut, ut
    n_at = (t[:, None] >= ut[None, :]).sum(0)
    d_at = ((t[:, None] == ut[None, :]) & e[:, None]).sum(0)
    return ut, np.cumprod(1 - d_at / n_at)


def rmst(t, e, tau):
    tt, S = km(t, e)
    k = tt <= tau
    tt, S = tt[k], S[k]
    return float(np.sum(np.r_[1, S] * np.diff(np.r_[0, tt, tau])))


def per_mouse(df, fn, mice, min_n=5):
    """Matrix mice x 4 conditions of fn(sub-dataframe)."""
    M = np.full((len(mice), 4), np.nan)
    for i, m in enumerate(mice):
        for c, (d, l) in enumerate(COND):
            s = df[(df.animal == m) & (df.dark == d) & (df.laser == l)]
            if len(s) >= min_n:
                M[i, c] = fn(s)
    return M


def paired(M):
    """Light and dark: (mean off, mean on, n mice up, signed-rank p)."""
    out = []
    for c0 in (0, 2):
        a, b = M[:, c0], M[:, c0 + 1]
        ok = ~np.isnan(a) & ~np.isnan(b)
        p = wilcoxon(b[ok], a[ok]).pvalue if ok.sum() > 1 else np.nan
        out.append((np.nanmean(a), np.nanmean(b), int((b[ok] > a[ok]).sum()), int(ok.sum()), p))
    return out


def permute_labels(trials, rng):
    """Shuffle trial laser labels within mouse x lighting; returns Series trial_id -> laser."""
    lp = trials.set_index("trial_id").laser.copy()
    for _, g in trials.groupby(["animal", "dark"]):
        ids = g.trial_id.to_numpy()
        lp.loc[ids] = rng.permutation(lp.loc[ids].to_numpy())
    return lp


def mh_logor(df, y="re", x="laser"):
    """Mantel-Haenszel common log odds ratio across mice (0.5 continuity)."""
    num = den = 0.0
    for _, s in df.groupby("animal"):
        n = len(s)
        a = ((s[x] == 1) & (s[y] == 1)).sum(); b = ((s[x] == 1) & (s[y] == 0)).sum()
        c = ((s[x] == 0) & (s[y] == 1)).sum(); d = ((s[x] == 0) & (s[y] == 0)).sum()
        num += a * d / n; den += b * c / n
    return np.log((num + 0.5) / (den + 0.5))


def paired_panel(ax, M, ylabel, title):
    xs = [(1, 2), (4, 5)]
    for k, c0 in enumerate((0, 2)):
        for i in range(M.shape[0]):
            ax.plot(xs[k], M[i, c0:c0 + 2], "-", color="0.6", lw=1.5)
        ax.plot([xs[k][0]] * M.shape[0], M[:, c0], "o", color=C_OFF, ms=9)
        ax.plot([xs[k][1]] * M.shape[0], M[:, c0 + 1], "o", color=C_ON, ms=9)
    ax.set_xticks([1, 2, 4, 5]); ax.set_xticklabels(["off", "on", "off", "on"]); ax.set_xlim(0.3, 5.7)
    yl = ax.get_ylim()
    for x, lab in ((1.5, "light"), (4.5, "dark")):
        ax.text(x, yl[1], lab, ha="center", va="top")
    ax.set_ylabel(ylabel); ax.set_title(title, fontsize=18)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True, help="folder with the exported CSVs")
    ap.add_argument("--out", default=None, help="output folder (default: --dir)")
    ap.add_argument("--nperm", type=int, default=1000)
    ap.add_argument("--nboot", type=int, default=1000)
    ap.add_argument("--tau_ep", type=float, default=5.0, help="RMST horizon for episodes (s)")
    ap.add_argument("--tau_bout", type=float, default=1.5, help="RMST horizon for bouts (s)")
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    out = a.out or a.dir
    os.makedirs(out, exist_ok=True)
    rng = np.random.default_rng(a.seed)
    written = []
    R = {}

    B = pd.read_csv(os.path.join(a.dir, "engaged_bouts_primary.csv"))
    E = pd.read_csv(os.path.join(a.dir, "engagement_episodes_primary.csv"))
    T = pd.read_csv(os.path.join(a.dir, "engaged_bout_trials_primary.csv"))
    B["ev_any"] = B.exit_type.isin([1, 2, 3]).astype(int)
    E["event"] = E["event"].astype(int)
    mice = sorted(E.animal.unique())
    print(f"{len(B)} bouts, {len(E)} episodes, {len(T)} trials, {len(mice)} mice")

    # --- episodes: descriptive
    print("\n=== Episodes ===")
    desc = E.groupby(["dark", "laser"]).agg(n=("dur_s", "size"), median_s=("dur_s", "median"),
                                             mean_s=("dur_s", "mean"), bouts_per_ep=("n_bouts", "mean"),
                                             p_disengage=("event", "mean")).round(3)
    print(desc)
    R["episode_descriptive"] = desc.reset_index().to_dict("records")

    # --- Cox: episodes (primary) and bouts with any exit (negative control)
    print("\n=== Cox, stratified by mouse (laser HR, Wald 95% CI, p) ===")
    obs = {}
    for label, D, t, ev in (("episodes", E, "dur_s", "event"), ("bouts_any_exit", B, "dur_s", "ev_any")):
        R[label] = {}
        for d in (0, 1):
            b, se, p = cox(D[D.dark == d], t, ev)
            R[label][LIGHT_NAMES[d]] = dict(HR=float(np.exp(b[0])), CI=[float(np.exp(b[0] - 1.96 * se[0])),
                                            float(np.exp(b[0] + 1.96 * se[0]))], p=float(p[0]))
            if label == "episodes":
                obs[d] = b[0]
            print(f"  {label:15s} {LIGHT_NAMES[d]}: HR {np.exp(b[0]):.2f} "
                  f"[{np.exp(b[0]-1.96*se[0]):.2f} {np.exp(b[0]+1.96*se[0]):.2f}], p = {p[0]:.3g}")
    Ei = E.assign(dark_x=E.dark.astype(float), lxd=E.laser * E.dark)
    b, se, p = cox(Ei, "dur_s", "event", ["laser", "dark_x", "lxd"])
    R["episodes"]["interaction_p"] = float(p[2])
    print(f"  episodes laser x dark interaction p = {p[2]:.3g}")

    # --- permutation (episodes)
    perm = np.zeros((a.nperm, 2))
    for k in range(a.nperm):
        lp = permute_labels(T, rng)
        Ep = E.assign(lp=lp.loc[E.trial_id].to_numpy())
        for d in (0, 1):
            perm[k, d] = cox(Ep[Ep.dark == d], "dur_s", "event", "lp")[0][0]
    pperm = [(1 + np.sum(np.abs(perm[:, d]) >= abs(obs[d]))) / (a.nperm + 1) for d in (0, 1)]
    R["episodes"]["perm_p"] = dict(zip(LIGHT_NAMES, map(float, pperm)))
    print(f"\n=== Episode permutation ({a.nperm}): light p = {pperm[0]:.3g}, dark p = {pperm[1]:.3g}")

    # --- mouse bootstrap (episodes)
    bb = np.full((a.nboot, 2), np.nan)
    for k in range(a.nboot):
        pick = rng.choice(mice, len(mice))
        S = pd.concat([E[E.animal == m].assign(animal=j) for j, m in enumerate(pick)])
        for d in (0, 1):
            s = S[S.dark == d]
            if s.laser.nunique() == 2:
                bb[k, d] = cox(s, "dur_s", "event")[0][0]
    ci = {LIGHT_NAMES[d]: list(map(float, np.exp(np.nanpercentile(bb[:, d], [2.5, 97.5])))) for d in (0, 1)}
    ci_ratio = list(map(float, np.exp(np.nanpercentile(bb[:, 1] - bb[:, 0], [2.5, 97.5]))))
    R["episodes"]["boot_CI"] = ci; R["episodes"]["boot_ratio_dark_light_CI"] = ci_ratio
    print(f"=== Episode bootstrap ({a.nboot}): light {ci['Light'][0]:.2f}-{ci['Light'][1]:.2f}, "
          f"dark {ci['Dark'][0]:.2f}-{ci['Dark'][1]:.2f}, ratio dark/light {ci_ratio[0]:.2f}-{ci_ratio[1]:.2f}")

    # --- per-mouse RMST
    rm_ep = per_mouse(E, lambda s: rmst(s.dur_s, s.event, a.tau_ep), mice)
    rm_bt = per_mouse(B, lambda s: rmst(s.dur_s, s.ev_any, a.tau_bout), mice)
    print("\n=== Per-mouse RMST, off -> on (mean), n up / n, signed-rank p ===")
    for label, M, tau in (("episodes", rm_ep, a.tau_ep), ("bouts_any_exit", rm_bt, a.tau_bout)):
        res = paired(M)
        R[label]["rmst"] = dict(tau=tau, light=res[0], dark=res[1])
        for d in (0, 1):
            mo, mn, up, n, p = res[d]
            print(f"  {label:15s} tau {tau}: {LIGHT_NAMES[d]} {mo:.3f} -> {mn:.3f}, {up}/{n} up, p = {p:.3g}")

    # --- episode onset rate (per s outside episodes)
    T["outside_s"] = T.inwindow_s - T.episode_s
    U = T[T.outside_s > 1]
    on = per_mouse(U, lambda s: s.n_episodes.sum() / s.outside_s.sum(), mice, min_n=1)
    res = paired(on); R["episode_onset"] = dict(light=res[0], dark=res[1])
    print("\n=== Episode onset rate (per s outside episodes) ===")
    for d in (0, 1):
        mo, mn, up, n, p = res[d]
        print(f"  {LIGHT_NAMES[d]} {mo:.3f} -> {mn:.3f}, {up}/{n} up, p = {p:.3g}")

    # --- bout exit composition
    comp = np.array([[np.mean(B[(B.dark == d) & (B.laser == l)].exit_type == e) for e in (1, 2, 3, 4)]
                     for d, l in COND])
    R["exit_composition"] = dict(zip(COND_NAMES, comp.round(4).tolist()))
    print("\n=== Bout exit composition ===")
    print(pd.DataFrame(comp, index=COND_NAMES, columns=EXIT_NAMES).round(3))

    # --- contact gaps (optional)
    gpath = os.path.join(a.dir, "reengage_gaps.csv")
    if os.path.exists(gpath):
        G = pd.read_csv(gpath)
        C = G[G.category == 1].copy(); C["re"] = (C.exit_type == 3).astype(int)
        R["contact"] = {}
        print(f"\n=== Contact gaps: {len(C)} of {len(G)} gaps ===")
        mix = pd.crosstab(G.exit_type, G.category, normalize="index").round(3)
        print("gap category (cols 1 contact .. 5 other) by exit type:\n", mix)
        obs_or = [mh_logor(C[C.dark == d]) for d in (0, 1)]
        permor = np.zeros((a.nperm, 2))
        for k in range(a.nperm):
            lp = permute_labels(T, rng)
            Cp = C.assign(laser=lp.reindex(C.trial_id).to_numpy())
            for d in (0, 1):
                permor[k, d] = mh_logor(Cp[Cp.dark == d])
        pmh = [(1 + np.sum(np.abs(permor[:, d]) >= abs(obs_or[d]))) / (a.nperm + 1) for d in (0, 1)]
        pre = per_mouse(C, lambda s: s.re.mean(), mice)
        lat = per_mouse(C[C.re == 1], lambda s: s.dur_s.median(), mice)
        for d in (0, 1):
            r1, r2 = paired(pre)[d], paired(lat)[d]
            R["contact"][LIGHT_NAMES[d]] = dict(MH_OR=float(np.exp(obs_or[d])), perm_p=float(pmh[d]),
                                                P_reengage=r1, latency_s=r2)
            print(f"  {LIGHT_NAMES[d]}: P(re-engage|contact) {r1[0]:.2f} -> {r1[1]:.2f} ({r1[2]}/{r1[3]} up, p = {r1[4]:.3g}); "
                  f"MH OR {np.exp(obs_or[d]):.2f}, perm p = {pmh[d]:.3g}; latency {r2[0]:.2f} -> {r2[1]:.2f} s "
                  f"({r2[2]}/{r2[3]} longer, p = {r2[4]:.3g})")
        tpath = os.path.join(a.dir, "reengage_gap_trials.csv")
        if os.path.exists(tpath):
            GT = pd.read_csv(tpath)
            GT = GT[GT.engaged_s > 0]
            cr = per_mouse(GT, lambda s: s.n_contact_gaps.sum() / s.engaged_s.sum(), mice, min_n=1)
            for d in (0, 1):
                r = paired(cr)[d]; R["contact"][LIGHT_NAMES[d]]["contact_rate"] = r
                print(f"  {LIGHT_NAMES[d]}: contacts / s engaged {r[0]:.3f} -> {r[1]:.3f} ({r[2]}/{r[3]} up, p = {r[4]:.3g})")
        print("  peak cricket speed in contact gaps (median, cm/s):\n",
              C.groupby(["dark", "laser"]).max_cricket_spd.median().round(1))

    # --- figure (2 x 4), mirrors the MATLAB figure
    fig, ax = plt.subplots(2, 4, figsize=(30, 14))
    ax = ax.ravel()
    for d in (0, 1):
        for l in (0, 1):
            s = E[(E.dark == d) & (E.laser == l)]
            tt, S = km(s.dur_s, s.event)
            ax[d].step(np.r_[0, tt], np.r_[1, S], where="post", color=C_ON if l else C_OFF, lw=2.5,
                       label=f"laser {'on' if l else 'off'} ({len(s)})")
        h = R["episodes"][LIGHT_NAMES[d]]["HR"]
        ax[d].set(xlim=(0, 10), ylim=(0, 1), xlabel="Episode duration (s)", ylabel="P(still engaged)")
        ax[d].set_title(f"{'AB'[d]}. {LIGHT_NAMES[d]} episodes: HR {h:.2f}, p = {pperm[d]:.2g}", fontsize=18)
        ax[d].legend(frameon=False)
    for d in (0, 1):
        for l in (0, 1):
            s = B[(B.dark == d) & (B.laser == l)]
            tt, S = km(s.dur_s, s.ev_any)
            ax[2].step(np.r_[0, tt], np.r_[1, S], where="post", ls="--" if d == 0 else "-",
                       color=C_ON if l else C_OFF, lw=2.5, label=COND_NAMES[2 * d + l])
    ax[2].set(xlim=(0, a.tau_bout), ylim=(0, 1), xlabel="Bout duration (s)", ylabel="P(still in bout)")
    ax[2].set_title(f"C. Continuous bouts: HR light {R['bouts_any_exit']['Light']['HR']:.2f}, "
                    f"dark {R['bouts_any_exit']['Dark']['HR']:.2f}", fontsize=18)
    ax[2].legend(frameon=False)
    y, yt, yl = 0, [], []
    for d in (0, 1):
        for label, mk in (("episodes", "o"), ("bouts_any_exit", "s")):
            r = R[label][LIGHT_NAMES[d]]; y += 1
            col = C_ON if d else C_OFF
            ax[3].plot(r["CI"], [y, y], "-", color=col, lw=2); ax[3].plot(r["HR"], y, mk, color=col, ms=10)
            yt.append(y); yl.append(f"{LIGHT_NAMES[d]} {label.split('_')[0]}")
        y += 1
    ax[3].axvline(1, ls="--", color="k"); ax[3].set_xscale("log")
    ax[3].xaxis.set_major_locator(matplotlib.ticker.FixedLocator([0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2]))
    ax[3].xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%.1f"))
    ax[3].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax[3].set_yticks(yt); ax[3].set_yticklabels(yl); ax[3].invert_yaxis()
    ax[3].set_xlabel("Laser hazard ratio (< 1 = longer)"); ax[3].set_title("D. Primary gap tolerance (Wald 95% CI)", fontsize=18)
    pe, pb, po = paired(rm_ep), paired(rm_bt), paired(on)
    paired_panel(ax[4], rm_ep, f"RMST to {a.tau_ep:g} s (s)", f"E. Per mouse, episodes: p = {pe[0][4]:.2g}, {pe[1][4]:.2g}")
    paired_panel(ax[5], rm_bt, f"RMST to {a.tau_bout:g} s (s)", f"F. Per mouse, bouts: p = {pb[0][4]:.2g}, {pb[1][4]:.2g}")
    paired_panel(ax[6], on, "Episodes / s outside episodes", f"G. Episode onsets: p = {po[0][4]:.2g}, {po[1][4]:.2g}")
    bottom = np.zeros(4)
    for e, col in enumerate([(0.8, 0.2, 0.2), (0.95, 0.6, 0.5), (0.6, 0.6, 0.6), (0.3, 0.6, 0.3)]):
        ax[7].bar(range(4), comp[:, e], bottom=bottom, color=col, edgecolor="k", label=EXIT_NAMES[e]); bottom += comp[:, e]
    ax[7].set_xticks(range(4)); ax[7].set_xticklabels(COND_NAMES, rotation=30)
    ax[7].set(ylim=(0, 1), ylabel="Fraction of bouts"); ax[7].set_title("H. How bouts end", fontsize=18)
    ax[7].legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.25))
    for x in ax[[0, 1, 2, 3, 7]]:
        for sp in ("top", "right"):
            x.spines[sp].set_visible(False)
    fig.suptitle("Engagement episodes (primary) and continuous bouts (Python replication). "
                 "p: A, B permutation; E-G light, dark signed-rank", fontsize=18)
    fig.tight_layout()
    f = os.path.join(out, "engagement_persistence_python.pdf"); fig.savefig(f); written.append(f)

    # --- save
    per = []
    for i, m in enumerate(mice):
        for c, (d, l) in enumerate(COND):
            per.append(dict(animal=m, dark=d, laser=l, rmst_episode=rm_ep[i, c], rmst_bout=rm_bt[i, c],
                            episode_onset_rate=on[i, c]))
    f = os.path.join(out, "engagement_persistence_python_per_mouse.csv"); pd.DataFrame(per).to_csv(f, index=False); written.append(f)
    f = os.path.join(out, "engagement_persistence_python_results.json")
    with open(f, "w") as fh:
        json.dump(R, fh, indent=2, default=float)
    written.append(f)
    print("\nFiles written:"); [print("  " + os.path.abspath(w)) for w in written]


if __name__ == "__main__":
    main()
