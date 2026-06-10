"""Local-anesthetic systemic-absorption PK + the double-uncertainty view — v0.6 LA0/LA1.

The one LA-specific fact that *is* the safety message: **systemic absorption is
site-driven, not milligram-driven.** The same dose produces wildly different peak
plasma concentrations depending on the vascularity of the injection site (the
well-documented rank order intercostal > caudal/epidural > brachial plexus >
subcutaneous infiltration). A single mg/kg ceiling is therefore unreliable on its
face — and saying so *is* the safety message this subsystem most needs to deliver.

**LA0** curated disposition + site absorption + binding only. **LA1** adds the
toxicity-threshold *ranges* and the headline **double-uncertainty view** (§6): the
predicted plasma-concentration trajectory shown against the *range* of published
toxicity thresholds, with the free-concentration trace beside the total — framed as
a research/education question (*"how do site, agent, and binding move the predicted
concentration relative to the range of published thresholds, and how uncertain is
all of it?"*), **never** "is this dose safe?". It computes no dose, no ceiling, no
margin-as-guarantee (v0.6 §7). The honest finding is usually that the threshold
uncertainty *dwarfs* the PK uncertainty — i.e., the answer is "the threshold itself
is too uncertain to draw a line," and that conclusion is the safety message.

Kernel. A site-selected first-order absorption feeds a one-compartment linear
disposition; the result is the classic **Bateman** (first-order in / first-order
out) closed form:

    C(t) = (F·Dose·ka) / (V·(ka − k10)) · (exp(−k10·t) − exp(−ka·t))

forward-only, exact, no dose ever computed (v0.6 §7). The disposition (V, CL) is
read from the curated model record; the site sets ``ka``. Everything is a function
of a *given* dose at a *given* site — Hypnos never inverts one.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np

from .load import Dataset
from .models import worst_tier


@dataclass
class SiteAbsorption:
    """Predicted systemic plasma-concentration trajectory of an LA at one site."""

    model_id: str
    drug: str
    site: str
    rank: int
    ka: Optional[float]                 # first-order absorption rate (1/min); None => rank-only
    dose_mg: float
    t: np.ndarray
    cp: np.ndarray                      # plasma concentration (ug/mL)
    cmax: float                         # peak plasma concentration (ug/mL)
    tmax_min: float                     # time of peak (min)
    tier: str
    warnings: List[str] = field(default_factory=list)


def absorption_pk(
    dose_mg: float, ka: float, V_L: float, k10: float, t: np.ndarray,
    bioavailability: float = 1.0,
) -> np.ndarray:
    """First-order-absorption one-compartment plasma concentration (Bateman), ug/mL.

    ``V_L`` in litres, rate constants in 1/min, ``dose_mg`` in mg -> ug/mL (= mg/L).
    The ``ka == k10`` degenerate case uses the analytic limit ``t·k·exp(−k·t)``.
    """
    t = np.asarray(t, dtype=float)
    coeff = bioavailability * dose_mg / V_L
    if abs(ka - k10) < 1e-9:
        return coeff * k10 * t * np.exp(-k10 * t)
    return coeff * ka / (ka - k10) * (np.exp(-k10 * t) - np.exp(-ka * t))


def _disposition(model) -> Dict[str, float]:
    """Read the curated one-compartment disposition (V litres, CL L/min -> k10)."""
    p = {x.symbol: x.central for x in model.parameters}
    V = p.get("V1") or p.get("V")
    Cl = p.get("Cl1") or p.get("CL")
    if V is None or Cl is None:
        raise ValueError(f"{model.id}: LA disposition needs V1 and Cl1 parameters")
    return {"V": float(V), "k10": float(Cl) / float(V)}


def _site_entry(model, site: str) -> Dict:
    block = model.absorption or {}
    for s in block.get("site_rates", []):
        if s.get("site") == site:
            return s
    raise ValueError(
        f"{model.id}: no curated absorption site '{site}' "
        f"(have: {[s.get('site') for s in block.get('site_rates', [])]})"
    )


def concentration_at_site(
    ds: Dataset, model_id: str, *, site: str, dose_mg: float,
    t_min: float = 90.0, n: int = 541,
) -> SiteAbsorption:
    """Forward systemic plasma-concentration trajectory for a *given* LA dose at a
    *given* site. Forward-only; computes no dose, no ceiling, no margin (v0.6 §7)."""
    model = ds[model_id]
    if model.subsystem != "local_anesthetics":
        raise ValueError(f"{model_id} is not a local_anesthetics model")
    disp = _disposition(model)
    entry = _site_entry(model, site)
    ka = entry.get("ka")
    warnings: List[str] = []
    if ka is None:
        raise ValueError(
            f"{model_id}: site '{site}' carries only a rank, no ka magnitude — "
            "the trajectory is not simulable (rank curated, rate not)"
        )
    absn = model.absorption or {}
    tier = worst_tier([model.tier, absn.get("tier", "D")])
    warnings.append(
        "site absorption is Tier-C in magnitude (the rank order is robust, the "
        "absolute ka is not); systemic concentration only — NOT block efficacy, "
        "and no toxicity threshold is drawn (v0.6 LA0)"
    )

    t = np.linspace(0.0, t_min, n)
    cp = absorption_pk(dose_mg, ka, disp["V"], disp["k10"], t,
                       bioavailability=entry.get("bioavailability") or 1.0)
    i = int(np.argmax(cp))
    return SiteAbsorption(
        model_id=model_id, drug=model.drug_name, site=site, rank=int(entry["rank"]),
        ka=float(ka), dose_mg=float(dose_mg), t=t, cp=cp,
        cmax=float(cp[i]), tmax_min=float(t[i]), tier=tier, warnings=warnings,
    )


def site_comparison(
    ds: Dataset, model_id: str, *, dose_mg: float, t_min: float = 90.0, n: int = 541,
) -> List[SiteAbsorption]:
    """The headline: the *same* dose of the *same* drug at every curated site,
    sorted by peak plasma concentration (highest/most-dangerous first). Makes the
    site-of-injection dominance visible — the mg/kg ceiling is the wrong model."""
    model = ds[model_id]
    block = model.absorption or {}
    out: List[SiteAbsorption] = []
    for s in block.get("site_rates", []):
        if s.get("ka") is None:
            continue   # rank curated but rate not — cannot simulate this site
        out.append(concentration_at_site(ds, model_id, site=s["site"],
                                          dose_mg=dose_mg, t_min=t_min, n=n))
    out.sort(key=lambda r: r.cmax, reverse=True)
    return out


# --------------------------------------------------------------------------- #
# v0.6 LA1 — free concentration + the double-uncertainty view
# --------------------------------------------------------------------------- #
def saturable_free_concentration(
    c_total: np.ndarray, fraction_bound: float, capacity_ug_ml: float,
) -> np.ndarray:
    """Capacity-limited (Langmuir, 1:1) bound→free transform — the v0.6 LA3 non-linear
    free-fraction model.

    A single saturable high-affinity binding site (α1-acid glycoprotein for the
    amide LAs) of capacity ``C`` (µg/mL drug-equivalent) and affinity ``Ka``:

        bound(free) = C · Ka·free / (1 + Ka·free),   total = free + bound

    The affinity is **pinned by the curated low-concentration fraction bound** ``fb0``
    (so no second number is invented): in the low-concentration limit ``C·Ka = fb0 /
    (1 − fb0)``. Inverting ``total → free`` is a quadratic with the physical (positive)
    root::

        Ka·free² + (1 + C·Ka − Ka·total)·free − total = 0

    As ``total`` rises toward ``C`` the site saturates and the free fraction climbs
    from ``fb0`` toward 1 — the documented failure mode a linear trace hides. The
    low-concentration limit reduces exactly to the linear ``free = total·(1 − fb0)``.
    """
    c_total = np.asarray(c_total, dtype=float)
    fb0 = float(fraction_bound)
    if not (0.0 < fb0 < 1.0) or capacity_ug_ml <= 0:
        return c_total * (1.0 - fb0)            # degenerate -> linear
    ka = (fb0 / (1.0 - fb0)) / capacity_ug_ml   # pinned by the low-conc fraction bound
    b = 1.0 + capacity_ug_ml * ka - ka * c_total
    free = (-b + np.sqrt(b * b + 4.0 * ka * c_total)) / (2.0 * ka)
    return np.clip(free, 0.0, None)


def free_concentration(c_total: np.ndarray, protein_binding: Optional[Dict]) -> "FreeConcentration":
    """Bound→free transform (v0.6 §5). Returns the *linear* free trace
    ``c_free = c_total·(1 − fraction_bound)`` and, when a ``free_fraction_model`` is
    curated (v0.6 LA3), the **non-linear** (capacity-limited) free trace beside it so
    the under-prediction gap is explicit. For a ``saturable`` drug with no curated
    non-linear model, the linear trace carries the saturation failure-mode caveat
    rather than a fabricated curve (never-invent)."""
    pb = protein_binding or {}
    fb = pb.get("fraction_bound")
    c_total = np.asarray(c_total, dtype=float)
    warnings: List[str] = []
    if fb is None:
        return FreeConcentration(c_free=None, c_free_linear=None, free_fraction=None,
                                 saturable=False, model="none",
                                 warnings=["no curated protein binding — free concentration not derivable"])
    free_fraction = 1.0 - float(fb)
    c_free_linear = c_total * free_fraction
    saturable = bool(pb.get("saturable"))
    ffm = pb.get("free_fraction_model") or {}
    if saturable and ffm.get("type") == "capacity_limited" and ffm.get("binding_capacity_ug_ml"):
        c_free = saturable_free_concentration(c_total, float(fb),
                                              float(ffm["binding_capacity_ug_ml"]))
        warnings.append(
            "binding is SATURABLE: the NON-LINEAR (capacity-limited, Tier-D illustrative) "
            "free trace rises ABOVE the linear one at high total concentration — total "
            "concentration under-predicts free-drug risk exactly when risk is highest "
            "(v0.6 LA3; the magnitude is representative, the rise is the documented fact).")
        return FreeConcentration(c_free=c_free, c_free_linear=c_free_linear,
                                 free_fraction=free_fraction, saturable=True,
                                 model="capacity_limited", warnings=warnings)
    if saturable:
        warnings.append(
            "binding is SATURABLE: this LINEAR free trace under-predicts the free "
            "(toxic) concentration at high total concentration — the free fraction "
            "rises non-linearly exactly when risk is highest (no non-linear model curated; "
            "the gap is surfaced, never hidden — v0.6 LA3).")
    return FreeConcentration(c_free=c_free_linear, c_free_linear=c_free_linear,
                             free_fraction=free_fraction, saturable=saturable,
                             model="linear", warnings=warnings)


@dataclass
class FreeConcentration:
    c_free: Optional[np.ndarray]                 # the best curated free trace (non-linear where modeled)
    free_fraction: Optional[float]               # low-concentration (linear) free fraction
    saturable: bool
    c_free_linear: Optional[np.ndarray] = None   # the linear free trace (for the under-prediction gap)
    model: str = "linear"                        # none | linear | capacity_limited
    warnings: List[str] = field(default_factory=list)


@dataclass
class EndpointReadout:
    """One toxicity endpoint placed beside the prediction — *not* a margin verdict.

    Reports where the predicted peak falls *relative to* the published threshold
    range, on the matching basis (total vs free), and the threshold's own relative
    width. Never 'safe'/'unsafe', never a ceiling (v0.6 §6/§7)."""

    endpoint: str
    basis: str
    low: float
    high: float
    units: str
    relative_width: Optional[float]          # the threshold's own uncertainty
    predicted_peak: Optional[float]          # peak on the matching basis (total or free)
    position: str                            # "below range" | "within range" | "above range" | "unknown"
    tier: str
    individual_variability: Optional[str]
    citation: Optional[str]


@dataclass
class DoubleUncertainty:
    """The headline LA1 instrument (v0.6 §6): the predicted concentration band shown
    against the *threshold band*, with the dominant uncertainty named — a research /
    education view, **never** a dosing or safety verdict."""

    model_id: str
    drug: str
    site: str
    dose_mg: float
    peak_total: float                        # ug/mL
    peak_free: Optional[float]               # ug/mL — the best curated free peak (non-linear where modeled)
    tmax_min: float
    site_cmax_spread: Optional[float]        # max/min Cmax across sites — the PK (site) spread
    endpoints: List[EndpointReadout] = field(default_factory=list)
    dominant_uncertainty: str = ""
    cardiotoxicity: Optional[Dict] = None    # the drug's cardiotoxicity_class (v0.6 LA2), or None
    peak_free_linear: Optional[float] = None  # the LINEAR free peak (v0.6 LA3 — the under-prediction baseline)
    free_model: str = "linear"               # none | linear | capacity_limited (v0.6 LA3)
    tier: str = "D"
    warnings: List[str] = field(default_factory=list)


def _position(peak: Optional[float], low: float, high: float) -> str:
    if peak is None:
        return "unknown"
    if peak < low:
        return "below range"
    if peak > high:
        return "above range"
    return "within range"


def double_uncertainty(
    ds: Dataset, model_id: str, *, site: str, dose_mg: float,
    t_min: float = 90.0, n: int = 541,
) -> DoubleUncertainty:
    """The double-uncertainty view (v0.6 §6) for a *given* dose at a *given* site.

    Overlays the predicted total-plasma peak (and the linear free-concentration peak)
    against each curated toxicity-threshold *range*, on the matching basis, and names
    which uncertainty dominates. It computes **no** dose, ceiling, margin, or
    safe/unsafe verdict (v0.6 §7) — the threshold is curated as a range precisely so
    no single line can be read off it.
    """
    model = ds[model_id]
    if not model.has_toxicity_thresholds:
        raise ValueError(
            f"{model_id}: no curated toxicity thresholds — the double-uncertainty view "
            "needs threshold ranges (a missing threshold is a stated gap, never a fabricated line)")
    abs_run = concentration_at_site(ds, model_id, site=site, dose_mg=dose_mg, t_min=t_min, n=n)
    pb = ds.drug(model.drug_name).get("protein_binding") if hasattr(ds, "drug") else None
    free = free_concentration(abs_run.cp, pb)
    peak_free = float(np.max(free.c_free)) if free.c_free is not None else None
    peak_free_linear = float(np.max(free.c_free_linear)) if free.c_free_linear is not None else None

    # PK (site) spread: the SAME dose's peak across every simulable site — the
    # education point that drives the comparison, and the PK-uncertainty proxy.
    sites = site_comparison(ds, model_id, dose_mg=dose_mg, t_min=t_min, n=n)
    cmaxes = [s.cmax for s in sites]
    site_spread = (max(cmaxes) / min(cmaxes)) if cmaxes and min(cmaxes) > 0 else None

    endpoints: List[EndpointReadout] = []
    fold_ranges: List[float] = []
    tiers = [abs_run.tier]
    for th in model.toxicity_thresholds:
        peak = peak_free if th.basis == "free_plasma" else abs_run.cmax
        if th.low and th.low > 0:
            fold_ranges.append(th.high / th.low)
        tiers.append(th.tier)
        endpoints.append(EndpointReadout(
            endpoint=th.endpoint, basis=th.basis, low=th.low, high=th.high, units=th.units,
            relative_width=th.relative_width, predicted_peak=peak,
            position=_position(peak, th.low, th.high), tier=th.tier,
            individual_variability=th.individual_variability, citation=th.primary_citation,
        ))
    endpoints.sort(key=lambda e: e.low)

    # Which uncertainty dominates? Compare on the SAME scale — a multiplicative
    # fold-range (high/low for the threshold band; max/min peak for the site-driven
    # PK spread). For LA the threshold band almost always wins, and saying so is the
    # safety message: the threshold itself is too uncertain to draw a line (v0.6 §6).
    max_fold = max(fold_ranges) if fold_ranges else None
    if max_fold is not None and site_spread is not None:
        if max_fold >= site_spread:
            dominant = (
                f"THRESHOLD uncertainty dominates: the widest published threshold band spans "
                f"a {max_fold:.1f}x fold-range vs a {site_spread:.1f}x site-driven PK spread — "
                "no single safe-concentration line is defensible. This is the answer, not a gap "
                "to be closed.")
        else:
            dominant = (
                f"SITE/PK spread ({site_spread:.1f}x across injection sites) exceeds the widest "
                f"threshold band ({max_fold:.1f}x fold-range) — both remain too wide for a single line.")
    else:
        dominant = "uncertainty undecomposable (a band is missing) — no line can be drawn"

    # the LA0 trajectory warning says "no toxicity threshold is drawn"; that is the
    # LA0 stance, stale inside the LA1 view (which DOES show thresholds, as ranges) —
    # drop it here so the warning set stays accurate to this context.
    warnings = [w for w in abs_run.warnings if "no toxicity threshold is drawn" not in w]
    warnings += list(free.warnings)
    cardiotox = (ds.drug(model.drug_name) or {}).get("cardiotoxicity_class") if hasattr(ds, "drug") else None
    if cardiotox and cardiotox.get("cns_to_cvs_margin") == "narrow":
        warnings.append(
            f"CARDIOTOXICITY ({cardiotox.get('rank')}): the CNS-to-CVS margin is NARROW — "
            "cardiovascular toxicity can occur with little or no CNS prodrome; the agent-choice "
            "comparison is `hypnos la --agents` (v0.6 LA2).")
    # v0.6 LA3: where a non-linear free model is curated, name the under-prediction
    # gap — how much higher the saturable free peak is than the linear estimate.
    if free.model == "capacity_limited" and peak_free and peak_free_linear and peak_free_linear > 0:
        gap = peak_free / peak_free_linear
        warnings.append(
            f"FREE-FRACTION SATURATION: the non-linear free peak ({peak_free:.3g} ug/mL) is "
            f"{gap:.2f}x the linear estimate ({peak_free_linear:.3g}) at this dose — the linear "
            "trace under-predicts free-drug risk (Tier-D illustrative magnitude; v0.6 LA3).")
    warnings.append(
        "RESEARCH/EDUCATION ONLY — this is a double-uncertainty view, NOT a dosing tool: "
        "no maximum dose, no margin-as-guarantee, no 'is this safe?' answer (v0.6 §7).")
    return DoubleUncertainty(
        model_id=model_id, drug=model.drug_name, site=site, dose_mg=float(dose_mg),
        peak_total=abs_run.cmax, peak_free=peak_free, tmax_min=abs_run.tmax_min,
        site_cmax_spread=site_spread, endpoints=endpoints, dominant_uncertainty=dominant,
        cardiotoxicity=cardiotox, peak_free_linear=peak_free_linear, free_model=free.model,
        tier=worst_tier(tiers), warnings=warnings,
    )


# --------------------------------------------------------------------------- #
# v0.6 LA2 — stereochemistry / cardiotoxicity differentiation (the agent-choice axis)
# --------------------------------------------------------------------------- #
_CARDIOTOX_RANK = {"high": 0, "intermediate": 1, "low": 2}   # most cardiotoxic first


@dataclass
class AgentCardiotoxicity:
    """One local anesthetic's cardiotoxicity standing (v0.6 §3.4 / LA2).

    The agent-choice axis: why two agents at a *similar CNS threshold* can have very
    different cardiovascular margins. ``cns_to_cvs_fold`` is the curated numeric
    margin (cardiovascular midpoint / CNS-first-symptoms midpoint, total basis) where
    both are curated — a higher fold = a wider, safer CNS-warning margin."""

    model_id: str
    drug: str
    rank: str                                # high | intermediate | low
    stereochemistry: str
    cns_to_cvs_margin: str                   # qualitative: narrow | moderate | wide
    cns_to_cvs_fold: Optional[float]         # numeric margin where curated
    mechanism: str
    citation: Optional[str]
    tier: str


def _total_midpoint(model, endpoint: str) -> Optional[float]:
    vals = [th.midpoint for th in model.toxicity_thresholds
            if th.endpoint == endpoint and th.basis == "total_plasma"]
    return vals[0] if vals else None


def cardiotoxicity_comparison(ds: Dataset) -> List[AgentCardiotoxicity]:
    """The LA2 agent-choice view (v0.6 §3.4 / §6): every local anesthetic carrying a
    curated ``cardiotoxicity_class``, sorted most-cardiotoxic first, with its numeric
    CNS-to-CVS margin where both thresholds are curated. Shows *why agent choice
    changes the cardiovascular margin* — the key educational comparison — and never
    ranks or recommends a dose (forward, comparative, research-only; v0.6 §7)."""
    out: List[AgentCardiotoxicity] = []
    for model in ds:
        if model.subsystem != "local_anesthetics":
            continue
        cc = (ds.drug(model.drug_name) or {}).get("cardiotoxicity_class")
        if not cc:
            continue
        cns = _total_midpoint(model, "cns_first_symptoms")
        cvs = _total_midpoint(model, "cardiovascular")
        fold = (cvs / cns) if (cns and cvs and cns > 0) else None
        out.append(AgentCardiotoxicity(
            model_id=model.id, drug=model.drug_name, rank=cc.get("rank", "?"),
            stereochemistry=cc.get("stereochemistry", "?"),
            cns_to_cvs_margin=cc.get("cns_to_cvs_margin", "?"), cns_to_cvs_fold=fold,
            mechanism=cc.get("mechanism", ""), citation=cc.get("citation"), tier=model.tier,
        ))
    # most-cardiotoxic first: by class rank, then by the numeric CNS-to-CVS margin
    # (a narrower margin = more dangerous = earlier); uncurated margins sort last.
    out.sort(key=lambda a: (_CARDIOTOX_RANK.get(a.rank, 9),
                            a.cns_to_cvs_fold if a.cns_to_cvs_fold is not None else float("inf"),
                            a.drug))
    return out
