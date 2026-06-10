"""Local-anesthetic systemic-absorption PK — v0.6 LA0 (the safety-first entry point).

The one LA-specific fact that *is* the safety message: **systemic absorption is
site-driven, not milligram-driven.** The same dose produces wildly different peak
plasma concentrations depending on the vascularity of the injection site (the
well-documented rank order intercostal > caudal/epidural > brachial plexus >
subcutaneous infiltration). A single mg/kg ceiling is therefore unreliable on its
face — and saying so *is* the safety message this subsystem most needs to deliver.

This module curates **disposition + site absorption only** — *no toxicity
thresholds* (those are deferred to v0.6 LA1, which needs its own safety framing).
So the subsystem's first release teaches the most important safety idea — that the
milligram ceiling is the wrong mental model — without yet drawing a single
threshold that could be misread.

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
