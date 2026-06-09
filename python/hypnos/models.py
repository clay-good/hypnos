"""Typed, read-only views over Hypnos dataset records.

These dataclasses are thin wrappers around the validated JSON. They add
attribute access, tier ordering helpers, and envelope/failure-mode evaluation,
but never mutate the underlying data — the dataset stays the single source of
truth.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

# Tier ordering: A is best, D is worst. "Worst input wins" => max by this rank.
TIER_RANK = {"A": 0, "B": 1, "C": 2, "D": 3}


def worst_tier(tiers: List[str]) -> str:
    """Return the worst (lowest-confidence) tier among inputs. D > C > B > A."""
    present = [t for t in tiers if t in TIER_RANK]
    if not present:
        return "D"
    return max(present, key=lambda t: TIER_RANK[t])


# Internal concentrations are always mg/L (== ug/mL). Drugs declare a conventional
# display unit (opioids and dexmedetomidine are reported in ng/mL); this maps that
# unit to the multiplier from the internal ug/mL value.
_CONC_FACTORS = {"ug/ml": 1.0, "µg/ml": 1.0, "mcg/ml": 1.0, "mg/l": 1.0, "ng/ml": 1000.0}


def concentration_factor(unit: Optional[str]) -> float:
    """Multiplier from the internal ug/mL concentration to the given display unit."""
    if not unit:
        return 1.0
    return _CONC_FACTORS.get(unit.strip().lower(), 1.0)


@dataclass(frozen=True)
class Range:
    min: Optional[float] = None
    max: Optional[float] = None

    def contains(self, x: float) -> bool:
        if x is None:
            return True
        if self.min is not None and x < self.min:
            return False
        if self.max is not None and x > self.max:
            return False
        return True


@dataclass(frozen=True)
class ParameterVariability:
    """Between-subject (and optional inter-occasion) variability for one parameter.

    ``omega2`` is the canonical eta-scale variance (the NONMEM ``$OMEGA`` diagonal
    entry); ``cv_percent`` is a derived convenience kept checked-consistent with it
    (v0.2 spec §4 Trap 1).
    """

    omega2: Optional[float] = None
    cv_percent: Optional[float] = None
    shrinkage_percent: Optional[float] = None
    distribution: str = "exponential"
    iov_omega2: Optional[float] = None
    tier: str = "D"
    primary_citation: Optional[str] = None
    extraction: Dict[str, Any] = field(default_factory=dict)

    @property
    def cv_from_omega2(self) -> Optional[float]:
        """Exact log-normal CV% recomputed from omega2 (the consistency reference)."""
        if self.omega2 is None or self.omega2 < 0:
            return None
        return 100.0 * math.sqrt(math.exp(self.omega2) - 1.0)


@dataclass(frozen=True)
class ResidualError:
    """The Sigma layer — residual unexplained variability (v0.2 spec §3.2)."""

    model: str
    proportional: Optional[Dict[str, Any]] = None
    additive: Optional[Dict[str, Any]] = None
    log: Optional[Dict[str, Any]] = None
    tier: str = "D"
    primary_citation: Optional[str] = None
    extraction: Dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class OmegaBlock:
    """Published off-diagonal Omega — correlated random effects (v0.2 spec §3.3)."""

    correlations: List[Dict[str, Any]] = field(default_factory=list)
    complete: bool = False
    tier: str = "D"
    primary_citation: Optional[str] = None
    extraction: Dict[str, Any] = field(default_factory=dict)


def _parse_param_variability(raw: Optional[Dict[str, Any]]) -> Optional[ParameterVariability]:
    if not raw:
        return None
    bsv = raw.get("bsv") or {}
    iov = raw.get("iov") or {}
    return ParameterVariability(
        omega2=bsv.get("omega2"),
        cv_percent=bsv.get("cv_percent"),
        shrinkage_percent=bsv.get("shrinkage_percent"),
        distribution=bsv.get("distribution", "exponential"),
        iov_omega2=iov.get("omega2"),
        tier=raw.get("tier", "D"),
        primary_citation=raw.get("primary_citation"),
        extraction=raw.get("extraction", {}),
    )


@dataclass(frozen=True)
class Parameter:
    symbol: str
    value: Dict[str, Any]
    tier: str
    covariate_model: Optional[str] = None
    label: Optional[str] = None
    primary_citation: Optional[str] = None
    extraction: Dict[str, Any] = field(default_factory=dict)
    variability: Optional[ParameterVariability] = None

    @property
    def central(self) -> Optional[float]:
        return self.value.get("central")

    @property
    def units(self) -> Optional[str]:
        return self.value.get("units")


@dataclass(frozen=True)
class FailureMode:
    condition: str
    behavior: str
    action: str
    predicate: Optional[str] = None
    citation: Optional[str] = None


@dataclass(frozen=True)
class Envelope:
    age_years: Range = field(default_factory=Range)
    weight_kg: Range = field(default_factory=Range)
    height_cm: Range = field(default_factory=Range)
    bmi_kg_m2: Range = field(default_factory=Range)
    populations: List[str] = field(default_factory=list)
    derivation_n: Optional[int] = None

    _COVARIATE_RANGES = ("age_years", "weight_kg", "height_cm", "bmi_kg_m2")

    def check(self, patient: Dict[str, Any]) -> List[str]:
        """Return a list of human-readable envelope-violation messages (empty == in envelope)."""
        violations: List[str] = []
        mapping = {
            "age_years": ("age", patient.get("age")),
            "weight_kg": ("weight", patient.get("weight")),
            "height_cm": ("height", patient.get("height")),
            "bmi_kg_m2": ("bmi", _bmi(patient)),
        }
        for attr in self._COVARIATE_RANGES:
            rng: Range = getattr(self, attr)
            label, val = mapping[attr]
            if val is None:
                continue
            if not rng.contains(val):
                violations.append(
                    f"{label}={val:g} outside derivation envelope "
                    f"[{_fmt(rng.min)}, {_fmt(rng.max)}]"
                )
        return violations


def _fmt(x: Optional[float]) -> str:
    return "-inf" if x is None else f"{x:g}"


def _bmi(patient: Dict[str, Any]) -> Optional[float]:
    w = patient.get("weight")
    h = patient.get("height")
    if w is None or h is None or h == 0:
        return None
    return w / ((h / 100.0) ** 2)


@dataclass(frozen=True)
class Model:
    """A single PK/PD model record."""

    raw: Dict[str, Any]

    # --- identity ---------------------------------------------------------
    @property
    def id(self) -> str:
        return self.raw["id"]

    @property
    def subsystem(self) -> str:
        return self.raw.get("subsystem", self.id.split(".")[0])

    @property
    def drug(self) -> Dict[str, Any]:
        return self.raw["drug"]

    @property
    def drug_name(self) -> str:
        return self.raw["drug"]["name"]

    @property
    def purpose(self) -> str:
        return self.raw["purpose"]

    @property
    def label(self) -> str:
        return self.raw.get("label", self.id)

    @property
    def tier(self) -> str:
        return self.raw["tier"]

    @property
    def primary_citation(self) -> str:
        return self.raw["primary_citation"]

    # --- structure --------------------------------------------------------
    @property
    def structure(self) -> Dict[str, Any]:
        return self.raw["structure"]

    @property
    def n_compartments(self) -> int:
        return self.structure["compartments"]

    @property
    def has_effect_compartment(self) -> bool:
        return bool(self.structure.get("effect_compartment", False))

    @property
    def kernel_implemented(self) -> bool:
        return bool(self.raw.get("kernel", {}).get("implemented", False))

    @property
    def kernel_function(self) -> Optional[str]:
        return self.raw.get("kernel", {}).get("function")

    # --- curation ---------------------------------------------------------
    @property
    def review_status(self) -> str:
        return self.raw["extraction"]["review_status"]

    @property
    def extraction(self) -> Dict[str, Any]:
        return self.raw["extraction"]

    @property
    def parameters(self) -> List[Parameter]:
        out = []
        for p in self.raw.get("parameters", []):
            out.append(
                Parameter(
                    symbol=p["symbol"],
                    value=p["value"],
                    tier=p["tier"],
                    covariate_model=p.get("covariate_model"),
                    label=p.get("label"),
                    primary_citation=p.get("primary_citation"),
                    extraction=p.get("extraction", {}),
                    variability=_parse_param_variability(p.get("variability")),
                )
            )
        return out

    def param(self, symbol: str) -> Parameter:
        for p in self.parameters:
            if p.symbol == symbol:
                return p
        raise KeyError(f"{self.id} has no parameter {symbol!r}")

    @property
    def covariates(self) -> Dict[str, Any]:
        return self.raw.get("covariates", {})

    @property
    def applicability_envelope(self) -> Envelope:
        env = self.raw.get("applicability_envelope", {})
        return Envelope(
            age_years=Range(**env.get("age_years", {})),
            weight_kg=Range(**env.get("weight_kg", {})),
            height_cm=Range(**env.get("height_cm", {})),
            bmi_kg_m2=Range(**env.get("bmi_kg_m2", {})),
            populations=env.get("populations", []),
            derivation_n=env.get("derivation_n"),
        )

    @property
    def known_failure_modes(self) -> List[FailureMode]:
        return [
            FailureMode(
                condition=f["condition"],
                behavior=f["behavior"],
                action=f["action"],
                predicate=f.get("predicate"),
                citation=f.get("citation"),
            )
            for f in self.raw.get("known_failure_modes", [])
        ]

    @property
    def predictive_performance(self) -> List[Dict[str, Any]]:
        return self.raw.get("predictive_performance", [])

    @property
    def predictive_mdape(self) -> List[Dict[str, Any]]:
        """Published MDAPE (inaccuracy) entries that apply *in envelope*.

        Returned as ``[{"value", "citation"}, ...]``. Out-of-envelope / failure-mode
        entries (e.g. the Minto James-LBM number measured in morbid obesity) are
        excluded: they characterize the model only where it would itself be greyed
        out, so they are not the accuracy a reader should attach to an *included*
        (in-envelope) model in the divergence view.
        """
        out: List[Dict[str, Any]] = []
        for pp in self.predictive_performance:
            if pp.get("metric") != "MDAPE":
                continue
            if "out of envelope" in (pp.get("population") or "").lower():
                continue
            out.append({"value": pp["value"], "citation": pp.get("citation")})
        return out

    # --- variability (v0.2 population-variability layer) -------------------
    @property
    def variability_status(self) -> str:
        """Rollup: 'none' | 'partial' | 'diagonal' | 'full' (default 'none')."""
        return self.raw.get("variability_status", "none")

    @property
    def has_published_variability(self) -> bool:
        """True when the model carries any curated random-effects structure.

        The never-synthesize rule (spec §5): a model with no published BSV draws
        no band — a missing band is a true statement, a borrowed one is a lie.
        """
        return self.variability_status not in (None, "none")

    @property
    def residual_error(self) -> Optional[ResidualError]:
        raw = self.raw.get("residual_error")
        if not raw:
            return None
        return ResidualError(
            model=raw["model"],
            proportional=raw.get("proportional"),
            additive=raw.get("additive"),
            log=raw.get("log"),
            tier=raw.get("tier", "D"),
            primary_citation=raw.get("primary_citation"),
            extraction=raw.get("extraction", {}),
        )

    @property
    def omega_block(self) -> Optional[OmegaBlock]:
        raw = self.raw.get("omega_block")
        if not raw:
            return None
        return OmegaBlock(
            correlations=raw.get("correlations", []),
            complete=raw.get("complete", False),
            tier=raw.get("tier", "D"),
            primary_citation=raw.get("primary_citation"),
            extraction=raw.get("extraction", {}),
        )

    def bsv_omegas(self) -> Dict[str, float]:
        """Map of structural-parameter symbol -> omega2 for parameters carrying BSV."""
        out: Dict[str, float] = {}
        for p in self.parameters:
            if p.variability and p.variability.omega2 is not None:
                out[p.symbol] = float(p.variability.omega2)
        return out

    @property
    def variability_tier(self) -> Optional[str]:
        """Worst tier among the curated variability components (None if uncurated).

        The variability layer is typically the least externally-validated component,
        so this usually drives ``band_tier`` below the point-estimate tier (spec §5).
        """
        tiers = [p.variability.tier for p in self.parameters
                 if p.variability and p.variability.omega2 is not None]
        if self.residual_error is not None:
            tiers.append(self.residual_error.tier)
        if self.omega_block is not None:
            tiers.append(self.omega_block.tier)
        return worst_tier(tiers) if tiers else None

    @property
    def band_tier(self) -> Optional[str]:
        """Tier of a prediction band: worst of the structural and variability tiers.

        The median line keeps its v0.1 tier; the band around it may be labeled
        lower (spec §5). ``None`` when the model publishes no variability (no band).
        """
        if not self.has_published_variability:
            return None
        vt = self.variability_tier
        return worst_tier([self.tier, vt]) if vt else self.tier

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        return f"<Model {self.id} tier={self.tier} review={self.review_status}>"
