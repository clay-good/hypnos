"""Typed, read-only views over Hypnos dataset records.

These dataclasses are thin wrappers around the validated JSON. They add
attribute access, tier ordering helpers, and envelope/failure-mode evaluation,
but never mutate the underlying data — the dataset stays the single source of
truth.
"""
from __future__ import annotations

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
class Parameter:
    symbol: str
    value: Dict[str, Any]
    tier: str
    covariate_model: Optional[str] = None
    label: Optional[str] = None
    primary_citation: Optional[str] = None
    extraction: Dict[str, Any] = field(default_factory=dict)

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

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        return f"<Model {self.id} tier={self.tier} review={self.review_status}>"
