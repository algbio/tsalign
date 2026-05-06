from __future__ import annotations

from dataclasses import dataclass
from typing import Union


@dataclass
class AlignmentRange:
    """Coordinate bounds for a pairwise alignment."""

    reference_start: int = 0
    reference_end: "int | None" = None
    query_start: int = 0
    query_end: "int | None" = None


@dataclass
class SimpleAlignmentOp:
    """A basic edit operation in the primary or secondary alignment track.

    The ``kind`` attribute is one of: ``PrimaryMatch``, ``PrimarySubstitution``,
    ``PrimaryInsertion``, ``PrimaryDeletion``, ``PrimaryFlankMatch``,
    ``PrimaryFlankSubstitution``, ``PrimaryFlankInsertion``, ``PrimaryFlankDeletion``,
    ``SecondaryMatch``, ``SecondarySubstitution``, ``SecondaryInsertion``,
    ``SecondaryDeletion``.
    """

    kind: str


@dataclass
class TemplateSwitchEntranceOp:
    """Entrance into a template switch region."""

    kind: str  # always "TemplateSwitchEntrance"
    first_offset: int
    primary: str    # "Reference" or "Query"
    secondary: str  # "Reference" or "Query"
    direction: str  # "Forward" or "Reverse"
    equal_cost_range: dict  # {min_start, max_start, min_end, max_end}


@dataclass
class TemplateSwitchExitOp:
    """Exit from a template switch region."""

    kind: str  # always "TemplateSwitchExit"
    anti_primary_gap: int


AlignmentOp = Union[SimpleAlignmentOp, TemplateSwitchEntranceOp, TemplateSwitchExitOp]


def _parse_op(raw: object) -> AlignmentOp:
    if isinstance(raw, str):
        return SimpleAlignmentOp(kind=raw)
    if isinstance(raw, dict):
        if "TemplateSwitchEntrance" in raw:
            d = raw["TemplateSwitchEntrance"]
            return TemplateSwitchEntranceOp(
                kind="TemplateSwitchEntrance",
                first_offset=d["first_offset"],
                primary=d["primary"],
                secondary=d["secondary"],
                direction=d["direction"],
                equal_cost_range=d["equal_cost_range"],
            )
        if "TemplateSwitchExit" in raw:
            return TemplateSwitchExitOp(
                kind="TemplateSwitchExit",
                anti_primary_gap=raw["TemplateSwitchExit"]["anti_primary_gap"],
            )
    raise ValueError(f"Unknown alignment op: {raw!r}")
