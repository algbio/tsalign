from __future__ import annotations

import pathlib
from typing import List, Optional, Tuple, Union

from tsalign._tsalign import Aligner as _Aligner  # type: ignore[import]
from tsalign._tsalign import TSPairwiseAlignment as _TSPairwiseAlignment  # type: ignore[import]
from tsalign._types import (
    AlignmentOp,
    AlignmentRange,
    SimpleAlignmentOp,
    TemplateSwitchEntranceOp,
    TemplateSwitchExitOp,
    _parse_op,
)

_ALIGNER_KWARG_NAMES = frozenset(
    {"no_ts", "min_length_strategy", "chaining_strategy", "total_length_strategy", "costs", "costs_file"}
)

__all__ = [
    "Aligner",
    "Alignment",
    "align",
    "AlignmentRange",
    "AlignmentOp",
    "SimpleAlignmentOp",
    "TemplateSwitchEntranceOp",
    "TemplateSwitchExitOp",
]


class Alignment:
    """Result of a pairwise alignment that may contain template switches.

    Obtain an instance via :meth:`Aligner.align` or the module-level :func:`align`.
    """

    def __init__(self, inner: _TSPairwiseAlignment) -> None:
        self._inner = inner

    def cigar(self) -> Optional[str]:
        """CIGAR string of the alignment.

        Template switch operations are encoded as extended CIGAR tokens.
        Returns ``None`` if no valid alignment target was found.
        """
        return self._inner.cigar()

    def stats(self) -> dict:
        """Dictionary of alignment statistics.

        Keys include: ``cost``, ``cost_per_base``, ``duration_seconds``,
        ``opened_nodes``, ``closed_nodes``, ``template_switch_amount``, and
        nested ``result`` and ``sequences`` dicts.
        """
        return self._inner.stats()

    def alignments(self) -> Optional[List[Tuple[int, AlignmentOp]]]:
        """Compact list of alignment operations.

        Each entry is ``(count, op)`` where ``count`` is the repetition count and
        ``op`` is one of:

        - :class:`SimpleAlignmentOp` — a basic edit (match, substitution,
          insertion, deletion) in the primary or secondary track.
        - :class:`TemplateSwitchEntranceOp` — start of a template switch.
        - :class:`TemplateSwitchExitOp` — end of a template switch.

        Returns ``None`` if the alignment has no valid target.
        """
        raw = self._inner.alignments()
        if raw is None:
            return None
        return [(count, _parse_op(op)) for count, op in raw]

    def viz_template_switches(self) -> None:
        """Print an ASCII visualisation of template switch jumps to stdout."""
        self._inner.viz_template_switches()


class Aligner:
    """Pairwise DNA sequence aligner with template switch detection.

    Template switches are short-range translocations where a query region
    aligns to a different part of the reference, possibly on the reverse
    complement strand.  The aligner uses A* search under a configurable
    gap-affine cost model.

    Parameters
    ----------
    no_ts : bool
        Disable template switch detection (plain gap-affine alignment).
        Default: ``False``.
    min_length_strategy : str
        Strategy for enforcing the minimum template switch length.
        One of ``"none"``, ``"lookahead"`` (default), ``"preprocess_price"``,
        ``"preprocess_filter"``, ``"preprocess_lookahead"``.
    chaining_strategy : str
        A* lower-bound chaining strategy.
        One of ``"none"`` (default), ``"lower_bound"``.
    total_length_strategy : str
        Total template switch length strategy.
        One of ``"none"``, ``"maximise"`` (default).
    costs : str, optional
        Cost configuration as a raw ``.tsa``-format string.
    costs_file : str or Path, optional
        Path to a ``.tsa`` cost configuration file.
        Mutually exclusive with ``costs``.
        See ``sample_tsa_config/config.tsa`` for a skeleton.

    Examples
    --------
    Default settings::

        aligner = tsalign.Aligner()

    Custom cost file::

        aligner = tsalign.Aligner(costs_file="sample_tsa_config/config.tsa")
    """

    def __init__(
        self,
        *,
        no_ts: bool = False,
        min_length_strategy: str = "lookahead",
        chaining_strategy: str = "none",
        total_length_strategy: str = "maximise",
        costs: Optional[str] = None,
        costs_file: Optional[Union[str, pathlib.Path]] = None,
    ) -> None:
        if costs is not None and costs_file is not None:
            raise ValueError("Provide at most one of 'costs' or 'costs_file'.")
        if costs_file is not None:
            costs = pathlib.Path(costs_file).read_text()
        kwargs: dict = {
            "no_ts": no_ts,
            "min_length_strategy": min_length_strategy,
            "chaining_strategy": chaining_strategy,
            "total_length_strategy": total_length_strategy,
        }
        if costs is not None:
            kwargs["costs"] = costs
        self._inner = _Aligner(**kwargs)

    def align(
        self,
        reference: object,
        query: object,
        *,
        reference_name: str = "reference",
        query_name: str = "query",
        range: Optional[AlignmentRange] = None,
        reference_start: Optional[int] = None,
        reference_limit: Optional[int] = None,
        query_start: Optional[int] = None,
        query_limit: Optional[int] = None,
        cost_limit: Optional[int] = None,
        memory_limit: Optional[int] = None,
    ) -> Optional[Alignment]:
        """Align two DNA sequences, accounting for template switches.

        Parameters
        ----------
        reference : str-like
            Reference sequence.  Any object with a string representation is
            accepted (e.g. ``Bio.Seq``).
        query : str-like
            Query sequence.
        reference_name : str
            Label for the reference sequence.  Default: ``"reference"``.
        query_name : str
            Label for the query sequence.  Default: ``"query"``.
        range : AlignmentRange, optional
            Coordinate window to align within.  When provided, overrides the
            individual ``reference_start`` / ``reference_limit`` /
            ``query_start`` / ``query_limit`` arguments.
        reference_start : int, optional
            Start position in the reference (inclusive).  Default: ``0``.
        reference_limit : int, optional
            End position in the reference (exclusive).  Default: full length.
        query_start : int, optional
            Start position in the query (inclusive).  Default: ``0``.
        query_limit : int, optional
            End position in the query (exclusive).  Default: full length.
        cost_limit : int, optional
            Abandon the search and return ``None`` if the alignment cost would
            exceed this value.
        memory_limit : int, optional
            Abandon the search and return ``None`` if the number of open A*
            nodes exceeds this count.

        Returns
        -------
        Alignment or None
            ``None`` if no valid alignment was found within the given limits.
        """
        if range is not None:
            reference_start = range.reference_start
            reference_limit = range.reference_end
            query_start = range.query_start
            query_limit = range.query_end

        inner = self._inner.align(
            reference,
            query,
            reference_name,
            query_name,
            reference_start,
            reference_limit,
            query_start,
            query_limit,
            cost_limit,
            memory_limit,
        )
        return Alignment(inner) if inner is not None else None


def align(
    reference: object,
    query: object,
    **kwargs: object,
) -> Optional[Alignment]:
    """Align two DNA sequences in a single call.

    A convenience wrapper that creates a temporary :class:`Aligner` and
    immediately calls :meth:`~Aligner.align`.

    Keyword arguments that belong to the aligner (``no_ts``,
    ``min_length_strategy``, ``chaining_strategy``, ``total_length_strategy``,
    ``costs``, ``costs_file``) are forwarded to the constructor; all remaining
    keyword arguments are forwarded to :meth:`~Aligner.align`.

    Returns
    -------
    Alignment or None
        ``None`` if no valid alignment was found within the given limits.

    Examples
    --------
    ::

        result = tsalign.align("ACGTACGT", "ACGACGT")
        print(result.cigar())
    """
    aligner_kwargs = {k: v for k, v in kwargs.items() if k in _ALIGNER_KWARG_NAMES}
    align_kwargs = {k: v for k, v in kwargs.items() if k not in _ALIGNER_KWARG_NAMES}
    return Aligner(**aligner_kwargs).align(reference, query, **align_kwargs)
