from __future__ import annotations

import os
import pathlib
from typing import Iterable, List, Optional, Tuple, Union

from tsalign._tsalign import Aligner as _Aligner  # type: ignore[import]
from tsalign._tsalign import ArenaFullError  # type: ignore[import]
from tsalign._tsalign import CHAIN_ALIGNER_DEFAULT_COSTS as _CHAIN_ALIGNER_DEFAULT_COSTS  # type: ignore[import]
from tsalign._tsalign import ChainAligner as _ChainAligner  # type: ignore[import]
from tsalign._tsalign import TSPairwiseAlignment as _TSPairwiseAlignment  # type: ignore[import]
from tsalign._types import (
    AlignmentOp,
    AlignmentRange,
    ArenaIndex,
    SimpleAlignmentOp,
    TemplateSwitchEntranceOp,
    TemplateSwitchExitOp,
    _parse_op,
)

_ALIGNER_KWARG_NAMES = frozenset(
    {
        "no_ts",
        "min_length_strategy",
        "chaining_strategy",
        "total_length_strategy",
        "costs",
        "costs_file",
        "arena_index",
    }
)

#: The default costs of :class:`ChainAligner` as a raw ``.tsa``-format string:
#: only reverse template switches.  Equal to
#: ``test_files/config/chainalign/config.tsa`` in the repository.
CHAIN_ALIGNER_DEFAULT_COSTS: str = _CHAIN_ALIGNER_DEFAULT_COSTS

_CHAIN_ALIGNER_KWARG_NAMES = frozenset(
    {
        "costs",
        "costs_file",
        "k",
        "max_anchor_mutations",
        "inexact_chaining_lower_bound",
        "max_chaining_successors",
        "max_exact_cost_function_cost",
        "chaining_closed_list",
        "chaining_open_list",
        "anchor_pruning",
        "inexact_pruning_rules",
        "exact_overlap_diagonal_distance",
        "cache_directory",
        "force_no_preprocessing",
        "arena_index",
    }
)

# The alignment methods, named as in the CLI's --alignment-method.
_METHOD_A_STAR_TEMPLATE_SWITCH = "a-star-template-switch"
_METHOD_A_STAR_CHAIN_TS = "a-star-chain-ts"

__all__ = [
    "Aligner",
    "ChainAligner",
    "CHAIN_ALIGNER_DEFAULT_COSTS",
    "Alignment",
    "align",
    "AlignmentRange",
    "AlignmentOp",
    "ArenaFullError",
    "ArenaIndex",
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

        Alignments computed by :class:`ChainAligner` additionally have the
        keys ``chain_cost`` (the cost of the chain the alignment was built
        from, which may be larger than ``cost`` for inexact anchors) and
        ``chain_statistics`` (see :meth:`chain_statistics`).
        """
        return self._inner.stats()

    def chain_statistics(self) -> Optional[dict]:
        """Statistics about the anchors and the chaining of tschainalign.

        Contains the anchor counts before and after pruning, the pruning mode
        and rules, the number of chaining rounds and the times of the phases.
        Returns ``None`` if the alignment was not computed by :class:`ChainAligner`.
        """
        return self.stats().get("chain_statistics")

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
        Only ``"none"`` (default) is supported.
        ``"lower_bound"`` is disabled and raises :class:`ValueError`
        (see https://github.com/algbio/tsalign/issues/248).
    total_length_strategy : str
        Total template switch length strategy.
        One of ``"none"``, ``"maximise"`` (default).
    costs : str, optional
        Cost configuration as a raw ``.tsa``-format string.
        Default: the built-in default costs, which are also the CLI's default
        for ``--alignment-method a-star-template-switch``.  They equal
        ``sample_tsa_config/config.tsa`` in the repository.
    costs_file : str or Path, optional
        Path to a ``.tsa`` cost configuration file.
        Mutually exclusive with ``costs``.
        See ``sample_tsa_config/config.tsa`` for a skeleton.
    arena_index : str
        The integer type of the node indices of the A* search: ``"u32"``
        (default) or ``"u64"``.  With ``"u32"``, the search can store at
        most 4294967295 nodes at the same time, which suffices unless it
        needs a few hundred GB of memory; if it needs more,
        :meth:`align` raises :class:`ArenaFullError`.  ``"u64"`` removes
        this limit, but needs more memory per node (about 3 % in total).
        It does not change the alignment.  Like the CLI's ``--arena-index``.

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
        arena_index: ArenaIndex = "u32",
    ) -> None:
        if costs is not None and costs_file is not None:
            raise ValueError("Provide at most one of 'costs' or 'costs_file'.")
        if chaining_strategy == "lower_bound":
            raise ValueError(
                "chaining_strategy 'lower_bound' is disabled because it is too expensive "
                "and did not help in practice (https://github.com/algbio/tsalign/issues/248). "
                "Use 'none'."
            )
        if costs_file is not None:
            costs = pathlib.Path(costs_file).read_text()
        kwargs: dict = {
            "no_ts": no_ts,
            "min_length_strategy": min_length_strategy,
            "chaining_strategy": chaining_strategy,
            "total_length_strategy": total_length_strategy,
            "arena_index": arena_index,
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
        gap_characters: str = "",
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
            All coordinates refer to the sequences after removing
            ``gap_characters``.
        reference_start : int, optional
            Start position in the reference (inclusive).  Default: ``0``.
        reference_limit : int, optional
            End position in the reference (exclusive).  Default: full length.
        query_start : int, optional
            Start position in the query (inclusive).  Default: ``0``.
        query_limit : int, optional
            End position in the query (exclusive).  Default: full length.
        gap_characters : str
            ASCII characters to remove from both sequences before aligning,
            e.g. ``"-"``.
            Default: ``""`` (remove nothing).
        cost_limit : int, optional
            Return ``None`` if there is no alignment with cost at most this
            value.  The limit is inclusive: an alignment whose cost equals the
            limit is returned.
        memory_limit : int, optional
            Approximate memory limit in bytes.  Abandon the search and return
            ``None`` if it is exceeded.

        Returns
        -------
        Alignment or None
            ``None`` if no valid alignment was found within the given limits.

        Raises
        ------
        ValueError
            If a sequence contains a character that is not part of the DNA
            alphabet (``A``, ``C``, ``G``, ``T``, ``N``) after removing
            ``gap_characters``, if ``gap_characters`` contains non-ASCII
            characters, or if the range does not fit the sequences without
            gap characters.
        ArenaFullError
            If the search needs more nodes at the same time than the
            ``arena_index`` of the aligner allows.  Create the aligner with
            ``arena_index="u64"`` to align these sequences.
        """
        if range is not None:
            reference_start = range.reference_start
            reference_limit = range.reference_end
            query_start = range.query_start
            query_limit = range.query_end

        inner = self._inner.align(
            reference,
            query,
            reference_name=reference_name,
            query_name=query_name,
            reference_start=reference_start,
            reference_limit=reference_limit,
            query_start=query_start,
            query_limit=query_limit,
            gap_characters=gap_characters,
            cost_limit=cost_limit,
            memory_limit=memory_limit,
        )
        return Alignment(inner) if inner is not None else None


def _read_costs_file(costs_file: Union[str, os.PathLike]) -> str:
    path = pathlib.Path(costs_file)
    if path.is_dir():
        # Like the CLI's --configuration-directory.
        path = path / "config.tsa"
    return path.read_text()


class ChainAligner:
    """Pairwise DNA sequence aligner with template switches, based on chaining anchors (tschainalign).

    This is the alignment method ``a-star-chain-ts`` of the CLI.  It chains
    k-mer anchors (exact, or with up to ``max_anchor_mutations`` mutations)
    with A* and aligns the gaps between them exactly.  It uses lower bounds of
    the gap costs that are expensive to compute (preprocessing), but depend
    only on the costs, ``k``, ``max_anchor_mutations``, the inexact lower
    bound and the maximum sequence length ``max_n`` (the longer sequence
    length, rounded up to a power of two).

    The aligner keeps all lower bounds it computed or loaded in memory, so
    repeated calls of :meth:`align` preprocess only once per ``k`` and
    ``max_n``.  With ``cache_directory``, it also loads them from and stores
    them into the same cache files as the CLI (``tsalign preprocess`` and
    ``tsalign align --cache-directory``), so Python and the CLI can share a
    cache directory.

    tschainalign supports only the characters ``A``, ``C``, ``G`` and ``T``
    in the sequences, and only cost configurations with uniform edit costs and
    zero-range jump costs, like its default costs :data:`CHAIN_ALIGNER_DEFAULT_COSTS`
    (``test_files/config/chainalign/config.tsa`` in the repository).

    All parameters are keyword-only, named and with the defaults as the
    options of ``tsalign align --alignment-method a-star-chain-ts``.

    Parameters
    ----------
    costs : str, optional
        Cost configuration as a raw ``.tsa``-format string.
    costs_file : str or Path, optional
        Path to a ``.tsa`` cost configuration file, or to a directory
        containing ``config.tsa`` (like the CLI's ``--configuration-directory``).
        At most one of ``costs`` and ``costs_file`` may be given.  Default:
        :data:`CHAIN_ALIGNER_DEFAULT_COSTS`, the built-in default costs of
        tschainalign, which are also the CLI's default for
        ``--alignment-method a-star-chain-ts``.  The classic default costs of
        :class:`Aligner` are not supported by tschainalign.
    k : int, optional
        The anchor size.  Default: inferred from the sequence lengths and
        ``max_anchor_mutations`` like in the CLI.  Fixing it saves
        preprocessing when aligning sequences of different lengths.
    max_anchor_mutations : int
        The maximum number of mutations in an anchor.  Default: ``0`` (exact anchors).
    inexact_chaining_lower_bound : str, optional
        The lower bound for inexact anchors.  Only ``"max-match-run"`` (default)
        is supported; ``"windowed-min-mutations"`` is disabled as in the CLI.
        Only for ``max_anchor_mutations >= 1``.
    max_chaining_successors : int
        The maximum number of successors to generate while processing a
        node during chaining.  Default: ``10``.
    max_exact_cost_function_cost : int
        The cost until which the chaining cost function is initialised
        exactly.  Default: ``0``.
    chaining_closed_list : str
        ``"special"`` (default) or ``"fx-hash-map"``.
    chaining_open_list : str
        ``"linear-heap"`` (default) or ``"std-heap"``.
    anchor_pruning : str, optional
        Which inexact anchors to remove before chaining: ``"exact"`` (default;
        only where this provably does not change the chain cost), ``"inexact"``
        (heuristic rules, may increase the cost) or ``"none"``.
        Only for ``max_anchor_mutations >= 1``.
    inexact_pruning_rules : str or iterable of str, optional
        The rules of ``anchor_pruning="inexact"``, by the same names as the
        CLI's ``--inexact-pruning-rules``: a comma-separated string such as
        ``"all,-secondary-exact-overlap"``, or a list of rule names and keywords
        (``all``, ``none``, ``primary``, ``secondary``).  Default: all rules.
    exact_overlap_diagonal_distance : int, optional
        The diagonal distance of the ``*-exact-overlap`` rules of
        ``anchor_pruning="inexact"``.  Default: ``1``.
    cache_directory : str or Path, optional
        Directory for the CLI-compatible cache files of the lower bounds.
        Default: no files, the lower bounds are kept in memory only.
    force_no_preprocessing : bool
        Never preprocess in :meth:`align`; raise :class:`FileNotFoundError`
        if the lower bounds are neither in memory (from :meth:`preprocess` or
        an earlier :meth:`align`) nor in the cache directory.  Default: ``False``.
    arena_index : str
        The integer type of the node indices of the chaining search:
        ``"u32"`` (default) or ``"u64"``.  With ``"u32"``, the chaining can
        store at most 4294967295 nodes at the same time (about 200 GB of
        memory); if it needs more, :meth:`align` raises
        :class:`ArenaFullError`.  ``"u64"`` removes this limit.  It does not
        change the alignment.  Like the CLI's ``--arena-index``.

    Raises
    ------
    ValueError
        If a parameter is invalid or not supported, if the costs cannot be
        parsed or are not supported by tschainalign (e.g. its length costs), or if
        options for inexact anchors are given with ``max_anchor_mutations=0``.

    Examples
    --------
    ::

        aligner = tsalign.ChainAligner(
            max_anchor_mutations=1,
            cache_directory="tsalign-cache",
        )
        result = aligner.align(reference, query)
        print(result.stats()["cost"], result.chain_statistics()["chaining_rounds"])
    """

    def __init__(
        self,
        *,
        costs: Optional[str] = None,
        costs_file: Optional[Union[str, os.PathLike]] = None,
        k: Optional[int] = None,
        max_anchor_mutations: int = 0,
        inexact_chaining_lower_bound: Optional[str] = None,
        max_chaining_successors: int = 10,
        max_exact_cost_function_cost: int = 0,
        chaining_closed_list: str = "special",
        chaining_open_list: str = "linear-heap",
        anchor_pruning: Optional[str] = None,
        inexact_pruning_rules: Optional[Union[str, Iterable[str]]] = None,
        exact_overlap_diagonal_distance: Optional[int] = None,
        cache_directory: Optional[Union[str, os.PathLike]] = None,
        force_no_preprocessing: bool = False,
        arena_index: ArenaIndex = "u32",
    ) -> None:
        if costs is not None and costs_file is not None:
            raise ValueError("Provide at most one of 'costs' or 'costs_file'.")
        if costs_file is not None:
            costs = _read_costs_file(costs_file)
        if inexact_pruning_rules is not None and not isinstance(inexact_pruning_rules, str):
            inexact_pruning_rules = ",".join(inexact_pruning_rules)
        self._inner = _ChainAligner(
            costs,
            k=k,
            max_anchor_mutations=max_anchor_mutations,
            inexact_chaining_lower_bound=inexact_chaining_lower_bound,
            max_chaining_successors=max_chaining_successors,
            max_exact_cost_function_cost=max_exact_cost_function_cost,
            chaining_closed_list=chaining_closed_list,
            chaining_open_list=chaining_open_list,
            anchor_pruning=anchor_pruning,
            inexact_pruning_rules=inexact_pruning_rules,
            exact_overlap_diagonal_distance=exact_overlap_diagonal_distance,
            cache_directory=None if cache_directory is None else os.fspath(cache_directory),
            force_no_preprocessing=force_no_preprocessing,
            arena_index=arena_index,
        )

    def preprocess(self, max_length: int) -> None:
        """Preprocess the lower bounds for sequences up to ``max_length``, like ``tsalign preprocess``.

        Computes the lower bounds for ``max_n`` = ``max_length`` rounded up to
        a power of two, and for all smaller powers of two, with the given
        ``k`` or with all ``k`` inferred for these lengths.  They are kept in
        memory and written to the cache directory, if any.  Lower bounds that
        are already in memory are not recomputed.

        Raises
        ------
        ValueError
            If ``max_length`` is less than 1, or an inferred ``k`` is not supported.
        """
        self._inner.preprocess(max_length)

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
        gap_characters: str = "",
        cost_limit: Optional[int] = None,
        iteration_limit: Optional[int] = None,
    ) -> Optional[Alignment]:
        """Align two DNA sequences with tschainalign.

        The parameters ``reference`` to ``gap_characters`` are the same as for
        :meth:`Aligner.align`.  Preprocesses the lower bounds if needed, see
        the class documentation.  The GIL is released while preprocessing and
        aligning.

        Parameters
        ----------
        cost_limit : int, optional
            Return ``None`` if the optimal chain costs more than this value
            (inclusive limit).  For inexact anchors, the limit applies to the
            chain cost, which may be larger than the alignment cost.
        iteration_limit : int, optional
            The maximum number of chaining rounds.  Return ``None`` if the
            optimal chain is not found within this many rounds.

        Returns
        -------
        Alignment or None
            ``None`` if a limit was reached.

        Raises
        ------
        ValueError
            If a sequence contains a character other than ``A``, ``C``, ``G``
            and ``T`` after removing ``gap_characters``, if both sequences are
            empty, if the range does not fit the sequences, if a limit is
            negative, or if the inferred ``k`` is not supported.
        FileNotFoundError
            If ``force_no_preprocessing`` is set and the lower bounds are not available.
        OSError
            If reading or writing a cache file fails.
        ArenaFullError
            If the chaining search needs more nodes at the same time than the
            ``arena_index`` of the aligner allows.  Create the aligner with
            ``arena_index="u64"`` to align these sequences.
        """
        if range is not None:
            reference_start = range.reference_start
            reference_limit = range.reference_end
            query_start = range.query_start
            query_limit = range.query_end

        inner = self._inner.align(
            reference,
            query,
            reference_name=reference_name,
            query_name=query_name,
            reference_start=reference_start,
            reference_limit=reference_limit,
            query_start=query_start,
            query_limit=query_limit,
            gap_characters=gap_characters,
            cost_limit=cost_limit,
            iteration_limit=iteration_limit,
        )
        return Alignment(inner) if inner is not None else None


def align(
    reference: object,
    query: object,
    *,
    method: str = _METHOD_A_STAR_TEMPLATE_SWITCH,
    **kwargs: object,
) -> Optional[Alignment]:
    """Align two DNA sequences in a single call.

    A convenience wrapper that creates a temporary aligner and immediately
    calls its ``align`` method.

    ``method`` chooses the aligner, named as the CLI's ``--alignment-method``:
    ``"a-star-template-switch"`` (default, :class:`Aligner`) or
    ``"a-star-chain-ts"`` (:class:`ChainAligner`).

    Keyword arguments that belong to the aligner's constructor (for
    :class:`Aligner`: ``no_ts``, ``min_length_strategy``,
    ``chaining_strategy``, ``total_length_strategy``, ``costs``,
    ``costs_file``, ``arena_index``; for :class:`ChainAligner`: all of its parameters) are
    forwarded to the constructor; all remaining keyword arguments are
    forwarded to ``align``.

    Since the aligner is temporary, :class:`ChainAligner` preprocesses on
    every call unless a ``cache_directory`` is given.  Create a
    :class:`ChainAligner` to align many sequences.

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
    if method == _METHOD_A_STAR_TEMPLATE_SWITCH:
        aligner_class, aligner_kwarg_names = Aligner, _ALIGNER_KWARG_NAMES
    elif method == _METHOD_A_STAR_CHAIN_TS:
        aligner_class, aligner_kwarg_names = ChainAligner, _CHAIN_ALIGNER_KWARG_NAMES
    else:
        raise ValueError(
            f"invalid method {method!r}, expected one of "
            f"{_METHOD_A_STAR_TEMPLATE_SWITCH!r}, {_METHOD_A_STAR_CHAIN_TS!r}"
        )
    aligner_kwargs = {k: v for k, v in kwargs.items() if k in aligner_kwarg_names}
    align_kwargs = {k: v for k, v in kwargs.items() if k not in aligner_kwarg_names}
    return aligner_class(**aligner_kwargs).align(reference, query, **align_kwargs)
