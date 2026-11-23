use std::fmt::Display;

use serde::{Deserialize, Serialize};

use crate::alignment_history::history_vec::{
    AlignmentHistory, unsigned_int_alignment_history_vec::UnsignedIntAlignmentHistoryVec,
};

/// Holds the parameters for the chaining algorithm.
///
/// Not every combination of parameters is supported, see [`Self::validate`].
#[derive(Clone, Serialize, Deserialize)]
pub struct ChainingLowerBoundConfig {
    pub anchor_k: u32,
    pub max_anchor_mutations: u8,
    pub inexact_lower_bound: InexactLowerBoundKind,
}

/// Chooses the kind of lower bound to use for chaining with inexact anchors.
#[derive(Clone, Serialize, Deserialize)]
pub enum InexactLowerBoundKind {
    /// Compute a lower bound assuming that in every window of size `anchor_k` there have to be at least `max_anchor_mutations + 1` mutations.
    WindowedMinMutations,
    /// Compute a lower bound assuming that there is no match run longer than `anchor_k - 1 - max_anchor_mutations`,
    /// except for alignments that are shorter than `anchor_k` on both sequences.
    MaxMatchRun,
}

/// A reason why a [`ChainingLowerBoundConfig`] is not supported.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InvalidChainingLowerBoundConfig {
    /// `max_anchor_mutations` is at least `anchor_k`.
    /// This also covers `anchor_k = 0`.
    TooManyAnchorMutations {
        anchor_k: u32,
        max_anchor_mutations: u8,
    },
    /// `anchor_k` is too large for the k-mer stores of exact anchors.
    ExactAnchorKTooLarge { anchor_k: u32 },
    /// `anchor_k` is too large for inexact anchors.
    InexactAnchorKTooLarge {
        anchor_k: u32,
        max_anchor_mutations: u8,
    },
    /// The alignment histories of [`InexactLowerBoundKind::WindowedMinMutations`] cannot store anchors of this size.
    WindowedMinMutationsCapacityExceeded {
        anchor_k: u32,
        max_anchor_mutations: u8,
    },
}

impl ChainingLowerBoundConfig {
    /// The largest supported `anchor_k` for exact anchors (`max_anchor_mutations = 0`).
    pub const MAX_EXACT_ANCHOR_K: u32 = 64;

    /// The largest supported `anchor_k` for inexact anchors (`max_anchor_mutations > 0`).
    pub const MAX_INEXACT_ANCHOR_K: u32 = u8::MAX as u32;

    /// Checks if this combination of parameters is supported by preprocessing and alignment.
    ///
    /// The supported parameters are:
    ///
    /// * `max_anchor_mutations < anchor_k` (so `anchor_k >= 1`);
    /// * `anchor_k <= 64` (see [`Self::MAX_EXACT_ANCHOR_K`]) for exact anchors (`max_anchor_mutations = 0`);
    /// * `anchor_k <= 255` (see [`Self::MAX_INEXACT_ANCHOR_K`]) for inexact anchors (`max_anchor_mutations > 0`);
    /// * for [`InexactLowerBoundKind::WindowedMinMutations`] with inexact anchors,
    ///   `anchor_k + ceil(max_anchor_mutations / 2)` must fit into the largest alignment history.
    pub fn validate(&self) -> Result<(), InvalidChainingLowerBoundConfig> {
        let anchor_k = self.anchor_k;
        let max_anchor_mutations = self.max_anchor_mutations;

        if anchor_k <= u32::from(max_anchor_mutations) {
            return Err(InvalidChainingLowerBoundConfig::TooManyAnchorMutations {
                anchor_k,
                max_anchor_mutations,
            });
        }

        if max_anchor_mutations == 0 {
            if anchor_k > Self::MAX_EXACT_ANCHOR_K {
                return Err(InvalidChainingLowerBoundConfig::ExactAnchorKTooLarge { anchor_k });
            }
        } else {
            let Ok(small_anchor_k) = u8::try_from(anchor_k) else {
                return Err(InvalidChainingLowerBoundConfig::InexactAnchorKTooLarge {
                    anchor_k,
                    max_anchor_mutations,
                });
            };

            if matches!(
                self.inexact_lower_bound,
                InexactLowerBoundKind::WindowedMinMutations
            ) && !UnsignedIntAlignmentHistoryVec::<u128>::has_enough_capacity(
                small_anchor_k,
                max_anchor_mutations,
            ) {
                return Err(
                    InvalidChainingLowerBoundConfig::WindowedMinMutationsCapacityExceeded {
                        anchor_k,
                        max_anchor_mutations,
                    },
                );
            }
        }

        Ok(())
    }
}

impl Display for InvalidChainingLowerBoundConfig {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::TooManyAnchorMutations {
                anchor_k,
                max_anchor_mutations,
            } => write!(
                f,
                "max_anchor_mutations = {max_anchor_mutations} must be less than k = {anchor_k}; \
                 use a larger k or fewer anchor mutations"
            ),
            Self::ExactAnchorKTooLarge { anchor_k } => write!(
                f,
                "k = {anchor_k} is too large for exact anchors (max_anchor_mutations = 0); \
                 the supported range is 1 <= k <= {}",
                ChainingLowerBoundConfig::MAX_EXACT_ANCHOR_K
            ),
            Self::InexactAnchorKTooLarge {
                anchor_k,
                max_anchor_mutations,
            } => write!(
                f,
                "k = {anchor_k} is too large for inexact anchors (max_anchor_mutations = {max_anchor_mutations}); \
                 the supported range is max_anchor_mutations < k <= {}",
                ChainingLowerBoundConfig::MAX_INEXACT_ANCHOR_K
            ),
            Self::WindowedMinMutationsCapacityExceeded {
                anchor_k,
                max_anchor_mutations,
            } => write!(
                f,
                "k = {anchor_k} with max_anchor_mutations = {max_anchor_mutations} is too large for the \
                 windowed-min-mutations lower bound; use a smaller k or the max-match-run lower bound"
            ),
        }
    }
}

impl std::error::Error for InvalidChainingLowerBoundConfig {}

#[cfg(test)]
mod tests {
    use super::{ChainingLowerBoundConfig, InexactLowerBoundKind, InvalidChainingLowerBoundConfig};

    fn config(
        anchor_k: u32,
        max_anchor_mutations: u8,
        inexact_lower_bound: InexactLowerBoundKind,
    ) -> ChainingLowerBoundConfig {
        ChainingLowerBoundConfig {
            anchor_k,
            max_anchor_mutations,
            inexact_lower_bound,
        }
    }

    fn validate(
        anchor_k: u32,
        max_anchor_mutations: u8,
    ) -> Result<(), InvalidChainingLowerBoundConfig> {
        config(
            anchor_k,
            max_anchor_mutations,
            InexactLowerBoundKind::MaxMatchRun,
        )
        .validate()
    }

    #[test]
    fn accepts_supported_parameters() {
        for (anchor_k, max_anchor_mutations) in [
            (1, 0),
            (9, 0),
            (64, 0),
            (2, 1),
            (6, 5),
            (9, 1),
            (255, 1),
            (255, 254),
        ] {
            assert_eq!(
                validate(anchor_k, max_anchor_mutations),
                Ok(()),
                "k = {anchor_k}, max_anchor_mutations = {max_anchor_mutations}"
            );
        }
    }

    #[test]
    fn rejects_too_many_anchor_mutations() {
        for (anchor_k, max_anchor_mutations) in [(0, 0), (1, 1), (6, 6), (6, 7), (255, 255)] {
            assert_eq!(
                validate(anchor_k, max_anchor_mutations),
                Err(InvalidChainingLowerBoundConfig::TooManyAnchorMutations {
                    anchor_k,
                    max_anchor_mutations
                }),
            );
        }
    }

    #[test]
    fn rejects_too_large_k() {
        assert_eq!(
            validate(65, 0),
            Err(InvalidChainingLowerBoundConfig::ExactAnchorKTooLarge { anchor_k: 65 })
        );
        assert_eq!(
            validate(256, 1),
            Err(InvalidChainingLowerBoundConfig::InexactAnchorKTooLarge {
                anchor_k: 256,
                max_anchor_mutations: 1
            })
        );
        assert_eq!(
            validate(u32::MAX, 200),
            Err(InvalidChainingLowerBoundConfig::InexactAnchorKTooLarge {
                anchor_k: u32::MAX,
                max_anchor_mutations: 200
            })
        );
    }

    #[test]
    fn windowed_min_mutations_checks_alignment_history_capacity() {
        let windowed = |anchor_k, max_anchor_mutations| {
            config(
                anchor_k,
                max_anchor_mutations,
                InexactLowerBoundKind::WindowedMinMutations,
            )
            .validate()
        };
        // The largest alignment history stores (128 - 6) / 2 = 61 operations, and k = 60 with one mutation needs 61.
        assert_eq!(windowed(60, 1), Ok(()));
        assert_eq!(
            windowed(61, 1),
            Err(
                InvalidChainingLowerBoundConfig::WindowedMinMutationsCapacityExceeded {
                    anchor_k: 61,
                    max_anchor_mutations: 1
                }
            )
        );
        // Exact anchors do not use alignment histories.
        assert_eq!(windowed(64, 0), Ok(()));
    }
}
