use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use compact_genome::interface::alphabet::Alphabet;
use generic_a_star::cost::{AStarCost, U32Cost};
use lib_ts_chainalign::costs::AlignmentCosts;
use lib_tsalign::config::TemplateSwitchConfig;
use log::{debug, info};
use sha::{
    sha1::Sha1,
    utils::{Digest, DigestExt},
};

pub fn load_tsa_config<AlphabetType: Alphabet, Cost: AStarCost>(
    configuration_directory: impl AsRef<Path>,
) -> Result<TemplateSwitchConfig<AlphabetType, Cost>> {
    let mut config_path = configuration_directory.as_ref().to_path_buf();
    info!("Loading alignment config directory {config_path:?}");

    config_path.push("config.tsa");
    let config_file = std::io::BufReader::new(
        std::fs::File::open(&config_path)
            .with_context(|| format!("Error opening config file {config_path:?}"))?,
    );
    let costs = TemplateSwitchConfig::read_plain(config_file)
        .with_context(|| format!("Error parsing template switch config from {config_path:?}"))?;
    Ok(costs)
}

pub fn infer_tschain_max_n(reference_length: usize, query_length: usize) -> usize {
    assert!(reference_length > 0 || query_length > 0);
    1 << (usize::BITS - (reference_length.max(query_length) - 1).leading_zeros())
}

pub fn infer_tschain_k(
    reference_length: usize,
    query_length: usize,
    max_anchor_mutations: u8,
) -> Result<u32> {
    // This evaluates to ceil(log_2(length_sum)).
    // The motivation is that there are length_sum k-mers,
    // so for each to be different, k needs to be at least ceil(log_4(length_sum)).
    // However, the birthday paradoxon states that for avoiding collisions,
    // the amount of possible k-mers needs to grow in the square of the amount available k-mers,
    // so we square that and arrive at ceil(log_2(length_sum)).
    let k = usize::BITS - ((reference_length + query_length) - 1).leading_zeros();
    // Decrease k a little, because we can hopefully afford a few more anchors.
    let k = k.saturating_sub(3).max(2);
    // Increase k by the maximum number of mutations in the anchors.
    let k = k.saturating_add(2 * u32::from(max_anchor_mutations));
    Ok(k)
}

pub fn tschain_preprocess_cache_file(
    alignment_costs: &AlignmentCosts<U32Cost>,
    cache_directory: impl AsRef<Path>,
    k: u32,
    max_anchor_mutations: u8,
    max_n: usize,
) -> Result<PathBuf> {
    let cost_hash = Sha1::default()
        .digest(
            &bincode::serde::encode_to_vec(alignment_costs, bincode::config::standard()).unwrap(),
        )
        .to_hex();
    debug!("Using cost_hash = {cost_hash}");

    Ok([
        cache_directory.as_ref().to_path_buf(),
        PathBuf::from(format!(
            "{cost_hash}-{k}-{max_anchor_mutations}-{max_n}.tsc"
        )),
    ]
    .iter()
    .collect())
}

#[cfg(test)]
mod tests {
    use crate::util::infer_tschain_max_n;

    #[test]
    fn fn_infer_tschain_max_n() {
        assert_eq!(infer_tschain_max_n(1, 1), 1);
        assert_eq!(infer_tschain_max_n(2, 2), 2);
        assert_eq!(infer_tschain_max_n(3, 3), 4);
        assert_eq!(infer_tschain_max_n(4, 4), 4);
        assert_eq!(infer_tschain_max_n(5, 5), 8);
        assert_eq!(infer_tschain_max_n(6, 6), 8);
        assert_eq!(infer_tschain_max_n(7, 7), 8);
        assert_eq!(infer_tschain_max_n(8, 8), 8);
        assert_eq!(infer_tschain_max_n(9, 9), 16);

        assert_eq!(infer_tschain_max_n(0, 1), 1);
        assert_eq!(infer_tschain_max_n(0, 2), 2);
        assert_eq!(infer_tschain_max_n(0, 3), 4);
        assert_eq!(infer_tschain_max_n(0, 4), 4);
        assert_eq!(infer_tschain_max_n(0, 5), 8);
        assert_eq!(infer_tschain_max_n(0, 6), 8);
        assert_eq!(infer_tschain_max_n(0, 7), 8);
        assert_eq!(infer_tschain_max_n(0, 8), 8);
        assert_eq!(infer_tschain_max_n(0, 9), 16);

        assert_eq!(infer_tschain_max_n(1, 0), 1);
        assert_eq!(infer_tschain_max_n(2, 0), 2);
        assert_eq!(infer_tschain_max_n(3, 0), 4);
        assert_eq!(infer_tschain_max_n(4, 0), 4);
        assert_eq!(infer_tschain_max_n(5, 0), 8);
        assert_eq!(infer_tschain_max_n(6, 0), 8);
        assert_eq!(infer_tschain_max_n(7, 0), 8);
        assert_eq!(infer_tschain_max_n(8, 0), 8);
        assert_eq!(infer_tschain_max_n(9, 0), 16);
    }
}
