use error::{Error, Result};
use lib_tsalign::{
    alignment_result::{AlignmentResult, AlignmentType},
    costs::U64Cost,
};
use resvg::usvg::Options;

pub mod error;
pub mod plain_text;
pub mod svg;
pub mod ts_arrangement;

/// Checks that the no-TS alignment, if given, aligns the same sequences as the TS alignment and has no template switches.
fn check_no_ts_result(
    result: &AlignmentResult<U64Cost>,
    no_ts_result: &Option<AlignmentResult<U64Cost>>,
) -> Result<()> {
    let Some(no_ts_result) = no_ts_result else {
        return Ok(());
    };

    let sequences = &result.statistics().sequences;
    let no_ts_sequences = &no_ts_result.statistics().sequences;
    if sequences.reference != no_ts_sequences.reference || sequences.query != no_ts_sequences.query
    {
        return Err(Error::NoTsAlignmentHasDifferentSequences);
    }

    if let AlignmentResult::WithTarget { alignment, .. } = no_ts_result
        && alignment.iter_compact().any(|(_, alignment_type)| {
            matches!(alignment_type, AlignmentType::TemplateSwitchEntrance { .. })
        })
    {
        return Err(Error::NoTsAlignmentHasTemplateSwitches);
    }

    Ok(())
}

pub fn svg_to_png(svg_in: &[u8], zoom: f32) -> Vec<u8> {
    log::info!("Converting SVG to PNG");

    let mut options = Options::default();
    // resvg ignores the fonts embedded into the SVG,
    // hence they are loaded from the same font files that were embedded.
    for font_file in svg::embedded_font_files() {
        options.fontdb_mut().load_font_data(font_file.to_vec());
    }
    options.fontdb_mut().load_system_fonts();
    let svg = resvg::usvg::Tree::from_data(svg_in, &options).unwrap();

    let raster_image_size = svg.size();
    let raster_image_width = (raster_image_size.width().ceil() * zoom) as u32;
    let raster_image_height = (raster_image_size.height().ceil() * zoom) as u32;
    log::info!("PNG size: {raster_image_width}x{raster_image_height}",);

    let mut raster_image =
        resvg::tiny_skia::Pixmap::new(raster_image_width, raster_image_height).unwrap();
    resvg::render(
        &svg,
        resvg::usvg::Transform::from_scale(zoom, zoom),
        &mut raster_image.as_mut(),
    );
    raster_image.encode_png().unwrap()
}
