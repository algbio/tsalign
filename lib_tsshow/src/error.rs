use thiserror::Error;

pub type Result<T> = std::result::Result<T, Error>;

#[derive(Debug, Error)]
pub enum Error {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("Alignment is incomplete, and hence cannot be rendered.")]
    AlignmentHasNoTarget,

    #[error("No-TS alignment is incomplete, and hence cannot be rendered.")]
    NoTsAlignmentHasNoTarget,

    #[error("The no-TS alignment contains template switches.")]
    NoTsAlignmentHasTemplateSwitches,

    #[error("The no-TS alignment aligns other sequences than the TS alignment.")]
    NoTsAlignmentHasDifferentSequences,

    #[error("A negative anti-descendant gap is not supported for SVG generation.")]
    SvgNegativeAntiDescendantGap,

    #[error("Forward TSes are not yet supported.")]
    ForwardTsNotSupported,

    #[error("Error reducing the font {0} to the characters used in the SVG: {1}")]
    FontSubsetting(String, String),
}
