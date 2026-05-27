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

    #[error("A negative anti-descendant gap is not supported for SVG generation.")]
    SvgNegativeAntiDescendantGap,

    #[error("Forward TSes are not yet supported.")]
    ForwardTsNotSupported,
}
