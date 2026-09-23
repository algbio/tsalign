pub struct PanicOnExtend;

impl<T> Extend<T> for PanicOnExtend {
    fn extend<Iter: IntoIterator<Item = T>>(&mut self, iter: Iter) {
        assert!(iter.into_iter().next().is_none());
    }
}

pub struct IgnoreExtend;

impl<T> Extend<T> for IgnoreExtend {
    fn extend<Iter: IntoIterator<Item = T>>(&mut self, _iter: Iter) {
        // Do nothing.
    }
}
