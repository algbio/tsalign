pub struct PanicOnExtend;

impl<T> Extend<T> for PanicOnExtend {
    fn extend<Iter: IntoIterator<Item = T>>(&mut self, iter: Iter) {
        assert!(iter.into_iter().next().is_none());
    }
}
