pub struct IndexedStr<'a> {
    index: Vec<usize>,
    string: &'a str,
}

impl<'a> IndexedStr<'a> {
    pub fn new(string: &'a str) -> Self {
        let index = string.char_indices().map(|(index, _)| index).collect();

        Self { index, string }
    }

    pub fn char_at(&self, index: usize) -> char {
        self.string[self.index[index]..].chars().next().unwrap()
    }
}

#[cfg(test)]
mod tests {
    use super::IndexedStr;

    #[test]
    fn char_at() {
        let string = IndexedStr::new("AcGT");
        assert_eq!(
            (0..4).map(|i| string.char_at(i)).collect::<String>(),
            "AcGT"
        );
    }

    #[test]
    fn char_at_multibyte() {
        let string = IndexedStr::new("aéb€c");
        assert_eq!(string.char_at(1), 'é');
        assert_eq!(string.char_at(2), 'b');
        assert_eq!(string.char_at(3), '€');
        assert_eq!(string.char_at(4), 'c');
    }

    #[test]
    #[should_panic]
    fn char_at_out_of_bounds() {
        IndexedStr::new("AC").char_at(2);
    }
}
