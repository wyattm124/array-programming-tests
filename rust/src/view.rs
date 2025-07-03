use std::ops;

// TODO : should compile time verify that we do not generate an index out of bounds error
// NOTE : indexer should be a FAST, simple, indexing function. 95% of the time it should be
//  a linear transformation.

pub trait Indexer: Fn(usize) -> usize + Copy + Clone {}
impl<G: Fn(usize) -> usize + Copy + Clone> Indexer for G {}

#[derive(Debug, Clone, Copy)]
pub struct View<'a, T, F: Indexer, const N: usize> {
    data: &'a [T],
    indexer: F,
}

#[derive(Debug)]
pub struct ViewMut<'a, T, F: Indexer, const N: usize> {
    data: &'a mut [T],
    indexer: F,
}

impl<'a, T, F: Indexer, const N: usize> View<'a, T, F, N> {
    pub fn new(data: &'a [T], indexer: impl Indexer) -> Self {
        Self { data, indexer }
    }
    pub fn ident_new(data: &'a [T]) -> Self {
        Self { data, indexer: |i| i }
    }
    pub fn subview<const M: usize>(&self, start: usize) -> View<'_, T, F, M> {
        View::new(&self.data[start..start + M], self.indexer.clone())
    }
    pub fn compose(self, indexer: impl Indexer) -> Self {
        View::new(self.data, move |i| (self.indexer)(indexer(i)))
    }
}

impl<'a, T, F: Indexer, const N: usize> ViewMut<'a, T, F, N> {
    pub fn new(data: &'a mut [T], indexer: impl Indexer) -> Self {
        Self { data, indexer }
    }
    pub fn ident_new(data: &'a mut [T]) -> Self {
        Self { data, indexer: |i| i }
    }
    pub fn as_view(&self) -> View<'_, T, F, N> {
        View::new(self.data, self.indexer.clone())
    }
}

impl<'a, T, F: Indexer, const N: usize> ops::Index<usize> for View<'a, T, F, N> {
    type Output = T;
    fn index(&self, index: usize) -> &Self::Output {
        &self.data[(self.indexer)(index)]
    }
}

impl<'a, T, F: Indexer, const N: usize> ops::Index<usize> for ViewMut<'a, T, F, N> {
    type Output = T;
    fn index(&self, index: usize) -> &Self::Output {
        &self.data[(self.indexer)(index)]
    }
}

impl<'a, T, F: Indexer, const N: usize> ops::IndexMut<usize> for ViewMut<'a, T, F, N> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.data[(self.indexer)(index)]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identity_view() {
        let data = [0, 1, 2, 3, 4];
        let view = View::<_, _, _, 5>::ident_new(&data);
        assert_eq!(view[0], 0);
        assert_eq!(view[1], 1);
        assert_eq!(view[2], 2);
        assert_eq!(view[3], 3);
        assert_eq!(view[4], 4);
    }

    #[test]
    fn permuted_view() {
        let data = [0, 1, 2, 3, 4];
        let view = View::<_, _, _, 5>::new(&data, |i| (4 - i));
        assert_eq!(view[0], 4);
        assert_eq!(view[1], 3);
        assert_eq!(view[2], 2);
        assert_eq!(view[3], 1);
        assert_eq!(view[4], 0);

        let view = View::<_, _, _, 5>::new(&data, |i| (2 * i) % 5);
        assert_eq!(view[0], 0);
        assert_eq!(view[1], 2);
        assert_eq!(view[2], 4);
        assert_eq!(view[3], 1);
        assert_eq!(view[4], 3);
    }

    #[test]
    fn identity_view_mut() {
        let mut data = [0, 1, 2, 3, 4];
        let mut view = ViewMut::<_, _, _, 5>::ident_new(&mut data);
        view[3] += 10;
        assert_eq!(view[0], 0);
        assert_eq!(view[1], 1);
        assert_eq!(view[2], 2);
        assert_eq!(view[3], 13);
        assert_eq!(view[4], 4);
    }

    #[test]
    fn permuted_view_mut() {
        let mut data = [0, 1, 2, 3, 4];
        let mut view = ViewMut::<_, _, _, 5>::new(&mut data, |i| (4 - i));
        view[0] = 10;
        assert_eq!(view[0], 10);
        assert_eq!(view[1], 3);
        assert_eq!(view[2], 2);
        assert_eq!(view[3], 1);
        assert_eq!(view[4], 0);

        let mut view = ViewMut::<_, _, _, 5>::new(&mut data, |i| (2 * i) % 5);
        view[3] *= 7;
        assert_eq!(view[0], 0);
        assert_eq!(view[1], 2);
        assert_eq!(view[2], 10);
        assert_eq!(view[3], 7);
        assert_eq!(view[4], 3);
    }
}