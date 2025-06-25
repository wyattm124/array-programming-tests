use std::ops;

// TODO : should compile time verify that we do not generate an index out of bounds error
// NOTE : indexer should be a FAST, simple, indexing function. 95% of the time it should be
//  a linear transformation.

trait Indexer: Fn(usize) -> usize + Copy + Clone {}
impl<F: Fn(usize) -> usize + Copy + Clone> Indexer for F {}

#[derive(Debug, Clone, Copy)]
pub struct View<'a, T, F: Indexer> {
    data: &'a [T],
    indexer: F,
    domain: usize,
}

#[derive(Debug)]
pub struct ViewMut<'a, T, F: Indexer> {
    data: &'a mut [T],
    indexer: F,
    domain: usize,
}

impl<'a, T, F: Indexer> View<'a, T, F> {
    pub fn new(data: &'a [T], indexer: F, domain: usize) -> Self {
        Self { data, indexer, domain }
    }
}

impl<'a, T, F: Indexer> ViewMut<'a, T, F> {
    pub fn new(data: &'a mut [T], indexer: F, domain: usize) -> Self {
        Self { data, indexer, domain }
    }

    pub fn as_view(&self) -> View<'_, T, F> {
        View::new(self.data, self.indexer.clone(), self.domain)
    }
}

impl<'a, T, F: Indexer> ops::Index<usize> for View<'a, T, F> {
    type Output = T;
    fn index(&self, index: usize) -> &Self::Output {
        &self.data[(self.indexer)(index)]
    }
}

impl<'a, T, F: Indexer> ops::Index<usize> for ViewMut<'a, T, F> {
    type Output = T;
    fn index(&self, index: usize) -> &Self::Output {
        &self.data[(self.indexer)(index)]
    }
}

impl<'a, T, F: Indexer> ops::IndexMut<usize> for ViewMut<'a, T, F> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.data[(self.indexer)(index)]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identity_view() {
        let data = vec![0, 1, 2, 3, 4];
        let view = View::new(&data, |i| i, data.len());
        assert_eq!(view[0], 0);
        assert_eq!(view[1], 1);
        assert_eq!(view[2], 2);
        assert_eq!(view[3], 3);
        assert_eq!(view[4], 4);
    }

    #[test]
    fn permuted_view() {
        let data = vec![0, 1, 2, 3, 4];
        let view = View::new(&data, |i| (4 - i), data.len());
        assert_eq!(view[0], 4);
        assert_eq!(view[1], 3);
        assert_eq!(view[2], 2);
        assert_eq!(view[3], 1);
        assert_eq!(view[4], 0);

        let view = View::new(&data, |i| (2 * i) % 5, data.len());
        assert_eq!(view[0], 0);
        assert_eq!(view[1], 2);
        assert_eq!(view[2], 4);
        assert_eq!(view[3], 1);
        assert_eq!(view[4], 3);
    }

    #[test]
    fn identity_view_mut() {
        let mut data = vec![0, 1, 2, 3, 4];
        let len = data.len();
        let mut view = ViewMut::new(&mut data, |i| i, len);
        view[3] += 10;
        assert_eq!(view[0], 0);
        assert_eq!(view[1], 1);
        assert_eq!(view[2], 2);
        assert_eq!(view[3], 13);
        assert_eq!(view[4], 4);
    }

    #[test]
    fn permuted_view_mut() {
        let mut data = vec![0, 1, 2, 3, 4];
        let len = data.len();
        let mut view = ViewMut::new(&mut data, |i| (len - 1) - i, len);
        view[0] = 10;
        assert_eq!(view[0], 10);
        assert_eq!(view[1], 3);
        assert_eq!(view[2], 2);
        assert_eq!(view[3], 1);
        assert_eq!(view[4], 0);

        let mut view = ViewMut::new(&mut data, |i| (2 * i) % len, len);
        view[3] *= 7;
        assert_eq!(view[0], 0);
        assert_eq!(view[1], 2);
        assert_eq!(view[2], 10);
        assert_eq!(view[3], 7);
        assert_eq!(view[4], 3);
    }
}

