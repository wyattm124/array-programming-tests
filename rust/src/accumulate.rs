use crate::view::View;
use crate::view::ViewMut;
use crate::view::Indexer;
use crate::complex::Num;

// Many Key components of efficient numerical algorithms come down to
//  a good compilation of these madd (multiply and add) loops
// NOTE : this is also a dot product
fn mac<'a, T: Num, F: Indexer, const N: usize>(
    a: View<'a, T, F, N>,
    b: View<'a, T, F, N>,
    c: &mut T
) -> () {
    for i in 0..N {
        *c += a[i] * b[i];
    }
}

fn correlate<'a, T: Num, F: Indexer, const M: usize, const N: usize, const O: usize>(
    a: View<'a, T, F, M>,
    b: View<'a, T, F, N>,
    c: &mut ViewMut<'a, T, F, O>,
) -> () {
    // Should have type system assert that M >= N
    // Should have type system assert that O == (M - N) + 1

    for i in 0..O {
        mac(a.subview::<N>(i), b, &mut c[i]);
    }
}

fn resample<'a, T: Num, F: Indexer, const M: usize, const N: usize, const O: usize>(
    a: View<'a, T, F, M>,
    b: View<'a, T, F, N>,
    c: &mut ViewMut<'a, T, F, O>,
    interp: usize,
    decim: usize,
) -> () {
    // Should have similar constraints to correlate, about accounting for interp
    //  and decim
    let a_p = a.compose(|i| i / interp);
    for i in (0..(M*interp)).step_by(decim) {
        mac(a_p.subview::<N>(i), b, &mut c[i]);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn simple_mac() {
        const FILTER_SIZE : usize = 5;
        const DATA_SIZE : usize = 16;
        let data : [i32; DATA_SIZE] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];
        let conv_start = View::<_, _, FILTER_SIZE>::ident_new(&data[0..FILTER_SIZE]);
        let filter = [0, 1, 2, 3, 4];
        let filter_flipper : fn(usize) -> usize = |i| (FILTER_SIZE-1)-i;
        let filter_view = View::<_, _, FILTER_SIZE>::new(&filter, filter_flipper);
        let corr_view = View::<_, _, FILTER_SIZE>::ident_new(&data[0..FILTER_SIZE]);
        let mut conv_ans : [i32; (DATA_SIZE - FILTER_SIZE) + 1] = [0; (DATA_SIZE - FILTER_SIZE) + 1];
       
        mac(conv_start, filter_view, &mut conv_ans[0]);
        assert_eq!(conv_ans[0], 10);
        conv_ans[0] = 0;
        mac(conv_start, corr_view, &mut conv_ans[0]);
        assert_eq!(conv_ans[0], 30);
    }

    #[test]
    fn simple_correlate() {
        const FILTER_SIZE : usize = 5;
        const DATA_SIZE : usize = 16;
        const ANS_SIZE : usize = (DATA_SIZE - FILTER_SIZE) + 1;
        let data : [i32; DATA_SIZE] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];
        let data_view = View::<_, _, DATA_SIZE>::ident_new(&data);
        let filter = [0, 1, 2, 3, 4];
        let filter_flipper : fn(usize) -> usize = |i| (FILTER_SIZE-1)-i;
        let conv_view = View::<_, _, FILTER_SIZE>::new(&filter, filter_flipper);
        let corr_view = View::<_, _, FILTER_SIZE>::ident_new(&filter);
        let mut corr_ans : [i32; ANS_SIZE] = [0; ANS_SIZE];
        let mut corr_ans_view: ViewMut<'_, _, _, ANS_SIZE> = ViewMut::ident_new(&mut corr_ans);
    
        // Try a simple correlation with filter
        correlate(data_view, corr_view, &mut corr_ans_view);
        for i in 0..ANS_SIZE {
            assert_eq!(corr_ans_view[i], (i as i32) * 10 + 30);
        }

        // Reset ans array
        for i in 0..ANS_SIZE {
            corr_ans_view[i] = 0;
        }

        // Try a simple convolution with filter flipped
        correlate(data_view, conv_view, &mut corr_ans_view);
        for i in 0..ANS_SIZE {
            assert_eq!(corr_ans_view[i], (i as i32) * 10 + 10);
        }
    }
}