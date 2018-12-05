#![feature(test)]

#[macro_use]
extern crate criterion;
//use criterion::black_box;
use criterion::Criterion;

extern crate test;
use self::test::black_box;

extern crate mersenne_ifma;
use mersenne_ifma::serial::*;

fn prime_field_add(c: &mut Criterion) {
    c.bench_function("F127 addition", |b| {
        let x = black_box(F127::from(2u128));
        let y = black_box(F127::from(1u128));

        b.iter(|| black_box(x + y));
    });
}

fn prime_field_sub(c: &mut Criterion) {
    c.bench_function("F127 subtraction", |b| {
        let x = F127::from(2u128);
        let y = F127::from(1u128);

        b.iter(|| black_box(black_box(x) - black_box(y)));
    });
}

fn prime_field_mul(c: &mut Criterion) {
    c.bench_function("F127 multiplication", |b| {
        let x = F127::from(2u128);
        let y = F127::from(1u128);

        b.iter(|| black_box(black_box(x) * black_box(y)));
    });
}

fn ext_field_add(c: &mut Criterion) {
    c.bench_function("F127Ext addition", |b| {
        let x = F127Ext::from((2u128, 9u128));
        let y = F127Ext::from((1u128, 8u128));

        b.iter(|| black_box(black_box(x) + black_box(y)));
    });
}

fn ext_field_sub(c: &mut Criterion) {
    c.bench_function("F127Ext subtraction", |b| {
        let x = F127Ext::from((2u128, 9u128));
        let y = F127Ext::from((1u128, 8u128));

        b.iter(|| black_box(black_box(x) - black_box(y)));
    });
}

fn ext_field_mul(c: &mut Criterion) {
    c.bench_function("F127Ext multiplication", |b| {
        let x = F127Ext::from((2u128, 9u128));
        let y = F127Ext::from((1u128, 8u128));

        b.iter(|| black_box(black_box(x) * black_box(y)));
    });
}

criterion_group!{
    name = prime_benches;
    config = Criterion::default();
    targets =
    prime_field_add,
    prime_field_sub,
    prime_field_mul,
}

criterion_group!{
    name = ext_benches;
    config = Criterion::default();
    targets =
    ext_field_add,
    ext_field_sub,
    ext_field_mul,
}

criterion_main!{
    prime_benches,
    ext_benches,
}
