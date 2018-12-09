An IFMA implementation of finite field arithmetic for the Mersenne field and a quadratic extension.

Eventually I plan to grow this into a FourQ implementation.

# TODO:

- split prime field type to handle reduction state (allow deferring reductions)
- add "compressed" field elements (fewer memory accesses for table lookups)
- rework extension field arithmetic once prime field can defer reductions
- benchmarks, etc.
