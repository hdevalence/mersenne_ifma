Vectorized arithmetic for the Mersenne field and its quadratic extension.

Implementation notes follow:

# Choice of radix

We want to use a binary radix \\(2^r\\) for some \\(r\\).  For
IFMA we need the limbs to fit within 52 bits, so we could choose:

* \\(r = 52\\) (saturated)
* \\(r = 51\\) (unsaturated, one carry bit)
* a much smaller radix.

We need at least 127 bits, so we need at least 3 limbs. Choosing
\\(r = 52\\) is not ideal because saturated arithmetic forces
sequential dependencies between instructions, making ILP more
difficult.  Choosing \\(r = 51\\) means that the limb boundaries
are not aligned with the bitsize of the prime, which means that
reducing high product terms crosses the boundaries of the lower
limbs.  For instance, using \\( r = 51 \\), we would compute a
52-bit product term \\(z_3 2^{153} = z_3 2^{127} 2^{26} = z_3
2^{26} \\). This term needs to be placed at position \\(2^{26}\\),
crossing a limb boundary.  The problem here is that weight of the
value is not distributed evenly across the limbs when using \\(r =
51\\).

The third option is to choose a much smaller radix, \\( r = 43 \\), so
that \\( 3r = 129 \\).  This means the weight is spread evenly across
the limbs and the limb boundaries are more closely aligned with the
bitsize of the prime.  With this choice, source limbs \\(x_i, y_j\\)
generate the product term
\\[
\begin{aligned}
x_i y_j 2^{43(i+j)} &= (\mathrm{lo}(x_i, y_j) + 2^{52}\mathrm{hi}(x_i, y_j)) 2^{43(i+j)} \\\\
&= \mathrm{lo}(x_i, y_j) 2^{43(i+j)} + 2^9 \mathrm{hi}(x_i, y_j) 2^{43(i+j+1)}.
\end{aligned}
\\]
Because IFMA accumulators have 12 bits of headroom, the high term
still fits into an accumulator, even if the \\(x_i, y_j\\) are
full-size.  If they are smaller (closer to \\(43\\) than \\(52\\)
bits), then there's even more room.  (More notes on tracing this
follows).

# Prime Field Multiplication

The inputs to a multiplication are two field elements
\\[
\begin{aligned}
x &= x_0 + x_1 2^{43} + x_2 2^{86} \\\\
y &= y_0 + y_1 2^{43} + y_2 2^{86},
\end{aligned}
\\]

If the \\(x_i, y_i\\) are used as inputs to IFMA operations, we have
\\( 9 = 52 -43 \\) bits of headroom for multiplications.  For now,
parameterize the bound on the input limbs by the _bit-excess_ \\(b\\),
so that \\(x_i, y_i < 2^{43+b}\\), and leave this variable free and
subject to further constraints.  The requirement that the limbs can be
used an inputs to IFMA operations forces \\(b \leq 9\\).

A schoolbook multiplication in product scanning form takes the form
\\[
\begin{aligned}
z_0 &= x_0 y_0 \\\\
z_1 &= x_1 y_0 + x_0 y_1 \\\\
z_2 &= x_2 y_0 + x_1 y_1 + x_0 y_2 \\\\
z_3 &=           x_2 y_1 + x_1 y_2 \\\\
z_4 &=                     x_2 y_2 \\\\
z_5 &= 0 \\\\
\end{aligned}
\\]
Writing 
\\[
x_i y_j 2^{43(i+j)} = \mathrm{lo}(x_i, y_j) 2^{43(i+j)} + 2^9 \mathrm{hi}(x_i, y_j) 2^{43(i+j+1)}
\\]
and substituting into the schoolbook multiplication, we get
\\[
\begin{aligned}
z_0 &= \mathrm{lo}(x_0, y_0) \\\\
&\qquad + 0 \\\\
z_1 &= \mathrm{lo}(x_1, y_0) + \mathrm{lo}(x_0, y_1) \\\\
&\qquad + 2^{9} \mathrm{hi}(x_0, y_0) \\\\
z_2 &= \mathrm{lo}(x_2, y_0) + \mathrm{lo}(x_1, y_1) + \mathrm{lo}(x_0, y_2) \\\\
&\qquad + 2^{9} \mathrm{hi}(x_1, y_0) + 2^{9} \mathrm{hi}(x_0, y_1) \\\\
z_3 &=                         \mathrm{lo}(x_2, y_1) + \mathrm{lo}(x_1, y_2) \\\\
&\qquad + 2^{9} \mathrm{hi}(x_2, y_0) + 2^{9} \mathrm{hi}(x_1, y_1) + 2^{9} \mathrm{hi}(x_0, y_2) \\\\
z_4 &=                                                 \mathrm{lo}(x_2, y_2) \\\\
&\qquad +                               2^{9} \mathrm{hi}(x_2, y_1) + 2^{9} \mathrm{hi}(x_1, y_2) \\\\
z_5 &= 0 \\\\
&\qquad +                                                             2^{9} \mathrm{hi}(x_2, y_2) \\\\
\end{aligned}
\\]

Since \\( p = 2^{127} - 1 \\), \\( 2^{3\cdot 43} = 2^{129} = 2^2 \pmod
p\\), so the high limbs can be reduced onto the low limbs by
multiplying by \\(4\\).  Combining terms, we get

\\[
\begin{aligned}
z_0 &= 2^0 \mathrm{lo}(x_0, y_0) + 2^2 \mathrm{lo}(x_2, y_1) + 2^2 \mathrm{lo}(x_1, y_2)
     + 2^{11} \mathrm{hi}(x_2, y_0) + 2^{11} \mathrm{hi}(x_1, y_1) + 2^{11} \mathrm{hi}(x_0, y_2) \\\\
z_1 &= 2^0 \mathrm{lo}(x_1, y_0) + 2^0 \mathrm{lo}(x_0, y_1) + 2^2 \mathrm{lo}(x_2, y_2)
     + 2^{9} \mathrm{hi}(x_0, y_0) + 2^{9} \mathrm{hi}(x_2, y_1) + 2^{9} \mathrm{hi}(x_1, y_2) \\\\
z_2 &= 2^0 \mathrm{lo}(x_2, y_0) + 2^0 \mathrm{lo}(x_1, y_1) + 2^0 \mathrm{lo}(x_0, y_2)
     + 2^{9} \mathrm{hi}(x_1, y_0) + 2^{9} \mathrm{hi}(x_0, y_1) + 2^{11} \mathrm{hi}(x_2, y_2) \\\\
\end{aligned}
\\]

Each of the \\(z_i\\) has six terms, which can be grouped into three
pairs, two with common coefficient \\(2^k\\) and one with different
coefficients.  The common-coefficient pairs can be computed by
accumulating into the same temporary register, and the
different-coefficient pair can do the higher term first, then shift
the result.  This gives \\(3 \times 3 = 9\\) independent chains of
IFMA operations.

What are the sizes of the product terms?  If \\(x_i, y_j <
2^{43+b}\\), then \\(\mathrm{lo}(x_i, y_j) < 2^{52}\\) and
\\(\mathrm{hi}(x_i, y_j) < 2^{(86 + 2b) - 52} = 2^{34 + 2b}\\).
Substituting these bounds gives

\\[
\begin{aligned}
z_0 &< 2^0 2^{52} + 2^2 2^{52} + 2^2 2^{52}
     + 2^{11} 2^{34+2b} + 2^{11} 2^{34+2b} + 2^{11} 2^{34+2b} \\\\
z_1 &< 2^0 2^{52} + 2^0 2^{52} + 2^2 2^{52}
     + 2^{9} 2^{34+2b} + 2^{9} 2^{34+2b} + 2^{9} 2^{34+2b} \\\\
z_2 &< 2^0 2^{52} + 2^0 2^{52} + 2^0 2^{52}
     + 2^{9} 2^{34+2b} + 2^{9} 2^{34+2b} + 2^{11} 2^{34+2b} \\\\
\end{aligned}
\\]

The largest is \\(z_0 < 2^{54.33} + 2^{46.59 + 2b}\\).  To ensure the
\\(z_i \\) fit in \\(64\\) bits, we can restrict \\(b < 8.5\\) so that
the \\(z_i < 2^{63.6}\\).

