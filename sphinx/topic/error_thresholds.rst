Error thresholds in RepLAB
==========================

We list here the uses of tolerances in RepLAB.

Detection of real representation type
-------------------------------------

In: `+replab.+irreducible.computeRealType`

To detect the real representation type, we take implicitly a Hermitian sample from the commutant of the complexification of the real representation. In practice, we sample a generic, non-symmetric sample from the commutant of the real representation, split it into symmetric and antisymmetric parts, and construct a complex Hermitian matrix $C$ from those parts.

We then sample a normally distributed complex vector $\vec{v}$ and compute the rank of the span of $\{ \vec{v}, C \vec{v}, C^2 \vec{v} \}$ (we use three vectors as a safety mechanism: the rank should always be one or two). The rank check uses the ``rank`` Matlab function, with a tolerance we prescribe (``doubleEigTol``).

However, we'd like to use a better bound here.

Suggested steps:

- Can we replace the ``rank`` check with Gram-Schmidt elimination? It's then easier to understand the meaning of the tolerance parameter.

- The distribution of the noise due to numerical error when the rank is one is probably difficult to estimate, as it depends on the quality of the basis of the subrepresentation, and we do not know yet how to quantify that part.

- The distribution of the residual norm when the rank if two is probably easier to compute.
