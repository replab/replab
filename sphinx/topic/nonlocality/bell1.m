% # Symmetries in Bell nonlocality, part 1
%
% For a compact introduction to representation theory, see the [book of Jean-Pierre Serre](https://link.springer.com/book/10.1007/978-1-4684-9458-7).
%
% Follow the [installation instructions](https://replab.github.io/replab/docs/installation.html).
%
% Then add the RepLAB path. The path below is the one for the machine that produced this notebook! Replace by your own. Or, better, move to the root RepLAB folder and just run the `replab_init` script.

addpath([pwd, '/../../..']);
replab_init('verbose', 0);

% ## Symmetries of the outcomes of a measurement
%
% We start by studying the symmetries of a single measurement on a quantum system. Let us say that the measurement has outcomes $a = 1, 2, \ldots, k$.
%
% ### Symmetry group
%
% We are interested in relabelings of these outcomes; those relabelings are represented by permutations. The group of all permutations of $k$ elements is the symmetric group of degree $k$, written $S_k$.
%
% For the motivation, see our [paper on the symmetries of nonlocality in the context of device-independence](https://arxiv.org/abs/1610.01833). See also the [description of the symmetric group on Wikipedia](https://en.wikipedia.org/wiki/Symmetric_group).
%
% In particular, there are different ways of writing elements of a symmetric group. RepLAB chooses to write down the permutation as the image row vector. For $\pi \in S_5$, we write $\pi = [\pi(1), \pi(2), \pi(3), \pi(4), \pi(5)]$ (using square brackets to distinguish this notation from the cycle notation), which corresponds to the *one-line notation* of https://en.wikipedia.org/wiki/Permutation#Notations .
%
% For the group axioms, Wikipedia has a [compact discussion](https://en.wikipedia.org/wiki/Group_(mathematics)) . When composing permutations, one needs to decide, for $\pi, \sigma \in S_k$, how to compute the product $\pi \cdot \sigma$. We make the choice of having $(\pi \cdot \sigma)(a) = \pi ( \sigma (a))$, which corresponds to a left action, i.e. $\sigma$ is applied first. See https://en.wikipedia.org/wiki/Permutation#Composition_of_permutations for details.
%
% This choice of composition corresponds to standard notation used in physics; however computational group theory uses a right action convention, most of the time using the *exponent notation*. Those notation choices should not concern the RepLAB user, as we made sure that the software uses consistent notation.
%
% Below, we construct the symmetric group $S_5$. The object `S5` we construct knows how to compose permutations, when they are written using row vectors of integers in Matlab/Octave. The composition example is the one from [Wikipedia](https://en.wikipedia.org/wiki/Symmetric_group#Multiplication), and we check that we recover the same results.

S5 = replab.S(5);
f = [3 2 1 5 4]
g = [2 5 4 3 1]
fg = S5.compose(f, g)

% Note that permutations carry the size of their domain with them: the permutation $[2, 3, 1, 4]$ is not the same as the permutation $[2, 3, 1]$ even if both have the same action on $1,2,3$ and the first one leaves $4$ invariant.

% **Exercice**: for fun, try to run `S5.order`, `S5.inverse(f)`, `S5.elements`, `S5.sample`. Construct `S50 = replab.S(50)`. Compute `S50.order`, and try `S50.elements`.

% ### Representation
%
% The distribution of outcomes of the measurement above can be represented by the distribution $P_\text{A}(a)$. We can also represent the distribution as a vector $\vec{P}_\text{A} \in V = \mathbb{R}^5$, with the following enumeration of coefficients: $\vec{P}_\text{A}^\top = \left ( P_\text{A}(1), P_\text{A}(2),  P_\text{A}(3),  P_\text{A}(4) ,  P_\text{A}(5) \right )$ . This action of elements of $S_5$ on $V = \mathbb{R}^5$ is encoded as a [group representation](https://en.wikipedia.org/wiki/Representation_theory#Definition); in RepLAB we call it the *defining representation* of $S_5$.

rep5 = S5.naturalRep
rep5.image(g)
rep5.image(f)*rep5.image(g) - rep5.image(S5.compose(f, g)) % check the representation composition axiom

% ### Invariant subspaces and subrepresentations
%
% Let $G$ be a group, and $\rho: G \to GL(V)$ a representation of $G$ on a vector space $V$. The space $GL(V)$ designates invertible linear maps from $V$ to $V$. In RepLAB, the vector space $V$ is defined either over the real or complex numbers.
%
% An invariant subspace of $V$ is a subspace $W \subseteq V$ that stays invariant under the action of $\rho$. Formally, we want for all $\vec{w} \in W$ and $g \in G$ that $\rho_g ~ \vec{w} \in W$.
%
% Any representation $V$ has two subrepresentations: the empty subspace $0$ and $V$ itself.
%
% For any permutation representation like `rep5` above, the subspace spanned by the vector of all ones $(1,1,1,1,1)^\top$ is a subrepresentation: any vector of the form $(x,x,x,x,x)^\top\in\mathbb{R}^5$ stays invariant under permutation of coefficients.
%
% The representations we construct are unitary, that is $\rho_g^+ = \rho_g^{-1} = \rho_{g^{-1}}$. Most of the literature assumes unitarity "without loss of generality", as a finite dimensional nonunitary representation for a compact group can always be made unitary after a change of basis.
%
% With the other assumptions we make, this means that the orthogonal complement to $(1,1,1,1,1)^\top$ is also an invariant subspace.
%
% By restricting a representation $\rho$ to one of its invariant subspaces, one then defines a subrepresentation.
%
% A representation that has no invariant subspaces except $0$ and the full space $V$ is *irreducible*. The decomposition of a representation into irreducible representations can be obtained by recursively decomposing subrepresentations until the components are irreducible; obtaining that decomposition is pretty tricky, but luckily RepLAB automates that.
%
% We obtain this decomposition and consider the first irreducible subrepresentation: the vector of all ones. Note that RepLAB recovers a "nice" algebraic expression for the basis of that subspace if we apply the `nice` method on the decomposition; this works in a few cases, the rest of the time that basis will be numerical/approximate to floating-point precision. We also check how the representation looks in that invariant subspace.

% The floating point approximate decomposition
D5 = rep5.decomposition

% Recovering nice coefficients for the basis
N5 = D5.nice

% Images in the subrepresentations
subrep5a = D5.irrep(1) % this is the first irreducible component
subrep5a.image(S5.generator(1))
subrep5a.image(S5.generator(2))

% We note that both generators of $S_5$ look like the identity matrix of dimension $1$ in that subrepresentation. Thus, any element of $S_5$ corresponds to the identity matrix for that subrepresentation. Such a representation with an identity image is the *trivial representation*.
%
% We move now to the other subrepresentation.

subrep5b = D5.irrep(2)
subrep5b.image(S5.generator(1))
subrep5b.image(S5.generator(2))

% Now the generators of $S_5$ look like they are doing something in that subrepresentation! That representation is actually the [standard representation](https://groupprops.subwiki.org/wiki/Standard_representation) of $S_5$.
%
% ### Physical interpretation
%
% The physical interpretation of invariant subspaces is that they correspond to properties that are preserved under symmetry.
%
% In the present case, the split into the two subrepresentations highlight a property of probability distributions invariant under relabeling of outcomes, the *normalization* of the distribution.
%
% ### The original representation in the basis that exhibits the decomposition
%
% The matrix $B$ is the change of basis that block diagonalizes the representation. Here, we can show that `rep5` has a 1x1 and 4x4 block structure in the relevant basis.

B = D5.basis
inv(B) * rep5.image(S5.generator(1)) * B
inv(B) * rep5.image(S5.generator(2)) * B

% ## Another example
%
% To continue this discussion, we consider a measurement with two outcomes $a=1,2$. The symmetry group is now $S_2$, and has two elements: identity $[1,2]$ and flip $[2,1]$.

S2 = replab.S(2)
S2.elements

% The defining representation has two images: the identity $\begin{pmatrix}1 & 0\\0 & 1\end{pmatrix}$ and the flip $\begin{pmatrix}0 & 1\\1 & 0 \end{pmatrix}$.

rep2 = S2.naturalRep
rep2.image([1,2])
rep2.image([2,1])

% Now, one can check that the two subrepresentations correspond to vector spaces spanned by $(1,1)^\top$ and $(1,-1)^\top$ respectively.

D2 = rep2.decomposition.nice
D2.irrep(1).basis
D2.irrep(2).basis

% Up to a multiplicative factor, these two vectors extract two quantities: $(1,1) \cdot \vec{P}$ computes the normalization, while $(1,-1) \cdot \vec{P}$ computes the correlator $P(1) - P(2)$ often used in quantum information for systems with binary outcomes.
%
% Note that the flip operation leaves $(1,1)$ invariant, while it changes the sign of $(1,-1)$. Accordingly, the subrepresentation in the subspace spanned by $(1,-1)$ is called the *sign representation*.

% ## Conditional probability distributions
%
% We now move to the study of conditional probability distributions $P_{A|X}(a|x)$. For simplicity, we consider binary outputs $a=1,2$ and binary inputs $x=1,2$. We enumerate the coefficients in a vector $\vec{P} = (P(1|1), P(2|1), P(1|2), P(2|2))^\top$.
%
% Different operations are available:
%
% - Permutation of inputs, which will permute coefficients of $\vec{P}$ according to $g_1 = [3,4,1,2]$.
%
% - Permutation of outputs of the first measurement, which will permute $\vec{P}$ according to $g_2 = [2,1,3,4]$.
%
% - The permutation of outputs of the second measurement is achieved by conjugating $g_2$ by $g_1$, i.e. by $g_1 ~ g_2 ~ g_1$, so we do not need to provide that group generator explicitly.
%
% We then construct the group $G=\left <g_1, g_2 \right>$ and its defining representation, which we decompose into irreducible subrepresentations.

g1 = [3,4,1,2];
g2 = [2,1,3,4];
generators = {g1 g2}; % arrays of things in Matlab are cell arrays
S4 = replab.S(4);
G = S4.subgroup(generators) % construction as a subgroup of a generic parent group
G.elements
Grep = G.naturalRep

% Let's decompose this representation into irreducible subrepresentations
Gdec = Grep.decomposition.nice
I1 = Gdec.irrep(1)
I2 = Gdec.irrep(2)
I3 = Gdec.irrep(3)

% What are the irreps we obtain? We have:
%
% - A first subspace encoding the overall normalization $(1,1,1,1) \cdot \vec{P} = 2$, by definition of $\sum_{a} P(a|x) = 1$.
%
% - A second subspace encoding the difference in normalization between the condition $x=1$ and the condition $x=2$: $(1,1,-1,-1) \cdot \vec{P} = 0$.
%
% - The third subspace, after a change of basis vectors, corresponds to the binary correlators $(1,-1,0,0) \cdot \vec{P} = P(1|1) - P(2|1)$ and $(0,0,1,-1) \cdot \vec{P} = P(1|2) - P(2|2)$.
%
% We leave, as an exercice, the construction of the symmetry group of the CHSH scenario using its action on the joint conditional distribution $P_{AB|XY}(ab|xy)$ -- one should recover the results of [our paper](https://arxiv.org/abs/1610.01833).

% ## A final word of warning
%
% The presentation above glossed over an important concept: the decomposition of a representation into irreducible subrepresentations is essentially unique, up to the degeneracy afforded by multiplicities. Indeed, when several copies of the same representation appear, the decomposition into irreducible subrepresentation is not exactly unique; refer to the end of Chapter 2 of the book by Jean-Pierre Serre cited in the preface.
%
% In particular, the method `dec.irrep(i)` needs to be replaced by `dec.irrep(i,j)` where $i$ is the index of the irreducible representation type, and $j$ is the copy index.
%
% Also, if one works using real representations, there could be a difference in the irreducible decomposition when the representation is treated as a complex one. This is because the real field is not algebraically closed. Over the reals, there are three kinds of irreducible representations: RepLAB identifies and puts them in a canonical representation. Not to worry: representation of symmetric groups and derivatives (direct products and the like) are decomposed the same way over the real and complex numbers.
