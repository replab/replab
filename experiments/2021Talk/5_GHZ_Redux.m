% # States with the GHZ symmetries
%
% Initialize the RepLAB toolbox (be in the /replab directory or use `run path/replab/replab_init.m`)

run ~/w/replab/replab_init
replab.globals.verbosity(0);
replab.globals.useReconstruction(1);

% In this notebook, we compute the family of states with the symmetry of the GHZ states, as in [C. Eltschka and J. Siewert, "Entanglement of Three-Qubit Greenberger-Horne-Zeilinger--Symmetric States"](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.108.020502).
%
%
% The GHZ state $\left|\text{GHZ}\right> = \left|000\right> + \left|111\right>$ is invariant (U |GHZ> = |GHZ>) under a family of matrices U
%
% \begin{equation}
% U = \begin{pmatrix} a_0 & 0 \\ 0 & a_1 \end{pmatrix} \otimes \begin{pmatrix} b_0 & 0 \\ 0 & b_1 \end{pmatrix} \otimes \begin{pmatrix} c_0 & 0 \\ 0 & c_1 \end{pmatrix} 
% \end{equation}
%
% who act on the state space $\mathbb{C}^8$. The eight coefficients, in order, are those of $\left|000\right>$,$\left|001\right>$,$\left|010\right>$,$\left|011\right>$,$\left|100\right>$,$\left|101\right>$,$\left|110\right>$,$\left|111\right>$.
%
% For the invariance to hold, the unit complex numbers $a_0$,$a_1$,$b_0$,$b_1$,$c_0$,$c_1$ obey the equation $a_0 b_0 c_0 = a_1 b_1 c_1 = 1$.
%
% The state is also invariant under permutation of subsystems (a group of order $3! = 6$), and permutation of the levels (a group of order $2$).
%
% Let us define the continuous connected group. $T(6)$ is the [torus group](https://en.wikipedia.org/wiki/Torus#n-dimensional_torus) with 6 elements, elements that we name according to our scenario.

T6 = replab.T(6);
T6 = T6.withNames({'a0' 'b0' 'c0' 'a1' 'b1' 'c1'});

% We then define the subgroup obeying the equation $a_0 b_0 c_0 = a_1 b_1 c_1 = 1$.

T = T6.subgroupWith('a0*b0*c0 = 1', 'a1*b1*c1 = 1')

% Now, how does that group act on the state space $\mathbb{C}^8$?
%
% We construct the representation whose image is the matrix $U$ above.

Trep = T.diagonalRepWith('a0 b0 c0', ...
                         'a0 b0 c1', ...
                         'a0 b1 c0', ...
                         'a0 b1 c1', ...
                         'a1 b0 c0', ...
                         'a1 b0 c1', ...
                         'a1 b1 c0', ...
                         'a1 b1 c1');

% We construct now the discrete part, by writing how the subsystem and level permutations affect the elements of the continuous connected part.
%
% This finite group permutes the three subsystems and the two levels, independently, so we write it as a direct product.

S3 = replab.S(3);
S2 = replab.S(2);
F = S3.directProduct(S2)

% Now, we write the action of the generators of this discrete group on the torus elements:

% Permutation of AB
gAB = {[2 1 3] [1 2]};
actAB = T.automorphism('b0', 'a0', 'c0', 'b1', 'a1', 'c1');
% Permutation of AC
gAC = {[3 2 1] [1 2]};
actAC = T.automorphism('c0', 'b0', 'a0', 'c1', 'b1', 'a1');
% Permutation of BC
gBC = {[1 3 2] [1 2]};
actBC = T.automorphism('a0', 'c0', 'b0', 'a1', 'c1', 'b1');
% Permutation of the two levels
gL = {[1 2 3] [2 1]};
actL = T.automorphism('a1', 'b1', 'c1', 'a0', 'b0', 'c0');

% We check that the generators generate the whole group.

assert(F.subgroup({gAB, gAC, gBC, gL}) == F);

% The [outer semidirect product construction](https://en.wikipedia.org/wiki/Semidirect_product#Inner_and_outer_semidirect_products) employs a morphism from a group $H$ to the automorphisms of another gorup $N$. This is what we do now.

G = T.semidirectProductByFiniteGroup(F, 'preimages', {gAB, gAC, gBC, gL}, 'images', {actAB, actAC, actBC, actL});

% ## How does this group act on the state space?
% Remember that the integers $1$, $2$, $3$, $4$, $5$, $6$, $7$, $8$ enumerate $\left|000\right>$, $\left|001\right>$, $\left|010\right>$, $\left|011\right>$, $\left|100\right>$, $\left|101\right>$, $\left|110\right>$, $\left|111\right>$.

imgAB = replab.Permutation.toMatrix([1 2 5 6 3 4 7 8]);
imgAC = replab.Permutation.toMatrix([1 5 3 7 2 6 4 8]);
imgBC = replab.Permutation.toMatrix([1 3 2 4 5 7 6 8]);
imgL = replab.Permutation.toMatrix([8 7 6 5 4 3 2 1]);

% We now construct the representation of the discrete group on the state space.

Frep = F.repByImages('C', 8, 'preimages', {gAB, gAC, gBC, gL}, 'images', {imgAB, imgAC, imgBC, imgL});

% We had the representation of the torus group `Trep`, and we now assemble the two.

rep = G.semidirectProductRep(Frep, Trep)

% ... and now it is time to decompose that decomposition.

rep.decomposition

% Thus we recover that the states invariant under this type of symmetries are the pure states $\left| \text{GHZ}_\pm \right> = \frac{\left| 000 \right> \pm \left| 111 \right>}{\sqrt{2}}$ and a mixed state containing all basis elements except $\left| 000 \right>$ and $\left| 111 \right>$ with equal weight $1/6$.


