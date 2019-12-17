%% Representations of the symmetric group S(n)
%
% The symmetric group $S(n)$ is one of the most fundamental group. This
% document presents various ways of defining representations of the
% symmetric group $S(n)$ with *RepLAB*.

%%
% Before trying any of the *RepLAB* commands, we must initialize the library:
replab_init

%% The Symmetric group
% Consider $n$ distinct elements, numbered $1,2,\ldots,n$. The symmetric
% group describes all possible ways in which these elements can be
% permuted.
% 
% With $n=3$, for instance, two possible arrangements are
[1 2 3]
%%
% and
[1 3 2]
%%
% In general there are $n!$ such arrangements (i.e. $6$ arrangements for
% $n=3$).
%
% In *RepLAB*, the symmetric group can be constructed with a simple command
n = 3;
S3 = replab.S(n);
%%
% The elements of the group can be easily listed
S3.elements
%%
% The number of elements is the order of the group
S3.order
%%
% As we see, permutations are represented by row-vectors with an
% arrangement of the elements 1, ..., n.
%
% One can use the group structure to compose elements of the group. For
% instance, permuting the first two elements in 1,2,3, followed by a
% permutaiton of the second and third elements results in a cyclic
% permutation:
element1 = [2 1 3];
element2 = [1 3 2];
S3.compose(element1, element2)

%% Representations of the symmetric group
% A finite group can always be represented by matrices of finite dimension.
% In this case, a d x d matrix is associated to each group element, called
% the image of the group element. The group composition then amounts to a
% matrix multiplication of these matrices.
%
% A group can have several matrix representations (with various dimension
% d). Here, we are going to explore some representations of our group S(3).
%
% When permuting n elements, a natural choice of representation consists in
% associating one dimension to each element, hence choosing d=n. Matrices
% then act on $\mathbb{R}^n$. This is called the *defining represention* of
% a group. It can be constructed easily:
defRep = S3.definingRep
%%
% The image of a group element in this representation can be extracted
defRep.image([2 1 3])
%%
% We see that the matrices here are permutation matrices which simply apply
% the permutation to the coordinates.
%
% Let us check that the representation is faithful, i.e. that is follows
% indeed the group structure
image1 = defRep.image(element1);
image2 = defRep.image(element2);
image1*image2
defRep.image(S3.compose(element1, element2))
%%
% The two images are indeed 

%
% The action of a group can always be represented by the product of d x d
% matrices
%
% When representing elements of a group as d x d matrices, the
% action of the group amounts to matrix multiplication. The 

% Elements of a group can always be represented by 
% The action of a group can be represented by the matrix product acting o.

%%
% Remark that permutations in *RepLAB* are simply row vectors of
% images. We use the convention that permutations act on indices
% $i=1,2,3,4$ *on the left*, thus the image of 3 under $g_2$
% is compatible with MATLAB syntax:
g2(3)
%%
% We construct the group by first accessing the group of
% permutations on 4 elements, then constructing a subgroup
S4 = replab.Permutations(4).subgroup({g1 g2})
%%
% noting that the generators of the group are named $a$, $b$,...
%
% We can perform a few operations on finite groups. Let $g$ be a
% random permutation, which we can obtain in (mostly) equivalent ways by either
g = randperm(n);
g = S4.sample;
g = S4.sampleUniformly;
%
% We can also compute the group order or enumerate the elements of
% the group
S4.order
%%
%
S4.elements
%%
% a technique that works even for big groups.
S30 = replab.Permutations(30).subgroup({[2:30 1] [2 1 3:30]})
S30.order
%%
%
S30.elements


%% Group representations
% The defining representation of $S_4$ simply permutes the
% coordinates of $R^4$:
rho = S4.definingRep
%%
% We take two elements of $S_4$:
g = [2 3 1 4]
h = [2 1 4 3]
%%
% and their composition
gh = g(h)
%%
% and verify that we have a representation of the group
rho.image(g) * rho.image(h)
rho.image(gh)
%%
% We can also define representations from their images. Here,
% we define the sign representation (using a permutation
% representation of it!). 
dim = 2;
isUnitary = true;
rho1 = S4.repByImages('R', dim, isUnitary, {[0 1; 1 0] [0 1; 1 0]})
rho2 = S4.permutationRep(2, {[2 1] [2 1]})
%%
%
rho1.image(g)
rho2.image(g)

%% Decomposing group representations
% *RepLAB* provides the irreducible decomposition of representations
% over the real numbers, identifying the representation type
% (real, complex or quaternionic):
I = rho.decomposition
%%
% We can get isotypic components and the copies of irreducible
% representations contained inside
I.component(1)
subrho1 = I.component(1).irrep(1)
%%
%
I.component(2)
subrho2 = I.component(2).irrep(1)
%%
% with their bases:
subrho1.U
subrho2.U

%% The commutant algebra
% The commutant algebra of $\rho$ is composed of all the matrices
% $M$ that commute with $\rho$, that is $M \rho_g = \rho_g M$ for
% all $g$ in the group.
%
% *RepLAB* gives an access to that algebra:
C = rho.commutant
%%
% and we can sample generic matrices from that algebra
%
C.sample
%%
% or perform an orthogonal projection of arbitrary matrices in that
% algebra
Mgen = rand(n, n)
%%
%
M = C.project(rand(n,n))
%%
%%
% Which is has a block diagonalization in the symmetry adapted basis:
U = I.U
%%
%
U*M*U'
