%% Representations of the symmetric group S(n)
%
% The symmetric group S(n) is group of primary importance. This
% document presents several representations of this group with **RepLAB**.

%%
% Before trying any of the *RepLAB* commands, we must initialize the library:
replab_init

%% The symmetric group
% Consider $n$ distinct elements, numbered $1,2,\ldots,n$. The symmetric
% group describes all possible ways in which these elements can be
% permuted.
% 
% With $n=3$, for instance, two possible arrangements are
[1 2 3];
%%
% and
[1 3 2];
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
% The group structure defines how elements of the group compose with each
% other. For instance, permuting the first two elements in 1,2,3, followed
% by a permutation of the second and third elements results in a cyclic
% permutation:
element1 = [2 1 3];
element2 = [1 3 2];
S3.compose(element1, element2)

%% The defining representation of S(n)
% A group can always be represented by matrices of finite dimension.
% In a $d$-dimensional matrix representation, each group element is
% associated with a $d \times d$ matrix called the image of the group
% element. The matrix multiplication rule then reflects the group
% composition rule.
%
% Several matrix representations are possible for the same group, possibly
% with various dimension $d$. Here, we are going to explore some
% representations of our group S(3).
%
% When seeing the group as a permutation within $n$ elements, a natural
% choice of representation arises in dimension $d=n$, where the permutation
% is applied to the coordinates of the vector space $R^n$. This is called
% the *defining representation* of a group. It can be constructed easily:
defRep = S3.definingRep
%%
% The image of a group element in this representation can be extracted
defRep.image([2 1 3])
%%
% We see that the matrices representing the group elements are $n \times n$
% permutation matrices. which simply apply the permutation to the
% coordinates of the vector space $R^n$.
%
% Let us check that the representation is valid, i.e. that is follows
% indeed the group structure
image1 = defRep.image(element1);
image2 = defRep.image(element2);
image1*image2
defRep.image(S3.compose(element1, element2))
%%
% We see that the image of the product of elements is indeed the product of
% the images of the respective elements, i.e. the algebra of matrices
% acting on the representation images reflects accurately the group
% algebra.
%
% It can be checked, moreover, that the defining representation is
% faithful, i.e. each group element has its own distinct image. This needs
% not always be the case as in the following two examples.

%% The trivial representation of S(n)

% Whether computed from, or from, we see that the matrix representing the permutation [2 3 1] ise

%
% The action of a group can always be represented by the product of d x d
% matrices
%
% When representing elements of a group as d x d matrices, the
% action of the group amounts to matrix multiplication. The 

% Elements of a group can always be represented by 
% The action of a group can be represented by the matrix product acting o.


%% Faithful representation

%% Trivial representation



