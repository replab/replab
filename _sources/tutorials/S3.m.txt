% # Tutorial for mathematicians
%
% The symmetric group $S_n$ is group of primary importance. Here, we consider this group for $n=3$ and discuss several of its representations with *RepLAB*.
%
% Before trying any of the *RepLAB* commands, we must initialize the library:

addpath([pwd, '/../..']);
replab_init('verbose', 0);

% ## The symmetric group
%
% Consider $n$ distinct elements, numbered $1,2,\ldots,n$. The symmetric group describes all possible ways in which these elements can be permuted.
%
% With $n=3$, for instance, two possible arrangements are

[1 2 3];

% and

[1 3 2];

% In general there are $n!$ such arrangements (i.e. $6$ arrangements for $n=3$).
%
% In *RepLAB*, the symmetric group can be constructed with a simple command

n = 3;
S3 = replab.S(n);

% The elements of the group can be easily listed

S3.elements

% The number of elements is the order of the group

S3.order

% As we see, permutations are represented by row-vectors with an arrangement of the elements $1,\ldots,n$.
%
% The group structure defines how elements of the group compose with each other. 
% For instance, permuting the first two elements in $1,2,3$, followed by a permutation of the second and third elements results in a cyclic permutation:

element1 = [2 1 3];
element2 = [1 3 2];
S3.compose(element1, element2)

% ## The natural representation of S(3)
%
% A group can always be represented by matrices of finite dimension. In a $d$-dimensional matrix representation, each group element is associated with a $d \times d$ matrix called the image of the group
% element. The action of matrix multiplication on the images then reflects the group composition rule.
%
% Several matrix representations are possible for the same group, possibly in various dimension $d$. Here, we are going to explore some representations of the group $S_3$
%
% When seeing the group as a permutation within $n$ elements, a natural choice of representation arises in dimension $d=n$, where the permutation
% is applied to the coordinates of the vector space $R^n$. This is called the *natural representation* of a group. It can be constructed easily:

natRep = S3.naturalRep;

% The image of a group element in this representation can be extracted

natRep.image([2 1 3])

% We see that the matrices representing the group elements are $n \times n$ permutation matrices, which simply apply the permutation to the coordinates of the vector space $R^n$.
%
% For a representation to be valid, the images must follow the structure of the group:

image1 = natRep.image(element1);
image2 = natRep.image(element2);
image1*image2
natRep.image(S3.compose(element1, element2))

% We see that the image of the product of elements is indeed the product of the images of the respective elements, i.e. the algebra of matrices
% acting on the representation images accurately reflects the group algebra.
%
% It can be checked, moreover, that the natural representation is faithful, i.e. each group element has its own distinct image. This needs
% not always be the case as in the following two examples. 
%
% ## The parity representation of $S_3$
%
% To construct an arbitrary group representation in *RepLAB*, it is sufficient to provide the image of the generators of the group. The group $S_3$ has two generators:

generators = S3.generators;
generators{:}

% Let us construct a representation of this group in dimension

d = 1;

% associating the each generator its parity (i.e. either $+1$ or $-1$). Since the parity of [2 3 1] is 1 and the parity of [2 1 3] is -1, this is achieved by instantiating the class *RepByImage* as follows:

parRep = replab.RepByImages(S3, 'R', d, {1 -1}, {1 -1});

% We can check that this representation is valid

image1 = parRep.image(element1)
image2 = parRep.image(element2)
image1*image2
parRep.image(S3.compose(element1, element2))

% The parity of each binary permutation ``[2 1 3]`` and ``[1 3 2]`` is $-1$, and their product, the cycle [2 3 1], has parity $1$.
%
% Clearly, this representation carries some information about the group, namely the parity, but some elements have the same image:

allElements = S3.elements.toCell;
for i = 1:length(allElements)
    disp(['Image of [', num2str(allElements{i}), '] : ', num2str(parRep.image(allElements{i}))]);
end

% Therefore, this representation is not *faithful*.
%
% ## The trivial representation of S(3)
% An even simpler representation of the group S(n) is one in which the image of all elements is set to $1$. This is called the *trivial representation*.
%
% It can be constructed similarly from the image of the generators:

trivRep = replab.RepByImages(S3, 'R', d, {1 1}, {1 1});

% This time, the group law is trivially satisfied since all images are indeed equal to $1$:

for i = 1:length(allElements)
    disp(['Image of [', num2str(allElements{i}), '] : ', num2str(trivRep.image(allElements{i}))]);
end

% ## Direct sum of representation of S(3)
%
% We have seen that a group can admit several representations. Any two representations of a group in dimension $d_1$ and $d_2$ can be joined
% together by taking their *direct sum* in order to create a new valid representation in dimension $d_1 + d_2$.
%
% For instance, we can construct a $4$-dimensional representation of S(3) by taking the direct sum of the natural representation with the parity
% representation

newRep = blkdiag(natRep, parRep)

% The images are then $4 \times 4$ matrices

newRep.image(element1)
newRep.image(element2)

% As expected, the images have a block diagonal structure. The properties
% of a direct sum representation are thus captured by the properties of
% the representation in each block. Such a representation is *reducible*.
%
% ## Irreducible representations of S(3)
%
% Given a representation, it is not always obvious whether it is the direct
% sum of more fundamental representations or if the representation is
% *irreducible*. For instance, the images of the natural representation
% introduced above are not obviously block-diagonal. Yet, this
% representation is also reducible. To see this, one needs to find a basis
% of the vector space $R^d$ in which the images have a block diagonal
% structure.
%
% This can be achieved in *RepLAB* by decomposing the representation

natDec = natRep.decomposition.nice;
natDec.nComponents

% This shows that the defininig representation of $S_3$ has two *irreducible components*, of dimension $1$ and $2$ respectively:

natDec.component(1).irrepDimension
natDec.component(2).irrepDimension

% The decomposition also provides the change of basis matrix which makes |this decomposition apparent

basis = natDec.basis

% The images of the two elements of $S_3$, ``[2 1 3]`` and ``[1 3 2]`` are indeed block-diagonal in this basis:

inv(basis)*natRep.image(element1)*basis
inv(basis)*natRep.image(element2)*basis

% ## Summary
%
% This shows that the natural representation of S(3) contains
%
% * One copy of the trivial representation of S(3)
% * One copy of the a faithful 2 x 2 representation of S(3)
