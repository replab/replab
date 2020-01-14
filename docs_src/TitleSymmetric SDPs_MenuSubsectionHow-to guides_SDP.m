%% Symmetric SDPs
%
% This document illustrated how *RepLAB* can be used to solve Semidefinite
% Programs (SDP) subject to symmetries. Several methods for constructing
% SDPs subject to symmetries are presented.


%% Preparation
% As always, before using *RepLAB* commands, initialize the library:
replab_init
%%
% Commands in this document makes use ot the YALMIP interface to solve
% convex optimization problems, see
% <../installation.html installation instructions>.


%% Introduction
% <https://en.wikipedia.org/wiki/Semidefinite_programming Semidefinite
% Programming> is a form of optimization that admits semidefinite
% constraints (as in the condition that all eigenvalues of a matrix must be
% positive). It naturally arises in various fields, including operational
% research and polynomial optimization. 
% 
% The ability to solve a semidefinite program depends heavily on:
%
% * the size of the positive semi-definite (PSD) blocks
% * the number of variables and constraints involved
%
% A positive semi-definite matrix that is invariant under some joint
% permutation of its lines and columns can be decomposed into a block
% diagonal form. This allows to:
%
% * decompose PSD blocks into smaller PSD blocks
% * set many variables to zero, hence reducing the number of variables in
% the problem
%
% As we shows below, *RepLAB* performs this simplification automatically.


%% Formulating a symmetric SDP
% To illustrate the usage of *RepLAB* for symmetric SDPs, we consider a
% simple example involving a 3x3 matrix $M$ with trace 1 that is symmetric
% under cyclic permutation of its indices, i.e. it satisfies $M([2\ 3\ 1], [2\ 3\ 1]) = M$.
% We ask what is the smallest value that the off-diagonal element
% $M(1,2)$ can take if this matrix $M$ has only positive eigenvalues.

%%
% <html>
% <h3>Symmetric formulation</h3>
% </html>
%%
% Using *RepLAB*, we can solve this problem as follows.
%
% We start by defining a matrix which satisfies the desired symmetry
permutation = [2 3 1];
MSym = replab.CommutantVar.fromPermutations({permutation}, 'symmetric', 'real');
%%
% We can then perform the optimization with:
constraintsSym = [trace(MSym) == 1, MSym >= 0];
diagnosticSym = optimize(constraintsSym, MSym(1,2), sdpsettings('verbose', 0));
MSymOpt = value(MSym)
%%
% We find the critical value of $-1/6$.
%
% At the end of this page, we discuss how this formulation is more
% efficient than a direct formulation which would not take advantage of the
% symmetry properties of the considered matrix.


%% Constructing a symmetric SDP matrix with additional structure
% An SDP matrix might sometimes be subject to more than only symmetry
% constraints. When these additional constraints take the form of equality
% between some of the matrix elements, it can be conveniently described by
% a matrix having at each element the index of the corresponding variable.
% All elements with identical index are then understood as being equal to
% each other.
%%
% For instance, imposing on a 3x3 matrix that elements (1,1) ans (2,1) are
% must be equal to each other can be described by the following index
% matrix:
indexMatrix = [1 2 3
               1 4 5
               6 7 8];
%%
% An SDP matrix satisfying this constraint is then obtained
cstrSdpMatrix = replab.CommutantVar.fromIndexMatrix(indexMatrix, {permutation}, 'symmetric', 'real')
%%
% In this trivial case, the SDP matrix is left to contain only one
% variable: the additional constraints collapsed the all group orbits
% together.


%% Imposing symmetry to an existing SDP matrix
% Symmetry constraints can also be straightforwardly imposed on existing
% SDP matrices with arbitrary structure.
%%
% For instance, consider the following special SDP matrix
x = sdpvar;
y = sdpvar
MSpecial = [1 x y
            x 1 y
            x y 1];
%%
% We can directly impose cyclic symmetry onto this matrix:
MSpecialSym = replab.CommutantVar.fromSdpMatrix(MSpecial, {[2 3 1]})
%%
% Requesting this matrix to be PSD now imposes both
%
% * Positivity of the 1x1 and 2x2 blocks
% * Equality with the imposed form
%
% as can be seen with
MSpecialSym >= 0


%% Block-diagonalizing a symmetric SDP matrix
% When an SDP matrix is invariant under the considered permutations,
% *RepLAB* can be used to block-diagonalize it. This allows imposing the
% positivity of the matrix through the positivity of small blocks. As an
% example, consider the following matrix
MInvariant = [x 1 y
              y x 1
              1 y x];
%%
% It is indeed invariant:
MInvariant - MInvariant(permutation,permutation)
%%
% But not transpose-invariant (which is necessay for PSD matrices), so we
% enforce it:
MInvariant = MInvariant + MInvariant';
MInvariant - MInvariant'
%%
% We can now block-diagonalize it by calling
MInvariantBlock = replab.CommutantVar.fromSymSdpMatrix(MInvariant, {[2 3 1]})
%%
% No new variable has been introduced in the new object, but the block
% structure has been found:
full(MInvariantBlock.blockMask)
%%
% The block structure is used when requesting this matrix to be PSD:
MInvariantBlock >= 0


%% Comparison with a direct formulation
%
% To conclude, let us show in more detail why the SDP formulation of a
% problem is more efficient if it takes advantage of the available symmetry
% properties. For this, we consider again the problem described at the
% beginning of this page. This problem can be solved directly as follows:
M = sdpvar(3);
constraints = [trace(M) == 1, M(permutation, permutation) == M, M >= 0];
diagnostic = optimize(constraints, M(1,2), sdpsettings('verbose', 0))
MOpt = value(M)
%%
% Again, we find that the lowest possible value of $M(1,2)$ which is
% compatible with the matrix $M$ having only positive eigenvalues is
% $-1/6$. However, this last SDP problem is more complex than the first one
% which takes into account symmetries.

%%
% <html>
% <h3>Complexity comparison</h3>
% </html>
%%
% The symmetric formulation of the above problem involves fewer variables
% and simpler constraints, as summarized in the following table:
%
% <html>
% <table align="center">
%  <tr>
%    <th>formulation</th>
%    <th>direct</th>
%    <th>symmetric</th>
%  </tr>
%  <tr>
%    <td># variables</td>
%    <td>6</td>
%    <td>2</td>
%  </tr>
%  <tr>
%    <td># equality constraints</td>
%    <td>10</td>
%    <td>1</td>
%  </tr>
%  <tr>
%    <td>size of PSD blocks</td>
%    <td>3</td>
%    <td>1 and 2</td>
%  </tr>
% </table> 
% </html>

%%
% This can be checked by examining the variables involved. In the non-symmetrized
% case, we have
M
%%
constraints
%%
% we see that it involves
%
% * 6 variables
% * 1+3x3=10 equality constraints
% * 1 PSD block of size 3x3
%
%%
% In the symmetrized case, we have
MSym
%%
% In other words, the matrix is made of two blocks of size 1x1 and 2x2, and
% involves altogether just 2 variables.
%%
% The constraints this time are
constraintsSym
%%
% This formulation thus involves :
%
% * 2 variables
% * 1 equality constraint
% * SDP blocks of size 1x1 and 2x2

