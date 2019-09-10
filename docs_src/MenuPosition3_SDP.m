%% Symmetric SDPs
%
% This document illustrated how *RepLAB* can be used to solve Semidefinite
% Programs (SDP) subject to symmetries.

%% Preparation
% As always, before using *RepLAB* commands, first add the paths:
replab_addpaths
%%
% In order to solve convex optimization problems, the
% YALMIP interface is needed, (see the
% <../installation.html installation instructions> for more details).

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

%% A simple example
% Consider a 3x3 matrix $M$ with trace 1 that is symmetric under cyclic
% permutation of its indices, i.e. it satisfies $M([2\ 3\ 1], [2\ 3\ 1]) = M$.
% In this example, we ask what is the smallest value that the off-diagonal
% element $M(1,2)$ can take if $M$ has only positive eigenvalues.

%%
% <html>
% <h3>Direct formulation</h3>
% </html>
%%
% Using YALMIP, this problem can be solved directly:
M = sdpvar(3);
permutation = [2 3 1];
constraints = [trace(M) == 1, M(permutation, permutation) == M, M >= 0];
diagnostic = optimize(constraints, M(1,2), sdpsettings('verbose', 0))
MOpt = value(M)
%%
% Hence the lowest possible value of $M(1,2)$ which is compatible
% with the matrix $M$ having only positive eigenvalues is $-1/6$.

%%
% <html>
% <h3>Symmetric formulation</h3>
% </html>
%%
% Using *RepLAB*, we can solve this problem while taking into account the
% structure of the matrix $M$.
%
% We start by defining a matrix which satisfies the desired symmetry
MSym = replab.CommutantVar.fromPermutations({permutation});
%%
% We can then perform the optimization with:
constraintsSym = [trace(MSym) == 1, MSym >= 0];
diagnosticSym = optimize(constraintsSym, MSym(1,2), sdpsettings('verbose', 0))
MSymOpt = value(MSym)
%%
% Again, we find the critical value of $-1/6$. This last formulation is
% however more concise as we now discuss.

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
% This can be checked by examining the variables involved. In the first
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
% In the second case, we have
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

%% Imposing symmetry to an existing SDP matrix
% Symmetry constraints can also be straightforwardly imposed on existing
% SDP matrices.
%%
% For instance, consider the special SDP matrix
vars = sdpvar(1,2);
MSpecial = [1 vars(1) vars(2)
            vars(1) 1 vars(2)
            vars(1) vars(2) 1];
%%
% We can directly impose cyclic symmetry onto this matrix:
MSpecialSym = replab.CommutantVar.fromSdpMatrix(MSpecial, {[2 3 1]})
%%
% Requesting this matrix to be PSD now imposes both
% * Positivity of the 1x1 and 2x2 blocks
% * Equality between the matrix and the imposed form
% as can be seen with
MSpecialSym >= 0
