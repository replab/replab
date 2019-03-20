%% Symmetric Semidefinite Programs
%
% This document illustrated how *RepLAB* can be used to solve Semidefinite
% Programs (SDP) subject to symmetries.

%% Preparation
% Before using *RepLAB* commands, we first add the paths:
replab_addpaths
%%
% In order to solve convex optimization problems, the
% YALMIP interface is needed, (see the
% <../installation.html#additional-resources additional resources> for more
% details).

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
% It turns out that a positive semi-definite matrix that is invariant under
% some joint permutation of its lines and columns can be decomposed into a
% block diagonal form. This allows to:
%
% * decompose PSD blocks into smaller PSD blocks
% * set many variables to zero, hence reducing the number of variables in
% the problem
%
% Here we show how *RepLAB* can be used to accomplish this simplification
% automatically.

%% A simple example
% Let us consider a 3x3 matrix $M$ with trace 1 that is symmetric under cyclic
% permutation of its indices, i.e. it satisfies $M([2\ 3\ 1], [2\ 3\ 1]) = M$.
% We are interested in some property of the matrix, for instance the
% smallest value of the off-diagonal element $M(1,2)$ for which the matrix
% must necessarily have a negative eigenvalue.

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
% This shows that the lowest possible value of $M(1,2)$ which is compatible
% with a matrix $M$ having only positive eigenvalues is $-1/6$.

%%
% <html>
% <h3>Symmetric formulation</h3>
% </html>
%%
% Using *RepLAB*, we can solve this problem while taking into account the
% structure of the matrix $M$.
%
% We start by defining a matrix which satisfies the desired symmetry
MSym = replab.Sdprep.fromGenerators({permutation});
%%
% We can then perform the optimization with:
constraintsSym = [MSym(1,1)+MSym(2,2)+MSym(3,3) == 1, MSym >= 0];
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
% and simpler constraints, as described in the following table:
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
% To see this, we examine the variables involved.
%
% In the first case, we have
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
% and notice that this variable is made up of two PSD blocks:
%
% * One block of size 1x1, which appears 1 time
% * One block of size 2x2, which appears 1 time
%
% Each block can be examined individually:
MSym.blocks{1}
%%
MSym.blocks{2}
%%
% and is found to contain exactly one variable.
%%
% The constraints this time are
constraintsSym
%%
% There are thus :
%
% * 2 variables
% * 1 equality constraint
% * SDP blocks of size 1x1 and 2x2
