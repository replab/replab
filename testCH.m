generators = {[1  4  7  2  5  8  3  6  9]
              [1 -2 -3 -4  5  6 -7  8  9]
              [1  3  2  4  6  5 -7 -9 -8]}';
G = replab.SignedSymmetricGroup(9).subgroup(generators);

% We construct the images of the generators in the Collins-Gisin
% picture. NOTE: Here, each "line" corresponds to the image of one
% basis element. The order of the basis elements is as follows:
% 1, PA(0|0), PA(0|1), PB(0|0), P(00|00), P(00|10), PB(0|1), P(00|01), P(00|11)
% First, permutation of parties:
image1 = [1  0  0  0  0  0  0  0  0
          0  0  0  1  0  0  0  0  0
          0  0  0  0  0  0  1  0  0
          0  1  0  0  0  0  0  0  0
          0  0  0  0  1  0  0  0  0
          0  0  0  0  0  0  0  1  0
          0  0  1  0  0  0  0  0  0
          0  0  0  0  0  1  0  0  0
          0  0  0  0  0  0  0  0  1];
% Second, permutation of all outcomes:
image2 = [1  0  0  0  0  0  0  0  0
          1 -1  0  0  0  0  0  0  0
          1  0 -1  0  0  0  0  0  0
          1  0  0 -1  0  0  0  0  0
          1 -1  0 -1  1  0  0  0  0
          1  0 -1 -1  0  1  0  0  0
          1  0  0  0  0  0 -1  0  0
          1 -1  0  0  0  0 -1  1  0
          1  0 -1  0  0  0 -1  0  1];
% Third, permutation of Alice's settings and Bob's outcome for his
% second setting:
image3 = [1  0  0  0  0  0  0  0  0
          0  0  1  0  0  0  0  0  0
          0  1  0  0  0  0  0  0  0
          0  0  0  1  0  0  0  0  0
          0  0  0  0  0  1  0  0  0
          0  0  0  0  1  0  0  0  0
          1  0  0  0  0  0 -1  0  0
          0  0  1  0  0  0  0  0 -1
          0  1  0  0  0  0  0  -1 0];
% line 7 PB(0|1) -> PB(1|1) = 1 - PB(0|1)
% line 8 P(00|01) -> P(01|11) = PA(0|1) - P(00|11)
% line 9 P(00|11) -> P(01|01) = PA(0|0) - P(00|01)
rep = G.repByImages('R', 9, 'images', {image1 image2 image3}, 'preimages', {image1 image2 image3});




% %%% THE FOLLOWING COMMAND FAILS %%%
D = rep.decomposition;




indexMatrix = [  1     2     3     6     7     8    11    12    13
                 2     2     4     7     7     9    12    12    14
                 3     4     3     8     9     8    13    14    13
                 6     7     8     6     7     8    16    17    18
                 7     7     9     7     7     9    17    17    19
                 8     9     8     8     9     8    18    20    18
                 11    12    13    16    17    18    11    12    13
                 12    12    14    17    17    20    12    12    14
                 13    14    13    18    19    18    13    14    13];

objective = [0 -1 0 -1 1 1 0 1 -1];
m = max(indexMatrix(:));
y = sdpvar(m, 1);
y(1) = 1;
A = cell(1, m);
Asym = cell(1, m);
X = zeros(9, 9);% We formulate the non-symmetrized SDP:
%vars = [0; sdpvar(max(max(indexMatrix)),1)];
%tmp = sparse(1:numel(indexMatrix), reshape(1+indexMatrix,1,numel(indexMatrix)), true);
%sdpMatrix = reshape(tmp*vars, size(indexMatrix));
%obj1 = (sqrt(2)-1)/2;
%evalc('solvesdp([sdpMatrix >= 0, sdpMatrix(1,1) == 1], -obj, sdpsettings(''verbose'', 0))');
%obj2 = value(obj);

Xsym = zeros(9, 9);
for i = 1:m
    A{i} = (indexMatrix == i);
    X = X + A{i} * y(i);
    % we use the embedding map, because we work with the Hermitian invariant,
    % not the commutant
    Asym{i} = D.E_internal * A{i} * D.E_internal';
    Asym{i} = (Asym{i} + Asym{i}')/2;
    Xsym = Xsym + Asym{i} * y(i);
end
mask = ones(9,9);
mask(1:2,1:2) = 0;
for i = 3:9
    mask(i, i) = 0;
end
% tricky part: we have lots of equality constraints, as we cannot use the orbit trick
% to reduce the number of scalar variables
%
% other possibility: project the variables "y" in the invariant subspace of the representation
% of the symmetry group (this is basically what we are doing here implicitly)
Csym = [Xsym >= 0
        Xsym(4,4) == Xsym(5,5)
        Xsym(6,6) == Xsym(7,7)
        Xsym(8,8) == Xsym(9,9)
        Xsym(logical(mask)) == 0];
obj = objective*X(:,1);
solvesdp(Csym, -obj)
double(obj)
%C = [X >= 0];
%solvesdp(C, -obj)
%double(obj)
