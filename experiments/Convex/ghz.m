n = 4; % number of subsystems
d = 2; % subsystem dimension
% the symmetry group of the GHZ states had three parts:
%
% 1. A continuous part given by phases (i.e. copies of U(1)). For the n-1 subsystems, each of the d
%    levels has a phase as an independent degree of freedom; the phase applied on the d-th level
%    of last n-th subsystem is equal to the inverse of the product of the phases of the d-th level
%    of the other subsystem.
%
%    We represent that symmetry using a torus group U(1)^((n-1)*d), as seen below. The torus
%    elements are stored by enumerating first i=1...n-1, and then j=1...d, so that, for a torus
%    element ``t``, the information can be accessed by performing ``reshape(t, [n-1 d])``
%
%    The torus group U(1)^N stores phases as elements of an additive group, with values in the
%    interval ``[0,1[``, and addition "modulo 1" as the binary operation.
%
% 2. A copy of S(n) that permutes the n subsystems.
%
% 3. A copy of S(d) that permutes the d levels.
%
% The discrete part acts on the torus; the action of S(d) is pretty easy to write and is a permutation
% action, while the action of S(n) involves the "virtual" phase for the n-th subsystem.

% The discrete part is a direct product of S(d) and S(n)
Sn = replab.S(n);
Sd = replab.S(d);
G = Sd.directProduct(Sn);

% We write the action of that finite group on torus elements
torusToFull = [eye(n-1); -ones(1, n-1)];
fullToTorus = [eye(n-1) zeros(n-1, 1)];
Sn_images = cell(1, Sn.nGenerators);
Sn_nr = Sn.naturalRep;
for i = 1:Sn.nGenerators
    Sn_images{i} = fullToTorus*Sn_nr.image(Sn.generator(i))*torusToFull;
end
Sn_rep = Sn.repByImages('R', n-1, 'preimages', Sn.generators, 'images', Sn_images);
% The action of a finite group on a torus can be written using a linear representation with
% integer coefficients
torusRep = G.tensorFactorRep('R', {Sd.naturalRep Sn_rep});
% We construct the GHZ group as a semidirect product
ghzGroup = replab.TorusGroup(zeros(0, torusRep.dimension)).semidirectProductFromRep(torusRep);

% Now we move to the action of the GHZ group on the GHZ state space

% First, the finite part acts as follows.
% * S(d) acts as a tensor product of dxd permutation matrices
% * S(n) acts by relabeling the subsystems.
% RepLAB has support for both representations.
stateFiniteRep = G.commutingFactorRepsRep('R', d^n, {Sd.naturalRep.tensorPower(n) Sn.indexRelabelingRep(d)});

% Now the action of the continuous part is trickier.
% Basically, we write a map between two torus groups.
%
% * The source torus has dimension d*(n-1), as stated above.
% * The target torus has dimension d^n, as it corresponds to the state space.
%
% Note that, as torus elements are written additively, morphisms between tori can be written
% using matrices with integer entries. We construct such a map below.
maps = cell(1, n);
lastmap = zeros(d, n-1, d);
for i = 1:n-1
    map = zeros(d, n-1, d);
    map(:,i,:) = eye(d);
    lastmap = lastmap - map;
    maps{i} = reshape(map, [d (n-1)*d]);
end
maps{n} = reshape(lastmap, [d (n-1)*d]);
subSysReps = cell(1, n);
for i = 1:n
    subSysReps{i} = ghzGroup.N.diagonalRep(maps{i});
end
contRep = kron(subSysReps{:});

% Variant:
% tm = zeros(d^n, n-1, d);
% for i = 1:n-1
%     tm = reshape(tm, [d^(n-i) d d^(i-1) n-1 d]);
%     for j = 1:d
%         tm(:,j,:,i,j) = tm(:,j,:,i,j) + 1;
%     end
% end
% tm = reshape(tm, [d d^(n-1) n-1 d]);
% for j = 1:d
%     tm(j,:,:,j) = tm(j,:,:,j) - 1;
% end
% tm = reshape(tm, [d^n (n-1)*d]);
% contRep = ghzGroup.N.mapRep(tm);

% finally, "concatenate" the compatible representations of the components of the semidirect product
rep = ghzGroup.semidirectProductRep(stateFiniteRep.complexification, contRep);
