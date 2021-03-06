n = 3; % number of subsystems
d = 3; % subsystem dimension
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
ghzGroup = replab.TorusGroup.semidirectProductFromRep(torusRep);

% Now we move to the action of the GHZ group on the GHZ state space

% First, the finite part acts as follows.
% * S(d) acts as a tensor product of dxd permutation matrices
% * S(n) acts by relabeling the subsystems.
% RepLAB has support for both representations.
stateFiniteRep = G.commutingProductFactorRep('R', d^n, {Sd.naturalRep.tensorPower(n), Sn.indexRelabelingRep(d)});

% Now the action of the continuous part is trickier.
% Basically, we write a map between two torus groups.
%
% * The source torus has dimension d*(n-1), as stated above.
% * The target torus has dimension d^n, as it corresponds to the state space.
%
% Note that, as torus elements are written additively, morphisms between tori can be written
% using matrices with integer entries. We construct such a map below.

% d^n is the target space, and we split the source space
tm = zeros(d^n, n-1, d);
% we iterate over the first n-1 subsystem whose phase is readily available
for i = 1:n-1
    % we split the target space into the ``d^(i-1), d^(n-i)`` that we leave aside, and the
    % d-dimensional subsystem that we target specifically
    tm = reshape(tm, [d^(i-1) d d^(n-i) n-1 d]);
    for j = 1:d
        % we apply each phase
        tm(:,j,:,i,j) = tm(:,j,:,i,j) + 1;
    end
end
% now, we address the last subsystem
tm = reshape(tm, [d^(n-1) d n-1 d]);
for j = 1:d
    tm(:,j,:,j) = tm(:,j,:,j) - 1;
end
% put it back as a torus morphism coefficient matrix
tm = reshape(tm, [d^n (n-1)*d]);
% construct the representation
contRep = ghzGroup.N.mapRep(tm);
% finally, "concatenate" the compatible representations of the components of the semidirect product
rep = ghzGroup.productRep(stateFiniteRep.complexification, contRep);