% # Symmetries in Bell nonlocality, part 2
%
% We describe the relabelings in the CHSH scenario, and decompose their action on probability distributions $P(ab|xy)$, now using an abstract group construction.
%
% We identify those relabelings with a [wreath product](https://en.wikipedia.org/wiki/Wreath_product) construction, which RepLAB supports.
% This topic guide can also be considered as a (small) introduction to the group constructions available in RepLAB.

% Before using *RepLAB* commands, we must first initialize the library:

run ../../../replab_init

% ## Only outputs
%
% Outputs are binary, so the relabelling of outputs is the symmetric group
% of domain size 2.

outputGroup = replab.S(2)

% Now, an element of `outputGroup` is just a permutation of the two outputs:

outputGroup.elements

% The representation of this relabelling on the probabilities $P(a)$ is then the defining representation

outputRep = outputGroup.naturalRep
outputRep.image([2 1])

% ## Inputs and outputs
%
% In the CHSH scenario, each party has two binary measurements. Therefore,
% a copy of $S_2$ acts on the outputs of the first measurement, and a second
% copy of $S_2$ acting on the outputs of the second measurement, while
% another copy of $S_2$ permutes the choice of measurement (input).
%
% This is described by the wreath product of $S_2$ by $S_2$.

inputGroup = replab.S(2)
ioGroup = inputGroup.wreathProduct(outputGroup)

% Elements of `ioGroup` are now composed of a permutation of inputs, and the conditional output permutation for $x=1,2$:

ioGroup.elements

% The representation on the conditional probility $P(a|x)$ is the imprimitive representation,
% given that we use the defining representation
% for the inner group

ioRep = ioGroup.imprimitiveRep(outputRep);

ioRep.image({[2 1] {[1 2] [1 2]}}) % permuting only the inputs

ioRep.image({[1 2] {[2 1] [1 2]}}) % permuting only the outputs, conditioned on $x=1$
