% # Symmetries in Bell nonlocality, part 3
%
% Let us now turn our attention to the symmetries of the full scenario. We first recall the definitions of part 1.

run ../../../replab_init
outputGroup = replab.S(2);
outputRep = outputGroup.naturalRep;
inputGroup = replab.S(2);
ioGroup = inputGroup.wreathProduct(outputGroup);
ioRep = ioGroup.imprimitiveRep(outputRep);

% The same story repeats for relabelings of parties: the scenario involves
% two homogeneous parties. We thus have two copies of the group relabeling
% inputs and/or outputs (one for Alice, one for Bob), and a copy of $S_2$ that
% permutes the parties. This group desribing all possible relabellings in
% this scenario is thus given by:

scenarioGroup = outputGroup.wreathProduct(ioGroup);

scenarioGroup.elements

% The representation on the behavior $P(ab|xy)$ is however a primitive
% representation, as $P(a|x; b|y)$ ressembles a tensor. Inside each party,
% we use the imprimitive representation constructed before.

probRep = scenarioGroup.primitiveRep(ioRep);

% Decomposition of the full probability space
%
% We can now compute the decomposition of this representation on $P(ab|xy)$:

dec = probRep.decomposition

% There are 6 irreducible components

dec.nComponents

% These components correspond to the following physical elements ( see [arXiv:1610.01833](https://arxiv.org/abs/1610.01833) ) for more details:
%
% - The global normalization of probabilities
% - The difference of normalization between $x=y$ and $x\ne y$
% - The difference of normalization between $x=0$ and $x=1$
% - The signaling between Alice and Bob
% - The marginal probabilities of Alice and Bob
% - The correlation between Alice and Bob
