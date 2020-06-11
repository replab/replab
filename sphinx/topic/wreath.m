% # Use of wreath product groups in quantum information
%
% We describe the output/input/party relabelings in the CHSH scenario, and decompose their action on probability distributions $P(ab|xy)$.
%
% We identify those relabelings with a [wreath product](https://en.wikipedia.org/wiki/Wreath_product) construction, which RepLAB supports. This topic guide can be considered as a (small) introduction to the group constructions available in RepLAB.

% Before using *RepLAB* commands, we must first initialize the library:

run ../../replab_init

% ## Relabelings the outputs
%
% Outputs are binary, so the relabelling of outputs is the symmetric group
% of domain size 2. The representation of this relabelling on the
% probabilities $P(a)$ is then the defining representation

outputGroup = replab.S(2)
outputRep = outputGroup.naturalRep

% ## Adding the relabelings of inputs
%
% In the CHSH scenario, each party has two binary measurements. Therefore,
% a copy of $S_2$ acts on the outputs of the first measurement, and a second
% copy of $S_2$ acting on the outputs of the second measurement, while
% another copy of $S_2$ permutes the choice of measurement (input).
%
% This is described by the wreath product of $S_2$ by $S_2$.

ioGroup = replab.S(2).wreathProduct(outputGroup)

% The representation on the conditional probility $P(a|x)$ is the
% imprimitive representation, given that we use the defining representation
% for the inner group

ioRep = ioGroup.imprimitiveRep(outputRep);

% ## Adding the relabelings of parties
%
% The same story repeats for relabelings of parties: the scenario involves
% two homogeneous parties. We thus have two copies of the group relabeling
% inputs and/or outputs (one for Alice, one for Bob), and a copy of $S_2$ that
% permutes the parties. This group desribing all possible relabellings in
% this scenario is thus given by:

scenarioGroup = outputGroup.wreathProduct(ioGroup);

% The representation on the behavior $P(ab|xy)$ is however a primitive
% representation, as $P(a|x; b|y)$ ressembles a tensor. Inside each party,
% we use the imprimitive representation constructed before.

probRep = scenarioGroup.primitiveRep(ioRep);

% ## Decomposition of the full probability space
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
