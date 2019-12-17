%% Relabelings in the CHSH scenario
%
% We describe the output/input/party relabelings in the CHSH scenario,
% and decompose their action on probability distributions P(ab|xy).
%

%% Defining the output relabelings
% Outputs are binary, so we use the symmetric group of domain size 2
% and the representation on P(a) is the defining representation
S2 = replab.S(2)
outputRep = S2.definingRep

%% Defining the party relabelings
% We have a copy of S2 acting on the outputs of the first measurement,
% and a second copy of S2 acting on the outputs of the second measurement,
% while another copy of S2 permutes the choice of measurement (input).
%
% This is described by the wreath product of S2 by S2.
%
% The representation on P(a|x) is the imprimitive representation,
% given that we use the defining representation for the inner group
party = S2.wreathProduct(S2)
partyRep = party.imprimitiveRep(outputRep)

%% Defining the scenario relabelings
% The same story repeats for scenario relabelings; we have two copies
% of the party relabeling group (one for Alice, one for Bob), and a copy
% of S2 that permutes the parties.
%
% The representation on P(ab|xy) is however a primitive representation,
% as P(a|x; b|y) ressembles a tensor. Inside each party, we use the
% imprimitive representation seen before.

scenario = S2.wreathProduct(party)
probRep = scenario.primitiveRep(partyRep);

%% We now compute the decomposition of the representation
probRep.decomposition
