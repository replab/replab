% ---
% jupyter:
%   jupytext:
%     formats: ipynb,m:light
%     text_representation:
%       extension: .m
%       format_name: light
%       format_version: '1.5'
%       jupytext_version: 1.3.2
%   kernelspec:
%     display_name: Octave
%     language: octave
%     name: octave
% ---

% # 4. A taste of composition with RepLAB
%
% This is part IV of the companion notebook to the RepLAB talk at the [Quantum Causal Structures](http://www.cs.ox.ac.uk/conferences/QCS2019/) workshop.

run ../../../replab_init.m % Init RepLAB library

% ## Party relabeling group
% We consider a party with $m$ measurement settings with $k$ outcomes.

m = 2; k = 2;
Goutcomes = replab.S(k);
Gsettings = replab.S(m);
Gparty = Gsettings.wreathProduct(Goutcomes)

% We now define the two canonical representations of a wreath product group. For party relabeling groups, they correspond to conditional probability distributions and to deterministic strategies.

probabilityRep = Gparty.imprimitiveRep(Goutcomes.naturalRep)

strategyRep = Gparty.primitiveRep(Goutcomes.naturalRep)

% We can decompose those representations. For $m=k=2$, the probability representation has invariant vectors `[1,1,1,1]` (corresponding to overall normalization), `[1,1,-1,-1]` (corresponding to equal normalization accross settings), and an additional orthogonal space corresponding to the correlations.

probabilityRep.decomposition.nice

probabilityRep.decomposition.nice.component(1).irrep(1)

probabilityRep.decomposition.nice.component(2).irrep(1)

% We could also examine the representation on deterministic strategies; this is left to the reader.

strategyRep.decomposition.nice

% ## Bell relabeling group
% We now consider a scenario of $n$ parties with $m$ settings and $k$ outcomes (as before), and construct the representation on joint conditional probabilities. For $n=m=k=2$, we reproduce the results of [our paper](https://iopscience.iop.org/article/10.1088/1751-8121/aa6f78).

n = 2;
m = 2; k = 2;
Gparties = replab.S(n);
Gscenario = Gparties.wreathProduct(Gparty)

probabilityScenarioRep = Gscenario.primitiveRep(Gparty.imprimitiveRep(Gsettings.naturalRep))

probabilityScenarioRep.decomposition


