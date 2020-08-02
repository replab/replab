function gens = permutationGeneratorsForRelators(nGenerators, relators)
% Computes a permutation realization of a finite group given by a presentation
%
% Args:
%   nGenerators (integer): Number of generators
%   relators (cell(1,\*) of integer(1,\*)): Relators given as words in letters
%
% Returns:
%   cell(1,\*) of permutation: Generators of a permutation group realizing the presentation
    ct = replab.fp.cosetEnumerationR(nGenerators, relators, {}, 2^50);
    gens = arrayfun(@(i) ct.C(:,nGenerators+i)', 1:nGenerators, 'uniform', 0);
end
