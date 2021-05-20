function [gens, domainSize, order] = permutationRealizationForRelators(nGenerators, relators)
% Computes a permutation realization of a finite group given by a presentation
%
% Args:
%   nGenerators (integer): Number of generators
%   relators (cell(1,\*) of integer(1,\*)): Relators given as words in letters
%
% Returns
% -------
%   gens: cell(1,\*) of permutation
%     Permutations realizing the presentation
%   domainSize: integer
%     Domain size
%   order: integer
%     Group order
    switch replab.globals.cosetEnumerationMethod
      case 'R'
        ct = replab.fp.CosetTable.cosetEnumerationR(nGenerators, relators, {});
      case 'C'
        ct = replab.fp.CosetTable.cosetEnumerationC(nGenerators, relators, {});
      otherwise
        error('Unknown coset enumeration method ''%s''', replab.globals.cosetEnumerationMethod);
    end
    gens = arrayfun(@(i) ct.C(:,nGenerators+i)', 1:nGenerators, 'uniform', 0);
    order = size(ct.C, 1);
    domainSize = size(ct.C, 1);
end
