function copyProperties(source, target)
% Copies the mutable properties of a representation to another
%
% Only the properties that can be "learned" are modified, so the object cannot become
% inconsistent; furthermore, we verify consistency with existing properties of the target.
%
% Args:
%   source (`+replab.Rep`): Source representation, is not modified
%   target (`+replab.Rep`): Target representation, mutable properties are updated
%
% Raises:
%   An error if inconsistencies are observed.
    target.isUnitary = replab.util.fusionOptional(source.isUnitary, target.isUnitary);
    target.trivialDimension = replab.util.fusionOptional(source.trivialDimension, target.trivialDimension);
    target.isIrreducible = replab.util.fusionOptional(source.isIrreducible, target.isIrreducible);
    target.frobeniusSchurIndicator = replab.util.fusionOptional(source.frobeniusSchurIndicator, target.frobeniusSchurIndicator);
    target.isDivisionAlgebraCanonical = replab.util.fusionOptional(source.isDivisionAlgebraCanonical, target.isDivisionAlgebraCanonical);
end
