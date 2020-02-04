function E = makeEquivariant(repR, repC, special)
% Returns the space of equivariant linear maps between two representations
%
% The equivariant vector space contains the matrices X such that
%
% ``repR.image(g) * X = X * repC.image(g)``
%
% Note: default implementations are registered in `+replab.dispatchDefaults`
%
% Args:
%   repR (`+replab.Rep`): Representation on the target/row space
%   repC (`+replab.Rep`): Representation on the source/column space
%   sepcial (charstring): Special structure see help on `+replab.Equivariant.special`
%
% Returns:
%   `+replab.Equivariant`: The equivariant vector space
    E = replab.dispatch('call', 'replab.makeEquivariant', repR, repC, special);
end
