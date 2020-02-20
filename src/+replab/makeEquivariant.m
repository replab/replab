function E = makeEquivariant(repC, repR, special)
% Returns the space of equivariant linear maps between two representations
%
% The equivariant vector space contains the matrices X such that
%
% ``repC.image(g) * X = X * repR.image(g)``
%
% Note: default implementations are registered in `+replab.dispatchDefaults`
%
% Args:
%   repC (`+replab.Rep`): Representation on the source/column space
%   repR (`+replab.Rep`): Representation on the target/row space
%   special (charstring): Special structure see help on `+replab.Equivariant.special`
%
% Returns:
%   `+replab.Equivariant`: The equivariant vector space
    E = replab.dispatch('call', 'replab.makeEquivariant', repC, repR, special);
end
