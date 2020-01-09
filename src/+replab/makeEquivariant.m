function E = makeEquivariant(repR, repC)
% Returns the space of equivariant linear maps between two representations
%
% The equivariant vector space contains the matrices X such that
%
% self.image(g) * X = X * repC.image(g)
%
% Note: default implementations are registered in `+replab.dispatchDefaults`
%
% Args:
%   repR (replab.Rep): Representation on the target/row space
%   repC (replab.Rep): Representation on the source/column space
%
% Returns:
%   replab.Equivariant: The equivariant vector space
    E = replab.dispatch('call', 'replab.makeEquivariant', repR, repC);
end
