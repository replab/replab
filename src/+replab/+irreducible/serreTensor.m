function [T, P, n] = serreTensor(rep, irrep, type)
% Returns, for a given representation, the tensor encoding the status of an irrep in a representation
%
% Args:
%   rep (`+replab.Rep`): Representation of dimension ``D``
%   irrep (`+replab.Rep`): Irreducible representation
%   type ('exact', 'double', 'double/sparse'): Type of the projector to compute
%
% Returns
% -------
%   T: double(D,D,d,d) or `+replab.cyclotomic`(D,D,d,d), may be sparse
%     Serre tensor
%   P: double(D,D) or `+replab.cyclotomic`(D,D), may be sparse
%     Projector into the isotypic component
%   n: integer
%     Returns the multiplicity of the given irreducible character
    assert(isa(rep, 'replab.Rep'));
    assert(rep.field == irrep.field);
    assert(rep.group == irrep.group);
    assert(irrep.isIrreducible);
    D = rep.dimension;
    d = irrep.dimension;
    t = rep.group.trivialRep(rep.field, d*D);
    K = kron(dual(irrep), rep);
    if nargin < 3 || isempty(type)
        type = 'double';
    end
    if strcmp(type, 'exact');
        T = K.equivariantFrom(t).project(replab.cyclotomic.eye(d*D), 'exact');
    else
        T = K.equivariantFrom(t).project(eye(d*D), 'double');
    end
    T = permute(reshape(T, [D d D d]), [1 3 2 4]);
    P = reshape(T, [D*D d*d]) * reshape(eye(d), [d*d 1]);
    P = reshape(P, [D D]);
    n = trace(P)/trace(P*P);
    T = T * n;
    P = P * n;
end
