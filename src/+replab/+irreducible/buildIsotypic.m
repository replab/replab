function iso = buildIsotypic(rep, samples, sub)
% Builds an isotypic canonical component from equivalent subrepresentations
%
% As part of the operation, we identify the division algebra type (if the representations are over R),
% put those division algebras in their canonical basis, and make the irreducible
% representations not only equivalent but identical.
%
% Args:
%   rep (replab.Rep): Real representation being decomposed
%   samples (replab.irreducible.OnDemandSamples): Lazy evaluation of various samples for `rep`
%   sub (row cell array of replab.SubRep): Equivalent irreducible subrepresentations of `rep`
%
% Returns:
%   replab.Isotypic: The isotypic component
    first = sub{1};
    if rep.overR
        % identifies the type of real representations and 
        % make the division algebra basis canonical        
        realType = replab.irreducible.computeRealType(rep, samples, first);
        switch realType
          case 'R'
            ii = replab.IrrepInfo(first.irrepInfo.label, realType, []);
            first = rep.subRepUnitary(first.U, first.niceBasis, ii);
          case 'C'
            W = replab.irreducible.enforceComplexEncoding(rep, samples, first);
            ii = replab.IrrepInfo(first.irrepInfo.label, realType, true);
            first = rep.subRepUnitary(W*first.U, [], ii);
          case 'H'
            W = replab.irreducible.enforceQuaternionEncoding(rep, samples, first);
            ii = replab.IrrepInfo(first.irrepInfo.label, realType, true);
            first = rep.subRepUnitary(W*first.U, [], ii);
        end
    end
    n = length(sub);
    copies = cell(1, n);
    copies{1} = first;
    for j = 2:n
        other = sub{j};
        W = replab.irreducible.findCommonBasis(rep, samples, first, other);
        copies{j} = rep.subRepUnitary(W * other.U, [], first.irrepInfo);
    end
    iso = replab.Isotypic(rep, copies);
end
