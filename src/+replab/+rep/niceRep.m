function sub1 = niceRep(sub)
% Tries to recover a nicer basis for the given subrepresentation
%
% Returns a new subrepresentation sub1 if successful,
% or [] in case of failure
%
% Returns a 'replab.Irrep' instance if 'sub' is a replab.Irrep, or only a
% 'replab.SubRep' if 'sub' is only a replab.SubRep.
%
% In particular, preserves the division algebra encoding of replab.Irrep.
    assert(isa(sub, 'replab.SubRep'));
    if isa(sub, 'replab.Irrep') && isequal(sub.field, 'R') && ~sub.realDivisionAlgebra.isReal
        % TODO: try to recover a nice basis even in those cases
        sub1 = [];
        return
    end
    U0rational = replab.rep.recoverRational(sub);
    if ~isequal(U0rational, [])
        if isa(sub, 'replab.Irrep')
            sub1 = replab.Irrep(sub.parent, U0rational, sub.realDivisionAlgebra);
        else
            sub1 = replab.SubRep(sub.parent, U0rational);
        end
        return
    end
    if isequal(sub.field, 'C')
        Ureal = replab.rep.recoverReal(sub);
        if ~isequal(Ureal, [])
            if isa(sub, 'replab.Irrep')
                sub1 = replab.Irrep(sub.parent, Ureal, sub.realDivisionAlgebra);
            else
                sub1 = replab.SubRep(sub.parent, Ureal);
            end
            return
        end
    end
    sub1 = []; % fallback, no transformation possible
end
