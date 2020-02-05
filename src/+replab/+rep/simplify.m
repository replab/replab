function rep1 = simplify(rep)
% Returns a representation equivalent to the given one, but whose computation is simpler
    newRep = [];
    switch class(rep)
      case 'replab.rep.DerivedRep'
        parent = rep.parent;
        switch class(parent)
          case 'replab.rep.DerivedRep'
            % collapse successive DerivedRep
            newRep = replab.rep.DerivedRep(parent.parent, ...
                                           xor(rep.conjugate, parent.conjugate), ...
                                           xor(rep.inverse, parent.inverse), ...
                                           xor(rep.transpose, parent.transpose));
          case 'replab.rep.ComplexifiedRep'
            % move DerivedRep inside ComplexifiedRep where it can be simplified
            newRep = replab.rep.ComplexifiedRep(...
                replab.rep.DerivedRep(parent.parent, rep.conjugate, rep.inverse, rep.transpose));
          otherwise
            if ~rep.conjugate && ~rep.inverse && ~rep.transpose
                % if DerivedRep has no effect, remove it
                newRep = parent;
            elseif isequal(rep.field, 'R') && rep.conjugate
                % Complex conjugation has no effect on reals
                newRep = replab.rep.DerivedRep(parent, false, rep.inverse, rep.transpose);
            elseif isequal(rep.isUnitary, true) && rep.inverse && rep.transpose
                % Replace dual by conjugate when unitary
                newRep = replab.rep.DerivedRep(parent, ~rep.conjugate, false, false);
            end
        end
      otherwise
    end
    if isempty(newRep)
        rep1 = rep;
    else
        % recurse
        rep1 = replab.rep.simplify(newRep);
    end
end
