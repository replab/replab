function res = collapse(rep)
% Simplifies the last level of a hierarchy of representation
%
% Works for the following types:
%
% +----------------------+-----------------------+-----------------------+
% |    ``class(rep)``    | ``class(rep.parent)`` | ``class(res)``        |
% +----------------------+-----------------------+-----------------------+
% | `+replab.SimilarRep` | `+replab.SimilarRep`  | `+replab.SimilarRep`  |
% +----------------------+-----------------------+-----------------------+
% |   `+replab.SubRep`   | `+replab.SimilarRep`  | `+replab.SubRep`      |
% +----------------------+-----------------------+-----------------------+
% | `+replab.SimilarRep` |   `+replab.SubRep`    | `+replab.SubRep`      |
% +----------------------+-----------------------+-----------------------+
% |   `+replab.SubRep`   |   `+replab.SubRep`    | `+replab.SubRep`      |
% +----------------------+-----------------------+-----------------------+
%
% Args:
%   rep (`+replab.Rep`): A subrepresentation with itself and ``rep.parent`` having supported types
%
% Returns:
%   `+replab.SubRep`: A subrepresentation of ``rep.parent.parent``
    switch class(rep)
      case 'replab.SubRep'
        parent = rep.parent;
        switch class(parent)
          case 'replab.SubRep'
            newB_internal = parent.B_internal * rep.B_internal;
            newE_internal = rep.E_internal * parent.E_internal;
            res = parent.parent.subRep(newB_internal, newE_internal);
            replab.rep.copyProperties(rep, res);
          case 'replab.SimilarRep'
            newB_internal = parent.Ainv_internal * rep.B_internal;
            newE_internal = rep.E_internal * parent.A_internal;
            res = parent.parent.subRep(newB_internal, newE_internal);
            replab.rep.copyProperties(rep, res);
          otherwise
            error('Not supported');
        end
      case 'replab.SimilarRep'
        parent = rep.parent;
        switch class(rep.parent)
          case 'replab.SubRep'
            newB_internal = parent.B_internal * rep.Ainv_internal;
            newE_internal = rep.A_internal * parent.E_internal;
            res = parent.parent.subRep(newB_internal, newE_internal);
            replab.rep.copyProperties(rep, res);
          case 'replab.SimilarRep'
            newAinv_internal = parent.Ainv_internal * rep.Ainv_internal;
            newA_internal = rep.A_internal * parent.A_internal;
            res = parent.parent.similarRep(newA_internal, newAinv_internal);
            replab.rep.copyProperties(rep, res);
          otherwise
            error('Not supported');
        end
      otherwise
        error('Not supported');
    end
end
