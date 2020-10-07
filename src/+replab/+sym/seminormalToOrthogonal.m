function diagonal = seminormalToOrthogonal(part,isCyclo)
% Find the diagonal of the diagonal change of basis vector from the seminormal irrep to the orthognal irrep
%
% Args:
%   part (integer(1,*\)): partition corresponding to irrep
%   isCyclo (boolean): If true return as replab.Cyclotomic, if false
%     return as double(1,*\)
%
% Returns:
%   diagonal (class is determined by inputs; see above):
%     diagonal of the diagonal change of basis matrix

    lat = replab.YoungLattice(part);
    diagonal = sqrt(lat.generateTabFun);
    if isCyclo
        diagonal = replab.Cyclotomic.fromDoubles(diagonal);
    end
end