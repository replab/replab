function prd = of(phi)
% Semidirect product construction by selecting the most refined type
    className = replab.directproduct.commonClass({phi.G phi.P});
    switch className
      case 'replab.NiceFiniteGroup'
        % we disabled the NiceFiniteGroup construction which is
        % buggy, see issue #132
        prd = replab.semidirectproduct.OfFiniteGroups(phi);
      case 'replab.FiniteGroup'
        prd = replab.semidirectproduct.OfFiniteGroups(phi);
      case 'replab.CompactGroup'
        prd = replab.semidirectproduct.OfCompactGroups(phi);
    end
end
