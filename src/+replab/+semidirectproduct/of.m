function prd = of(phi)
% Semidirect product construction by selecting the most refined type
    className = replab.directproduct.commonClass({phi.G phi.P});
    switch className
      case 'replab.NiceFiniteGroup'
        prd = replab.semidirectproduct.OfNiceFiniteGroups(phi);
      case 'replab.FiniteGroup'
        prd = replab.semidirectproduct.OfFiniteGroups(phi);
      case 'replab.CompactGroup'
        prd = replab.semidirectproduct.OfCompactGroups(phi);
    end
end
