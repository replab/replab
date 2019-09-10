function [G rep] = clifford_qubit
% Returns the Clifford group on a single qubit and its natural representation
%
% G is a permutation group isomorphic to the Clifford group
% rep is the group natural representation on qubits
    warning('Deprecated: use replab.quantum.clifford_qudit(2)');
    [G rep] = replab.quantum.clifford_qudit(2);
end
