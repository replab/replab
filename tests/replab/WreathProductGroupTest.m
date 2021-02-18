function test_suite = WreathProductGroupTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;

    S3 = replab.S(2);
    W = S3.wreathProduct(S3);
    test_suite = W.laws.addTestCases(test_suite);
end

function test_wreath_subgroup
    S2 = replab.S(2);
    W = S2.wreathProduct(S2);
    g = {[2 1] {[1 2] [2 1]}};
    assertEqual(W.subgroup({g}).order, vpi(4));
end

function mu = relabelings(n, m, k)
% Returns the permutation group that relabels a probability distribution/Bell inequality
%
% By convention, we assume that the permutations relabel a vectorized ``P(:)``, where ``P``
% is a multi-dimensional array ``P(a1,a2,...,an,x1,x2,...xn)``, where all ``ai = 1,...,k``
% and all ``xi = 1,..,m`.
%
% Args:
%   n (integer): Number of parties
%   m (integer): Number of inputs
%   k (integer): Number of outputs
%
% Returns:
%   `+replab.FiniteIsomorphism`: Isomorphism from a structured to a permutation group that represents all relabelings
    ind = reshape(1:(k*m)^n, [ones(1,n)*k ones(1,n)*m]);
    order = [n:-1:1
             n+(n:-1:1)];
    ind = permute(ind, order(:)');
    ind = ind(:)';
    d = length(ind);
    Sd = replab.S(d);
    generators = cell(1, 0);
    Sparties = replab.S(n);
    Sinputs = replab.S(m);
    Soutputs = replab.S(k);
    W1 = Sinputs.wreathProduct(Soutputs); % single party
    W = Sparties.wreathProduct(W1);
    mu = W.isomorphismByFunction(replab.S(d), [], @(w) Sd.leftConjugate(ind, W.primitivePermutation(w, @(w1) W1.imprimitivePermutation(w1))));
end

function test_relabelings
% Test the construction of relabelings in Bell scenarios, an application of primitive/imprimitive permutation representations
    n = 2;
    m = 2;
    k = 2;
    mu = relabelings(n, m, k);
    % We construct the relabeling action on P(a,b,x,y) by hand to verify the construction
    Iabxy = reshape(1:16, [2 2 2 2]);
    IPermParties = reshape(permute(Iabxy, [2 1 4 3]), 1, []);
    IPermAliceInput = Iabxy(:,:,[2 1],:);
    IPermAliceOutputForFirstInput = Iabxy;
    IPermAliceOutputForFirstInput(:,:,1,:) = IPermAliceOutputForFirstInput([2 1],:,1,:);

    % The corresponding elements of the double wreath product
    wPermParties = {[2, 1], { ... % Permutation of parties
        {[1, 2], ... % Permutation of Alice inputs
         {[1, 2], [1, 2]}}, ... % Permutation of Alice outputs for x=1 and x=2
        {[1, 2], ... % Permutation of Bob inputs
         {[1, 2], [1, 2]}} ... % Permutation of Bob outputs for y=1 and y=2
                   }};
    wPermAliceInput = {[1, 2], { ...
        {[2, 1], {[1, 2], [1, 2]}}, ...
        {[1, 2], {[1, 2], [1, 2]}}
                   }};
    wPermAliceOutputForFirstInput = {[1, 2], { ...
        {[1, 2], {[2, 1], [1, 2]}}, ...
        {[1, 2], {[1, 2], [1, 2]}}
                   }};
    % Verification
    assert(all(mu.imageElement(wPermParties) == reshape(IPermParties, 1, [])));
    assert(all(mu.imageElement(wPermAliceInput) == reshape(IPermAliceInput, 1, [])));
    assert(all(mu.imageElement(wPermAliceOutputForFirstInput) == reshape(IPermAliceOutputForFirstInput, 1, [])));
end
