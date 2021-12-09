function warmup
% Loads and uses a few key RepLAB classes to warmup Octave
%
% In Octave 6.1.0 and 6.3.0, we noticed that running ``replab_runtests`` right after ``replab_init``
% led to test failures (due to Octave mishandling ``Access=protected`` rules).
%
% However, loading and using the affected classes before running the test suite corrects the problem.
%
% This is what we do here.
    G = replab.PermutationGroup.of([2 1 3 4], [1 2 4 3]);
    G.abelianInvariants;
    h = [2 3 4 1];
    L = G.leftCoset(h);
    r = L.representative;
end
