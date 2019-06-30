function test_cd
% Test code for the Cayley-Dickson functions cd, dc and cdpolar.

% Copyright (c) 2008, 2010, 2013 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

disp('Testing Cayley-Dickson functions (cd, dc, cdpolar) ...')

% Fundamental tests on cd and dc.

check(dc(0,  0) == quaternion(0), 'dc function fails test 1');
check(dc(1,  0) == quaternion(1), 'dc function fails test 2');
check(dc(1i, 0) == qi,            'dc function fails test 3');
check(dc(0,  1) == qj,            'dc function fails test 4');
check(dc(0, 1i) == qk,            'dc function fails test 5');

[A, B] = cd(quaternion(0));
check(A ==  0 && B == 0, 'cd function fails test 6');

[A, B] = cd(quaternion(1));
check(A ==  1 && B == 0, 'cd function fails test 7');

[A, B] = cd(qi);
check(A == 1i && B == 0, 'cd function fails test 8');

[A, B] = cd(qj);
check(A ==  0 && B == 1, 'cd function fails test 9');

[A, B] = cd(qk);
check(A ==  0 && B == 1i,'cd function fails test 10');
 
% Test with random data. This is an exact test, because there is no
% approximation in decomposing or constructing quaternions in
% Cayley-Dickson form.

q = quaternion(randn(100), randn(100), randn(100), randn(100));

[A, B] = cd(q);

check(dc(A, B) == q, 'cd/dc functions fail test 11');

% Fundamental tests on cdpolar. These are not exact because of the logs and
% exponentials involved.

T = 1e-14;

% Test with single values 0, 1, qi, qj, qk.

cdtest(zerosq, T, '12');
cdtest(onesq,  T, '13');
cdtest(qi,     T, '14');
cdtest(qj,     T, '15');
cdtest(qk,     T, '16');

% Now test with combinations of values, and check that the reconstructed
% quaternion compares correctly with the original value.

cdtest(1 + qi + qj + qk, T, '17');
cdtest(    qi + qj + qk, T, '18');
cdtest(1 +      qj + qk, T, '19');
cdtest(1 + qi      + qk, T, '20');
cdtest(1 + qi + qj     , T, '21');
cdtest(1 + qi          , T, '22');
cdtest(1 +      qj     , T, '23');
cdtest(1 +           qk, T, '24');
cdtest(    qi      + qk, T, '25');
cdtest(         qj + qk, T, '26');
cdtest(    qi + qj     , T, '27');

% Test with random data on cdpolar. This test is not exact because it
% involves computing logs and exponentials.

Q = randq(50);

cdtest(Q,    T, '28');
cdtest(v(Q), T, '29'); % Check that it works with a pure quaternion array.

% Repeat test 28 with random sign changes in order to test the second
% (logical) parameter to cdpolar. (This causes some elements of A to be
% negated, but since the corresponding element of B will be changed to
% compensate, the result should still reconstruct correctly).

[A, B] = cdpolar(Q, randn(size(Q)) > 0); % Random logical array.
compare(Q, dc(A) .* exp(dc(B) .* qj), T, 'cdpolar function fails test 30');

% And again with false and true logical arrays ...

[A, B] = cdpolar(Q, false(size(Q))); % All false logical array.
compare(Q, dc(A) .* exp(dc(B) .* qj), T, 'cdpolar function fails test 31');

[A, B] = cdpolar(Q, true(size(Q))); % All true logical array.
compare(Q, dc(A) .* exp(dc(B) .* qj), T, 'cdpolar function fails test 32');

% Now do a few checks on octonions (more could be added).

check(dc(zerosq, zerosq) == octonion(0), 'dc function fails test 33');
check(dc(onesq,  zerosq) == octonion(1), 'dc function fails test 34');

o = rando(100);

[A, B] = cd(o);

check(dc(A, B) == o, 'cd/dc functions fail test 35');

disp('Passed');

% Function to decompose a given value and compare it with the original
% value.

function cdtest(q, T, E)

[A, B] = cdpolar(q);
compare(q, dc(A) .* exp(dc(B) .* qj), T, ['cdpolar function fails test ', E]);
end

end

% $Id: test_cd.m 1004 2017-11-15 17:14:09Z sangwine $
