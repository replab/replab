function test_exp
% Test code for the quaternion exponential function. This also tests the
% axis, modulus, and unit functions as well as many others.

% Copyright (c) 2005, 2013 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

disp('Testing exponential function ...')

T = 1e-12;

% Test 1. Real quaternion data.

q = randq(100) .* randn(100);

compare(q, abs(q) .* exp(axis(q) .* angle(q)), T, 'quaternion/exp failed test 1.');

% Test 2. Complex quaternion data.

b = complex(q, randq(100) .* randn(100));

compare(b, abs(b) .* exp(axis(b) .* angle(b)), T, 'quaternion/exp failed test 2.');

% Test 3. Imaginary quaternion data.

b = imag(b);

compare(b, abs(b) .* exp(axis(b) .* angle(b)), T, 'quaternion/exp failed test 3.');

% Test 4. Real octonion data.

q = rando(100) .* randn(100);

compare(q, abs(q) .* exp(axis(q) .* angle(q)), T, 'octonion/exp failed test 4.');

% Test 5. Complex octonion data.

b = complex(q, rando(100) .* randn(100));

compare(b, abs(b) .* exp(axis(b) .* angle(b)), T, 'octonion/exp failed test 5.');

% Test 6. Imaginary octonion data.

b = imag(b);

compare(b, abs(b) .* exp(axis(b) .* angle(b)), T, 'octonion/exp failed test 6.');

% Now test out the handling of nilpotents and idempotents, both alone and
% within an array of 'normal' values. The first two tests are done with the
% matrix exponential of an adjoint matrix to check the algorithm is OK.

nilq  = qi + qj .* 1i;

compare(exp(nilq), unadjoint(expm(adjoint(nilq, 'real')), 'real'), T, ...
                                          'quaternion/exp failed test 7.');
idem = (1/2) + (randv .* 1i) ./ 2;
                                      
compare(exp(idem), unadjoint(expm(adjoint(idem, 'real')), 'real'), T, ...
                                          'quaternion/exp failed test 9.');
                                      
nilo = octonion(0, 1, 0, 0, 1i, 0, 0, 0);

compare(exp(nilo), 1 + nilo, T, 'octonion/exp failed test 9.');

% Do not know how to do an idem test on octonions, as we can't use a
% matrix. TODO: devise one. It will be test 10, which is why we omit 10.

% Now test the handling of nilpotents and idempotents inside a larger
% array. The point here is to test the indexing.

q = randq(5);  o = rando(5);

q(3,4) = 0 + nilq; o(2,3) = nilo;
q(4,3) = idem;

expq = exp(q); expo = exp(o);

D = isdivisor(q); % Identify the divisors of zero.
E = isdivisor(o);

compare(expq(~D), exp(q(~D)), T, 'quaternion/exp failed test 11.');
compare(expo(~E), exp(o(~E)), T, 'octonion/exp failed test 12.');
compare(expq(D),  exp(q(D)),  T, 'quaternion/exp failed test 13.');
compare(expo(E),  exp(o(E)),  T, 'octonion/exp failed test 14.');

disp('Passed');

% $Id: test_exp.m 1031 2019-04-20 12:49:25Z sangwine $
