function test_octonion_fundamentals
% Test code for the fundamental octonion functions.

% Copyright (c) 2011 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

disp('Testing octonion fundamentals ...');

% Verify the multiplication table, in its entirety. This also tests that
% the seven named octonion 'constants' are correctly defined. Note that the
% table given here by T is not used to compute octonion products - times
% and mtimes are written to make use of quaternion products and conjugates,
% so the tests here verify the code against a multiplication table that is
% independently constructed.

e0 = 1; e1 = oi; e2 = oj; e3 = ok; e4 = ol; e5 = om; e6 = on; e7 = oo;

O = [ 1, oi, oj, ok, ol, om, on, oo];
T = [e0, e1, e2, e3, e4, e5, e6, e7; ...
     e1, -1, e3,-e2, e5,-e4,-e7, e6; ...
     e2,-e3, -1, e1, e6, e7,-e4,-e5; ...
     e3, e2,-e1, -1, e7,-e6, e5,-e4; ...
     e4,-e5,-e6,-e7, -1, e1, e2, e3; ...
     e5, e4,-e7, e6,-e1, -1,-e3, e2; ...
     e6, e7, e4,-e5,-e2, e3, -1,-e1; ...
     e7,-e6, e5, e4,-e3,-e2, e1, -1];
 
% First test the scalar multiplication function (times) by explicit calls,
% element by element of the table T above.

for r = 1:8
    for c = 1:8
        check(T(r, c) == O(r) .* O(c), ...
           ['Multiplication table incorrect at ', num2str(r), ',', num2str(c)]);
    end
end

% Now test the matrix multiplication (mtimes).

check(all(O.' * O == T), 'Multiplication table is incorrect.');

% Run a check on the basic multiplications, using explicit Cayley_Dickson
% products for comparison. Note the need for parentheses because of the lack
% of associativity.

A = rando(2,3); [a, b] = cd(A); Z = zerosq(2,3);
B = rando(2,3); [c, d] = cd(B);
compare(A .* B, dc(a .* c, Z) + ...
                dc(a, Z) .* (dc(d, Z) .* ol) + ...
               (dc(b, Z) .* ol) .* dc(c, Z) + ...
               (dc(b, Z) .* ol) .* (dc(d, Z) .* ol), ...
        1e-12, 'Scalar multiplication check fails.');

A = rando(3,3); [a, b] = cd(A); Z = zerosq(3,3);
B = rando(3,3); [c, d] = cd(B);
compare(A * B, dc(a * c, Z) + ...
               dc(a, Z) * (dc(d, Z) .* ol) + ...
              (dc(b, Z) .* ol) * dc(c, Z) + ...
              (dc(b, Z) .* ol) * (dc(d, Z) .* ol), ...
        1e-12, 'Matrix multiplication check 1 fails.');
    
% Another test of multiplication, using element by element computation.

V1 = [1,2,3,4,5,6,7,8];
V2 = [9,8,7,6,5,4,3,2];

O1 = O .* V1; % This makes two vectors of octonions consisting of scaled
O2 = O .* V2; % basis elements.

P =  times(sum(O1), sum(O2));
Q = mtimes(sum(O1), sum(O2));

% Now compute the product using the table in T.

R = octonion(0);
for r = 1:8
    for c = 1:8
        R = R + (V1(r) .* V2(c)) .* T(r, c);
    end
end
    
check(P == R, 'Multiplication test 1 by elements fails.');
check(Q == R, 'Multiplication test 2 by elements fails.');

% Check the matrix product against an explicit computation using the
% elementwise product to compute the inner products of rows and columns.

P = rando(4,5);
Q = rando(5,4);
S = octonion(zeros(4,4));
for r = 1:4
    for c = 1:4
        S(r, c) = sum(P(r, :) .* transpose(Q(:, c)));
    end
end
compare(P * Q, S, 1e-12, 'Matrix multiplication check 2 fails.');

% Check that multiplication of octonions and numerics is correct.

P = octonion(1,2,3,4,5,6,7,8);
Q = octonion(3,6,9,12,15,18,21,24);

check(P .* 3 == 3 .* P, 'Octonion/numeric multiplication test 1 fails.');
check(P  * 3 == 3  * P, 'Octonion/numeric multiplication test 2 fails.');
check(P .* 3 == Q,      'Octonion/numeric multiplication test 3 fails.');

P = rando(4);
Q = randn(4);

compare(P * Q, P * octonion(Q), 1e-12, ...
    'Octonion/numeric multiplication test 4 fails.');
compare(Q * P, octonion(Q) * P, 1e-12, ...
    'Octonion/numeric multiplication test 5 fails.');

% Now run some tests using the scalar and vector products.

a = rando(10) .* randn(10);
b = rando(10) .* randn(10);
c = rando(10) .* randn(10);

P = (a .*  b) .* c;
R =  a .* (b  .* c);

Q =   s(a) .* s(b) .* s(c) ...
    + v(a) .* s(b) .* s(c) ...
    + s(a) .* v(b) .* s(c) ...
    + s(a) .* s(b) .* v(c) ...
    - s(a) .* scalar_product(v(b), v(c)) ...
    - s(b) .* scalar_product(v(a), v(c)) ...
    - s(c) .* scalar_product(v(a), v(b)) ...
    + s(a) .* vector_product(v(b), v(c)) ...
    + s(b) .* vector_product(v(a), v(c)) ...
    + s(c) .* vector_product(v(a), v(b)) ...
    - v(c) .* scalar_product(v(a), v(b)) ...
    - scalar_product(vector_product(v(a), v(b)), v(c)) ...
    + vector_product(vector_product(v(a), v(b)), v(c));

S =   s(a) .* s(b) .* s(c) ...
    + v(a) .* s(b) .* s(c) ...
    + s(a) .* v(b) .* s(c) ...
    + s(a) .* s(b) .* v(c) ...
    - s(a) .* scalar_product(v(b), v(c)) ...
    - s(b) .* scalar_product(v(a), v(c)) ...
    - s(c) .* scalar_product(v(a), v(b)) ...
    + s(a) .* vector_product(v(b), v(c)) ...
    + s(b) .* vector_product(v(a), v(c)) ...
    + s(c) .* vector_product(v(a), v(b)) ...
    - v(a) .* scalar_product(v(b), v(c)) ...
    - scalar_product(v(a), vector_product(v(b), v(c))) ...
    + vector_product(v(a), vector_product(v(b), v(c)));

compare(P, Q, 1e-12, 'Octonion scalar/vector product test 1 fails')
compare(R, S, 1e-12, 'Octonion scalar/vector product test 2 fails')

disp('Passed');

end

% $Id: test_octonion_fundamentals.m 1004 2017-11-15 17:14:09Z sangwine $
