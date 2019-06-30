function c = vector_product(a, b, c)
% Vector (cross) product of two or three full/pure octonions.

% Copyright (c) 2012, 2013 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(2, 3), nargoutchk(0, 1)

if ~isa(a, 'octonion') || ~isa(b, 'octonion') || (nargin == 3 && ~isa(c, 'octonion'))
    error('Vector/cross product is not defined for an octonion and a non-octonion.')
end

% Restriction removed 16 January 2019, cf F Reese Harvey Definition 6.46.
% if nargin == 2 && (~ispure(a) || ~ispure(b))
%     error('A two-fold vector product is defined only for pure octonions.')
% end

% A three-fold cross product exists in eight dimensions (see the paper by
% Massey cited below), and hence for full octonions. If any of the three
% arguments is pure, the scalar part will be treated as zero, which does
% not invalidate the three-fold cross product.

% References:
%
% F. Reese Harvey, 'Spinors and Calibrations',
% (Perspectives In mathematics: vol. 9), Academic Press, 1990.
% [See Definition 6.46, p111, and Definition 6.54, p112]
%
% W. S. Massey,
% 'Cross Products of Vectors in Higher Dimensional Euclidean Spaces',
% The American Mathematical Monthly, 90(10) (December 1983), 697-701.
%
% Peter Zvengrowski, 'A 3-Fold Vector Product in R8',
% Commentarii Mathematici Helvetici, 40 (1965-1966), 149-152.
% DOI:10.5169/seals-30632
%
% Ronald Shaw (1987): Vector cross products in n dimensions,
% International Journal of Mathematical Education in Science and Technology,
% 18:6, 803-816, DOI:10.1080/0020739870180606
%
% Eugenio Calabi, Construction and Properties of Some 6-Dimensional Almost
% Complex Manifolds, Transactions of the American Mathematical Society,
% 87(2), March 1958, 407-438.

if nargin == 2

    % TODO One day, we could implement a detailed formula here like that
    % used in the quaternion vector product, to avoid the wasted time
    % computing the unused scalar part. See Calabi where there is an
    % explicit table for computing the product in Theorem 1 (equation 2.3)
    % on page 409.
    
    if ispure(a) && ispure(b)
        % F. Reese Harvey, Definition 6.50, p111. In the definition, the
        % inner product <y, x> will be equal in magnitude but opposite in
        % sign to the scalar part of the product. Here we just take the
        % vector part of the product, rather than computing the inner
        % product in order to cancel out the scalar part. This of course
        % also has the advantage that we return a pure octonion, rather
        % than a full octonion with zero scalar part.
        
        c = vector(a .* b);
    else
        % One or both arguments is a full octonion, and we need a formula
        % that correctly includes the scalar part.
        
        c = vector(conj(b) .* a); % F. Reese Harvey, Definition 6.47, p111.
    end
else
    % For the three-fold product we use Zvengrowski's formula, Theorem 2.1,
    % but change the sign to match the result given by Harvey's formula in
    % Definition 6.54:  (a .* (conj(b) .* c) - c .* (conj(b) .* a))./2.
    
    c = a .* (conj(b) .* c) - a .* scalar_product(b, c) ...
                            + b .* scalar_product(c, a) ...
                            - c .* scalar_product(a, b);
end

% TODO Check this code, and understand the results, in particular, how to
% construct a 7-dimensional orthogonal basis, with 7 mutually perpendicular
% unit pure octonions.

% TODO See: Wikipedia: Seven-dimensional cross product where the discussion
% under Generalizations suggests that there should be a six-vector product.
% This does not seem consistent with Harvey. A critical issue is how to
% construct a 7-dimensional orthonormal basis. Surely this function is the
% key, but how should it work?

% TODO See also:
%
% M. Zorn, 'The Automorphisms of Cayley's Non-Associative Algebra',
% Proc. Nat. Acad. Sci. (1935) v21, pp.355-358.
%
% This paper shows how a 7-vector basis can be generated from a 3-vector
% orthonormal basis by adding a 4th basis element orthonormal to the first
% three, AND VECTOR PRODUCTS alone. But he does not say what the definition
% of the vector product is.

% $Id: vector_product.m 1017 2019-02-12 17:14:06Z sangwine $
