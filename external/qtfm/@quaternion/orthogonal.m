function U = orthogonal(V, W)
% ORTHOGONAL(V) constructs unit pure quaternion(s) perpendicular to V, and
% perpendicular to W, if given. V (and W if given) must be pure quaternion
% arrays. Elements of W need not be perpendicular to elements of V, but
% they must not be parallel (that is the cross product of W and V must not
% be zero or close to it). This is an elementwise function: each element of
% U will be orthogonal to the correspond element of V.

% Copyright (c) 2005, 2010, 2017 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

% TODO This function can be generalized to work in 4-space using a very
% simple algorithm. Multiplying a quaternion q by an arbitrary unit pure
% quaternion gives a result p orthogonal to q. It is also possible that the
% result could be a set of three quaternions orthogonal to V and each
% other, but it may be better to do this as a separate function cf
% orthonormal_basis. To be considered.
%
% Reference:
%
% L. Meister and H. Schaeben, 'A concise quaternion geometry of rotations',
% Mathematical Methods in the Applied Sciences, 2005; 28(1), 101-126.
% DOI: 10.1002/mma.560.
% (See Proposition 1, page 104.)

narginchk(1, 2), nargoutchk(0, 1)

if ~isa(V, 'quaternion') || ~isempty(V.w)
    error('First argument must be a pure quaternion.')
end

if nargin == 2
    
    if ~isa(W, 'quaternion') || ~isempty(W.w)
        error('Second argument must be a pure quaternion.')
    end
    
    % The two arguments must be the same size, or one must be scalar.
    
    if ~(all(size(V) == size(W)) || isscalar(V) || isscalar(W))
        error('Parameters must be the same size, or one must be a scalar.');
    end

    if any(parallel(V(:), W(:)))
        error('Arguments must not have elements parallel to each other.');    
    end
end

% Now compute the required pure quaternion.  The method used was published
% in:
% 
% Ell, T. A. and Sangwine, S. J., 'Decomposition of 2D Hypercomplex Fourier
% Transforms into Pairs of Complex Fourier Transforms', in: Gabbouj, M. and
% Kuosmanen (eds), 'Signal Processing X, Theories and Applications',
% Proceedings of EUSIPCO 2000, Tenth European Signal Processing Conference,
% Tampere, Finland, 5-8 September 2000, Vol. II, 1061-1064, section 6.
% 
% In the paper, a specific case was shown. Here we have a general method.
%
% Subsequent to the coding of this function, and the publication of the
% paper above by Ell and Sangwine, it was discovered that the algorithm had
% been published in 1957 in section 5 of:
%
% J. K. Mackenzie and M. J. Thomson,  'Some Statistics Associated with the
% Random Disorientation of Cubes', Biometrika, 44(1-2), June 1957, 205-210.

if nargin == 1
        
    % W was not given, so we must choose an arbitrary vector not parallel
    % to any element of V.
    
    if isscalar(V)
        W = choose(V);
    else
        % V is an array, and the simplest way to handle this is to loop
        % over all the elements, choosing W for each case. This will result
        % in different choices for the various W only if V has elements
        % from the set {qi, qj, qk}, otherwise the choices will be
        % identical for all elements of V.
        
        W = zerosv(size(V));
        for j = 1:numel(V)
            % W(j) = choose(V(j));
            S = substruct('()', {j});
            W = subsasgn(W, S, choose(subsref(V, S)));
        end
    end
    
end

% Now compute the result.  We have to normalise this, because W may not be
% perpendicular to V, and therefore the result may not have unit modulus,
% even if W did. (This is why we do not normalise the values of W above:
% the normalization is done here to allow for the lack of perpendicularity,
% and it can therefore also adjust for the lack of unit modulus above.) It
% is possible that the call to unit( ) will fail if the cross product is a
% nilpotent biquaternion (whether this can happen is not known).

C = vector_product(V, W); % Formerly called cross, cuts out one call.
U = unit(C);

end

function W = choose(V)
% Given a SCALAR V, chooses a scalar W not parallel to V using a careful
% choice originally implemented for the non-vectorised version of the main
% function. We deal specially with the cases where elements of V are in the
% set {qi, qj, qk} in order to return a sensible choice. For example, if
% V == qi, we choose W so that the return result of the main function will
% be qj.

assert(isscalar(V));

L = parallel(V, [qi, qj, qk]);
if any(L)
    
    % V is parallel to one of the standard basis vectors. We choose one
    % of the other two, cyclically. The negative signs ensure that, for
    % example, orthogonal(qi) yields +qj, as you would expect.
    
    if L(1)
        W = -qk;
    elseif L(2)
        W = -qi;
    else
        W = -qj;
    end
else
    
    % V is not parallel to one of the standard basis vectors.
    
    W = quaternion(1, 1, 1); % The default, unless W == V.
    if parallel(V, W)
        % We have to make another choice in this case, somewhat
        % arbitrarily.
        W = quaternion(-1, -1, 0);
    end
end
end

function tf = parallel(x, y)
% Returns true if x and y are parallel or 'close' to it. x and/or y may be
% complex pure quaternions, which is why there are two calls on abs -- to
% ensure that the result is real before comparison with the threshold. Note
% that this test may give a falsely false result if x and/or y is a nilpotent
% biquaternion (that is has a vanishing modulus) - meaning that the return
% result of the function is false, even though x is not in fact parallel to
% y. This is because of the vanishing modulus of the nilpotent.
tf = abs(abs(vector_product(x, y))) < eps;
end

% Implementation notes: this function must be deterministic, so that it
% returns the same result for a given V every time it is called, otherwise
% the various fft functions that depend on this function would need to be
% changed to pass the orthornormal basis as a parameter. This means that we
% cannot choose W at random if only V is given.
%
% The case where V is a complex quaternion is treated the same way, but
% there can be problems with nilpotents. 

% $Id: orthogonal.m 1004 2017-11-15 17:14:09Z sangwine $
