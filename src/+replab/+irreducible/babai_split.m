function [allReps, epsilon] = babai_split(rep, rep_eps)
% Runs the Babai-Friedl split algorithm
%
%   -- WORK IN PROGRESS --
%
% c.f. "Approximate Representation Theory of Finite Groups" by L. Babai and
% K. Friedl.
%
% Args:
%   rep (replab.Rep, replab.SubRep): Approximate representation to split
%   rep_eps (double): Precision of the approximate representation rep
%
% Returns:
% --------
%   allReps: replab.Rep{1,*}
%     List of alpha-irreducible representations
%   epsilon: double
%     precision of the approximate decomposition
%
% Example:
%   >>> [allReps, epsilon] = replab.irreducible.babai_split(replab.S(6).definingRep, eps_gen);

if nargin < 2
    rep_eps = 0;
end

% Whether to display information along the way
verbose = true;

% The dimension of the considered representation
m = rep.dimension;

% We estimate the diameter of the group with what happens for S(n)
diameter = m^2;

% Choosing the alpha parameter. We need to choose it such that:
%   2*alpha*m^(5/2) + rep_eps < 1/(18*diameter*sqrt(m))
alphaMax = m^(-5)/36 - rep_eps*m^(-5/2)/2;
% Ideally, we would like alpha to be such that we can obtain some overall
% precision targetEpsilon
targetEpsilon = 1e-4;
expectedNbIrreps = 10; % wild guess of the expected number of irreducible blocks
alphaIdeal = (targetEpsilon/(12*expectedNbIrreps^2*diameter*m) - rep_eps)*m^(-5/2)/2;
if (alphaIdeal > alphaMax)
    alpha = alphaMax;
elseif (alphaIdeal < 0)
    alpha = min(1e-7, alphaMax);
else
    alpha = alphaIdeal;
end

% define further parameters
nbGen = length(rep.group.generators);
criticalGap = 1/(2*m^(5/2));
nbH = m^2;

% A basis for hermitian matrices
Hs = cell(1,m^2);
co = 1;
for i = 1:m
    Hii = sparse(i, i, 1, m, m);
    Hs{co} = Hii;
    co = co + 1;
    for j = i+1:m
        Hij = sparse([i j], [j i], [1 1], m, m);
        Hji = sparse([i j], [j i], [1i -1i], m, m);
        Hs{co} = Hij;
        Hs{co+1} = Hji;
        co = co + 2;
    end
end

% We precompute the image of the generators
imageGen = cell(1,nbGen);
for k = 1:nbGen
    gen = rep.group.generator(k);
    imageGen{k} = rep.image(gen);
end

% To keep track of the improvement
maxMaxNorms = zeros(1,10000);

% The main loop
for i = 1:10000
    if verbose
        disp(['Iteration number ', num2str(i)]);
    end
    
    % Compute the Dixon average
    for j = 1:nbH
        W = Hs{j};
        sigma = 0*W;
        for k = 1:nbGen
            sigma = sigma + imageGen{k}*W*imageGen{k}';
            sigma = sigma + imageGen{k}'*W*imageGen{k};
        end
        Hs{j} = sigma/(2*nbGen);
    end

    % We check the alpha-centralizing condition
    maxNorm = zeros(1,length(Hs));
    for j = 1:nbH
        for k = 1:nbGen
            maxNorm(j) = max(maxNorm(j), norm(Hs{j}*imageGen{k} - imageGen{k}*Hs{j}, 'fro'));
        end
    end
    if verbose
        disp(['Element of maxNorm are in the interval [', num2str(min(maxNorm)), ', ', num2str(max(maxNorm)), ']']);
    end
    centralizes = (maxNorm <= alpha);
    maxMaxNorms(i) = max(maxNorm);
    if sum(centralizes) == nbH
        if verbose
            disp(['All H are centralized after ', num2str(i), ' steps. "stably irreducible"']);
        end
        allReps = {rep};
        return;
    end
    which = find(centralizes);
    for j = which
        % We compute the eigenvalue gap
        vp = sort(real(eig(Hs{j})));
        if max(diff(vp)) >= criticalGap
            % We cut after the first gap goes too far
            [a, b] = eig(Hs{j});
            b = real(diag(b));
            [b, I] = sort(b);
            a = a(:,I);
            cut = find(diff(vp) >= criticalGap, 1, 'first');
            U1 = a(:,1:cut);
            U2 = a(:,cut+1:end);
            
            if verbose
                disp(['The space splits after ', num2str(i), ' steps. Overlaps:']);
                for k = 1:nbGen
                    U1'*imageGen{k}*U2
                end
            end
            
            rep1 = rep.subRepUnitary(U1');
            rep2 = rep.subRepUnitary(U2');
            
            % We further decompose
            allReps1 = babai_algo(rep1);
            allReps2 = babai_algo(rep2);
            allReps = cat(2, allReps1, allReps2);
            
            % bound the error
            delta = 2*alpha*m^(5/2);
            %[delta + rep_eps, 1/(18*diameter*sqrt(m))]
            epsilon = 12*length(allReps)^2*diameter*m*(rep_eps+delta);
            return;
        end
    end
end

if verbose
    % We did not converge, so we plot the convergence to see what happened
    plot(log10(maxMaxNorms))
end
error('Not enough steps')

