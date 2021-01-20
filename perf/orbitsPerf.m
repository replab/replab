% Here we do some speed performance analysis on several algorithms
% concerning the speed of orbits computation when having a large number of
% elements.
%
% The test case consists in the symmetrization of a N x N matrix under the
% joint action of S(N) on the lines and columns.


% The functions to be compared:
functions = {@(pairs) replab.graph.burningAlgorithm(pairs), ...
             @(pairs) replab.graph.burningAlgorithmFast(pairs)};


nbFunctions = length(functions);
functionsName = cell(1, nbFunctions);
for i = 1:nbFunctions
    functionsName{i} = char(functions{i});
end


log10Ns = 1:0.2:4;
Ns = round(10.^log10Ns);
nbNs = length(Ns);
timingsPreparation = [];
timings = zeros(1,nbNs);
co = 0;
for N = Ns
    co = co + 1;
    disp([num2str(co), '/', num2str(nbNs)]);
    
    % Preparing the data
    M = reshape(1:N^2, N, N);

    gens = {[2 1 3:N], [2:N 1]};

    % Further preparation
    tic;
    pairs = {};
    for i = 1:length(gens)
        pairs{i} = [[1:N^2].', reshape(M(gens{i}, gens{i}), N^2, 1)];
    end
    pairs = cat(1, pairs{:});

    %pairs = unique(pairs, 'rows');
    
    timingsPreparation(co) = toc;
    
    %% We test each function
    for i = 1:length(functions)
        timings(:,1:co)
        if (co == 1) || (max(timings(i, 1:co-1)) < 7)
            tic;
            subsets = functions{i}(pairs); 
            timings(i, co) = toc;
        end    
    end
end


%% Plot the comparison result
semilogy(Ns, timings);
xlabel('size N');
ylabel('time (s)');
legend(functionsName{:});
