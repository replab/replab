        function distances = sampleEquivariant(self, X, nSamples)
        % Samples random group elements, and computes the corresponding violation of equivariance
            distances = [];
            for i = 1:nSamples
                g = self.group.sample;
                gX = self.repC.action(g, self.repR.action(g, X)')';
                distances(i) = norm(X - gX, 'fro');
            end
        end

        
                function X1 = averageOver(self, X, elements)
        % Averages X over the given group elements (row cell vector)
            X1 = [];
            nElements = length(elements);
            for i = 1:nElements
                g = elements{i};
                gX = self.repC.action(g, self.repR.action(g, X)')';
                if i == 1
                    X1 = gX;
                else
                    X1 = X1 + gX;
                end
            end
            X1 = X1/nElements;
        end
        
        function X1 = randomAveraging(self, X, nSamples)
        % Computes the average of X using Monte Carlo sampling
            samples = arrayfun(@(i) self.group.sample, 1:nSamples, 'uniform', 0);
            X1 = self.averageOver(X, samples);
        end



                function X1 = averageOver(self, X, elements)
        % Averages X over the given group elements (row cell vector)
            X1 = [];
            nElements = length(elements);
            for i = 1:nElements
                g = elements{i};
                gX = self.rep.adjointAction(g, X);
                if i == 1
                    X1 = gX;
                else
                    X1 = X1 + gX;
                end
            end
            X1 = X1/nElements;
        end
        
        function X1 = randomAveraging(self, X, nSamples)
        % Computes the average of X using Monte Carlo sampling
            samples = arrayfun(@(i) self.group.sample, 1:nSamples, 'uniform', 0);
            X1 = self.averageOver(X, samples);
        end
        
        function distances = randomDistances(self, X, nSamples)
        % Samples random group elements, and computes the corresponding violation of commutativity
            distances = [];
            for i = 1:nSamples
                g = self.group.sample;
                gX = self.rep.adjointAction(g, X);
                distances(i) = norm(X - gX, 'fro');
            end
        end

