function sub1 = refine(sub, varargin)
% Refines a subrepresentation
%
% See `+replab.SubRep.refine`
%
% Returns:
%   `.SubRep`: Subrepresentation with refined subspace (i.e. injection/projection maps)
    if strcmp(sub.divisionAlgebraName, 'quaternion.rep')
        warning('Refine not implemented for quaternion-type representations');
        sub1 = sub;
        return
    end
    args = struct('numNonImproving', 20, 'largeScale', sub.parent.dimension > 1000, 'nSamples', 5, 'nInnerIterations', 10, 'maxIterations', 1000);
    args = replab.util.populateStruct(args, varargin);
    I = sub.injection;
    P = sub.projection;
    if sub.parent.knownUnitary
        if ~sub.mapsAreAdjoint
            [I, ~] = qr(I, 0);
        end
        if args.largeScale
            I1 = replab.rep.refine_unitaryLargeScale(sub.parent, I, sub.divisionAlgebraName, args.numNonImproving, args.nSamples, args.maxIterations, []);
        else
            I1 = replab.rep.refine_unitaryMediumScale(sub.parent, I, sub.divisionAlgebraName, args.nInnerIterations, args.maxIterations, []);
        end
        sub1 = sub.withUpdatedMaps(I1, I1', 'isUnitary', true);
    else
        if args.largeScale
            [I1, P1] = replab.rep.refine_nonUnitaryLargeScale(sub.parent, I, P, sub.divisionAlgebraName, args.numNonImproving, args.nSamples, args.maxIterations, [], []);
        else
            [I1, P1] = replab.rep.refine_nonUnitaryMediumScale(sub.parent, I, P, sub.divisionAlgebraName, args.nInnerIterations, args.maxIterations, [], []);
        end
        sub1 = sub.withUpdatedMaps(I1, P1);
    end
end
