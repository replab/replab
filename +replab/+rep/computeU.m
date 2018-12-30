function [U Uinv] = computeU(components)
    for i = 1:length(components)
        C = components{i};
        if isempty(C.U)
            U = [];
            Uinv = [];
            return
        elseif i == 1
            U = C.U;
            Uinv = C.Uinv;
        else
            U = [U C.U];
            Uinv = [Uinv; C.Uinv];
        end
    end
end
