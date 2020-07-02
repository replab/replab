function test_suite = CosetsTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_right_transversals
    n = 8;
    Sn = replab.S(n);
    G = replab.S(n).derivedSubgroup; % take the alternating group
    U = G;
    while U.order == G.order
        U = G.subgroup({G.sample G.sample});
    end
    TT = G.rightTransversals(U);
    UU = U.elements.toCell;
    els = zeros(n, 0);
    for t = 1:length(TT)
        for u = 1:length(UU)
            tel = TT{t};
            uel = UU{u};
            els(:,end+1) = Sn.compose(uel, tel);
        end
    end
    assert(size(unique(els', 'rows'), 1) == G.order);
end
function test_left_transversals
    n = 8;
    Sn = replab.S(n);
    G = replab.S(n).derivedSubgroup; % take the alternating group
    U = G;
    while U.order == G.order
        U = G.subgroup({G.sample G.sample});
    end
    TT = G.leftTransversals(U);
    UU = U.elements.toCell;
    els = zeros(n, 0);
    for t = 1:length(TT)
        for u = 1:length(UU)
            tel = TT{t};
            uel = UU{u};
            els(:,end+1) = Sn.compose(tel, uel);
        end
    end
    assert(size(unique(els', 'rows'), 1) == G.order);
end
