function test_suite = SignedPermTest()
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_group_laws()
    n = 100;
    cat = replab.cat.SignedPermAsGroup(n);
    randsp = @() randperm(n).*((randi(2, 1, n)-1)*2-1);
    for i = 1:100
        cat.verifyLaws(randsp);
    end
end

function test_domain_action_laws()
    n = 100;
    cat = replab.cat.SignedPermActingOnDomain(n);
    randsp = @() randperm(n).*((randi(2, 1, n)-1)*2-1);
    randel = @() randi(n).*((randi(2)-1)*2-1);    
    for i = 1:100
        cat.verifyLaws(randsp, randel);
    end
end
