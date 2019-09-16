function test_suite = CommutantVarTest()
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
    
    if ReplabTestParameters.onlyFastTests
        % We create just once two CommutantVar objects to be used by all
        % the tests of this class

        disp('Creating few simple CommutantVar objects');
        global matrix231 matrix23451 matrix23451H
        matrix231 = replab.CommutantVar.fromPermutations({[2 3 1]});
        matrix23451 = replab.CommutantVar.fromPermutations({[2 3 4 5 1]});
        matrix23451H = replab.CommutantVar.fromSdpMatrix(sdpvar(5,5,'hankel'), {[2 3 4 5 1]});
    end
end

function test_fromPermutations
    % We do a sanity check with one case
    global matrix231 matrix23451 matrix23451H
    generators = {[2 3 4 5 1]};
    matrix = matrix23451;
    fullMatrix = matrix.fullMatrix;
    for i = 1:length(generators)
        difference = fullMatrix - fullMatrix(generators{i}, generators{i});
        vars = getvariables(difference);
        for j = 1:length(vars)
            coeffs = getbasematrix(difference, vars(j));
            assert(norm(coeffs(:)) <= replab.Settings.doubleEigTol);
        end
    end
end

function test_fromSdpMatrix_SDP_CHSH
    indexMatrix = [  1   2   3   6   7   8  11  12  13
                     2   1   4   7   6   9  12  11  14
                     3   4   1   8   9   6  13  14  11
                     6   7   8   1   2   3  16  17  18
                     7   6   9   2   1   4  17  16  19
                     8   9   6   3   4   1  18  20  16
                    11  12  13  16  17  18   1   2   3
                    12  11  14  17  16  20   2   1   4
                    13  14  11  18  19  16   3   4   1];

    objective = [0 0 0 0 1 1 0 1 -1];

    generators = {[1  4  7  2  5  8  3  6  9]
                  [1 -2 -3 -4  5  6 -7  8  9]
                  [1  3  2  4  6  5 -7 -9 -8]}';
	
    % Check of the generators:
    isok = true;
    for i = 1:length(generators)
        if norm(reshape(sign(generators{i}).*objective(abs(generators{i}))-objective, 3, 3)) ~= 0
            isok = false;
        end
    end
    assert(isok, 'The symetry generators don''t leave the objective function invariant');

	% We formulate the non-symmetrized SDP:
    vars = [0; sdpvar(max(max(indexMatrix)),1)];
    tmp = sparse(1:numel(indexMatrix), reshape(1+indexMatrix,1,numel(indexMatrix)), true);
    sdpMatrix = reshape(tmp*vars, size(indexMatrix));
    obj = objective*sdpMatrix(:,1);
    if ReplabTestParameters.onlyFastTests
        obj1 = 2*sqrt(2);
    else
        solvesdp([sdpMatrix >= 0, sdpMatrix(1,1) == 1], -obj, sdpsettings('verbose', 0));
        obj1 = value(obj);
    end
    
    % We formulate the symmetrized SDP:
    symSdpMatrix = replab.CommutantVar.fromSdpMatrix(sdpMatrix,generators);
    solvesdp([symSdpMatrix >= 0, symSdpMatrix(1,1) == 1], -obj, sdpsettings('verbose', 0));
    obj2 = value(obj);
    
    % We compare the result:
    assert(abs(obj1 - obj2)/abs(obj1) < replab.Settings.doubleSdpTol, 'Symmetrized SDP doesn''t yield the same result as the non-symmetrized one');
end

function test_fromSdpMatrix_SDP_CHSH_FullProb
    if ReplabTestParameters.onlyFastTests
        return;
    end
    
    indexMatrix = [  1    2    4    3    5   14   15   17   16   18   40   41   43   42   44   27   28   30   29   31   53   54   56   55   57
                     2    2    0    6    7   15   15    0   19   20   41   41    0   45   46   28   28    0   32   33   54   54    0   58   59
                     4    0    4    9   11   17    0   17   22   24   43    0   43   48   50   30    0   30   35   37   56    0   56   61   63
                     3    6    9    3    0   16   19   22   16    0   42   45   48   42    0   29   32   35   29    0   55   58   61   55    0
                     5    7   11    0    5   18   20   24    0   18   44   46   50    0   44   31   33   37    0   31   57   59   63    0   57
                    14   15   17   16   18   14   15   17   16   18    0    0    0    0    0   66   67   69   68   70   79   80   82   81   83
                    15   15    0   19   20   15   15    0   19   20    0    0    0    0    0   67   67    0   71   72   80   80    0   84   85
                    17    0   17   22   24   17    0   17   22   24    0    0    0    0    0   69    0   69   75   76   82    0   82   88   89
                    16   19   22   16    0   16   19   22   16    0    0    0    0    0    0   68   73   74   68    0   81   86   87   81    0
                    18   20   24    0   18   18   20   24    0   18    0    0    0    0    0   70   77   78    0   70   83   90   91    0   83
                    40   41   43   42   44    0    0    0    0    0   40   41   43   42   44  105  106  108  107  109  131  132  134  133  135
                    41   41    0   45   46    0    0    0    0    0   41   41    0   45   46  106  106    0  112  116  132  132    0  136  137
                    43    0   43   48   50    0    0    0    0    0   43    0   43   48   50  108    0  108  113  117  134    0  134  140  141
                    42   45   48   42    0    0    0    0    0    0   42   45   48   42    0  107  110  114  107    0  133  138  139  133    0
                    44   46   50    0   44    0    0    0    0    0   44   46   50    0   44  109  111  115    0  109  135  142  143    0  135
                    27   28   30   29   31   66   67   69   68   70  105  106  108  107  109   27   28   30   29   31    0    0    0    0    0
                    28   28    0   32   33   67   67    0   73   77  106  106    0  110  111   28   28    0   32   33    0    0    0    0    0
                    30    0   30   35   37   69    0   69   74   78  108    0  108  114  115   30    0   30   35   37    0    0    0    0    0
                    29   32   35   29    0   68   71   75   68    0  107  112  113  107    0   29   32   35   29    0    0    0    0    0    0
                    31   33   37    0   31   70   72   76    0   70  109  116  117    0  109   31   33   37    0   31    0    0    0    0    0
                    53   54   56   55   57   79   80   82   81   83  131  132  134  133  135    0    0    0    0    0   53   54   56   55   57
                    54   54    0   58   59   80   80    0   86   90  132  132    0  138  142    0    0    0    0    0   54   54    0   58   59
                    56    0   56   61   63   82    0   82   87   91  134    0  134  139  143    0    0    0    0    0   56    0   56   61   63
                    55   58   61   55    0   81   84   88   81    0  133  136  140  133    0    0    0    0    0    0   55   58   61   55    0
                    57   59   63    0   57   83   85   89    0   83  135  137  141    0  135    0    0    0    0    0   57   59   63    0   57];

    objective = [0 0 0 0 0 0 1 -1 1 -1 0 -1 1 -1 1 0 1 -1 -1 1 0 -1 1 1 -1];

    generators = {[1   6  11  16  21   2   7  12  17  22   3   8  13  18  23   4   9  14  19  24   5  10  15  20  25]
                  [1   3   2   5   4  11  13  12  15  14   6   8   7  10   9  21  23  22  25  24  16  18  17  20  19]
                  [1   4   5   2   3   6   9  10   7   8  11  14  15  12  13  21  24  25  22  23  16  19  20  17  18]}';

    % Check of the generators:
    isok = true;
    for i = 1:length(generators)
        if norm(reshape(sign(generators{i}).*objective(abs(generators{i}))-objective, 5, 5)) ~= 0
            isok = false;
        end
    end
    assert(isok, 'The symetry generators don''t leave the objective function invariant');

    % non-symmetrized SDP solving:
    vars = [0; sdpvar(max(max(indexMatrix)),1)];
    tmp = sparse(1:numel(indexMatrix), reshape(1+indexMatrix,1,numel(indexMatrix)), true);
    sdpMatrix = reshape(tmp*vars, size(indexMatrix));
    obj = objective*sdpMatrix(:,1);
    solvesdp([sdpMatrix >= 0, sdpMatrix(1,1) == 1], -obj, sdpsettings('verbose', 0));
    obj1 = value(obj);
    
    % We formulate the symmetrized SDP:
    symSdpMatrix = replab.CommutantVar.fromSdpMatrix(sdpMatrix, generators);
    solvesdp([symSdpMatrix >= 0, symSdpMatrix(1,1) == 1], -obj, sdpsettings('verbose', 0));
    obj2 = value(obj);
    
    % We compare the result:
    assert(abs(obj1 - obj2)/abs(obj1) < replab.Settings.doubleSdpTol, 'Symmetrized SDP doesn''t yield the same result as the non-symmetrized one');
end

function test_fromSdpMatrix_SDP_CGLMP3_FullProb
    if ReplabTestParameters.onlyFastTests
        return;
    end
    
    indexMatrix = [   1    2    4    6    3    5    7   26   27   29   31   28   30   32   76   77   79   81   78   80   82  126  127  129  131  128  130  132   51   52   54   56   53   55   57  101  102  104  106  103  105  107  151  152  154  156  153  155  157
                      2    2    0    0    8    9   10   27   27    0    0   33   34   35   77   77    0    0   83   84   85  127  127    0    0  133  134  135   52   52    0    0   58   59   60  102  102    0    0  108  109  110  152  152    0    0  158  159  160
                      4    0    4    0   12   15   16   29    0   29    0   37   40   41   79    0   79    0   87   90   91  129    0  129    0  137  140  141   54    0   54    0   62   65   66  104    0  104    0  112  115  116  154    0  154    0  162  165  166
                      6    0    0    6   13   19   22   31    0    0   31   38   44   47   81    0    0   81   88   94   97  131    0    0  131  138  144  147   56    0    0   56   63   69   72  106    0    0  106  113  119  122  156    0    0  156  163  169  172
                      3    8   12   13    3    0    0   28   33   37   38   28    0    0   78   83   87   88   78    0    0  128  133  137  138  128    0    0   53   58   62   63   53    0    0  103  108  112  113  103    0    0  153  158  162  163  153    0    0
                      5    9   15   19    0    5    0   30   34   40   44    0   30    0   80   84   90   94    0   80    0  130  134  140  144    0  130    0   55   59   65   69    0   55    0  105  109  115  119    0  105    0  155  159  165  169    0  155    0
                      7   10   16   22    0    0    7   32   35   41   47    0    0   32   82   85   91   97    0    0   82  132  135  141  147    0    0  132   57   60   66   72    0    0   57  107  110  116  122    0    0  107  157  160  166  172    0    0  157
                     26   27   29   31   28   30   32   26   27   29   31   28   30   32    0    0    0    0    0    0    0    0    0    0    0    0    0    0  176  177  179  181  178  180  182  201  202  204  206  203  205  207  226  227  229  231  228  230  232
                     27   27    0    0   33   34   35   27   27    0    0   33   34   35    0    0    0    0    0    0    0    0    0    0    0    0    0    0  177  177    0    0  183  184  185  202  202    0    0  208  209  210  227  227    0    0  233  234  235
                     29    0   29    0   37   40   41   29    0   29    0   37   40   41    0    0    0    0    0    0    0    0    0    0    0    0    0    0  179    0  179    0  189  190  191  204    0  204    0  214  215  216  229    0  229    0  239  240  241
                     31    0    0   31   38   44   47   31    0    0   31   38   44   47    0    0    0    0    0    0    0    0    0    0    0    0    0    0  181    0    0  181  195  196  197  206    0    0  206  220  221  222  231    0    0  231  245  246  247
                     28   33   37   38   28    0    0   28   33   37   38   28    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0  178  186  187  188  178    0    0  203  211  212  213  203    0    0  228  236  237  238  228    0    0
                     30   34   40   44    0   30    0   30   34   40   44    0   30    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0  180  192  193  194    0  180    0  205  217  218  219    0  205    0  230  242  243  244    0  230    0
                     32   35   41   47    0    0   32   32   35   41   47    0    0   32    0    0    0    0    0    0    0    0    0    0    0    0    0    0  182  198  199  200    0    0  182  207  223  224  225    0    0  207  232  248  249  250    0    0  232
                     76   77   79   81   78   80   82    0    0    0    0    0    0    0   76   77   79   81   78   80   82    0    0    0    0    0    0    0  276  277  279  281  278  280  282  351  352  354  356  353  355  357  376  377  379  381  378  380  382
                     77   77    0    0   83   84   85    0    0    0    0    0    0    0   77   77    0    0   83   84   85    0    0    0    0    0    0    0  277  277    0    0  286  292  298  352  352    0    0  358  359  360  377  377    0    0  383  384  385
                     79    0   79    0   87   90   91    0    0    0    0    0    0    0   79    0   79    0   87   90   91    0    0    0    0    0    0    0  279    0  279    0  287  293  299  354    0  354    0  364  365  366  379    0  379    0  389  390  391
                     81    0    0   81   88   94   97    0    0    0    0    0    0    0   81    0    0   81   88   94   97    0    0    0    0    0    0    0  281    0    0  281  288  294  300  356    0    0  356  370  371  372  381    0    0  381  395  396  397
                     78   83   87   88   78    0    0    0    0    0    0    0    0    0   78   83   87   88   78    0    0    0    0    0    0    0    0    0  278  283  289  295  278    0    0  353  361  362  363  353    0    0  378  386  387  388  378    0    0
                     80   84   90   94    0   80    0    0    0    0    0    0    0    0   80   84   90   94    0   80    0    0    0    0    0    0    0    0  280  284  290  296    0  280    0  355  367  368  369    0  355    0  380  392  393  394    0  380    0
                     82   85   91   97    0    0   82    0    0    0    0    0    0    0   82   85   91   97    0    0   82    0    0    0    0    0    0    0  282  285  291  297    0    0  282  357  373  374  375    0    0  357  382  398  399  400    0    0  382
                    126  127  129  131  128  130  132    0    0    0    0    0    0    0    0    0    0    0    0    0    0  126  127  129  131  128  130  132  301  302  304  306  303  305  307  451  452  454  456  453  455  457  526  527  529  531  528  530  532
                    127  127    0    0  133  134  135    0    0    0    0    0    0    0    0    0    0    0    0    0    0  127  127    0    0  133  134  135  302  302    0    0  311  317  323  452  452    0    0  461  467  473  527  527    0    0  533  534  535
                    129    0  129    0  137  140  141    0    0    0    0    0    0    0    0    0    0    0    0    0    0  129    0  129    0  137  140  141  304    0  304    0  312  318  324  454    0  454    0  462  468  474  529    0  529    0  539  540  541
                    131    0    0  131  138  144  147    0    0    0    0    0    0    0    0    0    0    0    0    0    0  131    0    0  131  138  144  147  306    0    0  306  313  319  325  456    0    0  456  463  469  475  531    0    0  531  545  546  547
                    128  133  137  138  128    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0  128  133  137  138  128    0    0  303  308  314  320  303    0    0  453  458  464  470  453    0    0  528  536  537  538  528    0    0
                    130  134  140  144    0  130    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0  130  134  140  144    0  130    0  305  309  315  321    0  305    0  455  459  465  471    0  455    0  530  542  543  544    0  530    0
                    132  135  141  147    0    0  132    0    0    0    0    0    0    0    0    0    0    0    0    0    0  132  135  141  147    0    0  132  307  310  316  322    0    0  307  457  460  466  472    0    0  457  532  548  549  550    0    0  532
                     51   52   54   56   53   55   57  176  177  179  181  178  180  182  276  277  279  281  278  280  282  301  302  304  306  303  305  307   51   52   54   56   53   55   57    0    0    0    0    0    0    0    0    0    0    0    0    0    0
                     52   52    0    0   58   59   60  177  177    0    0  186  192  198  277  277    0    0  283  284  285  302  302    0    0  308  309  310   52   52    0    0   58   59   60    0    0    0    0    0    0    0    0    0    0    0    0    0    0
                     54    0   54    0   62   65   66  179    0  179    0  187  193  199  279    0  279    0  289  290  291  304    0  304    0  314  315  316   54    0   54    0   62   65   66    0    0    0    0    0    0    0    0    0    0    0    0    0    0
                     56    0    0   56   63   69   72  181    0    0  181  188  194  200  281    0    0  281  295  296  297  306    0    0  306  320  321  322   56    0    0   56   63   69   72    0    0    0    0    0    0    0    0    0    0    0    0    0    0
                     53   58   62   63   53    0    0  178  183  189  195  178    0    0  278  286  287  288  278    0    0  303  311  312  313  303    0    0   53   58   62   63   53    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
                     55   59   65   69    0   55    0  180  184  190  196    0  180    0  280  292  293  294    0  280    0  305  317  318  319    0  305    0   55   59   65   69    0   55    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
                     57   60   66   72    0    0   57  182  185  191  197    0    0  182  282  298  299  300    0    0  282  307  323  324  325    0    0  307   57   60   66   72    0    0   57    0    0    0    0    0    0    0    0    0    0    0    0    0    0
                    101  102  104  106  103  105  107  201  202  204  206  203  205  207  351  352  354  356  353  355  357  451  452  454  456  453  455  457    0    0    0    0    0    0    0  101  102  104  106  103  105  107    0    0    0    0    0    0    0
                    102  102    0    0  108  109  110  202  202    0    0  211  217  223  352  352    0    0  361  367  373  452  452    0    0  458  459  460    0    0    0    0    0    0    0  102  102    0    0  108  109  110    0    0    0    0    0    0    0
                    104    0  104    0  112  115  116  204    0  204    0  212  218  224  354    0  354    0  362  368  374  454    0  454    0  464  465  466    0    0    0    0    0    0    0  104    0  104    0  112  115  116    0    0    0    0    0    0    0
                    106    0    0  106  113  119  122  206    0    0  206  213  219  225  356    0    0  356  363  369  375  456    0    0  456  470  471  472    0    0    0    0    0    0    0  106    0    0  106  113  119  122    0    0    0    0    0    0    0
                    103  108  112  113  103    0    0  203  208  214  220  203    0    0  353  358  364  370  353    0    0  453  461  462  463  453    0    0    0    0    0    0    0    0    0  103  108  112  113  103    0    0    0    0    0    0    0    0    0
                    105  109  115  119    0  105    0  205  209  215  221    0  205    0  355  359  365  371    0  355    0  455  467  468  469    0  455    0    0    0    0    0    0    0    0  105  109  115  119    0  105    0    0    0    0    0    0    0    0
                    107  110  116  122    0    0  107  207  210  216  222    0    0  207  357  360  366  372    0    0  357  457  473  474  475    0    0  457    0    0    0    0    0    0    0  107  110  116  122    0    0  107    0    0    0    0    0    0    0
                    151  152  154  156  153  155  157  226  227  229  231  228  230  232  376  377  379  381  378  380  382  526  527  529  531  528  530  532    0    0    0    0    0    0    0    0    0    0    0    0    0    0  151  152  154  156  153  155  157
                    152  152    0    0  158  159  160  227  227    0    0  236  242  248  377  377    0    0  386  392  398  527  527    0    0  536  542  548    0    0    0    0    0    0    0    0    0    0    0    0    0    0  152  152    0    0  158  159  160
                    154    0  154    0  162  165  166  229    0  229    0  237  243  249  379    0  379    0  387  393  399  529    0  529    0  537  543  549    0    0    0    0    0    0    0    0    0    0    0    0    0    0  154    0  154    0  162  165  166
                    156    0    0  156  163  169  172  231    0    0  231  238  244  250  381    0    0  381  388  394  400  531    0    0  531  538  544  550    0    0    0    0    0    0    0    0    0    0    0    0    0    0  156    0    0  156  163  169  172
                    153  158  162  163  153    0    0  228  233  239  245  228    0    0  378  383  389  395  378    0    0  528  533  539  545  528    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0  153  158  162  163  153    0    0
                    155  159  165  169    0  155    0  230  234  240  246    0  230    0  380  384  390  396    0  380    0  530  534  540  546    0  530    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0  155  159  165  169    0  155    0
                    157  160  166  172    0    0  157  232  235  241  247    0    0  232  382  385  391  397    0    0  382  532  535  541  547    0    0  532    0    0    0    0    0    0    0    0    0    0    0    0    0    0  157  160  166  172    0    0  157];

    objective = [0  0  0  0  0  0  0  0 -1  0  1 -1  0  1  0  0  1 -1  1 -1  0  0  1 -1  0  0  1 -1  0 -1  1  0 -1  1  0  0  0 -1  1  1  0 -1  0  1  0 -1  0 -1  1];

    generators = {[1   8  15  22  29  36  43   2   9  16  23  30  37  44   3  10  17  24  31  38  45   4  11  18  25  32  39  46   5  12  19  26  33  40  47   6  13  20  27  34  41  48   7  14  21  28  35  42  49]
                  [1   3   4   2   7   5   6  22  24  25  23  28  26  27   8  10  11   9  14  12  13  15  17  18  16  21  19  20  36  38  39  37  42  40  41  43  45  46  44  49  47  48  29  31  32  30  35  33  34]
                  [1   5   6   7   2   3   4   8  12  13  14   9  10  11  22  26  27  28  23  24  25  15  19  20  21  16  17  18  29  33  34  35  30  31  32  43  47  48  49  44  45  46  36  40  41  42  37  38  39]}';

    % Check of the generators:
    isok = true;
    for i = 1:length(generators)
        if norm(reshape(sign(generators{i}).*objective(abs(generators{i}))-objective, 7, 7)) ~= 0
            isok = false;
        end
    end
    assert(isok, 'The symetry generators don''t leave the objective function invariant');

    % non-symmetrized SDP solving:
    vars = [0; sdpvar(max(max(indexMatrix)),1)];
    tmp = sparse(1:numel(indexMatrix), reshape(1+indexMatrix,1,numel(indexMatrix)), true);
    sdpMatrix = reshape(tmp*vars, size(indexMatrix));
    obj = objective*sdpMatrix(:,1);
    solvesdp([sdpMatrix >= 0, sdpMatrix(1,1) == 1], -obj, sdpsettings('verbose', 0));
    obj1 = value(obj);
    
    % We formulate the symmetrized SDP:
    symSdpMatrix = replab.CommutantVar.fromSdpMatrix(sdpMatrix, generators);
    solvesdp([symSdpMatrix >= 0, symSdpMatrix(1,1) == 1], -obj, sdpsettings('verbose', 0));
    obj2 = value(obj);
    
    % We compare the result:
    assert(abs(obj1 - obj2)/abs(obj1) < replab.Settings.doubleSdpTol, 'Symmetrized SDP doesn''t yield the same result as the non-symmetrized one');
end

function test_fromSymSdpMatrix
    indexMatrix = [  1   2   3   2   6   7   3   7  10
                     2   1   4   6   2   8   7   3  11
                     3   4   1   7   8   2  10  11   3
                     2   6   7   1   2   3   4   8  11
                     6   2   8   2   1   4   8   4  13
                     7   8   2   3   4   1  11  14   4
                     3   7  10   4   8  11   1   2   3
                     7   3  11   8   4  14   2   1   4
                    10  11   3  11  13   4   3   4   1];

    objective = [0 0 0 0 1 1 0 1 -1];

    generators = {[1  4  7  2  5  8  3  6  9]};

    % Non-block-diagonal SDP
    vars = [0; sdpvar(max(max(indexMatrix)),1)];
    tmp = sparse(1:numel(indexMatrix), reshape(1+indexMatrix,1,numel(indexMatrix)), true);
    sdpMatrix = reshape(tmp*vars, size(indexMatrix));
    obj = objective*sdpMatrix(:,1);
    if ReplabTestParameters.onlyFastTests
        obj1 = 2*sqrt(2);
    else
        solvesdp([sdpMatrix >= 0, sdpMatrix(1,1) == 1], -obj, sdpsettings('verbose', 0));
        obj1 = value(obj);
    end
    
    % We do a sanity check with one group
    blockSdpMatrix = replab.CommutantVar.fromSymSdpMatrix(sdpMatrix, generators);
    solvesdp([blockSdpMatrix >= 0, blockSdpMatrix(1,1) == 1], -obj, sdpsettings('verbose', 0));
    obj2 = value(obj);
    
    % We compare the result:
    assert(abs(obj1 - obj2)/abs(obj1) < replab.Settings.doubleSdpTol, 'Block-diagonalized SDP doesn''t yield the same result as the non-block-diagonalized one');
end

function test_inputs
    shouldProduceAnError(@(x) replab.CommutantVar.fromPermutations([]));
end
