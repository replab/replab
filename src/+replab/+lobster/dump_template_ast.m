function dump_template_ast(root)
% Dumps a template AST
    if isa(root, 'replab.lobster.Template')
        root = root.root;
    end
    replab.lobster.print_node(root, 0)
end

