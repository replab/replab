function print_node(node, indent)
    node_type = class(node);
    indent_str = repmat(' ', 1, indent);
    
    if strcmp(node_type, 'replab.lobster.Root')
        fprintf('%s[Root]\n', indent_str);
    elseif strcmp(node_type, 'replab.lobster.TextNode')
        fprintf('%s[Text Node] - <%s>\n', indent_str, ...
            strrep(node.text, sprintf('\n'), ''));
    elseif strcmp(node_type, 'replab.lobster.VarNode')
        fprintf('%s[Var Node] - <{{ %s }}>\n', indent_str, node.name);
    elseif strcmp(node_type, 'replab.lobster.IfNode')
        fprintf('%s[If Node] - <Expr: %s>\n', indent_str, node.expression);
    elseif strcmp(node_type, 'replab.lobster.ElseNode')
        fprintf('%s[Else Node]\n', indent_str);
    elseif strcmp(node_type, 'replab.lobster.ForNode')
        fprintf('%s[For Node] - Expr: %s in %s>\n', indent_str, node.lhs, node.rhs);
    elseif strcmp(node_type, 'replab.lobster.CallNode')
        fprintf('%s[Call Node] - <Expr: %s>\n', indent_str, node.expression);
    else
        fprintf('%s[Unknown Node Type] <%s>\n', indent_str, class(node));
    end
    
    for k = 1:length(node.children)
        replab.lobster.print_node(node.children{k}, indent+2);
    end
end
