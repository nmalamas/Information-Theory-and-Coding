function [code_book, code_lengths] = Huffman_codebook_generation(p, n)

    % 1) Compute probabilities for all n-bit symbols
    N = 2^n;
    symbols = dec2bin(0:N-1, n);
    probs = zeros(N,1);
    for i = 1:N
        count1 = sum(symbols(i,:) == '1');
        probs(i) = p^count1 * (1-p)^(n-count1);
    end

    % 2) Initialize leaf nodes with symbol index and probability
    nodes = struct('left', {}, 'right', {}, 'prob', {}, 'symbolIdx', {});
    for i = 1:N
        nodes(i).left      = [];
        nodes(i).right     = [];
        nodes(i).prob      = probs(i);
        nodes(i).symbolIdx = i;
    end

    % 3) Build the Huffman tree
    while numel(nodes) > 1
        [~, idx_sort] = sort([nodes.prob]);
        nodes         = nodes(idx_sort);
        newNode.left      = nodes(1);
        newNode.right     = nodes(2);
        newNode.prob      = nodes(1).prob + nodes(2).prob;
        newNode.symbolIdx = [];
        nodes = [nodes(3:end), newNode];
    end
    root = nodes;

    % 4) Generate codes via tree traversal (non-recursive)
    temp_book    = cell(N,1);
    code_lengths = zeros(N,1);
    stack = struct('node', root, 'code', '');
    while ~isempty(stack)
        item = stack(end);
        stack(end) = [];
        nd = item.node;
        cd = item.code;
        if isempty(nd.left) && isempty(nd.right)
            idx = nd.symbolIdx;
            temp_book{idx}    = cd;
            code_lengths(idx) = length(cd);
        else
            % push right then left so left gets '0' prefix first
            stack(end+1) = struct('node', nd.right, 'code', [cd '1']);
            stack(end+1) = struct('node', nd.left,  'code', [cd '0']);
        end
    end

    % 5) Convert temp_book to numeric array: rows = symbols, cols = max code length
    maxLen = max(code_lengths);
    %code_book = nan(N, maxLen);
    code_book = -1.*ones(N,maxLen);

    for i = 1:N
        cw = temp_book{i};           
        bits = cw - '0';              
        code_book(i, 1:length(bits)) = bits;
    end
end
