function [code_book, code_lengths] = Shannon_codebook_generation(p, n)
    N = 2^n;
    idx = (0:N-1)';     

    % Compute probabilities for each n-bit symbol
    w = zeros(N,1);
    for k = 1:n
        w = w + bitget(idx, k);
    end
    probs = p.^w .* (1-p).^(n - w);

    % Compute code lengths: l_i = ceil(-log2(p_i))
    code_lengths = ceil(-log2(probs));

    % Sort symbols by length then by symbol value
    [sorted_lengths, sort_idx] = sort(code_lengths);
    sorted_symbols = idx(sort_idx);

    % Prepare numeric array for code book: rows = symbols, cols = max code length
    maxLen = max(sorted_lengths);
    %code_book = nan(N, maxLen);
    code_book = -1.*ones(N,maxLen);


    next_code = 0;      
    prev_length = sorted_lengths(1);

    for i = 1:N
        len = sorted_lengths(i);
        symbol = sorted_symbols(i);

        % If the code length increases, shift the current code left
        if len > prev_length
            next_code = next_code * 2^(len - prev_length);
        end

        % Create the binary codeword string then convert to numeric bits
        cw = dec2bin(next_code, len);
        bits = cw - '0'; 

        % Assign bits into code_book numeric array
        code_book(symbol+1, 1:len) = bits;

        % Increment to the next code value
        next_code = next_code + 1;
        prev_length = len;
    end
end
