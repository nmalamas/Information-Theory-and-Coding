function Y = encoding(X, n, code_book, code_lengths)
    M = length(X);
    % Reshape into blocks
    X_blocks = reshape(X, n, []).';

    % Convert blocks to indices 0..2^n-1 (MSB first)
    weights = 2.^(n-1:-1:0);
    indices = X_blocks * weights.';  % column vector

    % Preallocate Y length
    total_bits = sum(code_lengths(indices + 1));
    Y = false(1, total_bits);

    ptr = 1;
    for i = 1:length(indices)
        idx = indices(i) + 1;       
        len = code_lengths(idx);
        % Extract numeric code bits
        cw = code_book(idx, 1:len);
        % Write to Y
        Y(ptr:ptr+len-1) = cw;
        ptr = ptr + len;
    end
end
