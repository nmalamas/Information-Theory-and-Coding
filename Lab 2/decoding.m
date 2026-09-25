function X_decoded = decoding(Y, n, code_book, code_lengths)
    N = size(code_book, 1);
    lenY = numel(Y);
    minL = min(code_lengths);
    maxSyms = floor(lenY / minL);

    % Preallocate output buffer
    Xbuf = false(1, maxSyms * n);

    % Precompute bit patterns for each symbol (0..N-1)
    sym_bits = false(N, n);
    for s = 0:N-1
        sym_bits(s+1, :) = bitget(s, n:-1:1);
    end

    % read pointer in Y
    idx = 1;  
    % write pointer in Xbuf
    ptr = 1;  

    % Decode loop: match each codeword
    while idx <= lenY
        matched = false;
        for sym = 1:N
            L = code_lengths(sym);
            if idx + L - 1 <= lenY
                % Compare next L bits of Y to code_book row (as row vectors)
                if all(Y(idx:idx+L-1) == code_book(sym, 1:L))
                    % Match found: write n bits
                    Xbuf(ptr:ptr+n-1) = sym_bits(sym, :);
                    ptr = ptr + n;
                    idx = idx + L;
                    matched = true;
                    break;
                end
            end
        end
        if ~matched
            error('Decoding failed at bit position %d', idx);
        end
    end

    % Return decoded bits
    X_decoded = double(Xbuf(1:ptr-1));
end
