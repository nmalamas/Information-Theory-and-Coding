function [recovered_symbols] = Huffman_decoding1(encoded_bits,n)
% null bit, used to inidcate the end of the encoding
% in the encoded bit stream
nul = -1;

recovered_symbols = zeros(1,n);

i = 1;
j = 1;
while encoded_bits(i) ~= nul
    if encoded_bits(i) == 1
        recovered_symbols(1,j) = 1;
        i = i + 1;
        j = j + 1;
    elseif encoded_bits(i) == 0 & encoded_bits(i+1) == 1
        recovered_symbols(1,j) = 2;
        i = i + 2;
        j = j + 1;
    elseif encoded_bits(i) == 0 & encoded_bits(i+1) == 0 & encoded_bits(i+2) == 0
        recovered_symbols(1,j) = 3;
        i = i + 3;
        j = j + 1;
    elseif encoded_bits(i) == 0 & encoded_bits(i+1) == 0 & encoded_bits(i+2) == 1
        recovered_symbols(1,j) = 4;
        i = i + 3;
        j = j + 1;
    end      
end

end