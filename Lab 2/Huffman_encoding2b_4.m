function [encoded_bits] = Huffman_encoding2b_4(x,n)
% null bit, used to inidcate the end of the encoding
% in the encoded bit stream
nul = -1;

% 3 because worst-case scenario is that we have to decode symbols
% that need 4 bits for their encoding(that is symbols 2 and 5)
encoded_bits = zeros(1,4*n); 

j=1;
for i=1:n
    if x(i) == 1
        encoded_bits(j) = 1;
        encoded_bits(j+1) = 1;
        encoded_bits(j+2) = 0;
        j = j + 3;
    elseif x(i) == 2
        encoded_bits(j) = 1;
        encoded_bits(j+1) = 1;
        encoded_bits(j+2) = 1;
        encoded_bits(j+3) = 1;
        j = j + 4;
    elseif x(i) == 3
        encoded_bits(j) = 0;
        j = j + 1;
    elseif x(i) == 4
        encoded_bits(j) = 1;
        encoded_bits(j+1) = 0;
        j = j + 2;
    elseif x(i) == 5
        encoded_bits(j) = 1;
        encoded_bits(j+1) = 1;
        encoded_bits(j+2) = 1;
        encoded_bits(j+3) = 0;
        j = j + 4;
    end
    
end
encoded_bits(j) = nul;

end