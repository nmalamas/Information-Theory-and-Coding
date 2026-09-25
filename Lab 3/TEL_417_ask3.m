%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Information Theory and Coding - TEL 417
% Exercise 3
% Patsea Ioanna-Georgia - 2020030033
% Malamas Nikolaos - 2020030180
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;close all; clear all;

% 1
p_err = 0.1;
Rates = [ 1 1/2 1/3 1/4 1/5 1/8 ];
no_bin_codes = 10^4;
no_codewords = 10;

figure; hold on;

for rIdx = 1:length(Rates)
    R = Rates(rIdx);
    Nlist = (1/R : 1/R : 10/R);
    BER = zeros(1,length(Nlist));

    for nIdx = 1:length(Nlist)
        N = Nlist(nIdx);
        k = R*N;  
        bit_errors = 0;
        total_bits = 0;

        for i = 1:no_bin_codes
            % create the codebook on the fly
            tmp_messages = round(rand(2^k,k));
            tmp_codebook = round(rand(2^k,N));

            for j = 1:no_codewords               
                % choose a random codeword from the
                % codebook and transmit it
                idx = randi(2^k);
                message = tmp_messages(idx, :);
                tmp_codeword = tmp_codebook(idx,:);

                % create the p_err codebook
                p_codeword = rand(1, N) < p_err;
                received_codeword = xor(tmp_codeword,p_codeword);

                % Calculate the Hamming distance for each codeword
                distances = abs(tmp_codebook - received_codeword);
                hammingDistances = sum(distances, 2);
                
                % Store the minimum distance found
                [minDistance, final_idx] = min(hammingDistances);
                
                % find the most likely sent codeword
                % based on the ML rule
                ML_word = tmp_messages(final_idx,:);
             
                bit_errors = bit_errors + sum(message ~= ML_word);
                total_bits = total_bits + k;
                
            end
        end

        BER(nIdx) = bit_errors / total_bits;
    end
    
    % plot the BER vs N curve for this rate
    semilogy( Nlist, BER, ...
          'LineWidth'  , 1.5, ...
          'DisplayName', sprintf('R=%.2f', R) );
       
end

% finish up
xlabel('Block length N');
ylabel('Average BER');
title(sprintf('BER vs N for various rates, p_{err}=%.2f', p_err));
legend('show','Location','southwest');
grid on;