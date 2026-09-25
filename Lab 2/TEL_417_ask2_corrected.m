%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Information Theory and Coding - TEL 417
% Exercise 2
% Patsea Ioanna-Georgia -  2020030033
% Malamas Nikolaos - 2020030180
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;close all;clear all;

%%
% -- Exercise 1 -- 
n = 100;
no_experiments = 1000;
% null bit, used to inidcate the end of the encoding
% in the encoded bit stream
nul = -1;

fprintf('========== Exercise 1 ==========\n');

% 1.1
% The calculation below is for verification, 
% compared to the theoretical one derived, in the report
H_X = -( (1/2)*log2(1/2)+(1/4)*log2(1/4)+(1/8)*log2(1/8)+(1/8)*log2(1/8) );

% 1.3
% 3 because worst-case scenario is that we have to decode symbols
% that need 3 bits for their encoding(that is symbols 3 and 4)
experiments_encoded_bits = zeros(no_experiments,3*n); 
experiments_initial_symbols = zeros(no_experiments,n);

% The Huffman encoding we will use is
% 1->1, 2->01, 3->000 and 4->001
for k=1:no_experiments
    x = generate_iid_symbols(n);
    experiments_initial_symbols(k,:) = x;

    experiments_encoded_bits(k,:) = Huffman_encoding1(x,n);
end

% 1.4
% The Huffman decoding we will use is
% 1->1, 01->2, 000->3 and 001->4
experiments_recovered_symbols = zeros(no_experiments,n);

for k=1:no_experiments
    experiments_recovered_symbols(k,:) = Huffman_decoding1(experiments_encoded_bits(k,:),n);
end

% Check whether or not decoding was successful
if experiments_recovered_symbols(:,:) == experiments_initial_symbols(:,:)
    fprintf('The recovered symbols, indeed match the initial ones!\n');
    fprintf('\n');
else
    fprintf('There are symbols NOT correctly recovered! Decoding process unsuccessful.\n');
    fprintf('\n');
end

% Now, calculate the numerical average over the
% average word length of every x
avg_word_lengths = zeros(1,no_experiments);

for k=1:no_experiments
    num_of_bits = 0;
    i = 1;
    while experiments_encoded_bits(k,i) ~= nul
        num_of_bits = num_of_bits + 1;
        i = i + 1;
    end
    avg_word_lengths(k) = num_of_bits/n;
end

num_avg_word_length = sum(avg_word_lengths)/no_experiments;

disp('Comparing the numerical average of encoded word length with the entropy of the source');
fprintf('Numerical average: %.4f bits/symb\n', num_avg_word_length);
fprintf('Source Entropy: %.4f bits/symb\n', H_X);
fprintf('\n');

% 1.5
d_5 = zeros(1,no_experiments);
d_6 = zeros(1,no_experiments);
d_56_pairs = zeros(2,no_experiments);

for k=1:no_experiments
    if experiments_encoded_bits(k,5)~=nul
        d_5(1,k) = experiments_encoded_bits(k,5);
    end
    if experiments_encoded_bits(k,6)~=nul
        d_6(1,k) = experiments_encoded_bits(k,6);
    end
    d_56_pairs(1,k) = d_5(1,k);
    d_56_pairs(2,k) = d_6(1,k);
end

disp('By counting the #of appearances of each bit we can obtain their distribution:');
% Find the marginal distributions of d_5 and d_6
d_5_eq0 = sum(d_5(1,:)==0);
d_5_eq1 = sum(d_5(1,:)==1);
fprintf('#(d_5 = 0): %d - #(d_5 = 1): %d\n', d_5_eq0, d_5_eq1);

d_6_eq0 = sum(d_6(1,:)==0);
d_6_eq1 = sum(d_6(1,:)==1);
fprintf('#(d_6 = 0): %d - #(d_6 = 1): %d\n', d_6_eq0, d_6_eq1);

% As well as the distribution of the (d_5,d_6) pairs
d_56_eq00 = sum(d_56_pairs(1,:)==0 & d_56_pairs(2,:)==0);
d_56_eq01 = sum(d_56_pairs(1,:)==1 & d_56_pairs(2,:)==0);
d_56_eq10 = sum(d_56_pairs(1,:)==0 & d_56_pairs(2,:)==1);
d_56_eq11 = sum(d_56_pairs(1,:)==1 & d_56_pairs(2,:)==1);
fprintf('#[(d_5,d_6) = 00]: %d\n', d_56_eq00);
fprintf('#[(d_5,d_6) = 01]: %d\n', d_56_eq01);
fprintf('#[(d_5,d_6) = 10]: %d\n', d_56_eq10);
fprintf('#[(d_5,d_6) = 11]: %d\n', d_56_eq11);
fprintf('\n');

%%
% -- Exercise 2 -- 
n = 100;
no_experiments = 1000;

fprintf('========== Exercise 2 ==========\n');

% 2.1
P = [1/2 1/8 1/8 1/8 1/8;
    1/4 1/8 1/16 1/16 1/2;
    1/4 1/8 1/8 1/4 1/4;
    1/8  0  1/2 1/4 1/8;
     0 1/2  1/4 1/4  0];

[Q,Lamda] = eig(P);
lamda_idx = find(Lamda(:,:)==1);
q_lamda1 = Q(:,mod(lamda_idx,5));

% Find the normalized eigenvector of P
% (that is, the convergence distribution)
p_inf = q_lamda1/(sum(q_lamda1(:)));

% Find the Entropy of the convergence distribution
H_cd = -sum(p_inf(:).*log2(p_inf(:)));

% Find the Rate of Entropy
H_P = zeros(1,length(P));
for i=1:length(P)
    tmp = 0;
    for j=1:length(P(:,i))
        
        if P(j,i) == 0
            continue;
        end
        tmp = tmp + P(j,i)*log2(P(j,i));
    end
    H_P(i) = -tmp;
end

R_entr = H_P*p_inf;

% 2.2
p_X0 = [0.1 0.25 0.35 0.15 0.15];
%simulate_MC(n,p_X0,P);

% 2.3
n = 100;
no_experiments = 1000;

% (a)
experiments_encoded_bits_a = zeros(1,3*n);
experiments_initial_symbols_a = zeros(no_experiments,n);

for k=1:no_experiments
    x = simulate_MC(n,p_inf,P);
    experiments_initial_symbols_a(k,:) = x;

    experiments_encoded_bits_a(k,:) = Huffman_encoding2a(x,n);
end

% (b)
% b.1
experiments_encoded_bits_b1 = zeros(1,3*n);
experiments_initial_symbols_b1 = zeros(no_experiments,n);

for k=1:no_experiments
    x = simulate_MC(n,p_inf,P);
    experiments_initial_symbols_b1(k,:) = x;

    experiments_encoded_bits_b1(k,:) = Huffman_encoding2b_1(x,n);
end


% b.2
experiments_encoded_bits_b2 = zeros(1,4*n);
experiments_initial_symbols_b2 = zeros(no_experiments,n);

for k=1:no_experiments
    x = simulate_MC(n,p_inf,P);
    experiments_initial_symbols_b2(k,:) = x;
    
    experiments_encoded_bits_b2(k,:) = Huffman_encoding2b_2(x,n);
end

% b.3
experiments_encoded_bits_b3 = zeros(1,3*n);
experiments_initial_symbols_b3 = zeros(no_experiments,n);

for k=1:no_experiments
    x = simulate_MC(n,p_inf,P);
    experiments_initial_symbols_b3(k,:) = x;
    
    experiments_encoded_bits_b3(k,:) = Huffman_encoding2b_3(x,n);
end

% b.4
experiments_encoded_bits_b4 = zeros(1,4*n);
experiments_initial_symbols_b4 = zeros(no_experiments,n);

for k=1:no_experiments
    x = simulate_MC(n,p_inf,P);
    experiments_initial_symbols_b4(k,:) = x;
    
    experiments_encoded_bits_b4(k,:) = Huffman_encoding2b_4(x,n);
end

% b.5
experiments_encoded_bits_b5 = zeros(1,4*n);
experiments_initial_symbols_b5 = zeros(no_experiments,n);

for k=1:no_experiments
    x = simulate_MC(n,p_inf,P);
    experiments_initial_symbols_b5(k,:) = x;
    
    experiments_encoded_bits_b5(k,:) = Huffman_encoding2b_5(x,n);
end

% 2.4
% (a)
% The Huffman decoding we will use is
% 10->1, 11->2, 000->3, 001->4 and 01->5
experiments_recovered_symbols_a = zeros(no_experiments,n);

for k=1:no_experiments
    experiments_recovered_symbols_a(k,:) = Huffman_decoding2a(experiments_encoded_bits_a(k,:),n);
end

fprintf('2.4(a): \n');

% Check whether or not decoding was successful
if experiments_recovered_symbols_a(:,:) == experiments_initial_symbols_a(:,:)
    fprintf('-The recovered symbols, indeed match the initial ones!\n');
    fprintf('\n');
else
    fprintf('-There are symbols NOT correctly recovered! Decoding process unsuccessful.\n');
    fprintf('\n');
end

% Now, calculate the numerical average over the
% average word length of every x
avg_word_lengths = zeros(1,no_experiments);

for k=1:no_experiments
    num_of_bits = 0;
    i = 1;
    while experiments_encoded_bits_a(k,i) ~= nul
        num_of_bits = num_of_bits + 1;
        i = i + 1;
    end
    avg_word_lengths(k) = num_of_bits/n;
end

num_avg_word_length = sum(avg_word_lengths)/no_experiments;

disp('-Comparing the numerical average of the encoded word average length with the entropy of the source');
fprintf('Numerical average: %.4f bits/symb\n', num_avg_word_length);
fprintf('Source Entropy: %.4f bits/symb\n', H_cd);
fprintf('Entropy Rate: %.4f bits/symb\n', R_entr);
fprintf('\n');


% (b)
% b.1
% The Huffman decoding we will use is
% 0->1, 100->2, 101->3, 110->4 and 111->5
experiments_recovered_symbols_b1 = zeros(no_experiments,n);

for k=1:no_experiments
    experiments_recovered_symbols_b1(k,:) = Huffman_decoding2b_1(experiments_encoded_bits_b1(k,:),n);
end

fprintf('2.4(b1): \n');

% Check whether or not decoding was successful
if experiments_recovered_symbols_b1(:,:) == experiments_initial_symbols_b1(:,:)
    fprintf('-The recovered symbols, indeed match the initial ones!\n');
    fprintf('\n');
else
    fprintf('-There are symbols NOT correctly recovered! Decoding process unsuccessful.\n');
    fprintf('\n');
end

% Now, calculate the numerical average over the
% average word length of every x
avg_word_lengths = zeros(1,no_experiments);

for k=1:no_experiments
    num_of_bits = 0;
    i = 1;
    while experiments_encoded_bits_b1(k,i) ~= nul
        num_of_bits = num_of_bits + 1;
        i = i + 1;
    end
    avg_word_lengths(k) = num_of_bits/n;
end

num_avg_word_length = sum(avg_word_lengths)/no_experiments;

disp('-Comparing the numerical average of the encoded word average length with the entropy of the source');
fprintf('Numerical average: %.4f bits/symb\n', num_avg_word_length);
fprintf('Source Entropy: %.4f bits/symb\n', H_cd);
fprintf('Entropy Rate: %.4f bits/symb\n', R_entr);
fprintf('\n');

% b.2
% The Huffman decoding we will use is
% 0->5, 10->1, 110->2, 1110->3 and 1111->4
experiments_recovered_symbols_b2 = zeros(no_experiments,n);

for k=1:no_experiments
    experiments_recovered_symbols_b2(k,:) = Huffman_decoding2b_2(experiments_encoded_bits_b2(k,:),n);
end

fprintf('2.4(b2): \n');

% Check whether or not decoding was successful
if experiments_recovered_symbols_b2(:,:) == experiments_initial_symbols_b2(:,:)
    fprintf('-The recovered symbols, indeed match the initial ones!\n');
    fprintf('\n');
else
    fprintf('-There are symbols NOT correctly recovered! Decoding process unsuccessful.\n');
    fprintf('\n');
end

% Now, calculate the numerical average over the
% average word length of every x
avg_word_lengths = zeros(1,no_experiments);

for k=1:no_experiments
    num_of_bits = 0;
    i = 1;
    while experiments_encoded_bits_b2(k,i) ~= nul
        num_of_bits = num_of_bits + 1;
        i = i + 1;
    end
    avg_word_lengths(k) = num_of_bits/n;
end

num_avg_word_length = sum(avg_word_lengths)/no_experiments;

disp('-Comparing the numerical average of the encoded word average length with the entropy of the source');
fprintf('Numerical average: %.4f bits/symb\n', num_avg_word_length);
fprintf('Source Entropy: %.4f bits/symb\n', H_cd);
fprintf('Entropy Rate: %.4f bits/symb\n', R_entr);
fprintf('\n');

% b.3
% The Huffman decoding we will use is
% 00->1, 110->2, 111->3, 01->4 and 10->5
experiments_recovered_symbols_b3 = zeros(no_experiments,n);

for k=1:no_experiments
    experiments_recovered_symbols_b3(k,:) = Huffman_decoding2b_3(experiments_encoded_bits_b3(k,:),n);
end

fprintf('2.4(b3): \n');

% Check whether or not decoding was successful
if experiments_recovered_symbols_b3(:,:) == experiments_initial_symbols_b3(:,:)
    fprintf('-The recovered symbols, indeed match the initial ones!\n');
    fprintf('\n');
else
    fprintf('-There are symbols NOT correctly recovered! Decoding process unsuccessful.\n');
    fprintf('\n');
end

% Now, calculate the numerical average over the
% average word length of every x
avg_word_lengths = zeros(1,no_experiments);

for k=1:no_experiments
    num_of_bits = 0;
    i = 1;
    while experiments_encoded_bits_b3(k,i) ~= nul
        num_of_bits = num_of_bits + 1;
        i = i + 1;
    end
    avg_word_lengths(k) = num_of_bits/n;
end

num_avg_word_length = sum(avg_word_lengths)/no_experiments;

disp('-Comparing the numerical average of the encoded word average length with the entropy of the source');
fprintf('Numerical average: %.4f bits/symb\n', num_avg_word_length);
fprintf('Source Entropy: %.4f bits/symb\n', H_cd);
fprintf('Entropy Rate: %.4f bits/symb\n', R_entr);
fprintf('\n');

% b.4
% The Huffman decoding we will use is
% 0->3, 10->4, 110->1, 1110->5 and 1111->2
experiments_recovered_symbols_b4 = zeros(no_experiments,n);

for k=1:no_experiments
    experiments_recovered_symbols_b4(k,:) = Huffman_decoding2b_4(experiments_encoded_bits_b4(k,:),n);
end

fprintf('2.4(b4): \n');

% Check whether or not decoding was successful
if experiments_recovered_symbols_b4(:,:) == experiments_initial_symbols_b4(:,:)
    fprintf('-The recovered symbols, indeed match the initial ones!\n');
    fprintf('\n');
else
    fprintf('-There are symbols NOT correctly recovered! Decoding process unsuccessful.\n');
    fprintf('\n');
end

% Now, calculate the numerical average over the
% average word length of every x
avg_word_lengths = zeros(1,no_experiments);

for k=1:no_experiments
    num_of_bits = 0;
    i = 1;
    while experiments_encoded_bits_b4(k,i) ~= nul
        num_of_bits = num_of_bits + 1;
        i = i + 1;
    end
    avg_word_lengths(k) = num_of_bits/n;
end

num_avg_word_length = sum(avg_word_lengths)/no_experiments;

disp('-Comparing the numerical average of the encoded word average length with the entropy of the source');
fprintf('Numerical average: %.4f bits/symb\n', num_avg_word_length);
fprintf('Source Entropy: %.4f bits/symb\n', H_cd);
fprintf('Entropy Rate: %.4f bits/symb\n', R_entr);
fprintf('\n');

% b.5
% The Huffman decoding we will use is
% 0->2, 10->3, 110->4, 1110->1 and 1111->5
experiments_recovered_symbols_b5 = zeros(no_experiments,n);

for k=1:no_experiments
    experiments_recovered_symbols_b5(k,:) = Huffman_decoding2b_5(experiments_encoded_bits_b5(k,:),n);
end

fprintf('2.4(b5): \n');

% Check whether or not decoding was successful
if experiments_recovered_symbols_b5(:,:) == experiments_initial_symbols_b5(:,:)
    fprintf('-The recovered symbols, indeed match the initial ones!\n');
    fprintf('\n');
else
    fprintf('-There are symbols NOT correctly recovered! Decoding process unsuccessful.\n');
    fprintf('\n');
end

% Now, calculate the numerical average over the
% average word length of every x
avg_word_lengths = zeros(1,no_experiments);

for k=1:no_experiments
    num_of_bits = 0;
    i = 1;
    while experiments_encoded_bits_b5(k,i) ~= nul
        num_of_bits = num_of_bits + 1;
        i = i + 1;
    end
    avg_word_lengths(k) = num_of_bits/n;
end

num_avg_word_length = sum(avg_word_lengths)/no_experiments;

disp('-Comparing the numerical average of the encoded word average length with the entropy of the source');
fprintf('Numerical average: %.4f bits/symb\n', num_avg_word_length);
fprintf('Source Entropy: %.4f bits/symb\n', H_cd);
fprintf('Entropy Rate: %.4f bits/symb\n', R_entr);
fprintf('\n');

% -- Exercise 3 -- 
clc;clear all; close all;

% 3.6
p = 0.1;
h_file_sizes = zeros(1,10);
sf_file_sizes = zeros(1,10);

fprintf('========== Exercise 3 ==========\n');

% Read the data from the given file
fid = fopen('original_data_5.bin', 'r');
X = fread(fid, 'ubit1');
fclose(fid);

for n=1:10
    disp(['n = ', num2str(n)]);

    % 3.6.1
    [huf_code_book,huf_code_lengths] = Huffman_codebook_generation(p,n);
    [sf_code_book,sf_code_lengths] = Shannon_Fano_codebook_generation(p,n);

    % 3.6.2
    Y_huf = encoding(X,n,huf_code_book,huf_code_lengths);
    Y_sf = encoding(X,n,sf_code_book,sf_code_lengths);

    % 3.6.3
    X_huf_dec = decoding(Y_huf,n,huf_code_book,huf_code_lengths);
    X_sf_dec  = decoding(Y_sf,n,sf_code_book,sf_code_lengths);

    % 3.6.4
    fname = strcat('H_compressed_data_5_',num2str(n),'.bin');
    fid_H = fopen(fname, 'w');
    fwrite(fid_H, Y_huf, 'ubit1');
    fclose(fid_H);
    file_info = dir(fname);
    h_file_sizes(1,n) = file_info.bytes;

    fname = strcat('SF_compressed_data_5_',num2str(n),'.bin');
    fid_SF = fopen(fname, 'w');
    fwrite(fid_SF, Y_sf, 'ubit1');
    fclose(fid_SF);
    file_info = dir(fname);
    sf_file_sizes(1,n) = file_info.bytes;

    % Check for the recovery of the data
    % Huffman
    if isequal(X(:,1)',X_huf_dec(1,:))
        disp('Successful data recovery!');
        disp('\n');
    else
        disp('fail');
    end

    % Shannon-Fano
    if isequal(X(:,1), X_sf_dec(1,:)')
        disp('Successful data recovery!');
        disp('\n');
    else
        disp('fail');
    end
end

%%
% 3.7
M = length(X);
n = [1:10];
LB = H(p).*ones(1,length(n)).*M;
UB = (H(p).*ones(1,length(n))+(1./n)).*M;

figure
plot(n,LB,'Color','k');
hold on;
plot(n,UB,'Color','k','LineStyle','--');
hold on;
plot(n,sf_file_sizes.*8,'Color','b');
hold on;
plot(n,h_file_sizes.*8,'Color','r');
hold off;
grid on;
xlabel('n','FontSize',14);
ylabel('Compressed file size (in bits)','FontSize',13);
axis([1 10 1.5*10^6 5.5*10^6]);
legend('LB','UB','SF','H');

%%
% 3.8

% X
zero_counts_X=0;
one_counts_X=0;
for i =1:M
    if X(i)==0
        zero_counts_X=zero_counts_X+1;
    else
        one_counts_X=one_counts_X+1;
    end
end

% Y_huf
zero_counts_Yhuf=0;
one_counts_Yhuf=0;
for i =1:length(Y_huf)
    if Y_huf(i)==0
        zero_counts_Yhuf=zero_counts_Yhuf+1;
    else
        one_counts_Yhuf=one_counts_Yhuf+1;
    end
end

% Y_sh
zero_counts_Ysf=0;
one_counts_Ysf=0;
for i =1:length(Y_sf)
    if Y_sf(i)==0
        zero_counts_Ysf=zero_counts_Ysf+1;
    else
        one_counts_Ysf=one_counts_Ysf+1;
    end
end

X_pmf=[zero_counts_X./M one_counts_X./M];
Y_sf_pmf=[zero_counts_Ysf./length(Y_sf) one_counts_Ysf./length(Y_sf)];
Y_huf_pmf=[zero_counts_Yhuf./length(Y_huf) one_counts_Yhuf./length(Y_huf)];


figure;
subplot(1,3,1)
stem([0,1],X_pmf,'Marker','none','LineWidth',2);
hold on;
title('X PMF');
ylabel('Probability');
legend('p_X(x)');
grid on;
hold off;
axis([-0.5 1.5 0 1]);

subplot(1,3,2)
stem([0,1],Y_sf_pmf,'Marker','none','LineWidth',2,'Color','r');
hold on;
title('Y_{sf} PMF');
ylabel('Probability');
legend('p_{Ysf}(y)');
grid on;
hold off;
axis([-0.5 1.5 0 1]);

subplot(1,3,3)
stem([0,1],Y_huf_pmf,'Marker','none','LineWidth',2,'Color','m');
hold on;
title('Y_{huf} PMF');
ylabel('Probability');
legend('p_{Yhuf}(y)');
grid on;
hold off;
axis([-0.5 1.5 0 1]);
