%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Information Theory and Coding - TEL 417
% Exercise 1
% Patsea Ioanna-Georgia -  2020030033
% Malamas Nikolaos - 2020030180
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;close all;clear all;

% 3.1
p=0.45;
n=15;
epsilon = 0.05;

no_seqs = 2^n;
P = zeros(1,no_seqs);

for i = 0:no_seqs-1
    bin_seq = dec2bin(i, n);
    ones_count= sum(bin_seq == '1');
    P(i+1) = (p^ones_count) * ( (1-p)^(n-ones_count) );
end

% Calculate the probability bounds
H1 = -p*log2(p)-(1-p)*log2(1-p);

LB = 2^(-n*(H1 + epsilon));
UB = 2^(-n*(H1 - epsilon));

sorted_probs = sort(P, 'descend');

figure;
semilogy(sorted_probs, 'o-');
hold on;
semilogy([1, no_seqs], [LB, LB], 'r--', 'LineWidth', 2);
semilogy([1, no_seqs], [UB, UB], 'g--', 'LineWidth', 2);
title("3.1: Sequence propabilities in decreasing order ($p=0.45 and n=15$)",'Interpreter','latex');
hold off;
axis([-0.25*10^4 3.5*10^4 0 0.00025]);
grid on;
xlabel('i');
ylabel('Probability');


% 3.2
p=0.1;
n=15;
epsilon = 0.05;

no_seqs = 2^n;
P = zeros(1,no_seqs);

for i = 0:no_seqs-1
    bin_seq = dec2bin(i, n);
    ones_count= sum(bin_seq == '1');
    P(i+1) = (p^ones_count) * ( (1-p)^(n-ones_count) );
end

% Calculate the probability bounds
H2 = -p*log2(p)-(1-p)*log2(1-p);

LB = 2^(-n*(H2 + epsilon));
UB = 2^(-n*(H2 - epsilon));

sorted_probs = sort(P, 'descend');

figure;
semilogy(sorted_probs, 'o-');
hold on;
semilogy([1, no_seqs], [LB, LB], 'r--', 'LineWidth', 2);
semilogy([1, no_seqs], [UB, UB], 'g--', 'LineWidth', 2);
title("3.2: Sequence propabilities in decreasing order ($p=0.1 and n=15$)",'Interpreter','latex');
hold off;
axis([-0.25*10^4 3.5*10^4 0 0.5]);
grid on;
xlabel('i');
ylabel('Probability');



n=50;
epsilon = 0.05;

% 3.3
p = 0.45;
P = [];

for i = 1:n-1
    bin_seq = dec2bin(2^(i)-1, n);
    ones_count = sum(bin_seq == '1');
    P_i = (p^ones_count) * ( (1-p)^(n-ones_count) );
    P = [P P_i P_i];
end

P_final = [(p^0) * ( (1-p)^(n-0) ) P (p^n) * ( (1-p)^(n-n) )];

% Calculate the probability bounds
H3 = -p*log2(p)-(1-p)*log2(1-p);

LB = 2^(-n*(H3 + epsilon));
UB = 2^(-n*(H3 - epsilon));

sorted_probs = sort(P_final, 'descend');
xaxis = [0];
for i=1:n-1
    cur_idx = nchoosek(n,i);
    xaxis = [xaxis xaxis(end) xaxis(end)+cur_idx];
end

xaxis = [xaxis xaxis(end)+nchoosek(n,n)];

figure;
semilogy(xaxis, sorted_probs, 'o-');
hold on;
yline(LB,'r--');
hold on;
yline(UB,'g--');
hold on;
axis([-10^14 12*10^14 10^(-18)/2 10^(-12)]);
grid on;
xlabel('i');
ylabel('Probability');
title('3.3: Sequence propabilities in decreasing order ($p=0.45$ and $n=50$)','Interpreter','latex');
hold off;

% 3.4
p = 0.1;

no_seqs = 2^n;
P = [];

for i = 1:n-1
    bin_seq = dec2bin(2^(i)-1, n);
    ones_count = sum(bin_seq == '1');
    P_i = (p^ones_count) * ( (1-p)^(n-ones_count) );
    P = [P P_i P_i];
end

P_final = [(p^0) * ( (1-p)^(n-0) ) P (p^n) * ( (1-p)^(n-n) )];

% Calculate the probability bounds
H4 = -p*log2(p)-(1-p)*log2(1-p);

LB = 2^(-n*(H3 + epsilon));
UB = 2^(-n*(H3 - epsilon));

sorted_probs = sort(P_final, 'descend');
xaxis = [0];
for i=1:n-1
    cur_idx = nchoosek(n,i);
    xaxis = [xaxis xaxis(end) xaxis(end)+cur_idx];
end

xaxis = [xaxis xaxis(end)+nchoosek(n,n)];

figure;
semilogy(xaxis, sorted_probs, 'o-');
hold on;
yline(LB,'r--');
hold on;
yline(UB,'g--');
hold on;
axis([-10^14 12*10^14 10^(-51) 10^(-1)]);
grid on;
xlabel('i');
ylabel('Probability');
title('3.4: Sequence propabilities in decreasing order ($p=0.1$ and $n=50$)','Interpreter','latex');
hold off;


n = [10 : 10 : 1000];

% 3.5
p1 = 0.45;
p2 = 0.1;

bin_seq = 2.^(n);
A_epsilon1= zeros(1,length(n));
A_epsilon2 = zeros(1,length(n));

P_typical1 = zeros(1,length(n));
P_typical2 = zeros(1,length(n));

for k=1:length(n)
    no_tseq1 = 0;
    tprob1 = 0;
    
    no_tseq2 = 0;
    tprob2 = 0;

    for m=0:n(k)
        if ( (p1^m)*((1-p1)^(n(k)-m)) <= 2^(-n(k)*(H1-epsilon)) ) && ( (p1^m)*((1-p1)^(n(k)-m)) >= 2^(-n(k)*(H1+epsilon)) )
            no_tseq1 = no_tseq1 + nchoosek(n(k),m);
            tprob1 = tprob1 + (p1^m)*((1-p1)^(n(k)-m))*nchoosek(n(k),m);
        end

        if ( (p2^m)*((1-p2)^(n(k)-m)) <= 2^(-n(k)*(H2-epsilon)) ) && ( (p2^m)*((1-p2)^(n(k)-m)) >= 2^(-n(k)*(H2+epsilon)) )
            no_tseq2 = no_tseq2 + nchoosek(n(k),m);
            tprob2 = tprob2 + (p2^m)*((1-p2)^(n(k)-m))*nchoosek(n(k),m);
        end
    end
    A_epsilon1(k)=no_tseq1;
    A_epsilon2(k)=no_tseq2;

    P_typical1(k) = tprob1;
    P_typical2(k) = tprob2;
end

LB_A1 = 2.^(n.*(H1+epsilon));
UB_A1 = (1-epsilon).*2.^(n.*(H1-epsilon));

LB_A2 = 2.^(n.*(H2+epsilon));
UB_A2 = (1-epsilon).*2.^(n.*(H2-epsilon));


% For p = 0.45
figure; 
plot(n,A_epsilon1,'LineWidth',1,'Color','b');
yscale("log");
grid on;
hold on;
plot(n, [LB_A1], 'r--', 'LineWidth', 2);
hold on;
plot(n, [UB_A1], 'g--', 'LineWidth', 2);
hold on;
plot(n,bin_seq,'LineStyle','--','Color','c');
hold on;
title("3.5 ($p=0.45$): $|A_{\epsilon}^{(n)}|$ and $\#Binary\ Sequences$ " ,'Interpreter','latex');
legend('#Typical Sequences','Lower Bound','Upper Bound','#Binary Sequences','Location','north');
xlabel('$n$','Interpreter','latex','FontSize',14);

% For p = 0.1
figure; 
plot(n,A_epsilon2,'LineWidth',1,'Color','b');
yscale("log");
grid on;
hold on;
plot(n, [LB_A2], 'r--', 'LineWidth', 2);
hold on;
plot(n, [UB_A2], 'g--', 'LineWidth', 2);
hold on;
plot(n,bin_seq,'LineStyle','--','Color','c');
hold on;
title("3.5 ($p=0.1$): $|A_{\epsilon}^{(n)}|$ and $\#Binary\ Sequences$ " ,'Interpreter','latex');
legend('#Typical Sequences','Lower Bound','Upper Bound','#Binary Sequences','Location','north');
xlabel('$n$','Interpreter','latex','FontSize',14);

% 3.6
% For p = 0.45
p = 0.45;
LB_epsilon = 1-epsilon;

figure; 
plot(n,P_typical1);
grid on;
hold on;
yline(LB_epsilon, 'r--', 'LineWidth', 2);
title("3.6 ($p=0.45$): $Pr$(Typical Sequence)",'Interpreter','latex');
xlabel('$n$','Interpreter','latex','FontSize',14);
ylabel('Probability');
legend('Pr(typ. seq)','1-ε');
hold off;

% For p = 0.1
p = 0.1;
LB_epsilon = 1-epsilon;

figure; 
plot(n,P_typical2);
grid on;
hold on;
yline(LB_epsilon, 'r--', 'LineWidth', 2);
title("3.6 ($p=0.1$): $Pr$(Typical Sequence)",'Interpreter','latex');
xlabel('$n$','Interpreter','latex','FontSize',14);
ylabel('Probability');
legend('Pr(typ. seq)','1-ε');
hold off;

% 3.7
% p = 0.45
no_bits1 = ( (ceil(log2(A_epsilon1)) + 1) .* P_typical1 + (n+1) .* (1-P_typical1))/n;
LB_bits1 = (H1 - epsilon);
UB_bits1 = (H1 + epsilon);

figure;
h1 = plot(n, no_bits1.*ones(length(n)), 'c', 'LineWidth', 2);  
hold on;
h2 = yline(LB_bits1, 'r--', 'LineWidth', 2); 
hold on;
h3 = yline(UB_bits1, 'g--', 'LineWidth', 2); 
hold on;
h4 = yline(H1, 'k--', 'LineWidth', 1);
hold on;
grid on;
axis([1 1000 0.9 1.1]);
set(gca, 'YScale', 'log');
title("3.7 ($p=0.45$): $\#$bits and $H(X)$", 'Interpreter', 'latex');
xlabel('$n$', 'Interpreter', 'latex', 'FontSize', 14);
legend([h1(1), h2(1), h3(1), h4(1)], '#Bits', 'Lower Bound', 'Upper Bound', 'H_1(X)', 'Location', 'best');
hold off;


% p = 0.1
no_bits2 = ( (ceil(log2(A_epsilon2)) + 1) .* P_typical2 + (n+1) .* (1-P_typical2))/n;
LB_bits2 = (H2 - epsilon);
UB_bits2 = (H2 + epsilon);

figure;
h1 = plot(n, no_bits2.*ones(length(n)), 'c', 'LineWidth', 2);  
hold on;
h2 = yline(LB_bits2, 'r--', 'LineWidth', 2); 
hold on;
h3 = yline(UB_bits2, 'g--', 'LineWidth', 2); 
hold on;
h4 = yline(H2, 'k--', 'LineWidth', 1);
hold on;
grid on;
axis([1 1000 0.4 0.7]);
set(gca, 'YScale', 'log');
title("3.7 ($p=0.1$): $\#$bits and $H(X)$", 'Interpreter', 'latex');
xlabel('$n$', 'Interpreter', 'latex', 'FontSize', 14);
legend([h1(1), h2(1), h3(1), h4(1)], '#Bits', 'Lower Bound', 'Upper Bound', 'H_2(X)', 'Location', 'best');
hold off;