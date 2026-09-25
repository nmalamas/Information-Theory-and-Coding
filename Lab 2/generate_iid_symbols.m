function X =  generate_iid_symbols(n)
X = zeros(1,n);    

    for i=1:n
        r = rand;
        if r < (1/2)
            X(i) = 1;
        elseif r < (3/4)
            X(i) = 2;
        elseif r < (7/8)
            X(i) = 3;
        else
            X(i) = 4;
        end
    end
end