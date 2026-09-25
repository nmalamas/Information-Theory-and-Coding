function symbol_states = simulate_MC(n,p_initial,transition_matrix)
    p = zeros(length(p_initial),n);
    p(:,1) = p_initial;
    P = transition_matrix;

    symbol_states = zeros(1,n);

    % Initialize the state of the MC
    % based on the initial probability distribution
    r = rand;
    if r < (p(1))
        init_state = 1;
    elseif r < (sum(p(1:2)))
        init_state = 2;
    elseif r < (sum(p(1:3)))
        init_state = 3;
    elseif r < (sum(p(1:4)))
        init_state = 4;
    else
        init_state = 5;
    end

    % Set the current MC strate to initial
    % and run the MC
    cur_state = init_state;
    symbol_states(1,1) = cur_state;

    for i=2:n
        p(:,i)=P*p(:,i-1);

        for j=1:length(P(:,1))
            r = rand;
            if r < (sum(p(1:j)))
                cur_state = j;
                break;
            end
        end

        symbol_states(1,i) = cur_state;
    end
end

% for j=1:length(P(:,1))
%             if cur_state == j
%                 r = rand;
%                 if r < (P(j,1))
%                     cur_state = 1;
%                     continue;
%                 elseif r < (sum(P(j,1:2)))
%                     cur_state = 2;
%                     continue;
%                 elseif r < (sum(P(j,1:3)))
%                     cur_state = 3;
%                     continue;
%                 elseif r < (sum(P(j,1:4)))
%                     cur_state = 4;
%                     continue;
%                 else
%                     cur_state = 5;
%                     continue;
%                 end
%             end
%             %x(i) = cur_state;
%         end