% pull in a sample sequences and the known states (for comparision after)
seq = dlmread('sample_sequence_1.txt');       
states = dlmread('sample_states.txt', '_'); 

n = length(seq);

a = log([0.9551  0.0449; 0.0880  0.9120]);     % state transition matrix

e = [1.8234  5.7812];                          % emission probabilities, where index = state


vprob = [1; 0];                                % assumes start state is 1 (p = 1)


% viterbi algorithm in log space (to avoid underflow)
for i = 2:(n + 1)
    for j = 1: 2
        vprob(j, i) = log( poisspdf(seq(i - 1), e(j)) ) + max( vprob(:, i - 1) + a(:, j) );
        [pi(j, i), ptr(j, i)] = max( vprob(:, i - 1) + a(:, j) );
    end
end

% traceback to find most likely path
[~, path(n + 1)] = max(pi(:, n + 1));
for i = n: -1: 2
    path(i) = ptr(path(i + 1), i + 1);
end

% compare states found by viterbi to the known states solution key
path = path(:,2:251);
states = states';

state_compare = [path; states(2, :)]          
        
