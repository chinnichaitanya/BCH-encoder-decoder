%% clear env %%
clear;
close all;
clc;

%% generate the GF(2^7) field %%
P = 2;
M = 7;
delta = 15;
erasure = '2';
N = P^M - 1;
PRIM_POLY = 'D^7 + D^3 + 1';
[FIELD, EXPFORM] = gftuple([-1:P^M-2]', PRIM_POLY, P);

%% find minimal polynomials %%
MINPOLS = gfminpol(EXPFORM, PRIM_POLY, P);
MINPOLS = log2(MINPOLS)/log2(P);

%% generate the narrow-sense BCH generator polynomial with design distance 14 %%
BCH_GEN_POLY = 0;
min_poly_vec = [1, 3, 5, 7, 9, 11, 13];
for i = min_poly_vec
    BCH_GEN_POLY = gfconv(BCH_GEN_POLY, MINPOLS(i+2, :), FIELD);
end

% calculate K %
K = N - (length(BCH_GEN_POLY)-1);

%% read the received messages from 'rx.txt' %%
rx_messages = cell(1, 0);
rx_msgfile = fopen('rx.txt');
fileline = fgetl(rx_msgfile);
while ischar(fileline)
    rx_messages{1, end+1} = fileline;
    fileline = fgetl(rx_msgfile);
end
fclose(rx_msgfile);

%% decode the received messages %%
est_codewords_vec = gf(zeros(2, N), M, PRIM_POLY);
num_errors = [];
for rx_msg = rx_messages
    for change = [0, 1]
        % replace erasures with 0's and 1's %
        rx = strrep(rx_msg{1, 1}, erasure, int2str(change));
        % generate the received vector vector form %
        Rx = [];
        for i = 1:length(rx)
            % convert string to number %
            Rx(i) = str2double(rx(i));
        end
        
        % flip the vector as `gf` function requires it in descending order of powers %
        Rx = fliplr(Rx);
        % generate the received vector polynomial in GF(P^M) %
        Rx_poly = gf(Rx, M, PRIM_POLY);
        
        % convert (alpha^1, ... alpha^14) to corresponding integer representation in GF-array %
        integer_representation = bi2de(FIELD(3:16, :), P, 'right-msb');
        % convert those integer represented numbers to GF-elements %
        bch_roots = gf(integer_representation, M, PRIM_POLY);
        % calculate the syndrome at the roots %
        syndrome = polyval(Rx_poly, bch_roots);
        
        % simplified BM-algorithm for finding error-locator-polynomial %
        
        % initializing the variables %
        % mu-vector %
        Mu = [-0.5, 0, 1, 2, 3, 4, 5, 6, 7];
        
        % discrepancy-vector %
        Discrepancy = gf([1, 0, 0, 0, 0, 0, 0, 0, 0], M, PRIM_POLY);
        
        % l-vector %
        L = [0, 0, 0, 0, 0, 0, 0, 0, 0];
        
        % difference (2*mu-rho) vector %
        Diff = [-1, 0, 0, 0, 0, 0, 0, 0, 0];
        
        % sigma polynomial initialization in decimal form %
        sigma_dec = zeros(9, 84);
        % constant term is always 1 %
        sigma_dec(:, 1) = 1;
        
        % convert sigma-polynomial to gf-field type %
        Sigma = gf(sigma_dec, M, PRIM_POLY);
        % update the sigma-polynomial for the first iteration %
        Sigma(3, 2) = syndrome(1);
        
        % update the discrepancy for the first iteration %
        Discrepancy(2) = syndrome(1);
        
        % iterate to find the error-locator-polynomial %
        for i = 1:6
            % find the degree of the current sigma polynomial %
            j = 84;
            while Sigma(i+2, j) == 0
                j = j-1;
            end
            % update the corresponding vectors %
            L(i+2) = j-1;
            Diff(i+2) = 2*Mu(i+2) - L(i+2);
            % calculate the discrepancy %
            Discrepancy(i+2) = 0;
            for k = 0:L(i+2) 
                Discrepancy(i+2)= Discrepancy(i+2) + Sigma(i+2, k+1)*syndrome(2*i+1-k);
            end
            
            % update the sigma polynomial depending on the discrepancy %
            if Discrepancy(i+2) == 0
                Sigma(i+3, :) = Sigma(i+2, :);
            else 
                % find the index of rho to maximize the difference %
                max = Diff(2);
                rho_index = 2;
                for t = 2:(i+1)
                    if Discrepancy(t) ~= 0
                        if max < Diff(t)
                            max = Diff(t);
                            rho_index = t;
                        end
                    end 
                end

                coeff = Discrepancy(i+2)/Discrepancy(rho_index);
                degree = 2*(Mu(i+2) - Mu(rho_index));
                add_poly = mul_poly(degree, Sigma(rho_index, :), L(rho_index), M, PRIM_POLY);
                Sigma(i+3, :) = Sigma(i+2, :) + coeff*add_poly;
            end
        end
        
        % find the degree of the final error-locator-polynomial %
        j = 84;
        while Sigma(end, j) == 0
            j = j-1;
        end
        num_errors(change+1) = j;
        
        % estimate the codeword from the error-locator-polynomial %
        % evaluate the error-locator-polynomial for each element in GF(P^M) %
        err_poly = gf(fliplr(Sigma(end, :)), M, PRIM_POLY);
        int_rep_all = bi2de(FIELD(2:(N+1), :), P, 'right-msb');
        ele_field = gf(int_rep_all, M, PRIM_POLY);
        % evaluated values for ELP at each element %
        err_vector = polyval(err_poly, ele_field);

        % initialize the estimated codeword in GF(P^M) %
        est_codeword = gf(zeros(1, N), M, PRIM_POLY);
        % update the estimated codeword depending on the evaluated vector %
        % update for 1 separately and for other elements separately %
        for h = 2:N
            if err_vector(h) == 0
                est_codeword(h-1) = 1 + Rx_poly(h-1);
            else
                est_codeword(h-1) = Rx_poly(h-1); 
            end
        end
        if err_vector(1) == 0
            est_codeword(N) = 1 + Rx_poly(N); 
        else
            est_codeword(N) = Rx_poly(N); 
        end

        % flip the estimated codeword to get in correct order %
        est_codeword = fliplr(est_codeword);
        % update the array %
        est_codewords_vec(change+1, :) = est_codeword;        
    end
    
    % estimate the original message based on the number of errors occurred %
    if num_errors(1) > num_errors(2)
        % most likely codeword is the second one %
        est_msg = est_codewords_vec(2, N-K+1:end);
    else
        % most likely code word is the first one %
        est_msg = est_codewords_vec(1, N-K+1:end);
    end
    
    % print the message and codeword to the file %
    
end