%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Chinni Chaitanya (EE13B072) and Prafullachandhra (EE16D402)
% Project-1: BCH-encoder-decoder
% EE5160: Error Control Coding
% Name: BCH_decoder.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% clear env %%
clear;
close all;
clc;

%% generate the GF(2^7) field %%
P = 2;
M = 7;
delta = 15;
T = (delta-1)/2;
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

%% clear the logfile for printing intermediate table values and syndrome %%
logfile = fopen('logfile.log', 'w');
fclose(logfile);

%% decode the received messages %%
all_est_codewords = cell(1, 0);
all_est_msgs = cell(1, 0);
est_codewords = gf(zeros(2, N), M, PRIM_POLY);
        
% open the logfile in append mode %
logfile = fopen('logfile.log', 'a');
% iterate over the messages %
for rx_msg = rx_messages
    % print the received codeword to logfile %
    fprintf(logfile, 'Received codeword: %s\n\n', rx_msg{1, 1});
    
    % initialize the falilure flag and message %
    failure_flag = false;
    failure_message = 'None';
    % flag to check if ELP has repeated roots or not %
    repeated_roots = [];
    % flag to check if all roots of ELP lie in GF(P^M) %
    roots_in_field = [];
    % variables to store number of roots of ELP %
    num_roots_err = [];
    num_roots_fun = [];
    
    % replcace erasures with 0 and 1 %
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

        % print the change with which erasure is replaced %
        fprintf(logfile, 'Replace erasure with: %d\n', change);
        % print syndrome to logfile %
        fprintf(logfile, 'Syndrome: ');
        fprintf(logfile, '%d\t', gf2exp(syndrome, FIELD, EXPFORM));
        fprintf(logfile, '\n');
        % print to logfile %
        fprintf(logfile, 'MU\t\t\t\t');
        fprintf(logfile, 'SIGMA(X)\t\t\t\t\t\t\t\t');
        fprintf(logfile, 'DISCREPANCY\t\t\t');
        fprintf(logfile, 'DEGREE\t\t\t');
        fprintf(logfile, 'DIFFERENCE\n');        
        
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

        % print the intermediate table values to logfile %
        for i = 1:2
            fprintf(logfile, '%.2f\t\t\t', Mu(i));

            % get degree of the sigma polynomial %
            j = degree_poly(Sigma(i, :))+1;
            sig = Sigma(i, 1:j);
            fprintf(logfile, '%d\t', gf2exp(sig, FIELD, EXPFORM));
            for tabIter = 1:(10-j)
                fprintf(logfile, '\t');            
            end       
            
            dis = Discrepancy(i);
            fprintf(logfile, '%d\t\t\t\t\t', gf2exp(dis, FIELD, EXPFORM));

            fprintf(logfile, '%d\t\t\t\t', L(i));
            fprintf(logfile, '%d\n', Diff(i));          
        end
        % iterate to find the error-locator-polynomial %
        for i = 1:6
            % find the degree of the current sigma polynomial %
            j = degree_poly(Sigma(i+2, :))+1;
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
            % print the intermediate table values to logfile %
            fprintf(logfile, '%.2f\t\t\t', Mu(i+2));

            % get degree of the sigma polynomial %
            j = degree_poly(Sigma(i+2, :))+1;            
            sig = Sigma(i+2, 1:j);
            fprintf(logfile, '%d\t', gf2exp(sig, FIELD, EXPFORM));
            for tabIter = 1:(10-j)
                fprintf(logfile, '\t');            
            end            
            dis = Discrepancy(i+2);
            fprintf(logfile, '%d\t\t\t\t\t', gf2exp(dis, FIELD, EXPFORM));

            fprintf(logfile, '%d\t\t\t\t', L(i+2));
            fprintf(logfile, '%d\n', Diff(i+2));             
        end
        % print the intermediate table values to logfile %
        fprintf(logfile, '%.2f\t\t\t', Mu(9));

        % get degree of the sigma polynomial %
        j = degree_poly(Sigma(9, :))+1;
        sig = Sigma(9, 1:j);
        fprintf(logfile, '%d\t', gf2exp(sig, FIELD, EXPFORM));
        for tabIter = 1:(10-j)
            fprintf(logfile, '\t');            
        end

        fprintf(logfile, '----\t\t\t\t');

        fprintf(logfile, '%d\t\t\t\t', j-1);
        fprintf(logfile, '%d\n', 2*Mu(9)-j+1);          
        
        fprintf(logfile, '\n');
        
        % estimate the codeword from the error-locator-polynomial %
        % evaluate the error-locator-polynomial for each element in GF(P^M) %
        err_poly = gf(fliplr(Sigma(9, :)), M, PRIM_POLY);
        int_rep_all = bi2de(FIELD(2:(N+1), :), P, 'right-msb');
        ele_field = gf(int_rep_all, M, PRIM_POLY);
        % evaluated values for ELP at each element %
        err_vector = polyval(err_poly, ele_field);
        % find the number of roots from error vector %
        num_roots_err(change+1) = length(find(err_vector == 0));
        % find the number of roots from roots function %
        num_roots_fun(change+1) = length(roots(fliplr(Sigma(9, :))));
        % check for repeated roots %
        if num_roots_err(change+1) < num_roots_fun(change+1)
            repeated_roots(change+1) = true;
        else
            repeated_roots(change+1) = false;
        end
        % check if all roots lie in the GF(P^M) %
        if num_roots_err(change+1) == degree_poly(Sigma(9, :))
            roots_in_field(change+1) = true;
        else
            roots_in_field(change+1) = false;
        end

        % initialize the estimated codeword in GF(P^M) %
        est_codeword = gf(zeros(1, N), M, PRIM_POLY);
        % update the estimated codeword depending on the evaluated vector %
        % update for first element separately and for other elements separately %
        for h = 2:N
            if err_vector(h) == 0
                % implies that element is a root of ELP %
                est_codeword(h-1) = 1 + Rx_poly(h-1);
            else
                % implies that element is NOT a root of ELP %
                est_codeword(h-1) = Rx_poly(h-1);
            end
        end
        % we update for the first element separately because of index notation %
        if err_vector(1) == 0
            est_codeword(N) = 1 + Rx_poly(N); 
        else
            est_codeword(N) = Rx_poly(N); 
        end

        % flip the estimated codeword to get in correct order %
        est_codeword = fliplr(est_codeword);
        % update the array %
        est_codewords(change+1, :) = est_codeword;
    end
    
    % check the syndrome of the two decoded codewords %
    est_code_poly_1 = gf(fliplr(est_codewords(1, :)), M, PRIM_POLY);
    est_code_poly_2 = gf(fliplr(est_codewords(2, :)), M, PRIM_POLY);
    syndrome_check_1 = polyval(est_code_poly_1, bch_roots);
    syndrome_check_2 = polyval(est_code_poly_2, bch_roots);
    
    % determine likely codeword depending on syndrome and examining roots of ELP %
    if syndrome_check_2 == 0 & syndrome_check_1 == 0
        if num_roots_err(1) < num_roots_err(2)
            if roots_in_field(1) == true
                if num_roots_err(1) <= T
                    % most likely code is the first one %
                    likely_codeword = est_codewords(1, :);
                    failure_flag = false;
                else
                    % no likely codeword %
                    % decoder failed %
                    % make them all 0 %
                    likely_codeword = 0*est_codewords(1, :);
                    failure_message = 'Decoder failure: Number of roots > t (=7)';
                    failure_flag = true;
                end
            else
                % no likely codeword %
                % decoder failed %
                % make them all 0 %
                likely_codeword = 0*est_codewords(1, :);
                failure_message = 'Decoder failure: Roots do not belong to GF[128]';
                failure_flag = true;
            end
        else
            if roots_in_field(2) == true
                if num_roots_err(2) <= T
                    % most likely code is the second one %
                    likely_codeword = est_codewords(2, :);
                    failure_flag = false;
                else
                    % no likely codeword %
                    % decoder failed %
                    % make them all 0 %
                    likely_codeword = 0*est_codewords(2, :);
                    failure_message = 'Decoder failure: Number of roots > t (=7)';
                    failure_flag = true;
                end
            else
                % no likely codeword %
                % decoder failed %
                % make them all 0 %
                likely_codeword = 0*est_codewords(1, :);
                failure_message = 'Decoder failure: Roots do not belong to GF[128]';
                failure_flag = true;
            end
        end 
    elseif syndrome_check_1 == 0 & syndrome_check_2 ~= 0
        if repeated_roots(1) == false 
            if num_roots_err(1) <= T 
                if roots_in_field(1) == true
                    % most likely code word is the first one %
                    likely_codeword = est_codewords(1, :);
                    failure_flag = false;
                else
                    % no likely codeword %
                    % decoder failed %
                    % make them all 0 %
                    likely_codeword = 0*est_codewords(1, :);
                    failure_message = 'Decoder failure: Roots do not belong to GF[128]';
                    failure_flag = true;    
                end
            else
                % no likely codeword %
                % decoder failed %
                % make them all 0 %
                likely_codeword = 0*est_codewords(1, :);
                failure_message = 'Decoder failure: Number of roots > t (=7)';
                failure_flag = true;
            end
        else
            % no likely codeword %
            % decoder failed %
            % make them all 0 %
            likely_codeword = 0*est_codewords(1, :);
            failure_message = 'Decoder failure: Error Locator Polynomial has repeated roots';
            failure_flag = true;
        end
    elseif syndrome_check_1 ~= 0 & syndrome_check_2 == 0
        if repeated_roots(2) == false 
            if num_roots_err(2) <= T 
                if roots_in_field(2) == true
                    % most likely code word is the first one %
                    likely_codeword = est_codewords(2, :);
                    failure_flag = false;
                else
                    % no likely codeword %
                    % decoder failed %
                    % make them all 0 %
                    likely_codeword = 0*est_codewords(2, :);
                    failure_message = 'Decoder failure: Roots do not belong to GF[128]';
                    failure_flag = true;    
                end
            else
                % no likely codeword %
                % decoder failed %
                % make them all 0 %
                likely_codeword = 0*est_codewords(2, :);
                failure_message = 'Decoder failure: Number of roots > t (=7)';
                failure_flag = true;
            end
        else
            % no likely codeword %
            % decoder failed %
            % make them all 0 %
            likely_codeword = 0*est_codewords(2, :);
            failure_message = 'Decoder failure: Error Locator Polynomial has repeated roots';
            failure_flag = true;
        end
    else
        % no likely codeword %
        % decoder failed %
        % make them all 0 %
        likely_codeword = 0*est_codewords(1, :);
        failure_message = 'Decoder failure: Number of roots > t (=7)';
        failure_flag = true;
    end
    % decode the original message from the most likely codeword %
    est_msg = likely_codeword(N-K+1:end);
    
    % add them to the vector %
    % convert the gf elements to double %
    likely_codeword_double = double(likely_codeword.x);
    est_msg_double = double(est_msg.x);
    % add them to the cell %
    code = cell(1, 1);
    msg = cell(1, 1);
    % check if atleast one syndrome is 0 %
    % condition for adding error message if decoder failed %
    if failure_flag == false
        for iter = 1:length(est_msg_double)
            msg{1, 1}(iter) = int2str(est_msg_double(iter));
        end
        for iter = 1:length(likely_codeword_double)
            code{1, 1}(iter) = int2str(likely_codeword_double(iter));
        end
    else
        msg{1, 1} = failure_message;
        code{1, 1} = failure_message;
    end
    all_est_msgs{1, end+1} = msg{1, 1};
    all_est_codewords{1, end+1} = code{1, 1};
    
    % print final syndromes for both cases %
    fprintf(logfile, 'Est. Syndrome for "0": ');
    fprintf(logfile, '%d\t', gf2exp(syndrome_check_1, FIELD, EXPFORM));
    fprintf(logfile, '\n');
    fprintf(logfile, 'Est. Syndrome for "1": ');
    fprintf(logfile, '%d\t', gf2exp(syndrome_check_2, FIELD, EXPFORM));
    fprintf(logfile, '\n');
    % print the mostlikely (estimated) codeword %
    fprintf(logfile, 'Est. codeword: %s\n', code{1, 1});
    % print to logfile %
    fprintf(logfile, '\n\n');
    for logIter_1 = 1:2
        for logIter_2 = 1:127
            fprintf(logfile, '=');
        end
        fprintf(logfile, '\n');
    end
    fprintf(logfile, '\n\n\n');
end
% close the logfile %
fclose(logfile);

%% print the message and codeword to the file %%
decodecodefile = fopen('decoderOut_codeword.txt', 'w');
decodemsgfile = fopen('decoderOut_msg.txt', 'w');
for i = 1:size(rx_messages, 2)
    codestr = all_est_codewords{1, i};
    msgstr = all_est_msgs{1, i};
    
    fprintf(decodecodefile, '%s\n', codestr);
    
    fprintf(decodemsgfile, '%s\n', msgstr);
end
fclose(decodecodefile);
fclose(decodemsgfile);
