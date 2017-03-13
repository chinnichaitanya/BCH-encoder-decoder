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
% [m, n] = size(MINPOLS);
% for i = 1:m
%     for j = 1:n
%         MINPOLS(i, j) = log2(MINPOLS(i, j))/log2(P);
%     end
% end

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
for rx_msg = rx_messages
    for change = [0, 1]
        % replace erasures with 0's and 1's %
        rx = strrep(rx_msg{1, 1}, erasure, int2str(change));
        % generate the received vector vector form %
        Rx = [];
        for i = 1:length(rx)
            % convert string to number %
            Rx(i) = str2double(rx(i));
            % convert to exponential form %
            % Rx(i) = log2(Rx(i))/log2(P);
        end
        
        % flip the vector as `gf` function requires it in descending order of powers %
        gfroots(Rx, PRIM_POLY, P)
        Rx = fliplr(Rx);
        % generate the received vector polynomial in GF(P^M) %
        Rx_poly = gf(Rx, M, PRIM_POLY);
        roots(Rx_poly)
        
        % calculate the syndrome %
        bch_roots = gf(min_poly_vec+1, M, PRIM_POLY);
        syndrome = polyval(Rx_poly, bch_roots);
    end
end
syndrome

%% alternate decoding %%
% for rx_msg = rx_messages
%     % replace the erasure first with 0 and then with 1 and error correct %
%     for change = [0, 1]
%         % replace erasures with 0's and 1's %
%         rx = strrep(rx_msg{1, 1}, erasure, int2str(change));
%         % generate the received vector in vector form %
%         Rx = [];
%         for i = 1:length(rx)
%             % convert string to number %
%             Rx(i) = str2double(rx(i));
%             % convert to exponential form %
%             % Rx(i) = log2(Rx(i))/log2(P);
%         end
%         
%         syndrome = [];
%         for i = min_poly_vec
%             temp_rx = Rx;
%             for j = 1:length(temp_rx)
%                 temp_rx(j) = mod(temp_rx(j)*i*j, N-1);
%             end
%             
%             sum = -Inf;
%             for j = temp_rx
%                 sum = gfadd(sum, EXPFORM(j+2), FIELD);
%             end
%             syndrome(end+1) = sum;
%         end
%     end
% end