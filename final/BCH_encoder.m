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
% convert to exponential form %
MINPOLS = log2(MINPOLS)/log2(P);

%% generate the narrow-sense BCH generator polynomial with design distance 15 %%
BCH_GEN_POLY = 0;
min_poly_vec = [1, 3, 5, 7, 9, 11, 13];
for i = min_poly_vec
    BCH_GEN_POLY = gfconv(BCH_GEN_POLY, MINPOLS(i+2, :), FIELD);
end

% calculate K %
K = N - (length(BCH_GEN_POLY)-1);

%% read the messages from 'msg.txt' %%
messages = cell(1, 0);
msgfile = fopen('msg.txt');
fileline = fgetl(msgfile);
while ischar(fileline)
    messages{1, end+1} = fileline;
    fileline = fgetl(msgfile);
end
fclose(msgfile);

%% encode the messages systematically %%
codewords = cell(1, 0);
for msgstr = messages
    mx = [];
    for i = 1:length(msgstr{1,1})
        % get individual message bits %
        mx(i) = str2double(msgstr{1,1}(i));
        
        % convert them to exponential form for computations in GF(2^7) %
        mx(i) = log2(mx(i))/log2(P);
    end
    % multiply mx with x^N-K to get Mx %
    Mx = padarray(mx, [0, N-K], -Inf, 'pre');
    % calculate the remainder to encode systematically %
    [Qx, Rx] = gfdeconv(Mx, BCH_GEN_POLY, FIELD);
    % compute the code word %
    Cx = gfsub(Mx, Rx, FIELD);
    % convert exponential to GF(2) numbers %
    for j = 1:length(Cx)
        if Cx(j) == -Inf
            Cx(j) = 0;
        elseif Cx(j) == 0
            Cx(j) = 1;
        end
    end

    % store the codeword in cell %
    code = cell(1, 1);
    for j = 1:length(Cx)
        code{1, 1}(j) = int2str(Cx(j));
    end
    % store each code in a cell %
    codewords{1, end+1} = code{1, 1};
end

%% store the encoded messages into 'codeword.txt' file %%
codefile = fopen('codeword.txt', 'w');
for codestr = codewords
   fprintf(codefile, '%s\n', codestr{1, 1});
end
fclose(codefile);

