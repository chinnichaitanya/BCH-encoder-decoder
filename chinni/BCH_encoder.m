%% clear env %%
clear;
close all;
clc;

%% generate the GF(2^7) field %%
P = 2;
M = 7;
N = P^M - 1;
PRIM_POLY = [1, 0, 0, 1, 0, 0, 0, 1];
[FIELD, EXPFORM] = gftuple([-1:P^M-2]', PRIM_POLY, P);

%% find minimal polynomials %%
MINPOLS = gfminpol(EXPFORM, PRIM_POLY, P);
[m, n] = size(MINPOLS);
for i = 1:m
    for j = 1:n
        if MINPOLS(i, j) == 0
            MINPOLS(i, j) = -Inf;
        else
            MINPOLS(i, j) = log2(MINPOLS(i, j))/log2(P);
        end
    end
end

%% generate the narrow-sense BCH generator polynomial with design distance 14 %%
BCH_GEN_POLY = 0;
for i = [1, 3, 5, 7, 9, 11, 13]
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
        if mx(i) == 0
            mx(i) = -Inf;
        else
            mx(i) = log2(mx(i))/log2(P);
        end
    end
    % multiply mx with x^N-K to get Mx %
    Mx = padarray(mx, [0, N-K], -Inf, 'pre');
    % calculate the remainder to encode systematically %
    [Qx, Rx] = gfdeconv(Mx, BCH_GEN_POLY, FIELD);
    % compute the code word %
    cx = gfsub(Mx, Rx, FIELD);
    % add zeros for higher order terms which don't exist to maintain uniform size %
    Cx = padarray(cx, [0, N-length(cx)], -Inf, 'post');
    % convert exponential to GF(2) numbers %
    for j = 1:length(Cx)
        if Cx(j) == -Inf
            Cx(j) = 0;
        else
            Cx(j) = P^Cx(j);
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
