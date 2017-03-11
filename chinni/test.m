%% clear window %%
clear;
close all;
clc;

%% testing GFADD %%
P = 2;
M = 7;
N = P^M - 1;
PRIM_POLY = [1, 0, 0, 1, 0, 0, 0, 1];
% FIELD = gftuple([-1:P^M-2]', M, P);
[FIELD, EXPFORM] = gftuple([-1:P^M-2]', PRIM_POLY, P);

%% finding minimal polynomials %%
MINPOLS = gfminpol(EXPFORM, PRIM_POLY, P);
[m, n] = size(MINPOLS);
for i = 1:m
    for j = 1:n
        if MINPOLS(i, j) == 0
            MINPOLS(i, j) = -Inf;
        elseif MINPOLS(i, j) == 1
            MINPOLS(i, j) = 0;
        end
    end
end

% generating the BCH code %
BCH_GEN_POLY = 0;
for i = [1, 3, 5, 7, 9, 11, 13]
    BCH_GEN_POLY = gfconv(BCH_GEN_POLY, MINPOLS(i+2, :), FIELD);
end

% calculating K %
K = N - (length(BCH_GEN_POLY)-1);

%% encoding using BCH code %
mx = [0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1];
for i = 1:length(mx)
    if mx(i) == 0
        mx(i) = -Inf;
    elseif mx(i) == 1
        mx(i) = 0;
    end
end
% multiplying mx with x^n-k to get Mx
Mx = padarray(mx, [0, N-K], -Inf, 'pre');
[Qx, Rx] = gfdeconv(Mx, BCH_GEN_POLY, FIELD);
% encoding systematically %
cx = gfsub(Mx, Rx, FIELD);
Cx = padarray(cx, [0, N-length(cx)], -Inf, 'post');
for i = 1:length(Cx)
    if Cx(i) == -Inf
        Cx(i) = 0;
    elseif Cx(i) == 0
        Cx(i) = 1;
    end
end
gfpretty(Cx);