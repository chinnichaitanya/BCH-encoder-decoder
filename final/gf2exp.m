%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Chinni Chaitanya (EE13B072) and Prafullachandhra (EE16D402)
% Project-1: BCH-encoder-decoder
% EE5160: Error Control Coding
% Name: gf2exp.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% convert gf array elements to exponential form %%
function exp_array = gf2exp(gf_array, FIELD, EXPFORM)
    field = bi2de(FIELD);
    exp_array = [];
    for i = 1:length(gf_array)
        index = find(field == gf_array(i));
        exp_array = [exp_array, EXPFORM(index(1))];
        
        % change -Inf to -1 %
        if exp_array(end) == -Inf
            exp_array(end) = -1;
        end
    end
end