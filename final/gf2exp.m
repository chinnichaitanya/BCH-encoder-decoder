%% converts gf array elements to exponential form %%
function exp_array = gf2exp(gf_array, FIELD, EXPFORM)
    field = bi2de(FIELD);
    exp_array = [];
    for i = 1:length(gf_array)
        index = find(field == gf_array(i));
        exp_array = [exp_array, EXPFORM(index(1))];
    end
end