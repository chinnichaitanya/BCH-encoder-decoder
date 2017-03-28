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
        s = polyval(Rx_poly, bch_roots);
        %simplified iterative algorithm for finding error location poly
        u = [-0.5,0,1,2,3,4,5,6,7];
        d = [1 0 0 0 0 0 0 0 0 ];
        l = [0 0 0 0 0 0 0 0 0 ];
        diff = [-1,0,0,0,0,0,0,0,0];
        a = zeros(9,84);
        a(:,:) = -inf; 
        a(:,1) = 1;
        a(3,2) = 1;
        d(2) = s(1);
        for i = 1:6
          j = 84;
          while a(i+2,j) == -inf
              j = j-1;
          end
          l(i+2) = j;
          diff(i+2) = 2(i)-l(i+2);
          d(i+2) = s(2i+1);
          for k = 1:j 
              d(i+2)= d(i+2)+a(i+2,k)s(2i+1-k);
          end
          if d(i+2) == 0
              a(i+3,:) = a(i+2,:);
          else 
              max = diff(2);
              for t = 2:(i+1)
                  if d(t) == 0
                     if max < diff(t)
                         max == diff(t);
                     v = t;
                     end
                  end 
              end
              b = d(i+2)/d(v);
              c = 2*(i - (v+2));
              w = fliplr(a(v,:));
              w = bitshift(w,c);
              w = fliplr(w);
              a(i+3,:) = a(i+2,:)+b*w;
          end
        end
        j = 84;
          while a(9,j) == -inf
              j = j-1;
          end
          
          if j > 8
             %print a message saying cannot correct
          else 
              %correcting each bit
            err_poly = a(9,2:8);
            int_rep_all = bi2de(FIELD(2:128, :), P, 'right-msb');
            err_vector = gf(int_rep_all, M, PRIM_POLY);
            corr_vector = zeros(127);
            for h = 1:127
                if err_vector(h) == 1
                    corr_vector(h) = 1+Rx_poly(h);
                    
                else
                    corr_vector(h) = Rx_poly(h);
                    
                end
                
            end
          end 
            
           
              
              
              
              
              
          
        
    end
end