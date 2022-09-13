%% Funzione ACS
% 0 -> Polarizzato
% 1 -> Non Polarizzato
function acf = acsF(vec, n_rit, type, shift)

N = length(vec);

% il vettore al primo posto ha ritardo nullo,
% poi ha n_rit positivi e successivamente n_rit negativi
% posizionamento scelto a lezione

acf(1) = sum(vec .* conj(vec))/N;

for i = 2 : (n_rit + 1) 

    vec_ = [zeros(1,i-1) vec];
    vec_j = [conj(vec) zeros(1,i-1)];
    
    
    acf(i) = sum(vec_ .* vec_j);
    
    
    if type == 1
        acf(i) = acf(i)/(N - i + 1);
    else
       acf(i) = acf(i)/N; 
    end
    
    acf(2*n_rit +3 -i) = conj(acf(i));
    % è necessario il + 3 in quanto io arrivo fino a n-rit +1
    
end
if shift == 1
    acf = fftshift(acf);
end
