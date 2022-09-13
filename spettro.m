function [spettro,f] = spettro(acs,i,nome)

% Queta funzione permette di calcolare lo Spettro di Potenza di un segnale
% applicando il metodo di Blakman-Tuckey. Esso si basa su metodo del
% Periodogramma (di Schuster) utilizzando direttamente la serie di Dati. 
% Il metodo del Periodogramma risulta essere uno stimatore polarizzato 
% dello spettro di potenza con polarizzazione data dalla finestra applicata
% alla serie dati. Il metodo si basa sul Teorema di Wiener-Khinchin che
% risulta essere un metodo indiretto nel calcolo della desnità spettrale di
% potenza. Il metodo di Blackman-Tuckey si basa sul calcolo della funzione
% di autocorellazione dei dati, successivamente si applica una finestra
% adatta (con basso valore di Lag) e infine si calcola la FFT della serie
% di dati. Se sequenze di dati corte si ha scarsa risoluzione, per
% diminuire varianza utilizzare altri metodi.

mat=zeros(length(acs),4);
mat(:,1)=rectwin(length(acs));
mat(:,2)=bartlett(length(acs));
mat(:,3)=hann(length(acs));
mat(:,4)=hamming(length(acs));

% provo differenti finestre per vedere l'effetto ottenuto

N=1024;

% può esser effettuato ZERO PADDING, avremmo lo stesso intervallo di 
% frequenza sia con che senza di esso, però quello effettuato con zero
% padding ha più punti spettrali

spettro=zeros(N,4);

close(figure(i));
figure(i);
hold on
for j=1:4   
    
    
    f_s=fft(acs(:,j),N)/sum(mat(:,j));
    % normalizzato, ma non essendo sempre la finestra rettangolare è meglio
    % normalizzare per la somma dell'intero vettore
    f_s_n=fftshift(f_s);

    ris_f=1.4/N; 
    
    %frequenza di campionamento ottenuta come nomero di battiti/tempo 
    % di campionamento: n. battiti/ tempo = 1.4
    
    if mod(N,2)==0
      f=(-N/2:N/2-1)*ris_f; %pari 
    else
     f=(-(N-1)/2:(N-1)/2)*ris_f; %dispari
    end
    
    
    plot(f,abs(f_s_n))
    title(sprintf('Spettro di potenza di %s',nome))
    xlabel('Frequenza [Hz]')
    ylabel('Potenza [ms^2]')
    
    spettro(:,j) = f_s_n;
end
xlim([0,0.7])
legend('f. rettangolare','f. Bartlett','f. Hann','f. Hamming')
hold off
end

