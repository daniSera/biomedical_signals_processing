function [f,f_s,ris_f] = discreteFt(s,Fs,Nfft,control)

% Script per lo svolgimento della TRASFORMATA DI FOURIER DISCRETA seguendo
% la fuzione FFT fornita da MATLAB
   
    ris_f = Fs/Nfft;
    f_s = fft(s-mean(s),Nfft)/length(s); 
    
    % FFT effettua ZERO-PADDING quindi devo dividere per la lunghezza del
    % vettore s in modo che il valore sia corretto. Va tolto il valore
    % medio in modo che sia normalizzato 
    
    
    if control == 1
        f_s = fftshift(f_s);
        
        % riordina le frequenza da 0 alla frequenza maggiore
        
        if mod(Nfft,2) == 0
            f = (-Nfft/2 : 1 : (Nfft/2)-1) * ris_f;
            
            % se NFFT è pari è necessario rappresentare dalla frequenza di
            % Nyquist fino ad essa meno un passo
            
        else
            f = (-(Nfft-1)/2 : 1 : (Nfft-1)/2) * ris_f;
            
            % se NFFT è dipari la frequenza di Nyquist non viene calcolata
        end
    else
        f = (0 : 1 : Nfft-1) * ris_f;
    end  
    
% La TRASFORMATA DI FOURIER permette di svolgere l'analisi in frequnza di
% segnali aperiodici considerando un periodo infinito. In particolare nel
% nostro caso però è necessario considerare una TRASFORMATA DISCRETA di 
% FOURIER che ci permetta di analizzare casi reali dove funzioni discrete
% in quanto campionate per una durata limitata vengao rappresentate nelle
% loro componenti frequenziali. 

% Da teoria una DFT applicata su una funzione discreta nel tempo,
% periodica con periodo N deve restituire una funzione discreta nelle 
% frequenze, periodica con periodo N. Nella realtà però il campionamento
% impone l'utilizzo di FINESTRATURE che attraverso la convoluzione
% influenzano la rappresentaione in frequenza tramite una distorsione, che
% può essere ridotta con un altra finestratura che smussi i dati all'inizio
% e fine di essa
    

end