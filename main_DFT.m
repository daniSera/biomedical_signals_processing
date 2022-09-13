%% Lezione 5 ESERCIZIO 1/3

% Script in cui vengono usate le nozioni di teoria riguardanti l' ANALISI
% IN FREQUENZA e la TRASFORMATA DI FOURIER, analisi della funzione FFT 
% rappresentante la TRASFORMATA DISCRETA DI FOURIER in MATLAB.

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

clear;
close all;
clc;

% Dati di Entrata

Fs = 10; % frequenza
Ts = 1/Fs;
t = -10:Ts:10;

Nfft = [199 200 201 256]; 

% numero di punti su cui calcolare la FFT, se PARI l'algoritmo fft calcola
% la componente alla frequenza di Nyquist (Fc/2) e la pone in uscita alla
% posizione 1+NFFT/2, mentre se DISPARI la componente non viene calcolata.
% La funzioen FFT risulta più performante se NFFT è una potenza di 2, viene
% effettuato in automatico il zero padding.


% Definizioone della Funzione

f0 = 2;
s_sin = sin(2*pi*f0*t);
s_cos = cos(2*pi*f0*t);

hold on;
figure(1);
plot(t,sin(2*pi*f0*t),'r');
plot(t,s_sin,'Marker','X');
title('Sin');
hold off;

hold on;
figure(2);
plot(t,cos(2*pi*f0*t),'r');
plot(t,s_cos,'Marker','X');
title('Cos');
hold off;

% Calcolo della Trasformata e Plot


for i = 1 : length(Nfft)
    f_ssin = fft(s_sin,Nfft(i))/ Nfft(i); 
    f_scos = fft(s_cos,Nfft(i))/ Nfft(i);
    % trasformata normalizzata, cioè divisa per il numero NFFT
    f_shiftsin = fftshift(f_ssin); 
    f_shiftcos = fftshift(f_scos);
    % Taratura Asse, sposta cioè la componente con frequenza nulla al
    % centro dello spettro
    
    ris_f = 1/(Nfft(i)*Ts);
    
    % risoluzione frequenziale, all'aumentare del numero di NFFT diminuisce
    % la risoluzione
   
    if mod(Nfft(i),2) == 0
        f = (-Nfft(i)/2 : 1 : (Nfft(i)/2)-1) * ris_f;
    else
        f = (-(Nfft(i)-1)/2 : 1 : (Nfft(i)-1)/2) * ris_f;
    end
    
    % questa parrte di codice mi permette di settare il limite su x alla 
    % frequenza di nyquist se pari o a più o meno la frequenza di nyquist 
    % se dipari    
    
    hold on;
    figure(i+2);
    plot(f,2*abs(f_shiftsin));    
    title(sprintf('Trasformata sin su %d punti', Nfft(i)))
    xlabel('f[Hz]')
    ylabel('f_s Normal')
    xlim([-(Nfft(i)/2)*ris_f (Nfft(i)/2)*ris_f])
    % cosi ho come limiti f Nyquist
    hold off;
    
    hold on;
    figure(i+2+length(Nfft));
    plot(f,2*abs(f_shiftcos)); 
    title(sprintf('Trasformata cos su %d punti', Nfft(i)))
    xlabel('f[Hz]')
    ylabel('f_s Normal')
    xlim([-(Nfft(i)/2)*ris_f (Nfft(i)/2)*ris_f])
    % cosi ho come limiti f Nyquist
    hold off;
end

%% Commento Generale

% Come si può notare dai grafici, impostando una finestratura su 200 punti
% si ha una corrisponenza della fuzione di trasferimento con i Delta di
% Dirac che rappresenterebbero la funzione periodica se si fosse fatto uso
% di una Trasformata di Furier teorica, questo perche si utilizza una
% frequenza a 10 Hz per il campionamento e la funzione ha la sua armonica
% principale a 2 Hz, quindi si ha una perfetta finestatura di 10 periodi.
% Già con 199 e 201 punti si può vedere una distorione nella funzione anche
% se non sono ancora presenti i lobi laterali, mentre con 256 si possono
% vedere lobi laterali.
 

%% ESERCIZIO 2

clear;
close all;
clc;

Fs = 1;
Ts = 1/Fs;
Nfft = 500;
ris_f = 1/(Nfft*Ts);
imp = [10 30];

% Campionamento

for i = 1 : length(imp)
    
    t = -(imp(i)-imp(i)/5) : Ts : 2*imp(i)-imp(i)/5; 
    s = square((t-1)*3/imp(i))*0.5 + 0.5; % per renderlo 0 e 1
    
    % queste due funzioni mi permettono di settare una finestra
    % rettangolare per il campionamento
    
    f_s = fft(s,Nfft);
    f_shift = fftshift(f_s);
    
     if mod(Nfft,2) == 0
        f = (-Nfft/2 : 1 : (Nfft/2)-1) * ris_f;
    else
        f = (-(Nfft-1)/2 : 1 : (Nfft-1)/2) * ris_f;
    end
    
    hold on;
    figure(i);
    subplot(2,1,1);
    plot(t,s,'r','Marker','X');
    xlabel('t[s]');
    ylabel('Impulse');
    title('Finestra di Campionamento');
    subplot(2,1,2);
    plot(f,abs(real(f_shift)));
    xlabel('f[Hz]');
    ylabel('f_s');
    hold off;

end

% ricordare che all'aumentare della lunghezza della finestra di
% campionamento (tempo dell' impulso) le gobbe laterali diventano sempre
% più piccole

% Maggiore è la lunghezza temporale della finestra, minore è la durata
% frequenziale, quindi all'aumentare di N i lobi diventano più piccoli
% in fase ma diventano più grandi in ampiezza in quanto l'area deve
% essere la medesima

%% ESERCIZIO 3

clear;
close all;
clc;

% Teorema Nyquist per freq. campionamento

% Il TEOREMA DI NYQUIST asserisce che la frequenza di campionamento
% necessaria affinche non ci sia ALIASING (sovrapposizione delle 
% ripetioni del segnale in frequnza) deve essere maggiore del doppio della
% frequenza massima che il segnale possiede

f1 = 20;
f2 = 70;
Fs = 3 * max(f1,f2); % > 2*fmax
Ts = 1/Fs;

t = -10 : Ts : 10;

% Segnale Composto

s_1 = 1.1*sin(2*pi*f1*t);
s_2 = 2.5*sin(2*pi*f2*t);

% nel DFT il picco più alto corrisponde a modulo maggiore

r = randn([1 length(t)]);

s = s_1 + s_2 + 1*r; % se rumore grande si perde la frequenza corretta

hold on;
figure(1);
plot(t,s);
title('Segnale composto');
xlabel('t[s]');
hold off;

% Trasformata Fourier, come numero opportuno di punti considero
% 2800 punti in quanto è multiplo di f1 e f2 

Nfft = 2800;
ris_f = 1/(Nfft * Ts);

f_s = fft(s,Nfft); % [0:(n/2)-1,0,(-n/2)+1:0]
f_sr = fft(r,Nfft);


% faccio la mia versione di ffshift

 if mod(Nfft,2) == 0
        f = (0 : 1 : (Nfft/2)-1) * ris_f; % calcola f_Nyquist in posizione
        f_scalata = f_s(1:(Nfft/2))/Nfft;
        f_srcalata = f_sr(1:(Nfft/2))/Nfft;

    else
        f = (0 : 1 : (Nfft-1)/2) * ris_f; % no f_Nyquist
        %(Nfft - 1) rende direttamente pari
        f_scalata = f_s(1:((Nfft-1)/2)+1)/Nfft;
        f_srcalata = f_sr(1:((Nfft-1)/2)+1)/Nfft;
 end
 
 hold on;
 figure(2);
 plot(f,2*abs(real(f_scalata)),f,abs(real(f_srcalata)),'y');
 %plot(f,abs(real(f_srcalata)),'y');
 title('DFT'),
 xlabel('f[Hz]');
 ylabel('f_scalata');
 legend('Signal','Noise');
 hold off;
 
 % Con il rndn aggiungimo del runore che si può notare anche nel DFT, in
 % particolare nel grafico cerco di avere sempre un pò di margine rispetto
 % alla frequenza maggiore del segnale


