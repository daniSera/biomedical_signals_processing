%% ESERCITAZIONE FILTRI

%% 1)Visualizzazione Segnali Grezzi

clear;
clc;
close all;

data = importdata('datiECG_PPG.mat');
d_ECG = data.ECG;
d_PPG = data.PPG;
t = data.t;
clear data;

F_s = 1000;

hold on;
figure(1);

subplot(2,1,1);
plot(t(1:F_s*10),d_ECG(1:F_s*10));
title('ECG');
xlabel('Tempo[s]');
ylabel('Segnale[V]');

subplot(2,1,2);
plot(t(1:F_s*10),d_PPG(1:F_s*10));
title('PPG');
xlabel('Tempo[s]');
ylabel('Segnale[V]');

hold off;

%% 2)Visualizzazione FFT segnali grezzi

[fECG,fftECG,risECG] = discreteFt(d_ECG,F_s,2^nextpow2(length(d_ECG)),0);
[fPPG,fftPPG,risPPG] = discreteFt(d_PPG,F_s,2^nextpow2(length(d_PPG)),0);

% La funzione discreteFt mi permette di valutare la FFT del vettore di dati
% inserito e mi riorganizza la rappresentazione in frequenza

lim = 2 : round(60/risECG); % voglio rappresentare fino 60 Hz
limzm = round(49/risECG) : round(51/risECG); % zoom tra 49 e 51 Hz
hold on;
figure(2);

subplot(2,1,1);
plot(fECG(lim),2*abs(fftECG(lim)));
% tolgo valore del segnale originario in 0 perchè a frequenza nulla
title('FFT-ECG');
xlabel('Frequency[Hz]');
ylabel('|Y(f)|');

subplot(2,1,2);
plot(fPPG(lim),2*abs(fftPPG(lim)));
% tolgo valore del segnale originario in 0 perchè a frequenza nulla
title('FFT-PPG');
xlabel('Frequency[Hz]');
ylabel('|Y(f)|');

hold off;

% Come si puo vedere nel grafico dell'ECG si vede un picco in
% corrispondenza dei 50 Hz, nota interferenza causata dalla linea di
% alimentazione elettrica

%% 3)ApplicazioneFiltro a Media Mobile

% E' utile ricordare che la media mobile induce un ritardo nel segale
% filtrato che deve essere considerato in casi specifici

hold on;
figure(3);

plot(t(F_s:F_s*2),d_ECG(F_s:F_s*2)); 
title('ECG');
xlabel('Tempo[s]');
ylabel('Segnale[V]');
legend('Dati Grezzi');

hold off;

hold on;
figure(4);

plot(t(F_s:F_s*2),d_PPG(F_s:F_s*2));
title('PPG');
xlabel('Tempo[s]');
ylabel('Segnale[V]');
legend('Dati Grezzi');

hold off;

span = [11 41 71]; % numero di campioni da mediare

for i = 1:length(span)
    
    d_ECG_sm = smooth(d_ECG,span(i));
    d_PPG_sm = smooth(d_PPG,span(i));
     
    % smooth restituisce il vettore di valori filtrato tramite filtro a
    % MEDIA MOBILE secondo la funzione
    % y[n] = (x[n] + x[n-1] + x[n-2] +...+ x[n-M+1])/M
    % dove M è il valore di span
    
    hold on;
    
    figure(3);
    plot(t(F_s:F_s*2),d_ECG_sm(F_s:F_s*2),'DisplayName',"Span = " + span(i));
    
    hold off;
    
    hold on;
    
    figure(4);
    plot(t(F_s:F_s*2),d_PPG_sm(F_s:F_s*2),'DisplayName',"Span = " + span(i));
        
    hold off;
end

% per gli ECG la scelta ottima risulta essere span = 41, mentre per gli PPG
% basta applicare smooth con uno dei qualsiasi span

%% 4)Visualizzazione della FFT con dati filtrati con Smooth

span = 41;

d_ECG_sm = smooth(d_ECG,span);
d_PPG_sm = smooth(d_PPG,span);

[fECG_sm,fftECG_sm,risECG_sm] = discreteFt(d_ECG_sm,F_s,2^nextpow2(length(d_ECG)),0);
[fPPG_sm,fftPPG_sm,risPPG_sm] = discreteFt(d_PPG_sm,F_s,2^nextpow2(length(d_PPG)),0);

figure(5);
subplot(2,1,1);
hold on;
plot(fECG(limzm),2*abs(fftECG(limzm)));
plot(fECG_sm(limzm),2*abs(fftECG_sm(limzm)));
hold off;
title('FFT-ECG');
xlabel('Frequency[Hz]');
ylabel('|Y(f)|');
legend({"Dati Grezzi","Span = " + span})


subplot(2,1,2);
hold on;
plot(fPPG(lim),2*abs(fftPPG(lim)));
plot(fPPG_sm(lim),2*abs(fftPPG_sm(lim)));
hold off;
title('FFT-PPG');
xlabel('Frequency[Hz]');
ylabel('|Y(f)|');
legend({"Dati Grezzi","Span = " + span})

% come si può vedere il valore a 50 Hz è stato eliminato senza eliminare
% altro

%% 5)Applicazione filtro Butterworth passabasso

% Il Filtro BUTTERWORTH è un filtro "maximally flat" che offre una risposta
% in frequenza il più possibile piatta nella banda passante, non presenta
% ripple ma ha una banda di transizione molto ampia quindi è necessario
% controllare di settare la frequenza di taglio giusta

% I filtri IIR, tipologia che viene utilizzata in questa sezione per 
% sviluppare il filtro voluto, sono caratterizzati da una risposta
% impulsiva infinita

hold on;
figure(6);

plot(t(1:F_s*1),d_ECG(1:F_s*1));
title('ECG');
xlabel('Tempo[s]');
ylabel('Segnale[V]');
legend("Dati Grezzi");

hold off;

hold on;
figure(7);

plot(t(1:F_s*1),d_PPG(1:F_s*1));
title('PPG');
xlabel('Tempo[s]');
ylabel('Segnale[V]');
legend("Dati Grezzi");

hold off;

%designfilt da usare per progettare il filtro

freq = [20 15 10]; % frequenze di taglio

for i = 1:length(freq)

    filter = designfilt('lowpassiir', 'FilterOrder', 4, 'HalfPowerFrequency', freq(i)/(F_s/2));
    
    d_ECG_fi = filtfilt(filter,d_ECG);
    d_PPG_fi = filtfilt(filter,d_PPG);
    
    hold on;
    
    figure(6);
    plot(t(1:F_s*1),d_ECG_fi(1:F_s*1),'DisplayName',"Freq = " + freq(i) + "Hz");
    
    hold off;
    
    hold on;
    
    figure(7);
    plot(t(1:F_s*1),d_PPG_fi(1:F_s*1),'DisplayName',"Freq = " + freq(i) + "Hz");
        
    hold off;
end

% fvtool(filter); % per vedere le caratteristiche del filtro

% come si puo vedere dai grafici la scelta ottima è quella filtrata con
% frequenza di taglio a 20 Hz

%% 6)Visualizzazione della FFT dei dati filtrati con filtro passabasso

% Filtro PASSABASSO, il risultato è il medesimo, avendo settato la
% frequenza di taglio a 20 Hz si può notare che il picco a 50 Hz non viene
% più plottato.

% La caratteristica di un filtro passa basso è quella di far passare 
% indisturbate le componenti a bassa frequenza del segnale, 
% eliminando le componenti a frequenza superiore della frequenza di taglio

filter1 = designfilt('lowpassiir', 'FilterOrder', 4, 'HalfPowerFrequency', freq(1)/(F_s/2)); 
d_ECG_fi = filtfilt(filter1,d_ECG);
filter2 = designfilt('lowpassiir', 'FilterOrder', 4, 'HalfPowerFrequency', freq(2)/(F_s/2));
d_PPG_fi = filtfilt(filter2,d_PPG);

% prendo le frequenze da sopra

[fECG_fi,fftECG_fi,risECG_fi] = discreteFt(d_ECG_fi,F_s,2^nextpow2(length(d_ECG_fi)),0);
[fPPG_fi,fftPPG_fi,risPPG_fi] = discreteFt(d_PPG_fi,F_s,2^nextpow2(length(d_PPG_fi)),0);

figure(8);
subplot(2,1,1);
hold on;
plot(fECG(limzm),2*abs(fftECG(limzm)));
plot(fECG_fi(limzm),2*abs(fftECG_fi(limzm)));
hold off;
title('FFT-ECG');
xlabel('Frequency[Hz]');
ylabel('|Y(f)|');
legend({"Dati Grezzi","Freq = " + freq(1)});


subplot(2,1,2);
hold on;
plot(fPPG(lim),2*abs(fftPPG(lim)));
plot(fPPG_fi(lim),2*abs(fftPPG_fi(lim)));
hold off;
title('FFT-PPG');
xlabel('Frequency[Hz]');
ylabel('|Y(f)|');
legend({"Dati Grezzi","Freq = " + freq(2)});

%% 7)Applicazione filtro Butterworth passa alto

hold on;
figure(9);

plot(t(1:F_s*1),d_ECG(1:F_s*1));
title('ECG');
xlabel('Tempo[s]');
ylabel('Segnale[V]');
legend("Dati Grezzi");

hold off;

hold on;
figure(10);

plot(t(1:F_s*1),d_PPG(1:F_s*1));
title('PPG');
xlabel('Tempo[s]');
ylabel('Segnale[V]');
legend("Dati Grezzi");

hold off;

%designfilt da usare per progettare il filtro

FREQ = [0.05 0.1 0.2];

for i = 1:length(freq)

    FILTER = designfilt('highpassiir', 'FilterOrder', 4, 'HalfPowerFrequency', FREQ(i)/(F_s/2));
    
    d_ECG_FI = filtfilt(FILTER,d_ECG);
    d_PPG_FI = filtfilt(FILTER,d_PPG);
    
    hold on;
    
    figure(9);
    plot(t(1:F_s*1),d_ECG_FI(1:F_s*1),'DisplayName',"Freq = " + FREQ(i) + "Hz");
    
    hold off;
    
    hold on;
    
    figure(10);
    plot(t(1:F_s*1),d_PPG_FI(1:F_s*1),'DisplayName',"Freq = " + FREQ(i) + "Hz");
        
    hold off;

end

% fvtool(FILTER); % per vedere le caratteristiche del filtro

%% 8)Visualizzazione della FFT dei dati filtrati con filtro passa alto

FILTER1 = designfilt('highpassiir','FilterOrder',4,'HalfPowerFrequency',FREQ(1)/(F_s/2));
d_ECG_FI = filtfilt(FILTER1,d_ECG);
FILTER2 = designfilt('highpassiir','FilterOrder',4,'HalfPowerFrequency',FREQ(1)/(F_s/2));
d_PPG_FI = filtfilt(FILTER2,d_PPG);

% prendo le frequenze da sopra

[fECG_FI,fftECG_FI,risECG_FI] = discreteFt(d_ECG_FI,F_s,2^nextpow2(length(d_ECG_FI)),0);
[fPPG_FI,fftPPG_FI,risPPG_FI] = discreteFt(d_PPG_FI,F_s,2^nextpow2(length(d_PPG_FI)),0);

figure(11);
subplot(2,1,1);
hold on;
plot(fECG(limzm),2*abs(fftECG(limzm)));
plot(fECG_FI(limzm),2*abs(fftECG_FI(limzm)));
hold off;
title('FFT-ECG');
xlabel('Frequency[Hz]');
ylabel('|Y(f)|');
legend({"Dati Grezzi","Freq = " + freq(1)});


subplot(2,1,2);
hold on;
plot(fPPG(lim),2*abs(fftPPG(lim)));
plot(fPPG_FI(lim),2*abs(fftPPG_FI(lim)));
hold off;
title('FFT-PPG');
xlabel('Frequency[Hz]');
ylabel('|Y(f)|');
legend({"Dati Grezzi","Freq = " + freq(1)});


%% 9)Applicazione di un filtro notch

%freqz
%mag2db
%utili

% Il filtro NOTCH è un filtro ARRESTA BANDA che filtra, lasciando 
% inalterata la risposta in frequenza di un segnale, le frequenze 
% all'interno della banda specificata nel design.

f_yu = [0 45 50 55 F_s/2]/(F_s/2);
m_yu = [1 0.8 0 0.8 1];
[b_yu,a_yu] = yulewalk(20,f_yu,m_yu);

% la fuzione yulewalk restituisce i coefficienti della funzioen di 
% trasferiemnto di un filtro IIR con ordine rappresentato dal primo
% coefficiente che segue le specifiche dettate dai seguenti coefficienti:
% il primo indica le frequenze che devono essere arrestate, il secondo la
% magnitudine con cui esse devono essere arrestate

d_ECG_yu = filtfilt(b_yu,a_yu,d_ECG);
d_PPG_yu = filtfilt(b_yu,a_yu,d_PPG);

figure(12);
hold on;
plot(t(1:F_s*1),d_ECG(1:F_s*1));
plot(t(1:F_s*1),d_ECG_yu(1:F_s*1));
hold off;
title('ECG');
xlabel('Tempo[s]');
ylabel('Segnale[V]');
legend("Dati Grezzi","Notch 50 Hz");

figure(13);
hold on;
plot(t(1:F_s*1),d_PPG(1:F_s*1));
plot(t(1:F_s*1),d_PPG_yu(1:F_s*1));
hold off;
title('PPG');
xlabel('Tempo[s]');
ylabel('Segnale[V]');
legend("Dati Grezzi","Notch 50 Hz");

% fvtool([b_yu,a_yu]);


%% 10 Visualizzazione della FFT dei dati filtrati con filtro notch

[fECG_yu,fftECG_yu,risECG_yu] = discreteFt(d_ECG_yu,F_s,2^nextpow2(length(d_ECG_yu)),0);
[fPPG_yu,fftPPG_yu,risPPG_yu] = discreteFt(d_PPG_yu,F_s,2^nextpow2(length(d_PPG_yu)),0);

figure(14);
subplot(2,1,1);
hold on;
plot(fECG(limzm),2*abs(fftECG(limzm)));
plot(fECG_yu(limzm),2*abs(fftECG_yu(limzm)));
hold off;
title('FFT-ECG');
xlabel('Frequency[Hz]');
ylabel('|Y(f)|');
legend("Dati Grezzi","Notch 50 Hz");


subplot(2,1,2);
hold on;
plot(fPPG(lim),2*abs(fftPPG(lim)));
plot(fPPG_yu(lim),2*abs(fftPPG_yu(lim)));
hold off;
title('FFT-PPG');
xlabel('Frequency[Hz]');
ylabel('|Y(f)|');
legend("Dati Grezzi","Notch 50 Hz");


%% 11)Confronto filtraggio con istruzione filtfilt e filter
% filtfilt

% purtoppo il comando filter risulta avere problemi.
% Anche seguendo gli script suggeriti da Matlab il comando filter sembra
% essere stato cambiato, anche utilizzando designfilt non riesco ad
% utilizzare il comando.

d = designfilt('lowpassiir','FilterOrder',4,'HalfPowerFrequency',10/(F_s/2),'DesignMethod','butter');

[b_bu,a_bu] = butter(4,10/(F_s/2),'low');

% il comando butter mi permette di costruire un filtro NUTTERWORTH IIR come
% richiesto da specifiche inserite nella funzione

d_ECG_fifi = filtfilt(b_bu,a_bu,d_ECG);
d_ECG_fibu = filtfilt(d,d_ECG);
d_PPG_fifi = filtfilt(b_bu,a_bu,d_PPG);
d_PPG_fibu = filtfilt(d,d_PPG);

figure(15);
hold on;
plot(t(1:F_s*1),d_ECG(1:F_s*1));
plot(t(1:F_s*1),d_ECG_fifi(1:F_s*1));
plot(t(1:F_s*1),d_ECG_fibu(1:F_s*1));
hold off;
title('ECG');
xlabel('Tempo[s]');
ylabel('Segnale[V]');
legend("Dati Grezzi","filtfilt","filtfilt design");

figure(16);

hold on;

plot(t(1:F_s*1),d_PPG(1:F_s*1));
plot(t(1:F_s*1),d_PPG_fifi(1:F_s*1));
plot(t(1:F_s*1),d_PPG_fibu(1:F_s*1));
title('PPG');
xlabel('Tempo[s]');
ylabel('Segnale[V]');
legend("Dati Grezzi","filtfilt","filtfilt design");

hold off;



