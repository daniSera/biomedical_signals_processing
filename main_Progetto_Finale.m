%% PROGETTO FINALE

data = importdata('DatiECGPPG _2minuti.csv');
col_h    = data.colheaders;
d_time   = data.data(:,1);
d_ECGone = data.data(:,2);
d_ECGtwo = data.data(:,3);
d_PPGsx  = -data.data(:,4);
d_PPGdx  = -data.data(:,5);

F = 1000;
tlim = 1:F*2;

%% Parte 1

%% DATI GREZZI

% Effettuo un controllo sui dati per vedere se ci sono anomalie visibili 
% ad occhio nudo. Passo al plot dei dati grezzi.

close(figure(1));
doubleplot(1,'PPG',col_h(1),col_h(4),col_h(5),d_time(tlim),d_PPGsx(tlim),d_PPGdx(tlim));
close(figure(2));
doubleplot(2,'ECG',col_h(1),col_h(2),col_h(3),d_time(tlim),d_ECGone(tlim),d_ECGtwo(tlim));

%% FFT

% Effettuo un analisi in frequenza per vedere se ci sono delle frequenze
% che prevalgono su altre, in particolare non si notano frequenze anomale
% quindi non è necessario filtrare il segnale in questa sezione.

[fPPGsx,fftPPGsx,risPPGsx] = discreteFt(d_PPGsx,F,2^nextpow2(length(d_PPGsx)),0);
[~,fftPPGdx,~] = discreteFt(d_PPGdx,F,2^nextpow2(length(d_PPGdx)),0);

plim = 2 : round(60/risPPGsx); % fino 60 Hz

close(figure(3));
doubleplot(3,'PPG','freq[Hz]',col_h(4),col_h(5),fPPGsx(plim),2*abs(fftPPGsx(plim)),2*abs(fftPPGdx(plim)));

[fECGone,fftECGone,risECGone] = discreteFt(d_ECGone,F,2^nextpow2(length(d_ECGone)),0);
[~,fftECGtwo,~] = discreteFt(d_ECGtwo,F,2^nextpow2(length(d_ECGtwo)),0);

elim = 2 : round(60/risECGone); % fino 60 Hz

close(figure(4));
doubleplot(4,'ECG','freq[Hz]',col_h(2),col_h(3),fECGone(elim),2*abs(fftECGone(elim)),2*abs(fftECGtwo(elim)));

%% ================ Algoritmo PAN_TOMPKIN segnale ECG =====================

% Passo all'Analisi del segnale con Algoritmo PAN_TOMPKIN.
% Il primo passo è quello di filtrare il segale con passa Banda a 5-15 Hz,
% per poi passare ad un filtro derivatore.
% E' utile dire che il campionamento è stato fatto con frequenza 1000 Hz 
% quindi noi potremmo andare a rappresentare fino a 500 Hz per evitare 
% il fenomeno dell'Aliasing.

% Il filtro di Pan Tompkins utilizza un derivatore ed un passabanda in 
% serie per la rilevazione dei complessi QRS di un segnale ECG.



%% DesignFilt PassaBanda

filt = designfilt('bandpassiir','FilterOrder',4,'HalfPowerFrequency1',0.01,'HalfPowerFrequency2',0.03);

% I filtri IIR, tipologia che viene utilizzata in questa sezione per 
% sviluppare il filtro voluto, sono caratterizzati da una risposta
% impulsiva infinita

% I filtri IIR sono tipicamente più efficienti sul piano matematico 
% rispetto ai filtri FIR e richiedono meno risorse hardware per 
% l’implementazione

% NORMALIZED FREQUENCY corrisponde ad una divisione per la frequenza di
% Nyquist rispetto alla frequenza di campionamento

% quindi 5 Hz e 15 Hz corrispondono a x / 500
% 0.01 e 0.03

d_ECGone_b = filtfilt(filt,d_ECGone);
d_ECGtwo_b = filtfilt(filt,d_ECGtwo);

close(figure(5));
doubleplot(5,'ECG Design',col_h(1),col_h(2),col_h(3),d_time(tlim),d_ECGone_b(tlim),d_ECGtwo_b(tlim));

% come mi aspettavo il DesignFilt è molto simile al Butter esendo entrambi
% basati sul Butterworth

fvtool(filt)

%% Filtro Derivatore

% Il FILTRO DERIVATORE

b_d = [-1 8 0 1 -8] * (F/12); %  [1 8 0 -1 -8]

d_ECGone_D = filtfilt(b_d,1,d_ECGone_b);
d_ECGtwo_D = filtfilt(b_d,1,d_ECGtwo_b);

close(figure(7));
doubleplot(7,'ECG Derivatore',col_h(1),col_h(2),col_h(3),d_time(tlim),d_ECGone_D(tlim),d_ECGtwo_D(tlim));

% in questo modo sono stati accentuati i picchi, cosi da poter selezionare
% meglio gli intervalli.

fvtool([b_d, 1])

%% Segnale alla seconda

% questo permette di far vedere ancora meglio i picchi dominanti

d_ECGone_D = d_ECGone_D .^2;
d_ECGtwo_D = d_ECGtwo_D .^2;

close(figure(9));
doubleplot(9,'ECG Squared',col_h(1),col_h(2),col_h(3),d_time(tlim),d_ECGone_D(tlim),d_ECGtwo_D(tlim));

%% Smooth

% Vado a provare differenti valori di span per vedere quali si addice
% meglio, va ricordato che smooth esegue una MEDIA MOBILE che
% intrinsecamente causa un ritardo nel segnale filtrato. Essendo necessario
% trovare i picchi del segnale originario bisogna fare attenzione a questo
% delay

% smooth restituisce il vettore di valori filtrato tramite filtro a
% MEDIA MOBILE secondo la funzione
% y[n] = (x[n] + x[n-1] + x[n-2] +...+ x[n-M+1])/M
% dove M è il valore di span

close(figure(10));
hold on;
figure(10);
plot(d_time(tlim),d_ECGone_D(tlim));
title('ECG');
xlabel(col_h(1));
ylabel('Segnale[V]');
legend(col_h(2));
hold off;

close(figure(11));
hold on;
figure(11);
plot(d_time(tlim),d_ECGtwo_D(tlim));
title('ECG');
xlabel(col_h(1));
ylabel('Segnale[V]');
legend(col_h(3));
hold off;

span = [10 50 100 150 500];

for i = 1:length(span)
    
    d_ECGone_D_sm = smooth(d_ECGone_D,span(i));
    d_ECGtwo_D_sm = smooth(d_ECGtwo_D,span(i));
    
    hold on;
    
    figure(10);
    plot(d_time(tlim),d_ECGone_D_sm(tlim),'DisplayName',"Span = " + span(i));
    
    hold off;
    
    hold on;
    
    figure(11);
    plot(d_time(tlim),d_ECGtwo_D_sm(tlim),'DisplayName',"Span = " + span(i));
        
    hold off;
end

% Come si puo vedere è il caso di prendere uno span pari a 100 o 150 per
% avere il segnale senza troppe oscillazioni. (accettabile per entrambi)

d_ECGone_D_sm = smooth(d_ECGone_D,100);
d_ECGtwo_D_sm = smooth(d_ECGtwo_D,100);
 
%% Media Mobile mediante Convoluzione
 
% E' possibile eseguire anche la media mobile mediante la convoluzione del
% seganle ed una finestra rettangolare, la prendo abbastanza grande in
% modo che si abbia una migliore rappresentazione del segnale.
 
dely = 35;
square = ones(1 ,dely)/(dely); % Normalizzato

d_ECGone_D_c = conv(d_ECGone_D,square,"same"); %zero Padding 
d_ECGtwo_D_c = conv(d_ECGtwo_D,square,"same"); %zero Padding 

close(figure(12));
hold on;
figure(12);
plot(d_time(tlim),d_ECGone_D_sm(tlim),d_time(tlim),d_ECGone_D_c(tlim));
legend('Smooth','Convoluzione');
hold off;

% Il risultato mediante Convoluzione sembra più accurato, proseguo con
% qiesto segnale filtrato

%% Rilevazione Picchi Segnale Filtrato

% Viene Utilizzata la Funzione findpeaks Matlab, come segnale utilizzo
% quello filtrato con la Convoluzione. Essendo stato inserito un ritardo, è
% necessario cercare in un intorno del picco filtrato, quindi utilizzo
% MINPEAKPROMINENCE come suggerito da Matlab

% The prominence of a peak measures how much the peak stands out due to its
% intrinsic height and its location relative to other peaks. A low isolated
% peak can be more prominent than one that is higher but is an otherwise 
% unremarkable member of a tall range.

[~,locs_one] = findpeaks(d_ECGone_D_c,d_time,'MinPeakProminence',30);
[~,locs_two] = findpeaks(d_ECGtwo_D_c,d_time,'MinPeakProminence',30);

close(figure(13));
hold on;
figure(13);
findpeaks(d_ECGone_D_c,d_time,'MinPeakProminence',30,'Annotate','peaks')
title(col_h(2));
xlabel(col_h(1));
ylabel('Segnale[V]');
legend('Signal','Peaks');
hold off;

% overkill

%% Rilevazione Picchi segnale originario

% Nel vettore locs sono note le posizioni dei massimi del segnale filtrato.
% Essendo la frequenza di campionamento F, i valori corrispondenti a locs
% nel segnale sono in posizione (locs[i] - t[0]) * F. Va preso come minimo
% un buffer di 40 elementi
% E' da ricordare che con la media mobile il segnale viene ritardato,
% utilizzo quindi il ritardo come finestra per visualizzare

loc_one = [];
loc_two = [];
val_s = 0;
for i = 1:length(locs_one)
 
    max = round((locs_one(i) - d_time(1)) * F);
    
    buff = d_ECGone( max-20 : max+20);
    
    [val,point]=findpeaks(buff,'MinPeakDistance',39);
    if val ~= val_s 
         loc_one = [loc_one (max - 20 + point -1)];
    end 
    val_s = val;
end


close(figure(15));
hold on;
figure(15);
plot(d_time,d_ECGone,d_time(loc_one),d_ECGone(loc_one),'or');
title('Peaks del segnale originario')
xlabel(col_h(1))
ylabel(col_h(2))
hold off;

val_s = 0;
for i = 1:length(locs_two)
   
    max = round((locs_two(i) - d_time(1)) * F); 
    
    buff = d_ECGtwo( max-20 : max+20);
    
    [val,point]=findpeaks(buff,'MinPeakDistance',39);
     if val ~= val_s
        loc_two = [loc_two (max - 20 + point -1)]; 
     end  
    val_s = val;
end

close(figure(16));
hold on;
figure(16);
plot(d_time,d_ECGtwo,d_time(loc_two),d_ECGtwo(loc_two),'or');
title('Peaks del segnale originario')
xlabel(col_h(1))
ylabel(col_h(3))
hold off;

%% ================= Algoritmo PAN_TOMPKIN segnale PPG ====================

%% Butterworth 

% Applico un filtro butterworth passalto in modo da lasciare intatte le
% onde PPG

% Il Filtro BUTTERWORTH è un filtro "maximally flat" che offre una risposta
% in frequenza il più possibile piatta nella banda passante, non presenta
% ripple ma ha una banda di transizione molto ampia quindi è necessario
% controllare di settare la frequenza di taglio giusta

f1 = 0.7;
f2 = 10;
F_t = [f1 f2] * 2 / F;

% Frequenze normalizzate necessarie per il butter. Utilizo ordine 4 come
% sviluppato a lezione

[b,a] = butter(4,F_t);
d_PPGsx_B = filtfilt(b,a,d_PPGsx);
d_PPGdx_B = filtfilt(b,a,d_PPGdx);

close(figure(17));
doubleplot(17,'PPG Butter',col_h(1),col_h(4),col_h(5),d_time(tlim),d_PPGsx_B(tlim),d_PPGdx_B(tlim));

%% Rilevazione Picchi segnale originario

[~,locs_sx] = findpeaks(d_PPGsx_B,d_time,'MinPeakProminence',0.001);
[~,locs_dx] = findpeaks(d_PPGdx_B,d_time,'MinPeakProminence',0.001);

close(figure(18));
hold on;
figure(18);
findpeaks(d_PPGsx_B,d_time,'MinPeakProminence',0.001,'Annotate','peaks')
title(col_h(4));
xlabel(col_h(1));
ylabel('Segnale[V]');
legend('Signal','Peaks');
hold off;
% sviluppo anche qua il caso precedente

loc_sx = [];
loc_dx = [];
val_s = 0;
for i = 1:length(locs_sx)
   
    max = round((locs_sx(i) - d_time(1)) * F);
    
    buff = d_PPGsx(max-20:max+20);
    
    [val,point]=findpeaks(buff,'MinPeakDistance',39);
     if val ~= val_s
        loc_sx = [loc_sx (max - 20 + point-1)]; 
     end  
    val_s = val;
end


val_s = 0;
for i = 1:length(locs_dx)
   
    max = round((locs_dx(i) - d_time(1)) * F);
    buff = d_PPGdx(max-20:max+20);
    
    [val,point]=findpeaks(buff,'MinPeakDistance',39);
     if val ~= val_s
        loc_dx = [loc_dx (max - 20 + point -1)]; 
     end  
    val_s = val;
end


close(figure(19));
figure(19);
plot(d_time,d_PPGsx,d_time(loc_sx),d_PPGsx(loc_sx),'or');
title('Peaks del segnale originario')
xlabel(col_h(1))
ylabel(col_h(4))

close(figure(20));
figure(20);
plot(d_time,d_PPGdx,d_time(loc_dx),d_PPGdx(loc_dx),'or');
title('Peaks del segnale originario')
xlabel(col_h(1))
ylabel(col_h(5))

%% ========== Tacogramma =================================================

% Creo un Verrote con gli R-R Intervals, il numero di battiti è dato dalla 
% lunghezza del vettore dei massimi

d_RRIone = [(d_time(loc_one(2)) - d_time(loc_one(1)))];
d_RRItwo = [(d_time(loc_two(2)) - d_time(loc_two(1)))];

for i = 2:length(loc_one)
    
    d_RRIone = [d_RRIone (d_time(loc_one(i)) - d_time(loc_one(i-1)))];    
    
end

for i = 2:length(loc_two)
    
    d_RRItwo = [d_RRItwo (d_time(loc_two(i)) - d_time(loc_two(i-1)))];    
    
end 

bat_one = (0:1:length(loc_one)-1);
bat_two = (0:1:length(loc_two)-1);

d_RRIone = d_RRIone*1000;
d_RRItwo = d_RRItwo*1000;

close(figure(21));
hold on;
figure(21);
plot(bat_one,d_RRIone,bat_two,d_RRItwo);
title('Serie temporali ECG')
legend('ECG 1','ECG 2');
xlabel('Battiti')
ylabel('RRI [ms]')
hold off;


%% ===== PPI ============================================================

d_PPIsx = [(d_time(loc_sx(2)) - d_time(loc_sx(1)))];
d_PPIdx = [(d_time(loc_dx(2)) - d_time(loc_dx(1)))];

for i = 2:length(loc_sx)
    
    d_PPIsx = [d_PPIsx (d_time(loc_sx(i)) - d_time(loc_sx(i-1)))];    
    
end

for i = 2:length(loc_dx)
    
    d_PPIdx = [d_PPIdx (d_time(loc_dx(i)) - d_time(loc_dx(i-1)))];    
    
end 

bat_sx = (0:1:length(loc_sx)-1);
bat_dx = (0:1:length(loc_dx)-1);

d_PPIsx = d_PPIsx*1000;
d_PPIdx = d_PPIdx*1000;


close(figure(22));
figure(22);
plot(bat_sx,d_PPIsx,bat_dx,d_PPIdx);
title('Serie temporali ECG')
legend('PPG sx','PPG dx');
xlabel('Battiti')
ylabel('PPI [ms]')

%% === Heart Rate Diagram =============================================

% noto il tempo tra ogni battito è possibile determinare la frequenza
% caridaca

freq_one = 60 * 1000  ./ d_RRIone;
freq_two = 60 * 1000  ./ d_RRItwo;
freq_sx = 60 * 1000  ./ d_PPIsx;
freq_dx = 60 * 1000  ./ d_PPIdx;


close(figure(23));
figure(23);
plot(d_time(loc_one),freq_one,d_time(loc_two),freq_two);
legend('HR RRI one','HR RRI two');
title('Frequenza Caridaca ECG')
xlabel('t[s]')
ylabel('Heart Rate [bpm]')

close(figure(24));
figure(24);
plot(d_time(loc_sx),freq_sx,d_time(loc_dx),freq_dx);
legend('HR PPI sx','HR PPI dx');
title('Frequenza Caridaca PPG')
xlabel('t[s]')
ylabel('Heart Rate [bpm]')

%% Considerazioni Statistiche RRI e PPI

statistics(d_RRIone,d_RRItwo,d_PPIsx,d_PPIdx,'Considerazioni statistiche tra RRI e PPI');

%% Considerazioni Statistiche Heart Rate

statistics(freq_one,freq_two,freq_sx,freq_dx,'Considerazioni statistiche Heart Rate segnali');
%% Parte 2


%% =========== Analisi Spettrale =======================================

% L'analisi spettrale è stata effettuata utilizzando 

% Nella scelta del lag ho selezionato un lag maggiore per i segnali ECG in
% modo da far risaltare di più i massimi che volevamo trovare nel grafico
% della potnza.

acs_RRIone = acsF((d_RRIone-mean(d_RRIone))/1000,40,0,1);
acs_RRItwo = acsF((d_RRItwo-mean(d_RRItwo))/1000,40,0,1);
acs_PPIsx = acsF((d_PPIsx-mean(d_PPIsx))/1000,35,0,1);
acs_PPIdx = acsF((d_PPIdx-mean(d_PPIdx))/1000,35,0,1);

% preferibile utilizzare lo stimatore polarizzato

% provo le varie finestre su d_RRIone

acs_RRIonefin = finestratura(acs_RRIone'); 
acs_RRItwofin = finestratura(acs_RRItwo');
acs_PPIsxfin = finestratura(acs_PPIsx');
acs_PPIdxfin = finestratura(acs_PPIdx');

% Avedo selezionato un lag pari a 40, le finestrature di Bartlett e Hann
% sono molto simili, come possiamo vedere nei seguenti grafici 
    
[spettro_RRIone,fone]=spettro(acs_RRIonefin,25,'RRI 1');
[spettro_RRItwo,ftwo]=spettro(acs_RRItwofin,26,'RRI 2');
[spettro_PPIsx,fsx]=spettro(acs_PPIsxfin,27,'PPI sx');
[spettro_PPIdx,fdx]=spettro(acs_PPIdxfin,28,'PPI dx');

plotarea(d_RRIone,spettro_RRIone(:,3),fone,29,'RRI 1');
plotarea(d_RRItwo,spettro_RRItwo(:,4),ftwo,30,'RRI 2');
plotarea(d_PPIsx,spettro_PPIsx(:,2),fsx,31,'PPI sx');
plotarea(d_PPIdx,spettro_PPIdx(:,2),fdx,32,'PPI dx');

% plotarea mi permette di ottenere la potenza del segnale in un range di
% frequenze come l'integrale dello spettro di potenza. In questa analisi
% sono utili le potenze da 0.003 a 0.04 Hz definite Very Low Frequency, da
% 0.04 a 0.15 Hz definite Low Frequency e da 0.15 a 0.4 Hz definite High 
% Frequency.  

