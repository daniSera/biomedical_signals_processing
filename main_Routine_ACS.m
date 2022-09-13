%% VERIFICA ROUTINE

% script riguardnte la teoria associata alla funzioen di AUTOCORRELAZIONE.
% Viene effettuato un confronto tra la funzione di MATLAB xcorr e una
% funzione eseguita da me per il calcolo dell'autocorrelazione. 

%% Parte 1

seq = (1:10);
lags = (-9:9);
x_funq = xcorr(seq,'unbiased'); % unbiased = non polarizzata

% xcorr con singolo vettore restituisce la funzione di autocorrelazione
% della sequenza di dati. L'AUTOCORRELAZIONE esprime la correlazione tra
% due istanti del processo in dipendenza del ritardo che li separa.

% L'Autocorrelazione fornisce informazioni sulla media del prodotto dei 
% campioni delle realizzazioni dello stesso processo in due istanti
% prefissati

acs_funq = acsF(seq,9,1,1);

% acsF è la funzione di Autocorrelazione sviluppata come esercizio e
% confrotata con la funzione a disposizione su MATLAB

close(figure(1));
figure(1);
hold on;
plot(lags,x_funq)
plot(lags,acs_funq,'or')
hold off;
title('Differenza tra Xcorr e ACS "unbiased"')
ylabel('Autocorrelazione [r_m]')
xlabel('Lags')
legend('xcorr','acs')

x_func = xcorr(seq,7,'unbiased');
acs_func = acsF(seq,7,1,1);
lags = (-7:7);

close(figure(2));
figure(2);
hold on;
plot(lags,x_func)
plot(lags,acs_func,'or')
hold off;
title('Differenza tra Xcorr e ACS "unbiased"')
ylabel('Autocorrelazione [r_m]')
xlabel('Lags')
legend('xcor','acs')


%% Parte 2

time = (0:1/100:5);
serie = sin(2*pi*15*time);

for i = 2 : 10
    
set = round(length(serie)/i);

x_func = xcorr(serie,set,'unbiased');
acs_func = acsF(serie,set,1,1);

close(figure(i+1));
figure(i+1);
hold on;
plot((-set:set),x_func);
plot((-set:set),acs_func,'or');
title(['Xcorr e ACS "unbiased" su segnale sinusoidale con ',num2str(set),' lags'])
ylabel('Autocorrelazione [r_m]')
xlabel('Lags')
legend('xcor','acsf')
hold off;
end

%% Parte 3 rndn signal

serie = randn(1,2000);
lag = 100;

x_func = xcorr(serie,lag,'unbiased');
acs_func = acsF(serie,lag,1,1);

close(figure(12));
figure(12);
hold on;
plot((-lag:lag),x_func);
plot((-lag:lag),acs_func,'+r');
title('Xcorr e ACS "unbiased" su segnale random')
ylabel('Autocorrelazione [r_m]')
xlabel('Lags')
legend('xcor','acs')
hold off;



%% Esercizio Finale, applicazione teorema Wiener-Khinchin Parte 1


F = 100;
time = (0:1/F:5);
serie = sin(2*pi*15*time);
lag = 50;

acs = acsF(serie,lag,0,1);

% il teorema predilige l'uso dello stimatore polarizzato

[freq,fft_acs,ris] = discreteFt(acs,F,2^nextpow2(length(serie)),0);
lim = 1 : round(F/(2*ris));

close(figure(13));
figure(13);
plot(freq(lim),abs(fft_acs(lim)));
title('Metodo indiretto su segnale deterministico')
ylabel('FFT absolute value')
xlabel('frequency (Hz)')

%cambio ritardo

for i = 2 : 5
   lags = round(lag/i); 
   acs = acsF(serie,lags,0,1);
   [freq,fft_acs,ris] = discreteFt(acs,F,2^nextpow2(length(serie)),0);
   lim = 1 : round(F/(2*ris));
   close(figure(12 + i));
    figure(12 + i);
    plot(freq(lim),abs(fft_acs(lim)));
    title(['Metodo indiretto su segnale deterministico con ',num2str(lags),' lags'])
    ylabel('FFT absolute value')
    xlabel('frequency (Hz)')
end

%% Parte 2

acs = acsF(serie,lag,0,1);

%  il teorema predilige l'uso dello stimatore polarizzato

acs_rect = acs.*rectwin(length(acs))'/mean(rectwin(length(acs)));
acs_bartlett = acs.*bartlett(length(acs))'/mean(bartlett(length(acs)));
acs_hann = acs.*hann(length(acs))'/mean(hann(length(acs)));
acs_hamming = acs.*hamming(length(acs))'/mean(hamming(length(acs)));

acs_mat = [acs_rect ; acs_bartlett ; acs_hann ; acs_hamming];

close(figure(18));
figure(18);
hold on;

for i = 1 : 4
    
    [acs_f,acs_fft,ris] = discreteFt(acs_mat(i,:),F,3^nextpow2(length(acs_mat(i,:))),0);
    lim = 1 : round(F/(2*ris));
    plot(acs_f(lim),abs(acs_fft(lim))); 
    
end

title('Metodo indiretto su segnale deterministico con finestre')
ylabel('FFT absolute value')
xlabel('frequency (Hz)')
legend('retwin','bartlett','hann','hamming');
hold off;

% Qua si può notare l'importanza della finestratura quando si svolgono le
% Trasformate di Fourier di segnali campionati e limitati. Una finestratura
% più dolce sui lati risulta in lobi laterali molto meno pronunciati, cioè
% in una distorsione molto meno significativa.

