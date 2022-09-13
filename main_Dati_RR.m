clc;
clear;
close;

%% Import Data

dataTot = importdata('Rest.txt');


%% Controllo Distribuzione


dataRR = sort(dataTot(:,3)); 
xdata = min(dataRR):10:max(dataRR);

hold on;
figure(1);
qqplot(dataRR)
hold off;

% qqplot restitutisce il grafico che rappresenta il confronto tra i
% quantili delle funzioni distribuzione di probabilità dei dati inseriti 
% e di una teroetica distribuzione normale che meglio fitta i dati forniti
% (che giustamente appare lineare)

pdfRR = normpdf(dataRR,mean(dataRR),std(dataRR));

% normpdf calcola la funzione densità di probabilità della distribuzione
% normale o gaussiana, 
% restituisce la funzione densità di probabilità calcolata sui
% valori assegnati e con parametri di media e deviazione standard forniti
% alla funzione

hold on;
figure(2);
Plot_obj = plot(dataRR, pdfRR,'ok');
title('Probability Density Function')
xlabel('dataRR')
obj = gca;
hold off;

%% Mean, Median, Mode

meanRR = mean(dataRR); 
% media (stima massima verosomiglianza)

medianRR = median(dataRR);
% mediana (termine preceduto e seguito dallo stesso numero di termini)

modeRR = mode(dataRR); 
% moda (valore più probabile o , meglio, più frequente)

%% Andamento con Indici Statistici

hold on;
figure(3);
copyobj(obj,figure(3));
xline(meanRR,'-r','Mean','linewidth',2);
xline(medianRR,'-r','Median','linewidth',2);
xline(modeRR,'-r','Mode','linewidth',2);
hold off;

%% Grafico Istogramma

hold on;
figure(4);
histfit(dataRR,length(dataRR)/20,'normal');
xline(meanRR,'-c','Mean','linewidth',2);
xline(medianRR,'-c','Median','linewidth',2);
xline(modeRR,'-c','Mode','linewidth',2);
title('Data RR');
legend('Data','Distribution');
obj1 = gca;
hold off;

% histfit restitutisce un istogramma rappresentante i dati inseriti e la
% rispettiva funzione densità di proabilità

% un istogramma rappresenta una distribuzione in frequenza di occorrenza 
% dei dati campionari 

%% Grafico Istogramma Completo

quart = prctile(dataRR ,[0 25 50 75 100]);

% prctile restituisce i percentili indicati nel vettore dei valori forniti
% alla funzione. I percentili rappresentano la probabilità associata a
% quella percentuale definita. I quartili sono dei specifici percentili

IQR = quart(4) - quart(2); 

% Range tra primo e terzo quartile, distanza interquantilica

out_up = quart(4) + 1.5 * IQR;
% valore più alto non outlier (valori "fuori istribuzione normale")

out_down = quart(2) - 1.5 * IQR;
% valore più basso non outlier

hold on;
figure(5);
copyobj(obj1,figure(5));
xline(out_down,'-g','Lower Limit','linewidth',2);
xline(quart(2),'-g','1Q','linewidth',2);
xline(quart(4),'-g','3Q','linewidth',2);
xline(out_up,'-g','Upper Limit','linewidth',2);
title('Data RR')
legend('Data','Distribution');
hold off;

%% Diversa Visualizzazione

dataRR_pure = dataTot(:,3);
figure(6);
plot(1:length(dataRR),dataRR_pure);
yline(out_down,'-g','Lower Limit','linewidth',2);
yline(quart(2),'-g','1Q','linewidth',2);
yline(quart(4),'-g','3Q','linewidth',2);
yline(out_up,'-g','Upper Limit','linewidth',2);
title('Data RR')

%% Grafico Box and Whisker

hold on;
figure(7);
boxplot(dataRR);
title('Box Plot Data RR');
hold off;

% oxplot restituisce il boxplot dei dati inseriti, la linea rossa centrale
% rappresenta la mediana mentre il box rappresenta il primo e terzo
% quartile, le "code" rappresentano i valori massimi non outliners mentre
% le croci rosse rappresentano i valori "esterni alla distribuzione normale"

%soglia superiore
%terzo quartile
%mediana
%primo quartile
%soglia inferiore








