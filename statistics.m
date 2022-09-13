function  statistics(ECG1, ECG2, PPG1, PPG2, titolo)
media = [mean(ECG1); mean(ECG2); mean(PPG1); mean(PPG2)];
mediana = [median(ECG1) ; median(ECG2); median(PPG1); median(PPG2)];
moda = [mode(ECG1); mode(ECG2); mode(PPG1); mode(PPG2)];
fprintf('%s\n',titolo);

fprintf('Media: ECG1=%d ECG2=%d PPG1=%d PPG2=%d\n',media(1),media(2),media(3),media(4)); 
fprintf('Moda: ECG1=%d ECG2=%d PPG1=%d PPG2=%d\n',moda(1),moda(2),moda(3),moda(4)); 
fprintf('Mediana: ECG1=%d ECG2=%d PPG1=%d PPG2=%d\n',mediana(1),mediana(2),mediana(3),mediana(4)); 
end

