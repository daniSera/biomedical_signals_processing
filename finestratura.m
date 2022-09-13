function finestratura = finestratura(acs)

%funzione che crea una matrice contenente diversi tipi di finestratura, è
%importante considerare differenti tipi di finestre in quanto la
%distorsione del segnale può essere differente

mat=zeros(length(acs),4);
mat(:,1)=rectwin(length(acs));
mat(:,2)=bartlett(length(acs));
mat(:,3)=hann(length(acs));
mat(:,4)=hamming(length(acs));

finestratura=zeros(length(acs),4);

for j=1:4
    finestratura(:,j)=acs.*(mat(:,j)); 
    %applico le varie finestrature ad acs
end
end

