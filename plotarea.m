function plotarea(battiti,spettro,f,i,nome)

close(figure(i));
figure (i);

subplot (2,1,1)

plot(battiti);
title(sprintf('%s',nome))
xlabel('battiti');
ylabel('[ms]');
xlim([0,154])



subplot (2,1,2) 
hold on

plot(f,abs(spettro))

in1 = find(f>=0 & f<0.04);
in2 = find(f>=0.04 & f<0.15);
in3 = find(f>=0.15 & f<0.4);


area(f(in1),abs(spettro(in1)),'FaceColor',[0.9290 0.6940 0.1250])
area(f(in2),abs(spettro(in2)),'FaceColor',[0.3010 0.7450 0.9330])
area(f(in3),abs(spettro(in3)),'FaceColor',[0.8500 0.3250 0.0980])

[sp,loc]=findpeaks(abs(spettro(in1)),f(in1));
v_area=trapz(f(in1),abs(spettro(in1)));

text(0.005,sp-0.00001,sprintf('%.1dW',v_area));            
text(loc+0.005,sp+0.00001,sprintf('%.3fHz', loc));  
plot(loc,sp,'rx');

[sp,loc]=findpeaks(abs(spettro(in2)),f(in2));
v_area=trapz(f(in2),abs(spettro(in2)));

text(loc-0.008,sp+0.00005,sprintf('%.1dW',v_area));
text(loc-0.008,sp+0.00003,sprintf('%.3fHz', loc));
plot(loc,sp,'ro');

[sp,loc]=findpeaks(abs(spettro(in3)),f(in3));
if length(sp) ~= 1
   
    val = max(sp);
    loc = loc(sp == val);
    
else
    
    val = sp;
    
end

v_area=trapz(f(in3),abs(spettro(in3)));

text(loc+0.008,val+0.00003,sprintf('%.1dW',v_area));
text(loc+0.008,val+0.00001,sprintf('%.3fHz', loc));
plot(loc,sp,'r+');

xlim([0,0.7])
ylim([0,1.5*10^(-4)])
xlabel('Frequenza [Hz]')
ylabel('Potenza [ms^2]')
legend ('Spettro','Potenza VLF','Potenza LF','Potenza HF','f_{VLF}','f_{LF}','f_{HF}')
hold off

end

