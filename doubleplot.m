function doubleplot(fig,tit,xl,ylone,yltwo,time,valone,valtwo)

hold on;
figure(fig);

subplot(2,1,1);
plot(time,valone);
title(tit);
xlabel(xl);
ylabel(ylone);

subplot(2,1,2);
plot(time,valtwo);
xlabel(xl);
ylabel(yltwo);

hold off;

end