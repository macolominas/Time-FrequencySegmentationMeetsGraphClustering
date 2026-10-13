figure;
subplot(441)
imagesc(abs(Fg))
set(gca,'ydir','normal')
set(gca,'xtick',[])
set(gca,'ytick',[])
title('STFT')
xlabel('time')
ylabel('frequency')
subplot(442)
plot(SNR,median(error_1_a,2),'-o')
hold on;
plot(SNR,median(error_1_b,2),'-d')
plot(SNR,median(error_1_c,2),'-s')
plot(SNR,median(error_1_d,2),'-^')
plot(SNR,median(error_1_e,2),'-<')
title('Median Error for linear chirp')
legend('A','B','C','D','E','location','northeast')
xlabel('SNR in')
subplot(443)
plot(SNR,median(error_2_a,2),'-o')
hold on;
plot(SNR,median(error_2_b,2),'-d')
plot(SNR,median(error_2_c,2),'-s')
plot(SNR,median(error_2_d,2),'-^')
plot(SNR,median(error_2_e,2),'-<')
title('Median Error for Hermite function')
%legend('A','B','C','D','E','F','G')
xlabel('SNR in')


subplot(4,4,6)
boxplot([error_1_a(1,:)' error_1_b(1,:)' error_1_c(1,:)' error_1_e(1,:)'])
set(gca,'xticklabel',{'A','B','C','E'})
ylim([3.5 7])
text(0.65,6.75,'5 dB')

subplot(4,4,10)
boxplot([error_1_a(4,:)' error_1_b(4,:)' error_1_c(4,:)' error_1_e(4,:)'])
set(gca,'xticklabel',{'A','B','C','E'})
ylim([0.5 2])
text(0.65,1.75,'20 dB')


subplot(4,4,14)
boxplot([error_1_b(8,:)' error_1_d(8,:)' error_1_e(8,:)'])
set(gca,'xticklabel',{'B','D','E'})
% ylim([0.05 2])
text(0.65,0.7,'40 dB')


subplot(4,4,7)
boxplot([error_2_a(1,:)' error_2_b(1,:)' error_2_c(1,:)' error_2_e(1,:)'])
set(gca,'xticklabel',{'A','B','C','E'})
ylim([2.8 6])
text(0.65,5.75,'5 dB')

subplot(4,4,11)
boxplot([error_2_a(4,:)' error_2_b(4,:)' error_2_c(4,:)' error_2_e(4,:)'])
set(gca,'xticklabel',{'A','B','C','E'})
ylim([0.5 1.4])
text(0.65,1.15,'20 dB')

subplot(4,4,15)
boxplot([error_2_b(8,:)' error_2_d(8,:)' error_2_e(8,:)'])
set(gca,'xticklabel',{'B','D','E'})
text(0.65,.45,'20 dB')

t = 0:1/N:1-1/N;
subplot(4,4,5)
plot(t,x)
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
xlabel('time')
text(0.05,3.5,'input signal (40 dB)')
ylim([-4 4])

subplot(4,4,9)
plot(t,sigma*2*real(sum(Fg.*S_b1)))
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
xlabel('time')
text(0.05,3.5,'first component (method B)')
ylim([-4 4])


subplot(4,4,13)
plot(t,sigma*2*real(sum(Fg.*S_b2)))
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
xlabel('time')
text(0.05,3.5,'second component (method B)')
ylim([-4 4])





