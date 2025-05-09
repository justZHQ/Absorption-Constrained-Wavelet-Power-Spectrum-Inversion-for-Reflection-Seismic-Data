clear all;
clc;
close all;
load R.mat;
load Q;

[n,m]=size(R);
dt=0.001;
t=0:dt:(n-1)*dt;
 

gcf1_Q=figure;
set(gcf1_Q,'position',[20 50 900 850]);
imagesc(1:m,t,Q);
colormap('jet');
colorbar;
caxis([-0 330]);
xlabel({'Trace Number'});
ylabel({'Time/s'});
set(gca,'FontName','Arial','FontSize',18,'linewidth',2);
set(gca,'TickLength',[0 0.001]);
set(gca,'position',[0.095 0.5 0.28*1.5 0.22*1.5]);

Q_eff=zeros(size(Q));
for i=1:m
    Q_eff(:,i)=Qeff(Q(:,i),t);
end

%%
fdom=40;%dominant frequency of minimum phase source waveform
fmax=80;
fa=0;
ricker_l=60;wavelet_type_method='Ricker_COM';
[tw,w,~,~]=Ricker_my(dt,ricker_l,fdom,fa,fmax);

load S.mat 
I=260;
gcf2_RS=figure;
set(gcf2_RS,'position',[20 50 850 750]);
subplot(1,2,1);
wiggle(1:m,t,R,'wiggle',0.05);
xlim([0 m]);
ylim([0 t(end)]);
xlabel({'Trace Number';'(a)'});
ylabel({'Time/s'});
set(gca,'FontName','Arial','FontSize',18,'linewidth',2);
set(gca,'TickLength',[0 0.001]);
% set(gcf,'unit','normalized','position',[0.1,0.1,0.25,0.5]);
set(gca,'position',[0.12 0.14 0.35 0.82]);

subplot(1,2,2);
imagesc(1:m,t,S);hold on;
colormap(seismic(1));
caxis([-0.02 0.02]);
colorbar('Ticks',[-0.02,0,0.02],'position',[0.885 0.14 0.02 0.82]);
stem(I,t(end),'-.k','linewidth',2,'Marker','none');
xlim([0 m]);
% ylim([0 t(end)]);
set(gca,'yticklabel',[]);
xlabel({'Trace Number';'(b)'});
% ylabel({'Time/s'});
set(gca,'FontName','Arial','FontSize',18,'linewidth',2);
set(gca,'TickLength',[0 0.001]);
set(gca,'position',[0.525 0.14 0.35 0.82]);

%%
GausiannVariance=10000000000;
step=200;L=step;
Winsize=200;


qmat=qmatrix(Q_eff(:,I),t,w,tw);%build the qmatrix
[qwavelet,qwaveletwin,signal_afterwin,N,Gausiann,nstart,nend]=Window_cutGausiann(GausiannVariance,Winsize,step,S(:,I),qmat);
gcf4_qmat=figure;
imagesc(t,t,qmat);
ylabel({'Time/s'});
xlabel({'Time/s'});
colormap(jet);
cb=colorbar;
hold on;
for i=1:N
   stem(t(L/2+(i-1)*L),t(end),'-.k','linewidth',2,'Marker','none');
end
set(gca,'FontName','Arial','FontSize',18,'linewidth',2);
set(gca,'TickLength',[0 0.001]);


 Swindows=signal_afterwin;
 Wwindows=qwaveletwin;
scal=1;z=1:N;

gcf5_SI=figure;
plot(t,S(:,I),'k','linewidth',1.2);
ylabel('Amplitude');
xlabel({'Time/s'});
set(gca,'FontName','Arial','FontSize',18,'linewidth',2);
set(gca,'TickLength',[0 0.001]);
xlim([0 2.03]);
gcf6_wavelet=figure;
wigb_wavelet(nstart,nend,qwavelet,scal,z,t);
ylabel('Window Number');
xlabel('Time/s');
set(gca,'FontName','Arial','FontSize',18,'linewidth',2);
set(gca,'TickLength',[0 0.001]);
xlim([0 2.03]);

gcf7_Swindows=figure;
wigb_wavelet(nstart,nend,-Swindows,scal,z,t);
ylabel('Window Number');
xlabel('Time/s');
set(gca,'FontName','Arial','FontSize',18,'linewidth',2);
set(gca,'TickLength',[0 0.001]);
xlim([0 2.03]);
Time=[];
for i=1:N-1
    Time=[Time;L*dt];
end

sigma = 0.01;
numBasis=100;
molecular=[10,20,30,40,50,60,70,80,90];
Signal=S(:,I);
[molecular_vecter,molecular_Signal]=POU_my(numBasis,molecular,sigma,dt,Signal);

gcf8_windows=figure;
wigb_my(-molecular_vecter',scal,z,t);
ylabel('Window Number');
xlabel('Time/s');
set(gca,'FontName','Arial','FontSize',18,'linewidth',2);
set(gca,'TickLength',[0 0.001]);
xlim([0 2.03]);

%%
fsample=1/dt;   % number of samples
df=fsample/L;
f=0:fsample/L:fsample/2;
Begf_high=15;
Endf_low=35;

number=1024;
Begfmin=0;
Endfmax=120;
% [~,~,~,~,W0_log_ampx,W0_ampx]=Effective_Amplitude_Spectrum_LengthFix(dt,W0,Begfmin,Endfmax,number);
[~,~,~,~,W_log_ampx,W_ampx]=Effective_Amplitude_Spectrum_LengthFix(dt,Wwindows,Begfmin,Endfmax,number);
[BegNum,Bandf,df,ff,log_ampx,ampx]=Effective_Amplitude_Spectrum_LengthFix(dt,Swindows,Begfmin,Endfmax,number);

gettrace=[2,3,5,6,8,9];
%%
k=0;
epsilon=0.01;
[Begf,Endf]=Set_wi(df,epsilon,W_ampx);
% Begf=zeros(N,1);
Begf=1*ones(N,1);Begf=Begf*df;Endf=Endf*df;
% Iterations=10;

% %% Q场参数
% p=0.26*ones(N,1);
lamda1=3000;
lamda2=3000;



%% 振幅谱参数

[log_ampx_Correction,ampx_Correction_SSAWPSI]=Fast_S(lamda1,lamda2,ampx,ff,k,Begf,Endf,Begf_high,Endf_low);
gcf_Wampx_SS=figure;
set(gcf_Wampx_SS,'position',[100 100 500 700]);
wigb_my((-W_ampx./max(W_ampx)),1,z,ff,1,gettrace);
title('Wavelet');
% ylabel('Window Number');
set(gca,'yticklabel',[]);
xlabel('Frequency/Hz');
set(gca,'FontName','Arial','FontSize',18,'linewidth',2);
set(gca,'TickLength',[0 0.001]);
xlim([-10 110]);
ylim([-0.5 11]);
for i=1:N-1
    BegNuma=floor(Begf(i)/df)+1; 
    BegNumb=floor(Begf(i+1)/df)+1; 
    EndNuma=floor(Endf(i)/df)+1;
    EndNumb=floor(Endf(i+1)/df)+1;
    
    Begf_highNuma=floor(Begf_high/df)+1;
    Begf_highNumb=floor(Begf_high/df)+1;
    Endf_lowNuma=floor(Endf_low/df)+1;
    Endf_lowNumb=floor(Endf_low/df)+1;
    
    hold on;
% 获取a和b点的时间值
time_a = z(i);
time_b = z(i+1);
% 获取a和b点的频率值
Begfreq_a = ff(BegNuma);
Begfreq_b = ff(BegNumb);
Endfreq_a = ff(EndNuma);
Endfreq_b = ff(EndNumb);
% 绘制红色线连接a和b点
% plot([Begfreq_a, Begfreq_b],[time_a, time_b], 'r', 'LineWidth', 2);
% plot([Endfreq_a, Endfreq_b],[time_a, time_b], 'r', 'LineWidth', 2);

[maxa,peak_a]=max(W_ampx(:,i));
[maxab,peak_b]=max(W_ampx(:,i+1));
peakfreq_a = ff(peak_a);
peakfreq_b = ff(peak_b);
plot([peakfreq_a, peakfreq_b],[time_a-W_ampx(peak_a,time_a)./max(W_ampx(:,time_a)), time_b-W_ampx(peak_b,time_b)./max(W_ampx(:,time_b))], '--oy', 'LineWidth', 2);
hold off;
end


gcf_ampxSSAWPSI=figure;
set(gcf_ampxSSAWPSI,'position',[100 100 500 700]);
title('AWSPI-SS');
wigb_my((-ampx_Correction_SSAWPSI./max(ampx_Correction_SSAWPSI)),1,z,ff,1,gettrace);
% ylabel('Window Number');
set(gca,'yticklabel',[]);
xlabel('Frequency/Hz');
set(gca,'FontName','Arial','FontSize',18,'linewidth',2);
set(gca,'TickLength',[0 0.001]);
xlim([-10 110]);
ylim([-0.5 11]);
for i=1:N-1
    BegNuma=floor(Begf(i)/df)+1; 
    BegNumb=floor(Begf(i+1)/df)+1; 
    EndNuma=floor(Endf(i)/df)+1;
    EndNumb=floor(Endf(i+1)/df)+1;
    
    Begf_highNuma=floor(Begf_high/df)+1;
    Begf_highNumb=floor(Begf_high/df)+1;
    Endf_lowNuma=floor(Endf_low/df)+1;
    Endf_lowNumb=floor(Endf_low/df)+1;
    
    hold on;
% 获取a和b点的时间值
time_a = z(i);
time_b = z(i+1);
% 获取a和b点的频率值
Begfreq_a = ff(BegNuma);
Begfreq_b = ff(BegNumb);
Endfreq_a = ff(EndNuma);
Endfreq_b = ff(EndNumb);

Begf_highfreq_a = ff(Begf_highNuma);
Begf_highfreq_b = ff(Begf_highNumb);
Endf_lowfreq_a = ff(Endf_lowNuma);
Endf_lowfreq_b = ff(Endf_lowNumb);
% 绘制红色线连接a和b点
plot([Endf_lowfreq_a, Endf_lowfreq_b],[time_a, time_b], 'g', 'LineWidth', 2);
plot([Begf_highfreq_a, Begf_highfreq_b],[time_a, time_b], 'g', 'LineWidth', 2);

[maxa,peak_a]=max(ampx_Correction_SSAWPSI(:,i));
[maxab,peak_b]=max(ampx_Correction_SSAWPSI(:,i+1));
peakfreq_a = ff(peak_a);
peakfreq_b = ff(peak_b);
plot([peakfreq_a, peakfreq_b],[time_a-W_ampx(peak_a,time_a)./max(W_ampx(:,time_a)), time_b-W_ampx(peak_b,time_b)./max(W_ampx(:,time_b))], '--oy', 'LineWidth', 2);

if i==N-1
text(Begf_highfreq_b+1,time_b +.3,'$$f^1$$','Interpreter', 'latex','FontSize',15,'Color','green');
text(Endf_lowfreq_b+1,time_b +.3,'$$f^2$$','Interpreter', 'latex','FontSize',15,'Color','green');
end

hold off;
end

[log_ampx_Correction0,ampx_Correction_SS]=Fast_S(0,0,ampx,ff,k,Begf,Endf,Begf_high,Endf_low);
gcf_SSampx=figure;
set(gcf_SSampx,'position',[100 100 500 700]);
title('SS');
wigb_my(-(ampx_Correction_SS./max(ampx_Correction_SS)),1,z,ff,1,gettrace);
% ylabel('Window Number');
set(gca,'yticklabel',[]);
xlabel('Frequency/Hz');
set(gca,'FontName','Arial','FontSize',18,'linewidth',2);
set(gca,'TickLength',[0 0.001]);
xlim([-10 110]);
ylim([-0.5 11]);
for i=1:N-1
    BegNuma=floor(Begf(i)/df)+1; 
    BegNumb=floor(Begf(i+1)/df)+1; 
    EndNuma=floor(Endf(i)/df)+1;
    EndNumb=floor(Endf(i+1)/df)+1;
    
    Begf_highNuma=floor(Begf_high/df)+1;
    Begf_highNumb=floor(Begf_high/df)+1;
    Endf_lowNuma=floor(Endf_low/df)+1;
    Endf_lowNumb=floor(Endf_low/df)+1;
    
    hold on;
% 获取a和b点的时间值
time_a = z(i);
time_b = z(i+1);
% 获取a和b点的频率值
Begfreq_a = ff(BegNuma);
Begfreq_b = ff(BegNumb);
Endfreq_a = ff(EndNuma);
Endfreq_b = ff(EndNumb);
% 绘制红色线连接a和b点
% plot([Begfreq_a, Begfreq_b],[time_a, time_b], 'r', 'LineWidth', 2);
% plot([Endfreq_a, Endfreq_b],[time_a, time_b], 'r', 'LineWidth', 2);


[maxa,peak_a]=max(ampx_Correction_SS(:,i));
[maxab,peak_b]=max(ampx_Correction_SS(:,i+1));
peakfreq_a = ff(peak_a);
peakfreq_b = ff(peak_b);
plot([peakfreq_a, peakfreq_b],[time_a-W_ampx(peak_a,time_a)./max(W_ampx(:,time_a)), time_b-W_ampx(peak_b,time_b)./max(W_ampx(:,time_b))], '--oy', 'LineWidth', 2);
hold off;
end

gcf_ampxIN_SS=figure;
set(gcf_ampxIN_SS,'position',[100 100 500 700]);
wigb_my(-ampx./max(ampx),1,z,ff,1,gettrace);
title('Seismic trace');
ylabel('Window Number');
xlabel('Frequence/Hz');
set(gca,'FontName','Arial','FontSize',18,'linewidth',2);
set(gca,'TickLength',[0 0.001]);
xlim([-10 110]);
ylim([-0.5 11]);
for i=1:N-1
    BegNuma=floor(Begf(i)/df)+1; 
    BegNumb=floor(Begf(i+1)/df)+1; 
    EndNuma=floor(Endf(i)/df)+1;
    EndNumb=floor(Endf(i+1)/df)+1;
    
    Begf_highNuma=floor(Begf_high/df)+1;
    Begf_highNumb=floor(Begf_high/df)+1;
    Endf_lowNuma=floor(Endf_low/df)+1;
    Endf_lowNumb=floor(Endf_low/df)+1;
    
    hold on;
% 获取a和b点的时间值
time_a = z(i);
time_b = z(i+1);
% 获取a和b点的频率值
Begfreq_a = ff(BegNuma);
Begfreq_b = ff(BegNumb);
Endfreq_a = ff(EndNuma);
Endfreq_b = ff(EndNumb);

Begf_highfreq_a = ff(Begf_highNuma);
Begf_highfreq_b = ff(Begf_highNumb);
Endf_lowfreq_a = ff(Endf_lowNuma);
Endf_lowfreq_b = ff(Endf_lowNumb);
% 绘制红色线连接a和b点
% plot([Begfreq_a, Begfreq_b],[time_a, time_b], 'r', 'LineWidth', 2);
% plot([Endfreq_a, Endfreq_b],[time_a, time_b], 'r', 'LineWidth', 2);

plot([Endf_lowfreq_a, Endf_lowfreq_b],[time_a, time_b], 'g', 'LineWidth', 2);
plot([Begf_highfreq_a, Begf_highfreq_b],[time_a, time_b], 'g', 'LineWidth', 2);
if i==N-1
text(Begf_highfreq_b+1,time_b +.3,'$$f^1$$','Interpreter', 'latex','FontSize',15,'Color','green');
text(Endf_lowfreq_b+1,time_b +.3,'$$f^2$$','Interpreter', 'latex','FontSize',15,'Color','green');
end

hold off;
end

%%
gcf_amp_SS=figure;
set(gcf_amp_SS,'position',[20 50 900 850]);
subplot(2,3,1);i=gettrace(1);
r=diag(corr(W_ampx(:,i)./max(W_ampx(:,i)),ampx_Correction_SS(:,i)./max(ampx_Correction_SS(:,i))));
rawpsi=diag(corr(W_ampx(:,i)./max(W_ampx(:,i)),ampx_Correction_SSAWPSI(:,i)./max(ampx_Correction_SSAWPSI(:,i))));
hold on;box on;
% plot(W0_ampx(:,i)./max(W0_ampx(:,i)),'y','linewidth',1.2);
plot(ff,W_ampx(:,i)./max(W_ampx(:,i)),'r','linewidth',1.2);
plot(ff,ampx(:,i)./max(ampx(:,i)),'--k','linewidth',1.2);
plot(ff,ampx_Correction_SS(:,i)./max(ampx_Correction_SS(:,i)),'-.g','linewidth',1.2);
plot(ff,ampx_Correction_SSAWPSI(:,i)./max(ampx_Correction_SSAWPSI(:,i)),'b','linewidth',1.2);
ylabel('Magnitude');
set(gca,'xticklabel',[]);
xlabel({'(a)'});
% legend('Wavelet','Seismic','COM','AWSPI-COM','location','NorthEast','FontSize',15,'linewidth',2);
set(gca,'FontName','Arial','FontSize',18,'linewidth',2);
set(gca,'TickLength',[0 0.001]);
str1 = append('SS:', string(round(r,3)));
str2 = append('AWPSI-SS:', string(round(rawpsi,3)));
str=[str1 str2];
annotation('textbox',[0.215 0.865 0.16 0.06],'String',str,'FitBoxToText','off','FontName','Arial','FontSize',12,'linewidth',1.5,'FontWeight','bold');
set(gca,'position',[0.095 0.555 0.29 0.38]);

subplot(2,3,2);i=gettrace(2);
r=diag(corr(W_ampx(:,i)./max(W_ampx(:,i)),ampx_Correction_SS(:,i)./max(ampx_Correction_SS(:,i))));
rawpsi=diag(corr(W_ampx(:,i)./max(W_ampx(:,i)),ampx_Correction_SSAWPSI(:,i)./max(ampx_Correction_SSAWPSI(:,i))));
hold on;box on;
% plot(W0_ampx(:,i)./max(W0_ampx(:,i)),'y','linewidth',1.2);
plot(ff,W_ampx(:,i)./max(W_ampx(:,i)),'r','linewidth',1.2);
plot(ff,ampx(:,i)./max(ampx(:,i)),'--k','linewidth',1.2);
plot(ff,ampx_Correction_SS(:,i)./max(ampx_Correction_SS(:,i)),'-.g','linewidth',1.2);
plot(ff,ampx_Correction_SSAWPSI(:,i)./max(ampx_Correction_SSAWPSI(:,i)),'b','linewidth',1.2);
set(gca,'yticklabel',[]);
set(gca,'xticklabel',[]);
xlabel({'(b)'});
% legend('Wavelet','AWSPI-COM','COM','Seismic','location','NorthEast','FontSize',15,'linewidth',2);
set(gca,'FontName','Arial','FontSize',18,'linewidth',2);
set(gca,'TickLength',[0 0.001]);
str1 = append('SS:', string(round(r,3)));
str2 = append('AWPSI-SS:', string(round(rawpsi,3)));
str=[str1 str2];
annotation('textbox',[0.515 0.865 0.16 0.06],'String',str,'FitBoxToText','off','FontName','Arial','FontSize',12,'linewidth',1.5,'FontWeight','bold');
set(gca,'position',[0.395 0.555 0.29 0.38]);

subplot(2,3,3);i=gettrace(3);
r=diag(corr(W_ampx(:,i)./max(W_ampx(:,i)),ampx_Correction_SS(:,i)./max(ampx_Correction_SS(:,i))));
rawpsi=diag(corr(W_ampx(:,i)./max(W_ampx(:,i)),ampx_Correction_SSAWPSI(:,i)./max(ampx_Correction_SSAWPSI(:,i))));
hold on;box on;
plot(ff,W_ampx(:,i)./max(W_ampx(:,i)),'r','linewidth',1.2);
plot(ff,ampx(:,i)./max(ampx(:,i)),'--k','linewidth',1.2);
plot(ff,ampx_Correction_SS(:,i)./max(ampx_Correction_SS(:,i)),'-.g','linewidth',1.2);
plot(ff,ampx_Correction_SSAWPSI(:,i)./max(ampx_Correction_SSAWPSI(:,i)),'b','linewidth',1.2);
set(gca,'yticklabel',[]);
set(gca,'xticklabel',[]);
xlabel({'(c)'});
% legend('Wavelet','AWSPI-COM','COM','Seismic','location','NorthEast','FontSize',15,'linewidth',2);
set(gca,'FontName','Arial','FontSize',18,'linewidth',2);
set(gca,'TickLength',[0 0.001]);
str1 = append('SS:', string(round(r,3)));
str2 = append('AWPSI-SS:', string(round(rawpsi,3)));
str=[str1 str2];
annotation('textbox',[0.815 0.865 0.16 0.06],'String',str,'FitBoxToText','off','FontName','Arial','FontSize',12,'linewidth',1.5,'FontWeight','bold');
set(gca,'position',[0.695 0.555 0.29 0.38]);

subplot(2,3,4);i=gettrace(4);
r=diag(corr(W_ampx(:,i)./max(W_ampx(:,i)),ampx_Correction_SS(:,i)./max(ampx_Correction_SS(:,i))));
rawpsi=diag(corr(W_ampx(:,i)./max(W_ampx(:,i)),ampx_Correction_SSAWPSI(:,i)./max(ampx_Correction_SSAWPSI(:,i))));
hold on;box on;
% plot(W0_ampx(:,i)./max(W0_ampx(:,i)),'y','linewidth',1.2);
plot(ff,W_ampx(:,i)./max(W_ampx(:,i)),'r','linewidth',1.2);
plot(ff,ampx(:,i)./max(ampx(:,i)),'--k','linewidth',1.2);
plot(ff,ampx_Correction_SS(:,i)./max(ampx_Correction_SS(:,i)),'-.g','linewidth',1.2);
plot(ff,ampx_Correction_SSAWPSI(:,i)./max(ampx_Correction_SSAWPSI(:,i)),'b','linewidth',1.2);
ylabel('Magnitude');
xlabel({'(d)';'Frequency/Hz'});
% legend('Wavelet','AWSPI-COM','COM','Seismic','location','NorthEast','FontSize',15,'linewidth',2);
set(gca,'FontName','Arial','FontSize',18,'linewidth',2);
set(gca,'TickLength',[0 0.001]);
str1 = append('SS:', string(round(r,3)));
str2 = append('AWPSI-SS:', string(round(rawpsi,3)));
str=[str1 str2];
annotation('textbox',[0.215 0.427 0.16 0.06],'String',str,'FitBoxToText','off','FontName','Arial','FontSize',12,'linewidth',1.5,'FontWeight','bold');
set(gca,'position',[0.095 0.12 0.29 0.38]);

subplot(2,3,5);i=gettrace(5);
r=diag(corr(W_ampx(:,i)./max(W_ampx(:,i)),ampx_Correction_SS(:,i)./max(ampx_Correction_SS(:,i))));
rawpsi=diag(corr(W_ampx(:,i)./max(W_ampx(:,i)),ampx_Correction_SSAWPSI(:,i)./max(ampx_Correction_SSAWPSI(:,i))));
hold on;box on;
% plot(W0_ampx(:,i)./max(W0_ampx(:,i)),'y','linewidth',1.2);
plot(ff,W_ampx(:,i)./max(W_ampx(:,i)),'r','linewidth',1.2);
plot(ff,ampx(:,i)./max(ampx(:,i)),'--k','linewidth',1.2);
plot(ff,ampx_Correction_SS(:,i)./max(ampx_Correction_SS(:,i)),'-.g','linewidth',1.2);
plot(ff,ampx_Correction_SSAWPSI(:,i)./max(ampx_Correction_SSAWPSI(:,i)),'b','linewidth',1.2);
set(gca,'yticklabel',[]);
xlabel({'(e)';'Frequency/Hz'});
% legend('Wavelet','AWSPI-COM','COM','Seismic','location','NorthEast','FontSize',15,'linewidth',2);
set(gca,'FontName','Arial','FontSize',18,'linewidth',2);
set(gca,'TickLength',[0 0.001]);
str1 = append('SS:', string(round(r,3)));
str2 = append('AWPSI-SS:', string(round(rawpsi,3)));
str=[str1 str2];
annotation('textbox',[0.515 0.427 0.16 0.06],'String',str,'FitBoxToText','off','FontName','Arial','FontSize',12,'linewidth',1.5,'FontWeight','bold');
set(gca,'position',[0.395 0.12 0.29 0.38]);
% 
subplot(2,3,6);i=gettrace(6);
r=diag(corr(W_ampx(:,i)./max(W_ampx(:,i)),ampx_Correction_SS(:,i)./max(ampx_Correction_SS(:,i))));
rawpsi=diag(corr(W_ampx(:,i)./max(W_ampx(:,i)),ampx_Correction_SSAWPSI(:,i)./max(ampx_Correction_SSAWPSI(:,i))));
hold on;box on;
plot(ff,W_ampx(:,i)./max(W_ampx(:,i)),'r','linewidth',1.2);
plot(ff,ampx(:,i)./max(ampx(:,i)),'--k','linewidth',1.2);
plot(ff,ampx_Correction_SS(:,i)./max(ampx_Correction_SS(:,i)),'-.g','linewidth',1.2);
plot(ff,ampx_Correction_SSAWPSI(:,i)./max(ampx_Correction_SSAWPSI(:,i)),'b','linewidth',1.2);
set(gca,'yticklabel',[]);
% set(gca,'xticklabel',[]);
% legend('Wavelet','AWSPI-COM','COM','Seismic','location','NorthEast','FontSize',15,'linewidth',2);
set(gca,'FontName','Arial','FontSize',18,'linewidth',2);
set(gca,'TickLength',[0 0.001]);
xlabel({'(f)';'Frequency/Hz'});
% legend('Ture','COM','AWSPI-COM','LSR','location','NorthEast','FontSize',15,'linewidth',2);
legend('Ture','Seismic','SS','AWPSI-SS','NumColumns',4,'location',[0.18 0.95 0.7 0.04],'FontName','Arial');
str1 = append('SS:', string(round(r,3)));
str2 = append('AWPSI-SS:', string(round(rawpsi,3)));
str=[str1 str2];
annotation('textbox',[0.815 0.427 0.16 0.06],'String',str,'FitBoxToText','off','FontName','Arial','FontSize',12,'linewidth',1.5,'FontWeight','bold');
set(gca,'position',[0.695 0.12 0.29 0.38]);
