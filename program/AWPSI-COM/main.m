clear all;
clc;
close all;
load R.mat;
load Q;

addpath(genpath(fullfile(pwd, 'inpaint_nans')));

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
S= awgn(S,1000000000000000000000);
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
% step=100;L=step;
% step=140;L=step;
step=200;L=step;
Winsize=200;

% step=140;L=step;
% Winsize=300;

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
numBasis=200;

molecular_start=[1,21,41,61,81,101,121,141,161,181];
molecular_end=[20,40,60,80,100,120,140,160,180,200];

Signal=S(:,I);
[molecular_vecter,molecular_Signal]=POU_my_startandend(numBasis,molecular_start,molecular_end,sigma,dt,Signal);

gcf8_windows=figure;
wigb_my(-molecular_vecter',scal*0.8,z,t);
ylabel('Window Number');
xlabel('Time/s');
set(gca,'FontName','Arial','FontSize',18,'linewidth',2);
set(gca,'TickLength',[0 0.001]);
xlim([0 2.03]);

figure;
wigb_my(-molecular_Signal',scal,z,t);
ylabel('Window Number');
xlabel('Time/s');
set(gca,'FontName','Arial','FontSize',18,'linewidth',2);
set(gca,'TickLength',[0 0.001]);
xlim([0 2.03]);
%%
fsample=1/dt;   % number of samples
df=fsample/L;
f=0:fsample/L:fsample/2;
Begf_high=5;
Endf_low=40;

number=1024;
Begfmin=0;
Endfmax=120;
% [~,~,~,~,W0_log_ampx,W0_ampx]=Effective_Amplitude_Spectrum_LengthFix(dt,W0,Begfmin,Endfmax,number);
[~,~,~,~,W_log_ampx,W_ampx]=Effective_Amplitude_Spectrum_LengthFix(dt,Wwindows,Begfmin,Endfmax,number);
[BegNum,Bandf,df,ff,log_ampx,ampx]=Effective_Amplitude_Spectrum_LengthFix(dt,Swindows,Begfmin,Endfmax,number);

gettracei=[2,8];
gettrace=[2,3,8,9];

%%
% k=1.5;
epsilon=0.01;
[Begf,Endf]=Set_wi(df,epsilon,W_ampx);
% Begf=zeros(N,1);
Begf=0*ones(N,1);Begf=Begf*df;Endf=Endf*df;
Iterations=10;

% %% Q场参数
p=0.26*ones(N,1);
lamda1=5000;
lamda2=5000;
GausiannVariance_AWPSI=4;

%% 振幅谱参数
[log_ampx_Correction,ampx_Correction_COMAWPSI]=Fast_COM(lamda1,lamda2,ampx,df,p,Iterations,Begf,Endf,Begf_high,Endf_low,GausiannVariance_AWPSI);
gcf10_Wampx=figure;
set(gcf10_Wampx,'position',[100 100 500 700]);
% wigb_my((-W_ampx./max(W_ampx)),1,z,ff,1,gettrace);
wigb_my((-W_ampx./max(W_ampx)),1,z,ff);
title('Wavelet');
% ylabel('Window Number');
set(gca,'yticklabel',[]);
xlabel('Frequency/Hz');
set(gca,'FontName','Arial','FontSize',18,'linewidth',2);
set(gca,'TickLength',[0 0.001]);
xlim([-10 110]);
ylim([-0.5 N+1]);
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
plot([Begfreq_a, Begfreq_b],[time_a, time_b], 'r', 'LineWidth', 2);
plot([Endfreq_a, Endfreq_b],[time_a, time_b], 'r', 'LineWidth', 2);

[maxa,peak_a]=max(W_ampx(:,i));
[maxab,peak_b]=max(W_ampx(:,i+1));
peakfreq_a = ff(peak_a);
peakfreq_b = ff(peak_b);
plot([peakfreq_a, peakfreq_b],[time_a-W_ampx(peak_a,time_a)./max(W_ampx(:,time_a)), time_b-W_ampx(peak_b,time_b)./max(W_ampx(:,time_b))], '--oy', 'LineWidth', 2);

txt1=strcat('$$f_{',num2str(i),'}^1$$');
txt2=strcat('$$f_{',num2str(i),'}^2$$');
text(Begfreq_a-9,time_a +.3,txt1,'Interpreter', 'latex','FontSize',15,'Color','red');
text(Endfreq_a+1,time_a +.3,txt2,'Interpreter', 'latex','FontSize',15,'Color','red');
if i==N-1
txt1=strcat('$$f_{',num2str(i+1),'}^1$$');
txt2=strcat('$$f_{',num2str(i+1),'}^2$$');
text(Begfreq_b-9,time_b +.3,txt1,'Interpreter', 'latex','FontSize',15,'Color','red');
text(Endfreq_b+1,time_b +.3,txt2,'Interpreter', 'latex','FontSize',15,'Color','red');
end
hold off;
end


gcf11_ampxCOMAWPSI=figure;
set(gcf11_ampxCOMAWPSI,'position',[100 100 500 700]);
title('AWPSI-COM');
% wigb_my((-ampx_Correction_COMAWPSI./max(ampx_Correction_COMAWPSI)),1,z,ff,1,gettrace);
wigb_my((-ampx_Correction_COMAWPSI./max(ampx_Correction_COMAWPSI)),1,z,ff);
% ylabel('Window Number');
set(gca,'yticklabel',[]);
xlabel('Frequency/Hz');
set(gca,'FontName','Arial','FontSize',18,'linewidth',2);
set(gca,'TickLength',[0 0.001]);
xlim([-10 110]);
ylim([-0.5 N+1]);
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
plot([Begfreq_a, Begfreq_b],[time_a, time_b], 'r', 'LineWidth', 2);
plot([Endfreq_a, Endfreq_b],[time_a, time_b], 'r', 'LineWidth', 2);

plot([Endf_lowfreq_a, Endf_lowfreq_b],[time_a, time_b], 'g', 'LineWidth', 2);
plot([Begf_highfreq_a, Begf_highfreq_b],[time_a, time_b], 'g', 'LineWidth', 2);

[maxa,peak_a]=max(ampx_Correction_COMAWPSI(:,i));
[maxab,peak_b]=max(ampx_Correction_COMAWPSI(:,i+1));
peakfreq_a = ff(peak_a);
peakfreq_b = ff(peak_b);
plot([peakfreq_a, peakfreq_b],[time_a-W_ampx(peak_a,time_a)./max(W_ampx(:,time_a)), time_b-W_ampx(peak_b,time_b)./max(W_ampx(:,time_b))], '--oy', 'LineWidth', 2);

txt1=strcat('$$f_{',num2str(i),'}^1$$');
txt2=strcat('$$f_{',num2str(i),'}^2$$');
text(Begfreq_a-9,time_a +.3,txt1,'Interpreter', 'latex','FontSize',15,'Color','red');
text(Endfreq_a+1,time_a +.3,txt2,'Interpreter', 'latex','FontSize',15,'Color','red');
if i==N-1
txt1=strcat('$$f_{',num2str(i+1),'}^1$$');
txt2=strcat('$$f_{',num2str(i+1),'}^2$$');

text(Begfreq_b-9,time_b +.3,txt1,'Interpreter', 'latex','FontSize',15,'Color','red');
text(Endfreq_b+1,time_b +.3,txt2,'Interpreter', 'latex','FontSize',15,'Color','red');

text(Begf_highfreq_b+1,time_b +.3,'$$f^1$$','Interpreter', 'latex','FontSize',15,'Color','green');
text(Endf_lowfreq_b+1,time_b +.3,'$$f^2$$','Interpreter', 'latex','FontSize',15,'Color','green');
end

hold off;
end


[log_ampx_Correction0,ampx_Correction_COM]=Fast_COM(0,0,ampx,df,p,Iterations,Begf,Endf,Begf_high,Endf_low,GausiannVariance_AWPSI);
gcf12_ampxCOM=figure;
set(gcf12_ampxCOM,'position',[100 100 500 700]);
title('COM');
% wigb_my(-(ampx_Correction_COM./max(ampx_Correction_COM)),1,z,ff,1,gettrace);
wigb_my(-(ampx_Correction_COM./max(ampx_Correction_COM)),1,z,ff);
% ylabel('Window Number');
set(gca,'yticklabel',[]);
xlabel('Frequency/Hz');
set(gca,'FontName','Arial','FontSize',18,'linewidth',2);
set(gca,'TickLength',[0 0.001]);
xlim([-10 110]);
ylim([-0.5 N+1]);
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
plot([Begfreq_a, Begfreq_b],[time_a, time_b], 'r', 'LineWidth', 2);
plot([Endfreq_a, Endfreq_b],[time_a, time_b], 'r', 'LineWidth', 2);

[maxa,peak_a]=max(ampx_Correction_COM(:,i));
[maxab,peak_b]=max(ampx_Correction_COM(:,i+1));
peakfreq_a = ff(peak_a);
peakfreq_b = ff(peak_b);
plot([peakfreq_a, peakfreq_b],[time_a-W_ampx(peak_a,time_a)./max(W_ampx(:,time_a)), time_b-W_ampx(peak_b,time_b)./max(W_ampx(:,time_b))], '--oy', 'LineWidth', 2);
txt1=strcat('$$f_{',num2str(i),'}^1$$');
txt2=strcat('$$f_{',num2str(i),'}^2$$');
text(Begfreq_a-9,time_a +.3,txt1,'Interpreter', 'latex','FontSize',15,'Color','red');
text(Endfreq_a+1,time_a +.3,txt2,'Interpreter', 'latex','FontSize',15,'Color','red');
if i==N-1
txt1=strcat('$$f_{',num2str(i+1),'}^1$$');
txt2=strcat('$$f_{',num2str(i+1),'}^2$$');
text(Begfreq_b-9,time_b +.3,txt1,'Interpreter', 'latex','FontSize',15,'Color','red');
text(Endfreq_b+1,time_b +.3,txt2,'Interpreter', 'latex','FontSize',15,'Color','red');
end
hold off;
end

gcf9_ampx=figure;
set(gcf9_ampx,'position',[100 100 500 700]);
% wigb_my(-ampx./max(ampx),1,z,ff,1,gettrace);
wigb_my(-ampx./max(ampx),1,z,ff);
title('Seismic trace');
ylabel('Window Number');
xlabel('Frequence/Hz');
set(gca,'FontName','Arial','FontSize',18,'linewidth',2);
set(gca,'TickLength',[0 0.001]);
xlim([-10 110]);
ylim([-0.5 N+1]);
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
plot([Begfreq_a, Begfreq_b],[time_a, time_b], 'r', 'LineWidth', 2);
plot([Endfreq_a, Endfreq_b],[time_a, time_b], 'r', 'LineWidth', 2);

plot([Endf_lowfreq_a, Endf_lowfreq_b],[time_a, time_b], 'g', 'LineWidth', 2);
plot([Begf_highfreq_a, Begf_highfreq_b],[time_a, time_b], 'g', 'LineWidth', 2);

txt1=strcat('$$f_{',num2str(i),'}^1$$');
txt2=strcat('$$f_{',num2str(i),'}^2$$');
text(Begfreq_a-9,time_a +.3,txt1,'Interpreter', 'latex','FontSize',15,'Color','red');
text(Endfreq_a+1,time_a +.3,txt2,'Interpreter', 'latex','FontSize',15,'Color','red');
if i==N-1
txt1=strcat('$$f_{',num2str(i+1),'}^1$$');
txt2=strcat('$$f_{',num2str(i+1),'}^2$$');

text(Begfreq_b-9,time_b +.3,txt1,'Interpreter', 'latex','FontSize',15,'Color','red');
text(Endfreq_b+1,time_b +.3,txt2,'Interpreter', 'latex','FontSize',15,'Color','red');

text(Begf_highfreq_b+1,time_b +.3,'$$f^1$$','Interpreter', 'latex','FontSize',15,'Color','green');
text(Endf_lowfreq_b+1,time_b +.3,'$$f^2$$','Interpreter', 'latex','FontSize',15,'Color','green');
end

hold off;
end

%%
gcf13_amplsr=figure;
set(gcf13_amplsr,'position',[20 50 900 850]);
subplot(2,3,1);i=gettrace(1);
r=diag(corr(W_ampx(:,i)./max(W_ampx(:,i)),ampx_Correction_COM(:,i)./max(ampx_Correction_COM(:,i))));
rawpsi=diag(corr(W_ampx(:,i)./max(W_ampx(:,i)),ampx_Correction_COMAWPSI(:,i)./max(ampx_Correction_COMAWPSI(:,i))));
hold on;box on;
% plot(W0_ampx(:,i)./max(W0_ampx(:,i)),'y','linewidth',1.2);
plot(ff,W_ampx(:,i)./max(W_ampx(:,i)),'r','linewidth',1.2);
plot(ff,ampx(:,i)./max(ampx(:,i)),'--k','linewidth',1.2);
plot(ff,ampx_Correction_COM(:,i)./max(ampx_Correction_COM(:,i)),'-.g','linewidth',1.2);
plot(ff,ampx_Correction_COMAWPSI(:,i)./max(ampx_Correction_COMAWPSI(:,i)),'b','linewidth',1.2);
ylabel('Magnitude');
set(gca,'xticklabel',[]);
xlabel({'(a)'});
% legend('Wavelet','Seismic','COM','AWSPI-COM','location','NorthEast','FontSize',15,'linewidth',2);
set(gca,'FontName','Arial','FontSize',18,'linewidth',2);
set(gca,'TickLength',[0 0.001]);
str1 = append('COM:', string(round(r,3)));
str2 = append('AWPSI-COM:', string(round(rawpsi,3)));
str=[str1 str2];
annotation('textbox',[0.205 0.865 0.17 0.06],'String',str,'FitBoxToText','off','FontName','Arial','FontSize',12,'linewidth',1.5,'FontWeight','bold');
set(gca,'position',[0.095 0.555 0.29 0.38]);

subplot(2,3,2);i=gettrace(2);
r=diag(corr(W_ampx(:,i)./max(W_ampx(:,i)),ampx_Correction_COM(:,i)./max(ampx_Correction_COM(:,i))));
rawpsi=diag(corr(W_ampx(:,i)./max(W_ampx(:,i)),ampx_Correction_COMAWPSI(:,i)./max(ampx_Correction_COMAWPSI(:,i))));
hold on;box on;
% plot(W0_ampx(:,i)./max(W0_ampx(:,i)),'y','linewidth',1.2);
plot(ff,W_ampx(:,i)./max(W_ampx(:,i)),'r','linewidth',1.2);
plot(ff,ampx(:,i)./max(ampx(:,i)),'--k','linewidth',1.2);
plot(ff,ampx_Correction_COM(:,i)./max(ampx_Correction_COM(:,i)),'-.g','linewidth',1.2);
plot(ff,ampx_Correction_COMAWPSI(:,i)./max(ampx_Correction_COMAWPSI(:,i)),'b','linewidth',1.2);
set(gca,'yticklabel',[]);
set(gca,'xticklabel',[]);
xlabel({'(b)'});
% legend('Wavelet','AWSPI-COM','COM','Seismic','location','NorthEast','FontSize',15,'linewidth',2);
set(gca,'FontName','Arial','FontSize',18,'linewidth',2);
set(gca,'TickLength',[0 0.001]);
str1 = append('COM:', string(round(r,3)));
str2 = append('AWPSI-COM:', string(round(rawpsi,3)));
str=[str1 str2];
annotation('textbox',[0.505 0.865 0.17 0.06],'String',str,'FitBoxToText','off','FontName','Arial','FontSize',12,'linewidth',1.5,'FontWeight','bold');
set(gca,'position',[0.395 0.555 0.29 0.38]);

subplot(2,3,4);i=gettrace(3);
r=diag(corr(W_ampx(:,i)./max(W_ampx(:,i)),ampx_Correction_COM(:,i)./max(ampx_Correction_COM(:,i))));
rawpsi=diag(corr(W_ampx(:,i)./max(W_ampx(:,i)),ampx_Correction_COMAWPSI(:,i)./max(ampx_Correction_COMAWPSI(:,i))));
hold on;box on;
% plot(W0_ampx(:,i)./max(W0_ampx(:,i)),'y','linewidth',1.2);
plot(ff,W_ampx(:,i)./max(W_ampx(:,i)),'r','linewidth',1.2);
plot(ff,ampx(:,i)./max(ampx(:,i)),'--k','linewidth',1.2);
plot(ff,ampx_Correction_COM(:,i)./max(ampx_Correction_COM(:,i)),'-.g','linewidth',1.2);
plot(ff,ampx_Correction_COMAWPSI(:,i)./max(ampx_Correction_COMAWPSI(:,i)),'b','linewidth',1.2);
ylabel('Magnitude');
xlabel({'(d)';'Frequency/Hz'});
% legend('Wavelet','AWSPI-COM','COM','Seismic','location','NorthEast','FontSize',15,'linewidth',2);
set(gca,'FontName','Arial','FontSize',18,'linewidth',2);
set(gca,'TickLength',[0 0.001]);
str1 = append('COM:', string(round(r,3)));
str2 = append('AWPSI-COM:', string(round(rawpsi,3)));
str=[str1 str2];
annotation('textbox',[0.205 0.427 0.17 0.06],'String',str,'FitBoxToText','off','FontName','Arial','FontSize',12,'linewidth',1.5,'FontWeight','bold');
set(gca,'position',[0.095 0.12 0.29 0.38]);

subplot(2,3,5);i=gettrace(4);
r=diag(corr(W_ampx(:,i)./max(W_ampx(:,i)),ampx_Correction_COM(:,i)./max(ampx_Correction_COM(:,i))));
rawpsi=diag(corr(W_ampx(:,i)./max(W_ampx(:,i)),ampx_Correction_COMAWPSI(:,i)./max(ampx_Correction_COMAWPSI(:,i))));
hold on;box on;
% plot(W0_ampx(:,i)./max(W0_ampx(:,i)),'y','linewidth',1.2);
plot(ff,W_ampx(:,i)./max(W_ampx(:,i)),'r','linewidth',1.2);
plot(ff,ampx(:,i)./max(ampx(:,i)),'--k','linewidth',1.2);
plot(ff,ampx_Correction_COM(:,i)./max(ampx_Correction_COM(:,i)),'-.g','linewidth',1.2);
plot(ff,ampx_Correction_COMAWPSI(:,i)./max(ampx_Correction_COMAWPSI(:,i)),'b','linewidth',1.2);
set(gca,'yticklabel',[]);
xlabel({'(e)';'Frequency/Hz'});
% legend('Wavelet','AWSPI-COM','COM','Seismic','location','NorthEast','FontSize',15,'linewidth',2);
set(gca,'FontName','Arial','FontSize',18,'linewidth',2);
set(gca,'TickLength',[0 0.001]);
str1 = append('COM:', string(round(r,3)));
str2 = append('AWPSI-COM:', string(round(rawpsi,3)));
str=[str1 str2];
annotation('textbox',[0.505 0.427 0.17 0.06],'String',str,'FitBoxToText','off','FontName','Arial','FontSize',12,'linewidth',1.5,'FontWeight','bold');
set(gca,'position',[0.395 0.12 0.29 0.38]);
%%
BegNum_high=floor(Begf_high/df)+1; 
EndNum_low=floor(Endf_low/df)+1;
BandfQ=((BegNum_high-1):(EndNum_low-1))*df;
W_ampx=W_ampx(BegNum_high:EndNum_low,:);
ampx=ampx(BegNum_high:EndNum_low,:);
ampx_Correction_COM=ampx_Correction_COM(BegNum_high:EndNum_low,:);
ampx_Correction_COMAWPSI=ampx_Correction_COMAWPSI(BegNum_high:EndNum_low,:);

type=3;
K=1;
Reference=W_ampx(:,1);
[Q_tureeff,Q_by_Qtureeff,SpectrumRatio_ture]=GetQeff(type,K,EndNum_low,BandfQ,df,W_ampx,Time,Reference);
[Qeff,Q_by_Qeff,SpectrumRatio]=GetQeff(type,K,EndNum_low,BandfQ,df,ampx,Time,Reference);
[Q_COMeff,Q_by_QCOMeff,SpectrumRatio_SCOM]=GetQeff(type,K,EndNum_low,BandfQ,df,ampx_Correction_COM,Time,Reference);
[Q_COMAWPSIeff,Q_by_QCOMAWPSIeff,SpectrumRatio_SAWPSI]=GetQeff(type,K,EndNum_low,BandfQ,df,ampx_Correction_COMAWPSI,Time,Reference);
Qt=nend(1:end-1);

%%
subplot(2,3,3);i=gettracei(1);
FitPar=polyfit(BandfQ,SpectrumRatio(:,i),1);
srm=BandfQ*FitPar(1)+FitPar(2);
hold on;box on;
% plot(W0_ampx(:,i)./max(W0_ampx(:,i)),'y','linewidth',1.2);
plot(BandfQ,SpectrumRatio_ture(:,i),'r','linewidth',1.2);
plot(BandfQ,SpectrumRatio(:,i),'--k','linewidth',1.2);
plot(BandfQ,SpectrumRatio_SCOM(:,i),'-.g','linewidth',1.2);
plot(BandfQ,SpectrumRatio_SAWPSI(:,i),'b','linewidth',1.2);
plot(BandfQ,srm,'m','linewidth',1.2);
% set(gca,'yticklabel',[]);
set(gca,'xticklabel',[]);
xlabel({'(c)'});
% legend('Ture','COM','AWSPI-COM','LSR','location','NorthEast','FontSize',15,'linewidth',2);
set(gca,'FontName','Arial','FontSize',18,'linewidth',2);
set(gca,'TickLength',[0 0.001]);
set(gca,'position',[0.725 0.555 0.26 0.38]);
yticks([-3 -2 -1 0]);

subplot(2,3,6);i=gettracei(2);
FitPar=polyfit(BandfQ,SpectrumRatio(:,i),1);
srm=BandfQ*FitPar(1)+FitPar(2);
hold on;box on;
% plot(W0_ampx(:,i)./max(W0_ampx(:,i)),'y','linewidth',1.2);
plot(BandfQ,SpectrumRatio_ture(:,i),'r','linewidth',1.2);
plot(BandfQ,SpectrumRatio(:,i),'--k','linewidth',1.2);
plot(BandfQ,SpectrumRatio_SCOM(:,i),'-.g','linewidth',1.2);
plot(BandfQ,SpectrumRatio_SAWPSI(:,i),'b','linewidth',1.2);
plot(BandfQ,srm,'m','linewidth',1.2);
% set(gca,'yticklabel',[]);
xlabel({'(f)';'Frequency/Hz'});
% legend('Ture','COM','AWSPI-COM','LSR','location','NorthEast','FontSize',15,'linewidth',2);
legend('Ture','Seismic','COM','AWPSI-COM','NumColumns',4,'location',[0.18 0.95 0.7 0.04],'FontName','Arial');
set(gca,'FontName','Arial','FontSize',18,'linewidth',2);
set(gca,'TickLength',[0 0.001]);
set(gca,'position',[0.725 0.12 0.26 0.38]);
yticks([-8 -4 0]);

%%
[Qeff2D,QCOMeff2D,QCOMAWPSIeff2D,Q_by_Qeff2D,Q_by_QCOMeff2D,Q_by_QCOMAWPSIeff2D]=QfilterCOM_2D(dt,L,GausiannVariance,Winsize,step,S,Begf_high,Endf_low,number,Begfmin,Endfmax,Begf,Endf,p,Iterations,lamda1,lamda2,Reference,GausiannVariance_AWPSI);

%%
window_size=1;
minQ=0;maxQ=150;
[QCOMAWPSIeff2D_1,QCOMAWPSIeff2D_2,QCOMAWPSIeff2D_3]=Qfit_interp(minQ,maxQ,window_size,QCOMAWPSIeff2D);
[QCOMeff2D_1,QCOMeff2D_2,QCOMeff2D_3]=Qfit_interp(minQ,maxQ,window_size,QCOMeff2D);
[QLSReff2D_1,QLSReff2D_2,QLSReff2D_3]=Qfit_interp(minQ,maxQ,window_size,Qeff2D);

%%
QLSReff2D_2=QLSReff2D_2;
QCOMeff2D_2=QCOMeff2D_2;
QCOMAWPSIeff2D_2=QCOMAWPSIeff2D_2;
gcf20_QLSReff2D_2=figure;
imagesc(1:m,t,QLSReff2D_2);hold on;
colormap('jet');
colorbar;
caxis([30 70]);
stem(I,t(end),'-.k','linewidth',2,'Marker','none');
xlim([0 m]);
ylim([0 t(end)]);
xlabel({'Trace Number'});
ylabel({'Time/s'});
set(gca,'FontName','Arial','FontSize',18,'linewidth',2);
set(gca,'TickLength',[0 0.001]);
annotation(gcf,'rectangle',...
    [0.37 0.17 0.13 0.60],...
    'Color',[0.65 0.65 0.65],...
    'LineWidth',2,...
    'LineStyle',':');

gcf21_QCOMeff2D_2=figure;
imagesc(1:m,t,QCOMeff2D_2);hold on;
colormap('jet');
colorbar;
caxis([30 70]);
stem(I,t(end),'-.k','linewidth',2,'Marker','none');
xlim([0 m]);
ylim([0 t(end)]);
xlabel({'Trace Number'});
ylabel({'Time/s'});
set(gca,'FontName','Arial','FontSize',18,'linewidth',2);
set(gca,'TickLength',[0 0.001]);
annotation(gcf,'rectangle',...
    [0.37 0.17 0.13 0.60],...
    'Color',[0.65 0.65 0.65],...
    'LineWidth',2,...
    'LineStyle',':');

gcf22_QCOMAWPSIeff2D_2=figure;
imagesc(1:m,t,QCOMAWPSIeff2D_2);hold on;
colormap('jet');
colorbar;
caxis([30 70]);
stem(I,t(end),'-.k','linewidth',2,'Marker','none');
xlim([0 m]);
ylim([0 t(end)]);
xlabel({'Trace Number'});
ylabel({'Time/s'});
set(gca,'FontName','Arial','FontSize',18,'linewidth',2);
set(gca,'TickLength',[0 0.001]);
annotation(gcf,'rectangle',...
    [0.37 0.17 0.13 0.60],...
    'Color',[0.65 0.65 0.65],...
    'LineWidth',2,...
    'LineStyle',':');


%%
QLSReff2D_3=QLSReff2D_3;
QCOMeff2D_3=QCOMeff2D_3;
QCOMAWPSIeff2D_3=QCOMAWPSIeff2D_3;

gcf17_IQeff=figure;
set(gcf17_IQeff,'position',[20 50 1050 300]);
hold on;box on;
plot(t,Q_eff(:,I),'r','linewidth',1.2);
plot(t,QLSReff2D_3(:,I),'--k','linewidth',1.2);
plot(t,QCOMeff2D_3(:,I),'-.g','linewidth',1.2);
plot(t,QCOMAWPSIeff2D_3(:,I),'b','linewidth',1.2);
ylabel('$$Q_{eff}$$','Interpreter','Latex');
xlabel({'(a)'});
set(gca,'xticklabel',[]);
legend('Ture','Seismic','COM','AWPSI-COM','location','SouthEast');
set(gca,'FontName','Arial','FontSize',18,'linewidth',2);
set(gca,'TickLength',[0 0.001]);
ylim([5 80]);
xlim([-0.1 2.1]);

%%
gcf20_Qget=figure;
set(gcf20_Qget,'position',[20 50 900 850]);
subplot(3,3,1);
imagesc(1:m,t,QLSReff2D_2);hold on;
colormap('jet');
% colorbar;
caxis([0 70]);
stem(I,t(end),'-.k','linewidth',3,'Marker','none');
% xlim([0 m]);
set(gca,'xticklabel',[]);
ylim([0 t(end)]);
xlabel({'(a)'});
ylabel({'Time/s'});
set(gca,'FontName','Arial','FontSize',18,'linewidth',2);
set(gca,'TickLength',[0 0.001]);
set(gca,'position',[0.095 0.7 0.28 0.22]);
% rectangle('Position',[200 0.3 75 1.7],'EdgeColor','w','LineWidth',2);

subplot(3,3,2);
imagesc(1:m,t,QCOMeff2D_2);hold on;
colormap('jet');
% colorbar;
caxis([0 70]);
stem(I,t(end),'-.k','linewidth',3,'Marker','none');
% xlim([0 m]);
set(gca,'xticklabel',[]);
ylim([0 t(end)]);
xlabel({'(b)'});
set(gca,'yticklabel',[]);
% ylabel({'Time/s'});
set(gca,'FontName','Arial','FontSize',18,'linewidth',2);
set(gca,'TickLength',[0 0.001]);
set(gca,'position',[0.385 0.7 0.28 0.22]);

subplot(3,3,3);
imagesc(1:m,t,QCOMAWPSIeff2D_2);hold on;
colormap('jet');
% colorbar;
caxis([0 70]);
stem(I,t(end),'-.k','linewidth',3,'Marker','none');
% xlim([0 m]);
set(gca,'xticklabel',[]);
ylim([0 t(end)]);
xlabel({'(c)'});
set(gca,'yticklabel',[]);
% ylabel({'Time/s'});
set(gca,'FontName','Arial','FontSize',18,'linewidth',2);
set(gca,'TickLength',[0 0.001]);
set(gca,'position',[0.675 0.7 0.28 0.22]);

subplot(3,3,4);
imagesc(1:m,t,QLSReff2D_3);hold on;
colormap('jet');
% colorbar;
caxis([0 70]);
% xlim([0 m]);
set(gca,'xticklabel',[]);
ylim([0 t(end)]);
xlabel({'(d)'});
ylabel({'Time/s'});
set(gca,'FontName','Arial','FontSize',18,'linewidth',2);
set(gca,'TickLength',[0 0.001]);
set(gca,'position',[0.095 0.43 0.28 0.22]);
rectangle('Position',[200 0.3 75 1.7],'EdgeColor','w','LineWidth',3);

subplot(3,3,5);
imagesc(1:m,t,QCOMeff2D_3);hold on;
colormap('jet');
% colorbar;
caxis([0 70]);
% stem(I,t(end),'-.k','linewidth',2,'Marker','none');
% xlim([0 m]);
set(gca,'xticklabel',[]);
% ylim([0 t(end)]);
xlabel({'(e)'});
set(gca,'yticklabel',[]);
% ylabel({'Time/s'});
set(gca,'FontName','Arial','FontSize',18,'linewidth',2);
set(gca,'TickLength',[0 0.001]);
set(gca,'position',[0.385 0.43 0.28 0.22]);
rectangle('Position',[200 0.3 75 1.7],'EdgeColor','w','LineWidth',3);

subplot(3,3,6);
imagesc(1:m,t,QCOMAWPSIeff2D_3);hold on;
colormap('jet');
% colorbar;
caxis([0 70]);
% stem(I,t(end),'-.k','linewidth',2,'Marker','none');
% xlim([0 m]);
set(gca,'xticklabel',[]);
% ylim([0 t(end)]);
xlabel({'(f)'});
set(gca,'yticklabel',[]);
% ylabel({'Time/s'});
set(gca,'FontName','Arial','FontSize',18,'linewidth',2);
set(gca,'TickLength',[0 0.001]);
set(gca,'position',[0.675 0.43 0.28 0.22]);
rectangle('Position',[200 0.3 75 1.7],'EdgeColor','w','LineWidth',3);

subplot(3,3,7);
imagesc(1:m,t,abs(QLSReff2D_3-Q_eff));hold on;
colormap('jet');
% colorbar;
caxis([0 70]);
xlim([0 m]);
% set(gca,'xticklabel',[]);
ylim([0 t(end)]);
xlabel({'(g)';'Trace Number'});
ylabel({'Time/s'});
set(gca,'FontName','Arial','FontSize',18,'linewidth',2);
set(gca,'TickLength',[0 0.001]);
set(gca,'position',[0.095 0.16 0.28 0.22]);
% rectangle('Position',[200 0.3 75 1.7],'EdgeColor','w','LineWidth',2);

subplot(3,3,8);
imagesc(1:m,t,abs(QCOMeff2D_3-Q_eff));hold on;
colormap('jet');
% colorbar;
caxis([0 70]);
% stem(I,t(end),'-.k','linewidth',2,'Marker','none');
xlim([0 m]);
% set(gca,'xticklabel',[]);
% ylim([0 t(end)]);
xlabel({'(h)';'Trace Number'});
set(gca,'yticklabel',[]);
% ylabel({'Time/s'});
set(gca,'FontName','Arial','FontSize',18,'linewidth',2);
set(gca,'TickLength',[0 0.001]);
set(gca,'position',[0.385 0.16 0.28 0.22]);

subplot(3,3,9);
imagesc(1:m,t,abs(QCOMAWPSIeff2D_3-Q_eff));hold on;
colormap('jet');
% colorbar;
caxis([0 70]);
% stem(I,t(end),'-.k','linewidth',2,'Marker','none');
xlim([0 m]);
% set(gca,'xticklabel',[]);
% ylim([0 t(end)]);
xlabel({'(i)';'Trace Number'});
set(gca,'yticklabel',[]);
% ylabel({'Time/s'});
set(gca,'FontName','Arial','FontSize',18,'linewidth',2);
set(gca,'TickLength',[0 0.001]);
set(gca,'position',[0.675 0.16 0.28 0.22]);
% cb = colorbar;
colorbar('southoutside','position',[0.095 0.94 0.86 0.02]);
% annotation({'subplot(3,3,9)'},'rectangle',...
%     [0.37 0.17 0.13 0.60],...
%     'Color',[0.65 0.65 0.65],...
%     'LineWidth',2,...
%     'LineStyle',':');


%%
gcf17_IQeffandSIcompare=figure;
set(gcf17_IQeffandSIcompare,'position',[20 50 1050 1000]);
subplot(4,1,1);hold on;box on;
y=[Q_tureeff',QLSReff2D_2(:,I),QCOMeff2D_2(:,I),QCOMAWPSIeff2D_2(:,I)];
b=bar(1:length(Q_tureeff),y);
b(1).FaceColor = [1 0 0];
b(2).FaceColor = [0 0 0];
b(3).FaceColor = [0 1 0];
b(4).FaceColor = [0 0 1];
h1=plot(Q_tureeff,'-or','linewidth',1.2);
h2=plot(QLSReff2D_2(:,I),'--+k','linewidth',1.2);
h3=plot(QCOMeff2D_2(:,I),'-.squareg','linewidth',1.2);
h4=plot(QCOMAWPSIeff2D_2(:,I),'-*b','linewidth',1.2);
h=[h1;h2;h3;h4];
ylabel('$$Q_{eff}$$','Interpreter','Latex');
% xlabel('Time/s');
xlabel({'(a)'});
set(gca,'xticklabel',[]);
set(gca,'FontName','Arial','FontSize',18,'linewidth',2);
set(gca,'TickLength',[0 0.001]);
set(gca,'position',[0.1 0.76 0.8 0.22]);
% legend(b,'Ture','Seismic','COM','AWPSI-COM','location','SouthEast');
% legend(h,'Ture','Seismic','COM','AWPSI-COM','location','SouthEast');
% set(gca,'position',[0.1 0.5 0.8 0.21]);
ax1 = axes( 'Position',get(gca,'Position'),'Visible','off');
legend([b(1) b(2)],'','','location',[0.4 0.9 0.1 0.05],'FontName','Arial','FontSize',18,'linewidth',2);
legend('boxoff');
ax2 = axes( 'Position',get(gca,'Position'),'Visible','off');
legend(ax2,[b(3) b(4)],'','','location',[0.59 0.9 0.1 0.05],'FontName','Arial','FontSize',18,'linewidth',2);
legend('boxoff');
ax3 = axes( 'Position',get(gca,'Position'),'Visible','off');
legend(ax3,[h1 h2],'Ture','Seismic','location',[0.47 0.9 0.1 0.05],'FontName','Arial','FontSize',18,'linewidth',2);
legend('boxoff');
ax4 = axes( 'Position',get(gca,'Position'),'Visible','off');
legend(ax4,[h3 h4],'COM','AWPSI-COM','location',[0.66 0.9 0.1 0.05],'FontName','Arial','FontSize',18,'linewidth',2);
legend('boxoff');

SI=conv(R(:,I),w,'same');
tol=0.0005;
QI=QLSReff2D_3(:,I);
QICOM=QCOMeff2D_3(:,I);
QICOMAWPSI=QCOMAWPSIeff2D_3(:,I);
% Q=60*ones(n,1);

qmatinv=invq(QI,t,tol);
SinvI=qmatinv*S(:,I);
crI=diag(corr(SI,SinvI));
subplot(4,1,2);
hold on;box on;
h1=plot(t,SI,'r','linewidth',1.2);
h2=plot(t,S(:,I),'b','linewidth',1.2);
h3=plot(t,SinvI,'--k','linewidth',1.2);
h=[h1;h2;h3];
ylabel('Amplitude');
% xlabel({'Time/s'});
set(gca,'xticklabel',[]);
set(gca,'FontName','Arial','FontSize',18,'linewidth',2);
set(gca,'TickLength',[0 0.001]);
str = strcat('Correlation coefficient:', string(crI));
annotation('textbox',[0.32 0.35 0.5 0.2],'String',str,'FitBoxToText','on','FontName','Arial','FontSize',18,'linewidth',2,'LineStyle','none');
xlim([-0.1 2.1]);
ylim([-0.5 0.5]);
yticks([-0.2, 0, 0.2]);
% box off;
% ax = gca;
% ax.YAxisLocation = 'left';
% ax.XAxis.Visible = 'off';
set(gca,'position',[0.1 0.5 0.8 0.21]);
ax2 = axes( 'Position',get(gca,'Position'),'Visible','off');
legend('boxoff')
legend(h(1:2),{'Ture','Nonstationary'},...
    'NumColumns',2,'location','north','FontName','Arial','FontSize',18,'linewidth',2);legend('boxoff')
legend(ax2,h(3),'Seismic','location',[0.7 0.7 0.1 0.1],'FontName','Arial','FontSize',18,'linewidth',2);legend('boxoff')

qmatinv=invq(QICOMAWPSI,t,tol);
SinvICOMAWPSI=qmatinv*S(:,I);
crICOMAWPSI=diag(corr(SI,SinvICOMAWPSI));
subplot(4,1,4);
hold on;box on;
plot(t,SI,'r','linewidth',1.2);
plot(t,S(:,I),'b','linewidth',1.2);
h=plot(t,SinvICOMAWPSI,'--k','linewidth',1.2);
ylabel('Amplitude');
xlabel({'Time/s';'(b)'});
set(gca,'FontName','Arial','FontSize',18,'linewidth',2);
set(gca,'TickLength',[0 0.001]);
legend(h,'AWPSI-COM','location',[0.6 0.08 0.3 0.1],'FontName','Arial');legend('boxoff')
str = strcat('Correlation coefficient:', string(crICOMAWPSI));
annotation('textbox',[0.32 0.05 0.5 0.1],'String',str,'FitBoxToText','on','FontName','Arial','FontSize',18,'linewidth',2,'LineStyle','none');
xlim([-0.1 2.1]);
ylim([-0.5 0.5]);
yticks([-0.2, 0, 0.2]);
% box on;
% ax = gca;
% ax.XAxisLocation = 'bottom';
% ax.YAxisLocation = 'left';
% ax.XRuler.Axle.Visible = 'off';
% ax.XAxis.Visible = 'off';
set(gca,'position',[0.1 0.1 0.8 0.21]);

qmatinv=invq(QICOM,t,tol);
SinvICOM=qmatinv*S(:,I);
crICOM=diag(corr(SI,SinvICOM));
subplot(4,1,3);
hold on;box on;
plot(t,SI,'r','linewidth',1.2);
plot(t,S(:,I),'b','linewidth',1.2);
h=plot(t,SinvICOM,'--k','linewidth',1.2);
ylabel('Amplitude');
% xlabel({'Time/s'});
set(gca,'xticklabel',[]);
set(gca,'FontName','Arial','FontSize',18,'linewidth',2);
set(gca,'TickLength',[0 0.001]);
legend(h,'COM','location',[0.6 0.3 0.3 0.1],'FontName','Arial');legend('boxoff')
str = strcat('Correlation coefficient:', string(crICOM));
annotation('textbox',[0.32 0.15 0.5 0.2],'String',str,'FitBoxToText','on','FontName','Arial','FontSize',18,'linewidth',2,'LineStyle','none');
xlim([-0.1 2.1]);
ylim([-0.5 0.5]);
yticks([-0.2, 0, 0.2]);
box on;
ax = gca;
ax.YAxisLocation = 'left';
ax.XAxis.Visible = 'off';
set(gca,'position',[0.1 0.3 0.8 0.21]);


