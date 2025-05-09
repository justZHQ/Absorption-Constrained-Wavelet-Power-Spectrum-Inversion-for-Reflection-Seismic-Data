function [Qeff2D,QCOMeff2D,QCOMAWPSIeff2D,Q_by_Qeff2D,Q_by_QCOMeff2D,Q_by_QCOMAWPSIeff2D]=QfilterCOM_2D(dt,L,GausiannVariance,Winsize,step,Data,Begf_high,Endf_low,number,Begfmin,Endfmax,Begf,Endf,p,Iterations,lamda1,lamda2,Reference,GausiannVariance_AWPSI)
SinvCOMAWPSI2D=[];
SinvCOM2D=[];
Sinv2D=[];
Q_by_Qeff2D=[];
Q_by_QCOMeff2D=[];
Q_by_QCOMAWPSIeff2D=[];
Qeff2D=[];
QCOMeff2D=[];
QCOMAWPSIeff2D=[];
parfor j=1:length(Data(1,:))
    Signal_i=Data(:,j);
    [signal_afterwin,N,~,~,nend]=Window_cutGausiann_filed(GausiannVariance,Winsize,step,Signal_i);
    S=signal_afterwin;
    %%    
%     fsample=1/dt;   % number of samples
%     df=fsample/L;
%     f=0:fsample/L:fsample/2;
%     Begf_high=10;
%     Endf_low=30;   
%     number=1024;
%     Begfmin=0;
%     Endfmax=120;
    
    [BegNum,Bandf,df,ff,log_ampx,ampx]=Effective_Amplitude_Spectrum_LengthFix(dt,S,Begfmin,Endfmax,number);
    %%
%     Begf=zeros(N,1);
%     Endf=[95,85,75,70,55,50];
%     p=0.5*ones(N,1);
    
%     Iterations=10;
%     lamda1=1000;
%     lamda2=10000;
%     [log_ampx_Correction,ampx_Correction_SSAWSPI]=Fast(lamda1,lamda2,ampx,df,p,Iterations,Begf,Endf,Begf_high,Endf_low);
    [log_ampx_Correction,ampx_Correction_COMAWSPI]=Fast_COM(lamda1,lamda2,ampx,df,p,Iterations,Begf,Endf,Begf_high,Endf_low,GausiannVariance_AWPSI);
%     [log_ampx_Correction0,ampx_Correction_SS]=Fast(lamda1,lamda2,ampx,df,p,Iterations,Begf,Endf,Begf_high,Endf_low);
        [log_ampx_Correction,ampx_Correction_COM]=Fast_COM(0,0,ampx,df,p,Iterations,Begf,Endf,Begf_high,Endf_low,GausiannVariance_AWPSI);
        
%     k=1.5;
% [log_ampx_Correction,ampx_Correction_COMAWSPI]=Fast_S(lamda1,lamda2,ampx,ff,k,Begf+2,Endf,Begf_high,Endf_low);
% [log_ampx_Correction,ampx_Correction_COM]=Fast_S(0,0,ampx,ff,k,Begf+2,Endf,Begf_high,Endf_low);
    %%
    Time=[];
    for i=1:N-1
        Time=[Time;L*dt];
    end
    
    BegNum_high=floor(Begf_high/df)+1;
    EndNum_low=floor(Endf_low/df)+1;
    BandfQ=((BegNum_high-1):(EndNum_low-1))*df;
    ampx=ampx(BegNum_high:EndNum_low,:);
    ampx_Correction_COM=ampx_Correction_COM(BegNum_high:EndNum_low,:);
    ampx_Correction_COMAWSPI=ampx_Correction_COMAWSPI(BegNum_high:EndNum_low,:);
    
    %%
    type=3;
    K=1;
    [Qeff,Q_by_Qeff,SpectrumRatio]=GetQeff(type,K,EndNum_low,BandfQ,df,ampx,Time,Reference);
    [Q_COMeff,Q_by_QCOMeff,SpectrumRatio_COM]=GetQeff(type,K,EndNum_low,BandfQ,df,ampx_Correction_COM,Time,Reference);
    [Q_COMAWPSIeff,Q_by_QCOMAWPSIeff,SpectrumRatio_COMAWPSI]=GetQeff(type,K,EndNum_low,BandfQ,df,ampx_Correction_COMAWSPI,Time,Reference);
    
    Q_by_Qeff2D=[Q_by_Qeff2D,Q_by_Qeff];
    Q_by_QCOMeff2D=[Q_by_QCOMeff2D,Q_by_QCOMeff];
    Q_by_QCOMAWPSIeff2D=[Q_by_QCOMAWPSIeff2D,Q_by_QCOMAWPSIeff];
    Qeff2D=[Qeff2D,Qeff'];
    QCOMeff2D=[QCOMeff2D,Q_COMeff'];
    QCOMAWPSIeff2D=[QCOMAWPSIeff2D,Q_COMAWPSIeff'];    
    
    %%
%     Qt=nend(1:end-1);
%     Q_SS=interp1(t(Qt),Q_SSeff,t,'nearest','extrap');
%     Q_SSAWPSI=interp1(t(Qt),Q_SSAWPSIeff,t,'nearest','extrap');
%     Q=interp1(t(Qt),Qeff,t,'nearest','extrap');
    
%     tol=0.2;
%     qmatinv=invq(Q,t,tol);
%     Sinv=qmatinv*Signal_i;
%     qmatinvSS=invq(Q_SS,t,tol);
%     SinvCOM=qmatinvSS*Signal_i;
%     qmatinvSSAWPSI=invq(Q_SSAWPSI,t,tol);
%     SinvSSAWPSI=qmatinvSSAWPSI*Signal_i;
%     
%     SinvSSAWPSI2D=[SinvSSAWPSI2D,SinvSSAWPSI];
%     SinvSS2D=[SinvSS2D,SinvCOM];
%     Sinv2D=[Sinv2D,Sinv];   
    %%
    j
end