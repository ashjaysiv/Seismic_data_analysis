function [Param] = EPSDParam(EWa,dt)
%dt = 0.005;
bb = length(EWa);  
fs = 1/dt;                                                         
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~EEMD~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
StandevData = std(EWa);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if StandevData > 50
    NstdSet = [0.1 0.2 0.3 0.5 0.7];
    VarP = 0;
    for  n = 1: length(NstdSet)
        Nstd = NstdSet(n);
        allmode = eemd(EWa,Nstd,3);
        varData=var(allmode(:,1));
        [lend wid]=size(allmode);
        for m=2:(wid-1)
            PerVar(m-1,:)=(var(allmode(:,m))./varData)*100;
        end
        VarTotal=sum(PerVar);
        if VarTotal <= 115 && VarTotal >= 85
            VarP = VarTotal;
            break
        end
    end
    if VarP == 0
        Nstd = 0.81;
        allmode = eemd(EWa, Nstd,3);
        varData = var(allmode(:,1));
        [lend, wid]=size(allmode);
        for m=2:(wid-1)
            PerVar(m-1,:)=(var(allmode(:,m))./varData)*100;
        end
        VarP=sum(PerVar);
    end
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if StandevData >= 10 && StandevData <= 50
    NstdSet = [0.01 0.05 0.09 0.1 0.2];
    VarP = 0;
    for  n = 1: length(NstdSet)
        Nstd = NstdSet(n);
        allmode = eemd(EWa,Nstd,3);
        varData=var(allmode(:,1));
        [lend wid]=size(allmode);
        for m=2:(wid-1)
            PerVar(m-1,:)=(var(allmode(:,m))./varData)*100;
        end
        VarTotal=sum(PerVar);
        if VarTotal <= 115 && VarTotal >= 85
            VarP = VarTotal;
            break
        end
    end
    if VarP == 0
        Nstd = 0.81;
        allmode = eemd(EWa, Nstd,3);
        varData = var(allmode(:,1));
        [lend, wid]=size(allmode);
        for m=2:(wid-1)
            PerVar(m-1,:)=(var(allmode(:,m))./varData)*100;
        end
        VarP=sum(PerVar);
    end
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if StandevData >= 1 && StandevData <= 10
    NstdSet = [0.005 0.009 0.01 0.05 0.1];
    VarP = 0;
    for  n = 1: length(NstdSet)
        Nstd = NstdSet(n);
        allmode = eemd(EWa,Nstd,3);
        varData=var(allmode(:,1));
        [lend wid]=size(allmode);
        for m=2:(wid-1)
            PerVar(m-1,:)=(var(allmode(:,m))./varData)*100;
        end
        VarTotal=sum(PerVar);
        if VarTotal <= 115 && VarTotal >= 85
            VarP = VarTotal;
            break
        end
    end
    if VarP == 0
        Nstd = 0.81;
        allmode = eemd(EWa, Nstd,3);
        varData = var(allmode(:,1));
        [lend, wid]=size(allmode);
        for m=2:(wid-1)
            PerVar(m-1,:)=(var(allmode(:,m))./varData)*100;
        end
        VarP=sum(PerVar);
    end
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if StandevData <= 1
    NstdSet = [0.001 0.003 0.005 0.009 0.01];
    VarP = 0;
    for  n = 1: length(NstdSet)
        Nstd = NstdSet(n);
        allmode = eemd(EWa,Nstd,3);
        varData=var(allmode(:,1));
        [lend wid]=size(allmode);
        for m=2:(wid-1)
            PerVar(m-1,:)=(var(allmode(:,m))./varData)*100;
        end
        VarT(n,1) = sum(PerVar);
        VarTotal=sum(PerVar);
        if VarTotal <= 115 && VarTotal >= 85
            VarP = VarTotal;
            break
        end
    end
    if VarP == 0
        Nstd = 0.81;
        allmode = eemd(EWa, Nstd,3);
        varData = var(allmode(:,1));
        [lend, wid]=size(allmode);
        for m=2:(wid-1)
            PerVar(m-1,:)=(var(allmode(:,m))./varData)*100;
        end
        VarP=sum(PerVar);
    end
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for k=1:(wid-2)
    IMF(:,k)=allmode(:,k+1);
end
clear allmode PerVar;
Time=bb*dt;
for h=1:(wid-2)
     AF(:,h)=hilbert(IMF(:,h));
     Cj(:,h)=abs(AF(:,h));
     Thetaj(:,h)=unwrap(angle(AF(:,h)));
     IF(:,h)= diff(Thetaj(:,h))./(2*pi*dt);
%      AF_1(:,h)=hilbert(IMF(20:end,h));
%      Cj_1(:,h)=abs(AF(:,h));
%      Thetaj_1(:,h)=unwrap(angle(AF(:,h)));
%       IF_1(:,h)= diff(Thetaj(:,h))./(dt);
%      Thetaj_1(:,h)=phase(AF(:,h));
%      IF_1(:,h)= diff(Thetaj_1(:,h))./(2*pi*dt);
%      Thetaj_1(:,h)= atand(AF(:,h)./abs(AF(:,h)));
%      IF_1(:,h)= diff(Thetaj_1(:,h))./(2*pi*dt);
end
IF(bb,:)=0; 
for k=1:(wid-2)
    for k1=1:bb
        if IF(k1,k)<0
            IF(k1,k)=0;
            %Cj(k1,k)=0;
        end
    end
end
% IF_1(bb,:)=0;
% for k=1:(wid-2)
%     for k1=1:bb
%         if IF_1(k1,k)<0
%             IF_1(k1,k)=0;
%             %Cj(k1,k)=0;
%         end
%     end
% end
% clear AF Thetaj IMF;
%=======================================================================
[npt,nimf]=size(Cj);
fw=fs/2;tw=Time;
tres=bb*0.1;
fres=fw*10;
nt=zeros(ceil(tres),fres);
dt=tw/tres;
dw=fw/fres;
%~~~~~~~~~~~~~~~~~~~~
%P=round(IF./dw); %%%just checking
P=round(IF*10)/10;
Tt=ceil((1:bb)*tres/bb);
Ws=(0:dw:fw)';
Ws1=round(Ws*10)/10;
for x=1:bb                                            %checking every FM values-loop A start  
    for imf=1:nimf                                     %checking the FM values for every IMF-loop B start 
        freqid=P(x,imf);	                            % use P as the FM position index ,called freqidx
        freqidx=(find(Ws1==freqid));
        if (freqidx >= 1 && freqidx <= fres)              %checking the position is 'inside' or 'outside' the grid 
            tx = Tt(x);                                   %tx is the final position of frequency on time axis 
            nt(tx,freqidx)=nt(tx,freqidx)+((Cj(x,imf)^2).*0.5);   %put energy(AM*AM) in its position
        end
    end                                                %checking the FM values for every IMF-loop B end
end   
clear Ws Ws1;
Ws=linspace(0,fw,fres)';
Ts=linspace(0,tw,ceil(tres))';
GTW=flipud(rot90((nt)));
%========================================================================
% % % Pa = sum(GTW)*dw; % Instantaneous average power
% % % Fc1 = sum((Ws')*GTW)*dw/(length(Ws)*Ws(end)); 
% % % Fc = Fc1./Pa; % Central frequency
% % % Fb1 = sum((Ws'.^2)*GTW)*dw/(length(Ws)*Ws(end)*Ws(end));
% % % Fb2 = Fb1./Pa;
% % % Fb = sqrt(abs(Fb2-Fc.^2)); % Frequency Bandwidth
% % % Moments = [Pa; Fc; Fb]; 
% % %1.Total Energy--------------------------------------------------------
% % %Eacc=trapz(Ws,trapz(Ts,GTW));
% % Eacc2=sum(sum(GTW).*dw)*dt;
% % EccK=Eacc2/(981*981);
% % %AI2D(i,1)=((Eacc2*pi./(2*981)));
% % %2.Spectral centroid-------------------------------------------------------
% % sp1=sum(sum((Ws')*GTW).*dw)*dt;
% % Ew=sp1/Eacc2;
% % EwK=Ew;
% % %3.Spectral standard deviation---------------------------------------------
% % for r=1:fres
% %     Wss(r,1)=(Ws(r,1)-Ew)^2;
% % end
% % sp2=sum(sum((Wss')*GTW).*dw)*dt;
% % Sw=sp2/Eacc2;
% % SwK=sqrt(Sw);
% % %4.Temporal centroid-------------------------------------------------------
% % tem1=sum(sum(GTW*Ts).*dw)*dt;
% % Et=tem1/Eacc2;
% % EtK=Et;
% % %5.Temporal standard deviation---------------------------------------------
% % for p=1:ceil(tres)
% %     Tss(p,1)=(Ts(p,1)-Et)^2;
% % end
% % tem2=sum(sum(GTW*Tss).*dw)*dt;
% % St=tem2/Eacc2;
% % StK=sqrt(St);
% % %6.Correlation of time and frequency
% % for w=1:ceil(tres)
% %     Tr(w,1)=(Ts(w,1)-Et);
% % end
% % for y=1:fres
% %     Wr(y,1)=(Ws(y,1)-Ew);
% % end
% % Corr1=sum(sum(((Wr)'*GTW)*Tr))*dw*dt;
% % Corr=Corr1/(Eacc2*sqrt(St)*sqrt(Sw));
% % CorrK=Corr;
%Param = [EccK EwK SwK EtK StK CorrK VarP StandevData Nstd];
Param.GTW = GTW;
Param.Ws = Ws;
Param.Ts = Ts;
end
% figure()
% surf(Ts,Ws,GTW)
% shading interp
% colormap(jet)
