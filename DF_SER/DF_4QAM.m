clc;clear;
close all;

tStart=tic;

N=1e+9;
Nb=1;
SNRdb=0:2:40;
dsd=1;dsr=1;drd=1;
M=4;
b_qam=3/(M-1);
K=1-1/sqrt(M);
SER_exact=zeros(1,length(SNRdb));
SER_up=zeros(1,length(SNRdb));
SER_tight=zeros(1,length(SNRdb));
simber_DF=zeros(1,length(SNRdb));
for i=1:length(SNRdb)
    tSNRStart=tic;
    r_p=10^(SNRdb(i)/10);
    AQ=(M-1)/2/M+K.^2/pi;
    BQ=3*(M-1)/8/M+K.^2/pi;
%%%%%%%%功率分配%%%%%
   % r_1=(sqrt(dsr)+sqrt(dsr+8*(AQ.^2/BQ)*drd))/(3*sqrt(dsr)+sqrt(dsr+8*(AQ.^2/BQ)*drd))*r_p;%式（25）、（26）
   % r_2=2*sqrt(dsr)/(3*sqrt(dsr)+sqrt(dsr+8*(AQ.^2/BQ)*drd))*r_p;
     r_1=r_p/2;
     r_2=r_p/2;
%%%%%%%%%%%%%仿真%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 Nerror=0;
    for ip=1:N
        yuanma=1/sqrt(2)*[-1-1i,-1+1i,1-1i,1+1i ];%QPSK/4QAM
        bit_data=randi([0,1],1,2*Nb);
        so=[2,1]*reshape(bit_data,2,Nb)+1;
        s=yuanma(so);
        hsr=sqrt(dsr)*sqrt(1/2)*(randn+1i*randn);
        noise1=sqrt(1/2)*(randn+1i*randn);
        hrd=sqrt(drd)*sqrt(1/2)*(randn+1i*randn);
        noise2=sqrt(1/2)*(randn+1i*randn);        
        hsd=sqrt(dsd)*sqrt(1/2)*(randn+1i*randn);
        noise3=sqrt(1/2)*(randn+1i*randn);         
        rsr=sqrt(r_1)*hsr*s+noise1;
        sr_panjue=rsr*sqrt(r_1)*conj(hsr)/(r_1*abs(hsr).^2);
        min1=abs(bsxfun(@minus,yuanma,sr_panjue.')).^2;
        [min_x1,index1]=min(min1,[],2);
        if index1==so
        rrd=sqrt(r_2)*hrd*yuanma(index1)+noise2;
        rsd=sqrt(r_1)*hsd*s+noise3;        
        r2_panjue=(sqrt(r_2)*conj(hrd)*rrd+sqrt(r_1)*conj(hsd)*rsd)/(r_2*abs(hrd).^2+r_1*abs(hsd).^2);
        min2=abs(bsxfun(@minus,yuanma,r2_panjue.')).^2;
        [min_x2,index2]=min(min2,[],2);
        Nerror=Nerror+sum(index2~=so);
        end
        if index1~=so
        rsd=sqrt(r_1)*hsd*s+noise3;        
        r3_panjue=sqrt(r_1)*conj(hsd)*rsd/(r_1*abs(hsd).^2);
        min3=abs(bsxfun(@minus,yuanma,r3_panjue.')).^2;
        [min_x3,index3]=min(min3,[],2);
        Nerror=Nerror+sum(index3~=so);
        end
    end
    simber_DF(i)=Nerror/Nb/N;
    
%%%%%%%%%%%%%%%%%%%精确,N0=1%%%%%%%%%%%%%%%%%%%%%%%%%%
    fsd=@(x) 1./(1+b_qam*dsd*r_1./sin(x).^2/2);
    Fsd=4*K/pi*integral(fsd,0,1/2*pi)-4*K.^2/pi*integral(fsd,0,1/4*pi);
    fsr=@(x) 1./(1+b_qam*dsr*r_1./sin(x).^2/2);
    Fsr=4*K/pi*integral(fsr,0,1/2*pi)-4*K.^2/pi*integral(fsr,0,1/4*pi);
    fsdrd=@(x) 1./((1+b_qam*dsd*r_1./sin(x).^2/2).*(1+b_qam*drd*r_2./sin(x).^2/2));
    Fsdrd=4*K/pi*integral(fsdrd,0,1/2*pi)-4*K.^2/pi*integral(fsdrd,0,1/4*pi);
    SER_exact(i)=Fsd*Fsr+Fsdrd*(1-Fsr);
    %%%%%%%%上限%%%%%%%%%%%%%
    p1=M*(b_qam/2)*r_1*dsd+(M-1)*(b_qam/2)*r_2*drd+(2*M-1);
    p2=(1+(b_qam/2)*r_1*dsd)*(1+(b_qam/2)*r_1*dsr)*(1+(b_qam/2)*r_2*drd);
    SER_up(i)=((M-1)/M.^2).*p1./p2;
    %%%%%%%%%渐近%%%%%%%%%%%%%%%
    SER_tight(i)=AQ.^2/dsr/dsd/(b_qam/2).^2/r_1.^2+BQ/(b_qam/2).^2/drd/dsd/r_1/r_2; 
    tSNREnd=toc(tSNRStart);
    disp(['i=',num2str(i),', takes ',num2str(tSNREnd),' s']);
end

figure(1);
semilogy(SNRdb,simber_DF,'r*',SNRdb,SER_exact,'gd-',SNRdb,SER_up,'k:',SNRdb,SER_tight,'bo-');
grid on
xlabel('SNR(dB)');
ylabel('SER');
legend('Monte Carlo','SER Exact','Upper bound','Tight Appro')
axis([0 40 1e-7 1e+1]);

tEnd=toc(tStart);
disp(['The total time is ',num2str(tEnd),' seconds']);

save data.mat
print('-f1','-dpdf','savepic1.pdf'); 
