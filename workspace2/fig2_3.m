close all
clear
clc

tic
N_sim = 1000;       %代表信道数量，最后对一千个信道的BER取平均
%N_data = 10000;    %代表每个信道传输的bit总数
SNR = 5:1:40;         %在不同的信噪比下进行仿真
noise_dB = -134;    %这个高斯白噪声功率是实际测量得到的，这里取多少并不重要，因为BER只与信噪比有关，只要保证信号功率和噪声功率的比是个定值就可以了
noise = 10^(noise_dB/10);     %注意噪声和信号功率都要转化为真值

a = sin(pi/8);   %画出8PSK的星座图
b = cos(pi/8);
yuanma_8psk = [b+a*1i, a+b*1i, -b+a*1i, -a+b*1i, b-a*1i, a-b*1i, -b-a*1i, -a-b*1i];   %8PSK原码，格雷码
%yuanma_16qam = 1/sqrt(10)*[1+1i, 1+3i, 3+1i, 3+3i, 1-1i, 1-3i, 3-1i, 3-3i, -1+1i, -1+3i, -3+1i, -3+3i, -1-1i, -1-3i, -3-1i, -3-3i];%格雷码
yuanma_16qam = 1/sqrt(10)*[1+1i, 1+3i, 3+3i, 3+1i, 3-1i, 3-3i, 1-3i, 1-1i, -1-1i, -1-3i, -3-3i, -3-1i, -3+1i, -3+3i, -1+3i, -1+1i];%自然码

%% 3bits/s/Hz；one receive antenna
[one_ats_BER,two_ats_BER,three_ats_BER,four_ats_BER] = deal(zeros(N_sim,length(SNR)));      %给BER变量预分配空间，方便后面赋值
[one_ats_SER,two_ats_SER,three_ats_SER,four_ats_SER] = deal(zeros(N_sim,length(SNR)));

for n_sim = 1:N_sim
    fprintf('n_sim = %d,', n_sim);
    tic
    for j = 1:length(SNR)
        Pt = noise*10^(SNR(j)/10);
        A1 = sqrt(Pt);
        A2 = sqrt(Pt/2);      %使用G2码单根天线的功率
        A3 = sqrt(4*Pt/9);    %使用H3码单根天线的功率
        A4 = sqrt(Pt/3);      %使用H4码单根天线的功率
        %fprintf('n_sim = %d, j_Pt = %d\n', n_sim, j);
        
        %% uncoded/8PSK
        one_Enabled = 1;
        if one_Enabled ==1
            N_data = 10000;
            N_data_8psk = N_data-mod(N_data,3);  %因为8psk每个符号对应3bit，这里把比特数改写为3的倍数
            s0 = randi([0,1],1,N_data_8psk);   %随机产生【0  1】伪随机序列
            L = N_data_8psk/3;  %代表发射的符号数量，等于比特数除以log2(8)
            St = A1*yuanma_8psk([4 2 1]*reshape(s0,3,L)+1);  %这里去了解一下reshape函数，[4 2 1]分别代表8psk符号每个比特的权值，比如101对应5，因为数组标号是从1开始的，而格雷码是从0开始的，所以需要加1
            
            h = sqrt(1/2)*(randn(1,L)+1i*randn(1,L));  %产生瑞利信道，我们这里不考虑大尺度衰落，瑞利信道方差为1
            n = 1/sqrt(2)*(wgn(1,L,noise_dB)+1i*wgn(1,L,noise_dB)); %产生噪声，同样是个复信号
            
            r = h.*St+n;  %天线处的接收信号，因为这里St 和  R都是矩阵，乘法除法都需要加一个  .  号
            S_e = conj(h).*r/A1; %检测信号
            minx = abs(bsxfun(@minus,yuanma_8psk,S_e.')).^2;  %了解一下bsxfun函数，这里就是求检测信号与每个原码的距离
            %注意  .'代表转置，   ‘代表共轭转置，         区分使用
            [ms, index] = min(minx,[],2);  %找到差距最小的那个原码的标号
            Sr_d = reshape((double(dec2bin(index'-1))-48)',1,N_data_8psk); %把那个标号变回原来的比特流，了解一下dec2bin和double函数
            [number_BER, ratio_BER] =symerr(s0, Sr_d);  %了解一下symerr函数，就是求BER的
            one_ats_BER(n_sim, j) = ratio_BER;
            
            Sr_c = A1*yuanma_8psk(index);
            [number_SER, ratio_SER] =symerr(St, Sr_c);
            one_ats_SER(n_sim, j) = ratio_SER; %误符号率
        end
        
        %% 2 antennas/8PSK/G2
        one_Enabled = 1;
        if one_Enabled ==1
            if j<15
                N_data = 1000;
            elseif j<26
                N_data = 10000;
            else
                N_data = 2*100000;
            end
            N_data_8psk = N_data-mod(N_data,6);
            s0 = randi([0,1],1,N_data_8psk);
            L = N_data_8psk/3;
            St = A2*yuanma_8psk([4 2 1]*reshape(s0,3,L)+1);
            St_re = reshape(St,2,L/2);
            St1 = St_re(1,:);
            St2 = St_re(2,:);
            
            h1 = sqrt(1/2)*(randn(1,L/2)+1i*randn(1,L/2));
            h2 = sqrt(1/2)*(randn(1,L/2)+1i*randn(1,L/2));
            n1 = 1/sqrt(2)*(wgn(1,L/2,noise_dB)+1i*wgn(1,L/2,noise_dB));
            n2 = 1/sqrt(2)*(wgn(1,L/2,noise_dB)+1i*wgn(1,L/2,noise_dB));
            
            r1 = h1.*St1 + h2.*St2 + n1;
            r2 = - h1.*conj(St2) + h2.*conj(St1) + n2;
            
            S_e1 = conj(h1).*r1+h2.*conj(r2);
            S_e2 = conj(h2).*r1-h1.*conj(r2);
            minx1 = abs(bsxfun(@minus,yuanma_8psk,S_e1.'/A2)).^2;
            minx2 = abs(bsxfun(@minus,yuanma_8psk,S_e2.'/A2)).^2;
            [ms1, index1] = min(minx1,[],2);
            [ms2, index2] = min(minx2,[],2);
            Sr_d1 = (double(dec2bin(index1.'-1,3))-48).';  %把那个标号变回原来的比特流，dec2bin(D)把十进制数D转换成二进制形式，double函数
            Sr_d2 = (double(dec2bin(index2.'-1,3))-48).';
            Sr_d = reshape([Sr_d1; Sr_d2],1,N_data_8psk);
            [number_BER, ratio_BER] =symerr(s0, Sr_d);
            two_ats_BER(n_sim, j) = ratio_BER; %误比特率
            
            Sr_c1 = A2*yuanma_8psk(index1);
            Sr_c2 = A2*yuanma_8psk(index2);
            [number_SER1, ratio_SER1] =symerr(St1, Sr_c1);
            [number_SER2, ratio_SER2] =symerr(St2, Sr_c2);
            two_ats_SER(n_sim, j) = (ratio_SER1+ratio_SER2)/2; %误符号率
        end
        
        %% 3 antennas/16qam/H3
        one_Enabled = 1;
        if one_Enabled ==1
            if j<18
                N_data = 10000;
            elseif j<32
                N_data = 2*100000;
            else
                N_data = 20;
            end
            N_data_16qam = N_data-mod(N_data,12);
            L = N_data_16qam/4;
            s0 = randi([0,1],1,N_data_16qam); 
            St = A3*yuanma_16qam([8 4 2 1]*reshape(s0,4,L)+1);
            St_re = reshape(St,3,L/3);
            St1 = St_re(1,:);
            St2 = St_re(2,:);
            St3 = St_re(3,:);
            
            h1 = sqrt(1/2)*(randn(1,L/3)+1i*randn(1,L/3));
            h2 = sqrt(1/2)*(randn(1,L/3)+1i*randn(1,L/3));
            h3 = sqrt(1/2)*(randn(1,L/3)+1i*randn(1,L/3));
            n1 = 1/sqrt(2)*(wgn(1,L/3,noise_dB)+1i*wgn(1,L/3,noise_dB));
            n2 = 1/sqrt(2)*(wgn(1,L/3,noise_dB)+1i*wgn(1,L/3,noise_dB));
            n3 = 1/sqrt(2)*(wgn(1,L/3,noise_dB)+1i*wgn(1,L/3,noise_dB));
            n4 = 1/sqrt(2)*(wgn(1,L/3,noise_dB)+1i*wgn(1,L/3,noise_dB));
            
            r1 = h1.*St1+h2.*St2+h3.*St3/sqrt(2)+n1;
            r2 = -h1.*conj(St2)+h2.*conj(St1)+h3.*St3/sqrt(2)+n2;
            r3 = h1.*conj(St3)/sqrt(2)+h2.*conj(St3)/sqrt(2)+h3.*(-St1-conj(St1)+St2-conj(St2))/2+n3;
            r4 = h1.*conj(St3)/sqrt(2)-h2.*conj(St3)/sqrt(2)+h3.*(St2+conj(St2)+St1-conj(St1))/2+n4;
            
            S_e1 = r1.*conj(h1)+conj(r2).*h2+(r4-r3).*conj(h3)/2-conj(r3+r4).*h3/2;
            S_e2 = r1.*conj(h2)-conj(r2).*h1+(r4+r3).*conj(h3)/2+conj(-r3+r4).*h3/2;
            S_e3 = (r1+r2).*conj(h3)/sqrt(2)+conj(r3).*(h1+h2)/sqrt(2)+conj(r4).*(h1-h2)/sqrt(2);
            K = ((abs(h1).^2+abs(h2).^2+abs(h3).^2).'-1)*abs(yuanma_16qam).^2;
            minx1 = abs(bsxfun(@minus,yuanma_16qam,S_e1.'/A3)).^2+K;
            minx2 = abs(bsxfun(@minus,yuanma_16qam,S_e2.'/A3)).^2+K;
            minx3 = abs(bsxfun(@minus,yuanma_16qam,S_e3.'/A3)).^2+K;
            [ms1, index1] = min(minx1,[],2);
            [ms2, index2] = min(minx2,[],2);
            [ms3, index3] = min(minx3,[],2);
            Sr_d1 = (double(dec2bin(index1.'-1,4))-48).';
            Sr_d2 = (double(dec2bin(index2.'-1,4))-48).';
            Sr_d3 = (double(dec2bin(index3.'-1,4))-48).';
            Sr_d = reshape([Sr_d1; Sr_d2; Sr_d3],1,N_data_16qam);
            [number_BER, ratio_BER] = symerr(s0, Sr_d);
            three_ats_BER(n_sim, j) = ratio_BER; 
            
            Sr_c1 = A3*yuanma_16qam(index1);
            Sr_c2 = A3*yuanma_16qam(index2);
            Sr_c3 = A3*yuanma_16qam(index3);
            [number_SER1, ratio_SER1] =symerr(St1, Sr_c1);
            [number_SER2, ratio_SER2] =symerr(St2, Sr_c2);
            [number_SER3, ratio_SER3] =symerr(St3, Sr_c3);
            three_ats_SER(n_sim, j) = (ratio_SER1+ratio_SER2+ratio_SER3)/3; 
        end
        
        %% 4 antennas/16QAM/H4
        one_Enabled = 1;
        if one_Enabled ==1
            if j<17
                N_data = 10000;
            elseif j<26
                N_data = 100000;
            else
                N_data = 20;
            end
            N_data_16qam = N_data-mod(N_data,12);
            L = N_data_16qam/4;
            s0 = randi([0,1],1,N_data_16qam); 
            St = A4*yuanma_16qam([8 4 2 1]*reshape(s0,4,L)+1);
            St_re = reshape(St,3,L/3);
            St1 = St_re(1,:);
            St2 = St_re(2,:);
            St3 = St_re(3,:);
            
            h1 = sqrt(1/2)*(randn(1,L/3)+1i*randn(1,L/3));
            h2 = sqrt(1/2)*(randn(1,L/3)+1i*randn(1,L/3));
            h3 = sqrt(1/2)*(randn(1,L/3)+1i*randn(1,L/3));
            h4 = sqrt(1/2)*(randn(1,L/3)+1i*randn(1,L/3));
            n1 = 1/sqrt(2)*(wgn(1,L/3,noise_dB)+1i*wgn(1,L/3,noise_dB));
            n2 = 1/sqrt(2)*(wgn(1,L/3,noise_dB)+1i*wgn(1,L/3,noise_dB));
            n3 = 1/sqrt(2)*(wgn(1,L/3,noise_dB)+1i*wgn(1,L/3,noise_dB));
            n4 = 1/sqrt(2)*(wgn(1,L/3,noise_dB)+1i*wgn(1,L/3,noise_dB));
            
            r1 = h1.*St1+h2.*St2+h3.*St3/sqrt(2)+h4.*St3/sqrt(2)+n1;
            r2 = -h1.*conj(St2)+h2.*conj(St1)+h3.*St3/sqrt(2)-h4.*St3/sqrt(2)+n2;
            r3 = h1.*conj(St3)/sqrt(2)+h2.*conj(St3)/sqrt(2)+h3.*(-St1-conj(St1)+St2-conj(St2))/2+h4.*(-St2-conj(St2)+St1-conj(St1))/2+n3;
            r4 = h1.*conj(St3)/sqrt(2)-h2.*conj(St3)/sqrt(2)+h3.*(St2+conj(St2)+St1-conj(St1))/2-h4.*(St1+conj(St1)+St2-conj(St2))/2+n4;
            
            S_e1 = r1.*conj(h1)+conj(r2).*h2+(r4-r3).*(conj(h3)-conj(h4))/2-conj(r3+r4).*(h3+h4)/2;
            S_e2 = r1.*conj(h2)-conj(r2).*h1+(r4+r3).*(conj(h3)-conj(h4))/2+conj(-r3+r4).*(h3+h4)/2;
            S_e3 =(r1+r2).*conj(h3)/sqrt(2)+(r1-r2).*conj(h4)/sqrt(2)+conj(r3).*(h1+h2)/sqrt(2)+conj(r4).*(h1-h2)/sqrt(2);
            K = ((abs(h1).^2+abs(h2).^2+abs(h3).^2+abs(h4).^2).'-1)*abs(yuanma_16qam).^2;
            minx1 = abs(bsxfun(@minus,yuanma_16qam,S_e1.'/A4)).^2+K;
            minx2 = abs(bsxfun(@minus,yuanma_16qam,S_e2.'/A4)).^2+K;
            minx3 = abs(bsxfun(@minus,yuanma_16qam,S_e3.'/A4)).^2+K;
            [ms1, index1] = min(minx1,[],2);
            [ms2, index2] = min(minx2,[],2);
            [ms3, index3] = min(minx3,[],2);
            Sr_d1 = (double(dec2bin(index1.'-1,4))-48).';
            Sr_d2 = (double(dec2bin(index2.'-1,4))-48).';
            Sr_d3 = (double(dec2bin(index3.'-1,4))-48).';
            Sr_d = reshape([Sr_d1; Sr_d2; Sr_d3],1,N_data_16qam);
            [number_BER, ratio_BER] = symerr(s0, Sr_d);
            four_ats_BER(n_sim, j) = ratio_BER; 
            
            Sr_c1 = A4*yuanma_16qam(index1);
            Sr_c2 = A4*yuanma_16qam(index2);
            Sr_c3 = A4*yuanma_16qam(index3);
            [number_SER1, ratio_SER1] =symerr(St1, Sr_c1);
            [number_SER2, ratio_SER2] =symerr(St2, Sr_c2);
            [number_SER3, ratio_SER3] =symerr(St3, Sr_c3);
            four_ats_SER(n_sim, j) = (ratio_SER1+ratio_SER2+ratio_SER3)/3; 
        end
    end
    toc
end

one_ats_mean = mean(one_ats_BER);
two_ats_mean = mean(two_ats_BER);
three_ats_mean = mean(three_ats_BER);
four_ats_mean = mean(four_ats_BER);

one_ats_sym_mean = mean(one_ats_SER);
two_ats_sym_mean = mean(two_ats_SER);
three_ats_sym_mean = mean(three_ats_SER);
four_ats_sym_mean = mean(four_ats_SER);   

figure(1)
semilogy(SNR,one_ats_mean,'m-+','LineWidth',1);
hold on;
semilogy(SNR,two_ats_mean,'r-*','LineWidth',1);
hold on;
semilogy(SNR,three_ats_mean,'g-x','LineWidth',1);
hold on;
semilogy(SNR,four_ats_mean,'c-o','LineWidth',1);
hold on;
legend('uncoded','2 antennas','3 antennas','4 antennas');
axis([5 39 1e-6 1e-1]);
xlabel('SNR(dB)'); 
ylabel('BER');
title('3bits/sec/Hz');
grid on;

figure(2);
semilogy(SNR,one_ats_sym_mean,'m-+','LineWidth',1);
hold on;
semilogy(SNR,two_ats_sym_mean,'r-*','LineWidth',1);
hold on;
semilogy(SNR,three_ats_sym_mean,'g-x','LineWidth',1);
hold on;
semilogy(SNR,four_ats_sym_mean,'c-o','LineWidth',1);
hold on;
legend('uncoded','2 antennas','3 antennas','4 antennas');
axis([5 39 1e-6 1]);
xlabel('SNR(dB)'); 
ylabel('SER');
title('3bits/sec/Hz');
grid on;
        
toc

save data.mat;
print('-f1','-dpng','pic1.png');
print('-f2','-dpng','pic2.png');
