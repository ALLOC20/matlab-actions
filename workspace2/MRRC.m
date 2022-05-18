close all
clear
clc

N_sim = 1000;%样本数
N_data = 100000;%信号比特数
SNR = 1:2:50;%信噪比
noise_dB = -134;%噪声功率dBw
noise = 10^(noise_dB/10);%噪声功率W

yuanma_bpsk = [-1 1];%BPSK码本

%% 1bits/s/Hz??one receive antenna
[T1R1_BER, T1R2_BER, T1R4_BER, T2R1_BER, T2R2_BER] = deal(zeros(N_sim,length(SNR)));%定义为N_sim*length(SNR)的矩阵

for n_sim = 1:N_sim
    fprintf('n_sim = %d;', n_sim);
    tic
    for j = 1:length(SNR)
        Pt = noise*10^(SNR(j)/10);%噪声功率乘信噪比为信号功率
        A1 = sqrt(Pt);
        A2 = sqrt(Pt/2);
        %fprintf('n_sim = %d, j_Pt = %d\n', n_sim, j);
        %% uncoded/BPSK 1TX,1RX
        T1R1_Enabled = 1;
        if T1R1_Enabled == 1
            s0 = 2*randi([0,1],1,N_data)-ones(1,N_data);%产生BPSK信号[-1 1]
            St = A1*s0;
            R = sqrt(1/2)*(randn(1,N_data)+1i*randn(1,N_data));%瑞利信道
            Noise = sqrt(1/2)*(wgn(1,N_data,noise_dB)+1i*wgn(1,N_data,noise_dB));%复高斯白噪声
            
            Sr = St.*R+Noise;%接收天线收到的信号（不考虑信道衰减）
            S_e = conj(R).*Sr/A1;%解码
            minx = abs(bsxfun(@minus,yuanma_bpsk,S_e.')).^2;%最大似然判决
            [ms, index] = min(minx,[],2);
            Sr_d = 2*(index-1)-1;
            [number, ratio] =symerr(s0, Sr_d');
            T1R1_BER(n_sim, j) = ratio;
        end
        %% 1TX,2RX
        T1R2_Enabled = 1;
        if T1R2_Enabled ==1
            s0 = 2*randi([0,1],1,N_data)-ones(1,N_data);%产生BPSK信号[-1 1]
            St = A1*s0;
            
            h0 = sqrt(1/2)*(randn(1,N_data)+1i*randn(1,N_data));
            h1 = sqrt(1/2)*(randn(1,N_data)+1i*randn(1,N_data));
            
            n0 = 1/sqrt(2)*(wgn(1,N_data,noise_dB)+1i*wgn(1,N_data,noise_dB));
            n1 = 1/sqrt(2)*(wgn(1,N_data,noise_dB)+1i*wgn(1,N_data,noise_dB));

            r0 = h0.*St + n0;
            r1 = h1.*St + n1;

            Se = conj(h0).*r0 + conj(h1).*r1;

            minx = abs(bsxfun(@minus,yuanma_bpsk,Se.'/A1)).^2;%最大似然判决
            [ms, index] = min(minx,[],2);%取出最小值及其索引
            
            Sd = 2*(index-1)-1;
            [number, ratio] =symerr(s0, Sd');%计算误码率
            
            T1R2_BER(n_sim, j) = ratio;
        end
        %% 1TX,4RX
        T1R4_Enabled = 1;
        if T1R4_Enabled ==1
            s0 = 2*randi([0,1],1,N_data)-ones(1,N_data);%产生BPSK信号[-1 1]
            St = A1*s0;
            
            h0 = sqrt(1/2)*(randn(1,N_data)+1i*randn(1,N_data));
            h1 = sqrt(1/2)*(randn(1,N_data)+1i*randn(1,N_data));
            h2 = sqrt(1/2)*(randn(1,N_data)+1i*randn(1,N_data));
            h3 = sqrt(1/2)*(randn(1,N_data)+1i*randn(1,N_data));
            
            n0 = 1/sqrt(2)*(wgn(1,N_data,noise_dB)+1i*wgn(1,N_data,noise_dB));
            n1 = 1/sqrt(2)*(wgn(1,N_data,noise_dB)+1i*wgn(1,N_data,noise_dB));
            n2 = 1/sqrt(2)*(wgn(1,N_data,noise_dB)+1i*wgn(1,N_data,noise_dB));
            n3 = 1/sqrt(2)*(wgn(1,N_data,noise_dB)+1i*wgn(1,N_data,noise_dB));

            r0 = h0.*St + n0;
            r1 = h1.*St + n1;
            r2 = h2.*St + n2;
            r3 = h3.*St + n3;

            Se = conj(h0).*r0 + conj(h1).*r1 + conj(h2).*r2 + conj(h3).*r3;

            minx = abs(bsxfun(@minus,yuanma_bpsk,Se.'/A1)).^2;%最大似然判决
            [ms, index] = min(minx,[],2);%取出最小值及其索引
            
            Sd = 2*(index-1)-1;
            [number, ratio] =symerr(s0, Sd');%计算误码率
            
            T1R4_BER(n_sim, j) = ratio;
        end
        %% 2TX,1RX
        T2R1_Enabled = 1;
        if T2R1_Enabled ==1
            s01 = 2*randi([0,1],1,N_data/2)-ones(1,N_data/2);%交错发送
            s02 = 2*randi([0,1],1,N_data/2)-ones(1,N_data/2);
            St01 = A2*reshape([s01;-conj(s02)],1,N_data);%发射天线1的信号，已经交错重整
            St02 = A2*reshape([s02;conj(s01)],1,N_data);%发射天线2的信号，已经交错重整
            
            r01 = sqrt(1/2)*(randn(1,N_data/2)+1i*randn(1,N_data/2));
            r02 = sqrt(1/2)*(randn(1,N_data/2)+1i*randn(1,N_data/2));
            R01 = reshape([r01;r01],1,N_data);
            R02 = reshape([r02;r02],1,N_data);
            Noise = 1/sqrt(2)*(wgn(1,N_data,noise_dB)+1i*wgn(1,N_data,noise_dB));
            
            Sr = R01.*St01+R02.*St02+Noise;%接收天线收到的信号（不考虑信道衰减）
            Sr_re = reshape(Sr,2,N_data/2);%将信号交错拆开，重整为2*5000
            Sr1 = Sr_re(1,:);%取第一行
            Sr2 = Sr_re(2,:);%取第二行
            
            S_e1 = conj(r01).*Sr1+r02.*conj(Sr2);%文章中公式(12)
            S_e2 = conj(r02).*Sr1-r01.*conj(Sr2);
            minx1 = abs(bsxfun(@minus,yuanma_bpsk,S_e1.'/A2)).^2;%最大似然判决
            minx2 = abs(bsxfun(@minus,yuanma_bpsk,S_e2.'/A2)).^2;
            [ms1, index1] = min(minx1,[],2);%取出最小值及其索引
            [ms2, index2] = min(minx2,[],2);
            Sr_d1 = 2*(index1-1)-1;
            Sr_d2 = 2*(index2-1)-1;
            [number1, ratio1] =symerr(s01, Sr_d1');%计算误码率
            [number2, ratio2] =symerr(s02, Sr_d2');
            T2R1_BER(n_sim, j) = (ratio1+ratio2)/2;
        end
        %% 2TX,2RX
        T2R2_Enabled = 1;
        if T2R2_Enabled ==1
            s0 = 2*randi([0,1],1,N_data/2)-ones(1,N_data/2);%交错发送
            s1 = 2*randi([0,1],1,N_data/2)-ones(1,N_data/2);
            St0 = A2*s0;
            St1 = A2*s1;
            
            h0 = sqrt(1/2)*(randn(1,N_data/2)+1i*randn(1,N_data/2));
            h1 = sqrt(1/2)*(randn(1,N_data/2)+1i*randn(1,N_data/2));
            h2 = sqrt(1/2)*(randn(1,N_data/2)+1i*randn(1,N_data/2));
            h3 = sqrt(1/2)*(randn(1,N_data/2)+1i*randn(1,N_data/2));
            
            n0 = 1/sqrt(2)*(wgn(1,N_data/2,noise_dB)+1i*wgn(1,N_data/2,noise_dB));
            n1 = 1/sqrt(2)*(wgn(1,N_data/2,noise_dB)+1i*wgn(1,N_data/2,noise_dB));
            n2 = 1/sqrt(2)*(wgn(1,N_data/2,noise_dB)+1i*wgn(1,N_data/2,noise_dB));
            n3 = 1/sqrt(2)*(wgn(1,N_data/2,noise_dB)+1i*wgn(1,N_data/2,noise_dB));
            
            r0 = h0.*St0 + h1.*St1 + n0;
            r1 = -h0.*conj(St1) + h1.*conj(St0) + n1;
            r2 = h2.*St0 + h3.*St1 + n2;
            r3 = -h2.*conj(St1) + h3.*conj(St0) + n3;
            
            S_e0 = conj(h0).*r0+h1.*conj(r1)+conj(h2).*r2+h3.*conj(r3);%文章中公式(12)
            S_e1 = conj(h1).*r0-h0.*conj(r1)+conj(h3).*r2-h2.*conj(r3);
            
            minx0 = abs(bsxfun(@minus,yuanma_bpsk,S_e0.'/A2)).^2;%最大似然判决
            minx1 = abs(bsxfun(@minus,yuanma_bpsk,S_e1.'/A2)).^2;
            [ms0, index0] = min(minx0,[],2);%取出最小值及其索引
            [ms1, index1] = min(minx1,[],2);
            Sr_d0 = 2*(index0-1)-1;
            Sr_d1 = 2*(index1-1)-1;
            [number0, ratio0] =symerr(s0, Sr_d0');%计算误码率
            [number1, ratio1] =symerr(s1, Sr_d1');
            T2R2_BER(n_sim, j) = (ratio0+ratio1)/2;
        end
    end
    toc
end

T1R1_mean = mean(T1R1_BER);
T1R2_mean = mean(T1R2_BER);
T1R4_mean = mean(T1R4_BER);
T2R1_mean = mean(T2R1_BER);
T2R2_mean = mean(T2R2_BER);

figure(1);
semilogy(SNR,T1R1_mean,'m-o','LineWidth',1.5);
hold on;
semilogy(SNR,T1R2_mean,'r-v','LineWidth',1.5);
hold on;
semilogy(SNR,T1R4_mean,'b-s','LineWidth',1.5);
hold on;
semilogy(SNR,T2R1_mean,'g-d','LineWidth',1.5);
hold on;
semilogy(SNR,T2R2_mean,'k-^','LineWidth',1.5);
hold on;
grid on;
legend('no diversity(1Tx,1Rx)','MRRC(1Tx,2Rx)','MRRC(1Tx,4Rx)','new scheme(2Tx,1Rx)','new scheme(2Tx,2Rx)');
xlabel('SNR'); 
ylabel('BER');

save data.mat
print('-f1','-dpng','pic1.png')