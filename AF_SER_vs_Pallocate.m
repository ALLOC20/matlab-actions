close all;
clc; clear;

tStart = tic;

Nb = 1;
SNRdB = [10 20 30];
rate = 0:0.1:1;
dsd = [0.1 1 10]; dsr = 1; drd = 1;
global M b_psk
M = 4;
b_psk = sin(pi / M).^2;
K = 1 - 1 / sqrt(M);

for j = 1:length(dsd)
    SER_rate = zeros(length(SNRdB), length(rate));

    for i = 1:length(SNRdB)
        tSNRStart = tic;
        P_all = 10^(SNRdB(i) / 10);
        AQ = (M - 1) / 2 / M + K.^2 / pi;
        BQ = 3 * (M - 1) / 8 / M + K.^2 / pi;
        %%%%%%%%功率分配%%%%%
        for ii = 1:length(rate)
            P_1 = rate(ii)*P_all;
            P_2 = P_all - P_1;
            %%%%%%%%%%%%%%%%%%%精确,N0=1%%%%%%%%%%%%%%%%%%%%%%%%%%\
            tIntStart = tic;
            SER_rate(i,ii) = Exact_SER(P_1, P_2, dsd(j), dsr, drd);
            tIntEnd = toc(tIntStart);
            disp(['rate = ', num2str(rate(ii)), '积分耗时 ', num2str(tIntEnd), ' 秒']);
        end

        tSNREnd = toc(tSNRStart);
        disp(['i=', num2str(i), '，耗时 ', num2str(tSNREnd), ' 秒']);
    end

    figure(1);
    subplot(130+j);
    semilogy(rate, SER_rate(1,:), 'r-.', rate, SER_rate(2,:), 'k--', rate, SER_rate(3,:), 'b-');
    hold on
    plot([0.6270,0.6270],ylim,'m:');
    hold on
    grid on
    title(['\delta_{r,d}^2 = ',num2str(dsd(j))]);
    xlabel('P1/P (%)');
    ylabel('SER');
    legend('P=10dBw', 'P=20dBw', 'P=30dBw', 'Location', 'south')
    axis([0 1 10^(-6.5) 1e0]);
end
tEnd = toc(tStart);
disp(['总耗时 ', num2str(tEnd), ' 秒']);

