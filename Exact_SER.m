function f = Exact_SER(SNR1, SNR2, Deta_sd, Deta_sr, Deta_rd)
global M b_psk

syms x;
sj = sqrt((1./SNR1/Deta_sr-1./SNR2/Deta_rd).^2+2*(1./SNR1/Deta_sr+1./SNR2/Deta_rd)*b_psk/(sin(x)^2)+b_psk^2/(sin(x)^4));
y = 1./(1+b_psk*SNR1*Deta_sd./(sin(x)^2)).*(((1./SNR1/Deta_sr-1./SNR2/Deta_rd).^2+(1./SNR1/Deta_sr+1./SNR2/Deta_rd)*b_psk/(sin(x)^2))./(sj.^2)+...
    +2*b_psk.*log(SNR1.*SNR2*Deta_sr*Deta_rd.*(1./SNR1/Deta_sr+1./SNR2/Deta_rd+b_psk/sin(x)^2+sj).^2/4)./SNR1./SNR2./(sj.^3)/Deta_sr/Deta_rd/sin(x)^2);
f = 1/pi*int(y, x, 0, (M-1)*pi/M);

end