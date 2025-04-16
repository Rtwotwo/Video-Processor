function [STO_est, Mag] = STO_by_correlation(y, Nfft, Ng, com_delay)
% STO estimation by maximizing the correlation between CP and rear part of OFDM symbol
% estimates STO by maximizing the correlation between CP (cyclic prefix)  
%     and rear part of OFDM symbol
% Input:  y         = Received OFDM signal including CP
%         Ng        = Number of samples in Guard Interval (CP)
%         com_delay = Common delay
% Output: STO_est   = STO estimate
%         Mag       = Correlation function trajectory varying with time

%MIMO-OFDM Wireless Communications with MATLAB㈢   Yong Soo Cho, Jaekwon Kim, Won Young Yang and Chung G. Kang
%?2010 John Wiley & Sons (Asia) Pte Ltd

N_ofdm = Nfft + Ng; 
if nargin < 4, com_delay = N_ofdm / 2; end
%maximum=1e-8; 
nn = 0 : Ng-1; 
yy = y(nn + com_delay) * y(nn+com_delay+Nfft)'; % 计算两个窗中采样的相关值的和，从com_delay往后Ng长度
maximum = abs(yy);
STO_est = N_ofdm - com_delay - 1; 
Mag(1) = maximum;
for n = 1 : (N_ofdm - 1)
   n1 = n - 1;
   yy1 = y(n1 + com_delay) * y(n1 + com_delay + Nfft)';  % 当前窗的第一个采样
   yy2 = y(n1 + com_delay + Ng) * y(n1 + com_delay + Nfft + Ng)'; % 当前窗后面第一个采样
   yy = yy - yy1 + yy2; % 减掉第一个采样的相关值，加上创后面的第一个采样值，相当于窗向后划了一个采样   
   Mag(n + 1) = abs(yy); % Eq.(5.13)
   if (Mag(n + 1) > maximum)
     maximum = Mag(n + 1); 
     STO_est = N_ofdm - com_delay - n1; 
   end
end

[Mag_cor_max, ind_max] = max(Mag); 
%if (STO_est>=Ng/2), STO_est= Ng/2-1;
% elseif (STO_est<-Ng/2), STO_est= -Ng/2;
%end    
% [a b] =max(Mag)
% figure(3), plot(Mag), pause