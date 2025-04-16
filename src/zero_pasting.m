function y=zero_pasting(x)
% Padd zeros at the center of the input sequence x

%MIMO-OFDM Wireless Communications with MATLAB¢ç   Yong Soo Cho, Jaekwon Kim, Won Young Yang and Chung G. Kang
%?2010 John Wiley & Sons (Asia) Pte Ltd

N=length(x); M=ceil(N/4);
y = [x(1:M) zeros(1,N/2) x(N-M+1:N)];
