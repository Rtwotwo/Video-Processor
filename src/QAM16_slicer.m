function [X_hat] = QAM16_slicer(X,N)
%MIMO-OFDM Wireless Communications with MATLAB¢ç   Yong Soo Cho, Jaekwon Kim, Won Young Yang and Chung G. Kang
%2010 John Wiley & Sons (Asia) Pte Ltd
sq10=sqrt(10); b = [-2 0 2]/sq10; c = [-3 -1 1 3]/sq10;
[Nr,Nc]=size(X);  M=Nr*Nc;
X_super = X(:);
complex_X = 0;
if norm(imag(X))>1e-5
  complex_X=1;  X_super=[real(X_super); imag(X_super)]; 
end
X_hat(find(X_super<b(1))) = c(1);  
X_hat(find(b(1)<=X_super&X_super<b(2))) = c(2);  
X_hat(find(b(2)<=X_super&X_super<b(3))) = c(3);  
X_hat(find(b(3)<=X_super)) = c(4);  
if complex_X==1,  X_hat = X_hat(1:M)+j*X_hat(M+1:2*M);  end
X_hat = reshape(X_hat,Nr,Nc);

%if nargin<2,  N = length(X);  end
%sq10=sqrt(10); b = [-2 0 2]/sq10; c = [-3 -1 1 3]/sq10;
%Xr = real(X);  Xi = imag(X);
%for i=1:N
%   R(find(Xr<b(1))) = c(1);  I(find(Xi<b(1))) = c(1);
%   R(find(b(1)<=Xr&Xr<b(2))) = c(2);  I(find(b(1)<=Xi&Xi<b(2))) = c(2);
%   R(find(b(2)<=Xr&Xr<b(3))) = c(3);  I(find(b(2)<=Xi&Xi<b(3))) = c(3);
%   R(find(b(3)<=Xr)) = c(4);  I(find(b(3)<=Xi)) = c(4);
%end
%X_hat = R + j*I;