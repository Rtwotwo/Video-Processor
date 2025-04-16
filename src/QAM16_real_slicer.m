function [X_hat] = QAM16_real_slicer(X,N)
sq10=sqrt(10); b = [-2 0 2]/sq10; c = [-3 -1 1 3]/sq10;
X = real(X);  % Against the case where X happens to be complex 
X_hat(find(X<b(1))) = c(1);  
X_hat(find(b(1)<=&X<b(2))) = c(2);  
X_hat(find(b(2)<=X&X<b(3))) = c(3);  
X_hat(find(b(3)<=X)) = c(4);  