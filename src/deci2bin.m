function y=deci2bin(x,b,KC)
%converts given decimal numbers into binary numbers of b bits
y=[];
if nargin<2
  xmax=max(x); b=1; 
  while xmax>1, xmax=floor(xmax/2); b=b+1; end
end
for n=1:length(x)
  xn=x(n); 
  for i=b:-1:1
    xn_1= floor(xn/2);
    yn(i)= xn-2*xn_1; xn= xn_1;
  end
  if xn>0, fprintf('\n The %dth input number exceeding %d bits in deci2bin()!!\n',n,b); end
  y= [y yn];
end
if nargin<3, KC=0; end
if KC>0 
  for i=1:length(y)
     if y(1)==0, y=y(2:end); end
  end    
end