function [BER,Neb,Ntb]=ber(x,x_ref,Nbps)
%   Bit Error Rate computes
%   1) ber2(): internal buffer clear
%   2) BER=ber2(x,x_ref,Nbps)
%   Inputs:
%       x       : input bits(integer type)
%       x_ref   : reference bits
%       Nbps     : bits per symbol
%   outputs:
%       BER   : BER 
%
%   Original Author : Toystar, CAU (toystar@dreamwiz.com)
%   Date            : 2004-02-04

persistent N_error_bits
persistent N_total_bits

if nargin==0
    N_error_bits=0;  N_total_bits=0;  BER=0;
 elseif nargin==3
  for i=1:length(x)
     N_error_bits=N_error_bits+error_check(bitxor(x(i),x_ref(i)),Nbps);
  end
  N_total_bits= N_total_bits + length(x)*Nbps;
  BER= N_error_bits/N_total_bits;
 else error('Dont match with # of arguments'); 
end
Neb=N_error_bits;  Ntb=N_total_bits;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ecount=error_check(x_xor,Nbps)
%   Integer bit check
%
ecount=0;
for i=0:Nbps-1
   ecount=ecount+bitand(bitshift(x_xor,-i),1);
end