% function [modulated_symbols,Mod] = mapper(b,N)
% % If N is given, it generates a block of N random 2^b-PSK/QAM modulated symbols.
% % Otherwise, it generates a block of 2^b-PSK/QAM modulated symbols for [0:2^b-1].
% 
% %MIMO-OFDM Wireless Communications with MATLAB¢ç   Yong Soo Cho, Jaekwon Kim, Won Young Yang and Chung G. Kang
% %?2010 John Wiley & Sons (Asia) Pte Ltd
% 
% M=2^b; % Modulation order or Alphabet (Symbol) size
% if b==1, Mod='BPSK'; A=1; mod_object=pskmod('M',M);
%  elseif b==2, Mod='QPSK';  A=1;
%       mod_object=pskmod('M',M,'PhaseOffset',pi/4);
%  else Mod=[num2str(2^b) 'QAM']; Es=1; A=sqrt(3/2/(M-1)*Es); 
%       mod_object=qammod('M',M,'SymbolOrder','gray');   
% end
% if nargin==2 % generates a block of N random 2^b-PSK/QAM modulated symbols 
% %   modulated_symbols = A*modulate(mod_object,randint(1,N,M));
%   modulated_symbols = A*modulate(mod_object,randi(M,1,N)-1);
%  else
%   modulated_symbols = A*modulate(mod_object,[0:M-1]);
% end
function [modulated_symbols, Mod] = mapper(b, N)
    % If N is given, it generates a block of N random 2^b-PSK/QAM modulated symbols.
    % Otherwise, it generates a block of 2^b-PSK/QAM modulated symbols for [0:2^b-1].
    
    M = 2^b; % Modulation order or Alphabet (Symbol) size
    
    if b == 1
        Mod = 'BPSK';
        A = 1;
        data = randi([0 M-1], N, 1); % Generate random data for modulation
        modulated_symbols = pskmod(data, M);
    elseif b == 2
        Mod = 'QPSK';
        A = 1;
        data = randi([0 M-1], N, 1); % Generate random data for modulation
        modulated_symbols = pskmod(data, M, pi/4); % Phase offset for QPSK
    else
        Mod = [num2str(2^b) 'QAM'];
        Es = 1;
        A = sqrt(3/2/(M-1)*Es);
        if nargin == 2 % Generates a block of N random 2^b-QAM modulated symbols
            data = randi([0 M-1], N, 1); % Generate random data for modulation
        else
            data = [0:M-1]; % Use all possible symbols
        end
        modulated_symbols = qammod(data, M, 'gray');
    end
    
    % Apply amplitude scaling to the modulated symbols
    modulated_symbols = A * modulated_symbols;
end