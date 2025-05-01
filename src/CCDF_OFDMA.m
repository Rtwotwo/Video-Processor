function CCDF = CCDF_OFDMA(N, Nos, b, dBs, Nblk)
% CCDF of OFDM signal with no PAPR reduction technique.
% N    : Number of total subcarriers
% Nos  : Oversampling factor
% b    : Bits per symbol (e.g., 2 for QPSK, 4 for 16QAM)
% dBs  : PAPR threshold vector in dB
% Nblk : Number of blocks

dBcs = dBs + (dBs(2)-dBs(1))/2;
NNos = N * Nos; 
M = 2^b;

% 正确归一化因子：使得 E{|x|^2} = 1
A = modnorm(qammod(0:M-1, M), 'avpow', 1); % 归一化因子

for nblk = 1:Nblk
    % 随机生成二进制消息 -> 符号索引 -> QAM 映射
    bits = randi([0 1], b, N);
    sym_idx = bi2de(bits', 'left-msb');     % 二进制转十进制索引
    mod_sym = qammod(sym_idx, M) * A;        % QAM调制 + 功率归一化
    
    % Oversampling using zero-padding in frequency domain
    mod_sym = mod_sym(:)';
    X_pad = [mod_sym, zeros(1, N*(Nos - 1))]; % 插零实现过采样
    ifft_sym = ifft(ifftshift(X_pad));       % IFFT
    
    % 计算PAPR
    power = abs(ifft_sym).^2;
    mean_pow(nblk) = mean(power);
    max_pow(nblk) = max(power);
end

PAPRdB = 10*log10(max_pow ./ mean_pow); % 单位：dB

% 计算 CCDF
CCDF = zeros(size(dBs));
for i = 1:length(dBs)
    CCDF(i) = sum(PAPRdB > dBs(i)) / Nblk;
end

end