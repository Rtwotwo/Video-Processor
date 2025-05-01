function CCDF = CCDF_PTS(N, Nos, Nsb, b, dBs, Nblk)
% CCDF of PTS (Partial Transmit Sequence) technique.
% N    : Number of Subcarriers 
% Nos  : Oversampling factor
% Nsb  : Number of subblocks
% b    : Number of bits per QAM symbol
% dBs  : dB vector
% Nblk : Number of OFDM blocks for iteration

%MIMO-OFDM Wireless Communications with MATLAB㈢   Yong Soo Cho, Jaekwon Kim, Won Young Yang and Chung G. Kang
%?2010 John Wiley & Sons (Asia) Pte Ltd

NNos = N * Nos;                       % FFT size
M = 2^b;                              % Alphabet size
Es = 1;                               % Energy per symbol
A = sqrt(3/2/(M-1)*Es);               % Normalization factor for M-QAM

% 手动构造 Gray 映射的星座图
tmp = qammod(0:M-1, M);               % 自然编码的星座点
[~, idx] = bitrevorder(tmp);          % 得到 Gray 编码的排列顺序
gray_constellation = tmp(idx);        % Gray 映射星座点

PAPRs = zeros(1, Nblk);               % 初始化 PAPR 数组

for iter = 1:Nblk
    data_symbols = randi([0 M-1], 1, N);        % 随机生成符号索引
    mod_sym = A * gray_constellation(data_symbols + 1);  % 星座点映射 + 功率归一化
    
    zero_pad_sym = zeros(1, NNos);              % 初始化过采样信号
    step = Nos;
    indices = 1:step:NNos;
    zero_pad_sym(indices) = mod_sym;            % 插零实现频域扩展

    sub_block = zeros(Nsb, NNos);               % 初始化子块矩阵
    for k = 1:Nsb                               % 子块划分
        kk = (k-1)*NNos/Nsb+1 : k*NNos/Nsb;
        sub_block(k, kk) = zero_pad_sym(kk);
    end

    ifft_sym = ifft(sub_block, [], 2);          % 对每一行做 IFFT

    % Phase Factor Optimization
    w = ones(Nsb, 1);                           % 相位因子初始化为列向量
    x = sum(bsxfun(@times, ifft_sym, w), 1);    % 加权和
    sym_pow = abs(x).^2;
    PAPR = max(sym_pow)/mean(sym_pow);

    for m = 1:Nsb
        w(m) = -w(m);                           % 反转当前相位因子
        x_new = sum(bsxfun(@times, ifft_sym, w), 1); % 新的加权和
        sym_pow_new = abs(x_new).^2;
        PAPR_new = max(sym_pow_new)/mean(sym_pow_new);

        if PAPR_new < PAPR                       % 如果新PAPR更小，则更新
            PAPR = PAPR_new;
        else
            w(m) = -w(m);                        % 否则恢复原来的相位因子
        end
    end

    PAPRs(iter) = PAPR;                          % 保存最小PAPR
end

PAPRdBs = 10*log10(PAPRs);                   % 将PAPR转换为dB

% 计算 CCDF
CCDF = zeros(size(dBs));
for i = 1:length(dBs)
    CCDF(i) = sum(PAPRdBs > dBs(i)) / Nblk;
end

plot_or_not = 0;
if plot_or_not > 0
    figure(1), clf
    semilogy(dBs, CCDF, '-s'); axis([dBs([1 end]) 1e-4 1]); grid on; hold on
    title(['16 QAM CCDF of OFDMA PAPR, ', num2str(N), '-point ', num2str(Nblk), '-blocks']);
    xlabel('PAPR_0 [dB]'); ylabel('Pr(PAPR>PAPR_0)');
end