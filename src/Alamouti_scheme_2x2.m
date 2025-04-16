% Alamouti_scheme_2x2.m
clear; %clf;
L_frame=130; N_Packets=1e5; % Number of frames/packet and Number of packets 
NT=2; NR=2; b=2; M=2^b;
h_mod = modem.pskmod('M',M, 'PhaseOffset',pi/4, 'InputType','Bit', 'SymbolOrder','gray'); % for QPSK modulation
h_dem = modem.pskdemod(h_mod); % for QPSK demodulation
EbN0dBs=[0:5:30];  sq_NT=sqrt(NT); sq2=sqrt(2); sqb=sqrt(b);
sym_tab = modulate(h_mod,de2bi([0:M-1],b,'left-msb')'); % for QPSK modulation
for i_EbN0=1:length(EbN0dBs)
   EbN0dB=EbN0dBs(i_EbN0);  
   sigma=sqrt(b/2/(10^(EbN0dB/10)));
   H = (randn(NR,NT)+j*randn(NR,NT))/sq2;
   Habs = sum(sum(abs(H).^2));  nobe = 0;
   for i_packet=1:N_Packets
      msg_bits=randint(b*NT,L_frame);
      X = [];
      for nT=1:NT
        %tmp1 = sym_tab(b2d*reshape(msg_bits(:,nT),b,L_frame)+1); 
        X = [X; modulate(h_mod,msg_bits(nT*b+[-b+1:0],:))];
      end
      X1=X;  X2=[-conj(X(2,:)); conj(X(1,:))];
      R1 = H/sq_NT*X1 + sigma*(randn(NR,L_frame)+j*randn(NR,L_frame));
      R2 = H/sq_NT*X2 + sigma*(randn(NR,L_frame)+j*randn(NR,L_frame));
      Z(1,:) = H(:,1)'*R1 + H(:,2).'*conj(R2);
      Z(2,:) = H(:,2)'*R1 - H(:,1).'*conj(R2);
      for m=1:M
        tmp = (-1+sum(Habs,2))*abs(sym_tab(m))^2;
        d1(m,:) = abs(Z(1,:)-sym_tab(m)).^2 + tmp;
        d2(m,:) = abs(Z(2,:)-sym_tab(m)).^2 + tmp;
      end
      [y1,i1]=min(d1);  S1d=sym_tab(i1); % clear d1
      [y2,i2]=min(d2);  S2d=sym_tab(i2); % clear d2
      msg_hat_bits = [demodulate(h_dem,S1d); demodulate(h_dem,S2d)];
      nobe = nobe + sum(sum(msg_bits~=msg_hat_bits)); 
      if nobe>100,  break;  end
   end % End of FOR loop for i_packet
   BER(i_EbN0) = nobe/(i_packet*L_frame*b*NT);
end    % End of FOR loop for i_EbN0
figure(2), clf
semilogy(EbN0dBs,BER,'r'), axis([EbN0dBs([1 end]) 1e-6 1e0]); 
grid on;  xlabel('EbN0[dB]'); ylabel('BER');  hold on