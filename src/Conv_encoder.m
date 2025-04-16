% 	Convolutional encoding with binary data
function coded_bits = Conv_encoder(in_bits)

%Constraint Length K=7 (Memory size=6 => State=64)
ConvCodeGenPoly=[1 0 1 1 0 1 1 ;...
                              1 1 1 1 0 0 1 ];
number_rows = size(ConvCodeGenPoly, 1);
number_bits = size(ConvCodeGenPoly,2)+length(in_bits)-1;

uncoded_bits = zeros(number_rows, number_bits);

for row=1:number_rows
    uncoded_bits(row,1:number_bits)=rem(conv(in_bits, ConvCodeGenPoly(row,:)),2);
end

coded_bits=uncoded_bits(:);
