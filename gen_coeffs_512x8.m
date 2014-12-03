
% will run inside octave
%%%  -- Grant's logic

clear;

 N=512; 
 D = 8;
 cof_bit = 12;

 b = fir1(N*D-1,1/N,'low',kaiser(N*D,5));    % FIR filter

 b = (1-2^-(cof_bit-1))*b/max(b);            % scale to less than 1 (signed)
 ROM_cof12b_512x8 = round(b*2^(cof_bit-1));  % 2009 Coeff

 save MwaPfbProtoFilterCoeff2009_512x8.dat ROM_cof12b_512x8 -ascii

%%% -- from Grant