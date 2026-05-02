% Note: This code expects that MCXLAB (link below) in MATLAB's path
% http://mcx.space/wiki/index.cgi?Doc/MCXLAB

addpath("deps/");
setALLdefault2LaTex;

A1_MCXLAB_Homo_sen;
A2_MCXLAB_Homo_sen_SD0;
B1_pullPckNoi;
B2_calcEta;
C1_makeMap_simSig_SNRvZ;
D1_realVsSimSig;
D2_SNRvsNonCanNoi;

rmpath("deps/");