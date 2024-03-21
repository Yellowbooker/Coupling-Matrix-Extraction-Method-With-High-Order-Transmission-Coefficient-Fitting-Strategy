% =========================================================================

% Import the s-parameters from the .s2p data file

% Output:
% freq: real sampled frequency points
% S11, S12, S22: sampled S-parameters

% Data: 2016/05/01

% =========================================================================


function [freq,S11,S12,S22]=Reads2p()
FILE_NAME=input('File name (.s2p): ','s');
S= sparameters(FILE_NAME);
freq = S.Frequencies;
S11 = rfparam(S,1,1);
S12 = rfparam(S,1,2);
S22 = rfparam(S,2,2);
end
