# LA_d11B
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dB_DR_v6: Data reduction script for processing laser ablation boron
% isotope data. The methodology is outlined in:
%   Evans, Gerdes, Coenen, Marschall, & Mueller (2021) Accurate correction 
%   for the matrix interference on laser ablation MC-ICPMS boron isotope
%   measurements in CaCO3 and silicate matrices. JAAS 36:1607.
%
% If this script is useful for your research, please cite the above paper.
%
% Derive carbonate standardised d11B values. Briefly:
% 1) Remove outliers and baseline subtract sample/standard analyses
% 2) Identify the location of the NIST analyses to be used for mass bias
%   and drift correction, and apply these corrections
% 3) Calculate 'raw' (NIST-standardised) d11B values
% 4) Calculate Dd11B for three carbonate standards, and derive the
%   power relationship between ~B/Ca (11/10.035) and Dd11B
% 5) Correct sample data using the derived relationship between Dd11B and
%   the 11/10.035 ratio
%
% This script is designed to be used in conjuction with a log file from the
% laser ablation system; functionality is greatly limited without.
