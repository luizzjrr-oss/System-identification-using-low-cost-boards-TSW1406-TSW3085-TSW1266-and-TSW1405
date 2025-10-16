%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculates NMSE between two data vectors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [nmse_dB]=test_nmse(sig1, sig2)
r_meas=real(sig1);
i_meas=imag(sig1);
r_model=real(sig2);
i_model=imag(sig2);
nmse_dB=10*log10( eps + sum( (r_meas-r_model).^2 + (i_meas-i_model).^2 )/sum( r_meas.^2+i_meas.^2 ) );
% disp( sprintf('NMSE (dB): %2.3f',nmse_dB) )
