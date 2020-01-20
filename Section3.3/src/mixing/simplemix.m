function vnew = simplemix(vin, vout, betamix)
%
% usage: [vnew,dvmat,vmat] = simplemix(vin, vout, betamix);
%
% purpose: simple mixing for accelerating the SCF iteration
%
% arguements:
%  input: 
%     vin (3-D array) --- input potential
%    vout (3-D array) --- output potential produced from eigenvalue 
%                         calculation
%    betamix (real scalar)  --- mixing parameter
%  output:
%    vnew (3-D array) --- new potential to be used for the next SCF iteration
%  

vnew = vin + betamix*(vout-vin);
