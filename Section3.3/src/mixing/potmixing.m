function [vin, dfmat, dvmat, cdfmat] = potmixing(mol, vin, vout, iterscf, ...
                                                 mixtype, betamix, dfmat, ...
                                                 dvmat, cdfmat, mixdim, brank);
%
% Purpose: entry point for selecting different potential mixing algorithms.
% 
% Arguements:
%  input: 
%     vin (3-D array)  --- input potential
%    vout (3-D array)  --- output potential produced from eigenvalue 
%                          calculation
%    iterscf (integer) --- current SCF iteration number
%    betamix (scalar)  --- mixing parameter between 0 and 1
%    df   (matrix)     --- work array that stores the changes in function
%                          evaluations.  The dimension of df is n1*n2*n3 by mixdim,
%                          where [n1,n2,n3]=size(vin) if Anderson or Pulay mixing is used.
%                          If Broyden mixing is used, the second dimension is
%                          mixdim*brank.
%    dv   (matrix)     --- work array that stores the changes in input 
%                          potentials. Its dimension is the same as that of df.
%    mixdim(integer)   --- the number of columns in df and dv. This is the 
%                          the mixing memory. Only information from the
%                          previous mixdim iterations are used.
%    brank (integer)   --- The number of secand equations satisfied by
%                          a Broyden update. This is also the rank of
%                          one Broyden update.
%  output:
%    vnew (3-D array) --- new potential to be used for the next SCF iteration
%    df   (matrix)    --- updated work array that keeps the changes in 
%                         function evaluations
%    dv   (matrix)    --- updated work array that keeps the changes in 
%                           input potentials

if (betamix == 0)
   betamix = 0.8;
end

switch mixtype
    case {'pulay','Pulay'}
      [vin,dfmat,dvmat] = pulaymix(vin, vout, betamix, dfmat, dvmat, iterscf, mixdim);
    case {'anderson','Anderson'}
      [vin,dfmat,dvmat] = andersonmix(vin,vout,betamix,dfmat,dvmat,...,
                                      iterscf, mixdim);
    case {'broyden1','Broyden1'}
      % the type I (or what is known as the good) Broyden 
      [vin,dfmat,dvmat,cdfmat] = broydenmix1(vin,vout,betamix,dfmat,dvmat,...
                                             cdfmat,iterscf,mixdim,brank);
    case {'broyden2','Broyden2'}
      % the type II (or what is known as the bad) Broyden 
      [vin,dfmat,dvmat] = broydenmix(vin,vout,betamix,dfmat,dvmat,iterscf,...
                                     mixdim, brank);
    case {'broyden','Broyden'}
      [vin,dfmat,dvmat] = broydenmix(vin,vout,betamix,dfmat,dvmat,iterscf,...
                                     mixdim, brank);
    case {'johnson','Johnson'}
      [vin,dfmat,dvmat] = johnsonmix(vin,vout,dfmat,dvmat,iterscf,...
                                     mixdim, brank);
    case {'kerker','Kerker'}
       vin = kerkmix(mol,vin,vout);
    case {'pulay+kerker','Pulay+Kerker'}
      [vin,dfmat,dvmat] = pulaymix(vin, vout, dfmat, dvmat, iterscf, mixdim);
       vin = kerkmix(mol,vin,vout);
    case {'simple','Simple'}
       vin = simplemix(vin, vout, betamix);
    case {'off','OFF'}
       % no mixing at all 
       vin = vout;
    otherwise
       fprintf('Invalid mixing type %s\n', mixtype);
       fprintf('Using Anderson mixing by default with beta=0.8 ...\n');
       betamix = 0.8;
       [vin,dfmat,dvmat] = andersonmix(vin,vout,betamix,dfmat,dvmat,...,
                                       iterscf, mixdim);
end  
