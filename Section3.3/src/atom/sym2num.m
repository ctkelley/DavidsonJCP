function anum = sym2num(symbol)
% SYM2NUM  returns the atom number and mass of the given symbol.
%   anum = SYM2NUM(symbol) returns the atom number of the given symbol.
%
%   Example:
%
%       anum = sym2num('Li')
%
%       anum =
% 
%            3
% 
%   See also num2sym, atominfo, GlobalPeriodicTable.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.
symbol = strtrim(symbol);
switch symbol
    case {'H','h'}
        anum = 1;
    case {'He','he'}
        anum = 2;
    case {'Li','li'}
        anum = 3;
    case {'Be','be'}
        anum = 4;
    case {'B','b'}
        anum = 5;
    case {'C','c'}
        anum = 6;
    case {'N','n'}
        anum = 7;
    case {'O','o'}
        anum = 8;
    case {'F','f'}
        anum = 9;
    case {'Ne','ne'}
        anum = 10;
    case {'Na','na'}
        anum = 11;
    case {'Mg','mg'}
        anum = 12;
    case {'Al','al'}
        anum = 13;
    case {'Si','si'}
        anum = 14;
    case {'P','p'}
        anum = 15;
    case {'S','s'}
        anum = 16;
    case {'Cl','cl'}
        anum = 17;
    case {'Ar','ar'}
        anum = 18;
    case {'K','k'}
        anum = 19;
    case {'Ca','ca'}
        anum = 20;
    case {'Sc','sc'}
        anum = 21;
    case {'Ti','ti'}
        anum = 22;
    case {'V','v'}
        anum = 23;
    case {'Cr','cr'}
        anum = 24;
    case {'Mn','mn'}
        anum = 25;
    case {'Fe','fe'}
        anum = 26;
    case {'Co','co'}
        anum = 27;
    case {'Ni','ni'}
        anum = 28;
    case {'Cu','cu'}
        anum = 29;
    case {'Zn','zn'}
        anum = 30;
    case {'Ga','ga'}
        anum = 31;
    case {'Ge','ge'}
        anum = 32;
    case {'As','as'}
        anum = 33;
    case {'Se','se'}
        anum = 34;
    case {'Br','br'}
        anum = 35;
    case {'Kr','kr'}
        anum = 36;
    case {'Rb','rb'}
        anum = 37;
    case {'Sr','sr'}
        anum = 38;
    case {'Y','y'}
        anum = 39;
    case {'Zr','zr'}
        anum = 40;
    case {'Nb','nb'}
        anum = 41;
    case {'Mo','mo'}
        anum = 42;
    case {'Tc','tc'}
        anum = 43;
    case {'Ru','ru'}
        anum = 44;
    case {'Rh','rh'}
        anum = 45;
    case {'Pd','pd'}
        anum = 46;
    case {'Ag','ag'}
        anum = 47;
    case {'Cd','cd'}
        anum = 48;
    case {'In','in'}
        anum = 49;
    case {'Sn','sn'}
        anum = 50;
    case {'Sb','sb'}
        anum = 51;
    case {'Te','te'}
        anum = 52;
    case {'I','i'}
        anum = 53;
    case {'Xe','xe'}
        anum = 54;
    case {'Cs','cs'}
        anum = 55;
    case {'Ba','ba'}
        anum = 56;
    case {'La','la'}
        anum = 57;
    case {'Ce','ce'}
        anum = 58;
    case {'Pr','pr'}
        anum = 59;
    case {'Nd','nd'}
        anum = 60;
    case {'Pm','pm'}
        anum = 61;
    case {'Sm','sm'}
        anum = 62;
    case {'Eu','eu'}
        anum = 63;
    case {'Gd','gd'}
        anum = 64;
    case {'Tb','tb'}
        anum = 65;
    case {'Dy','dy'}
        anum = 66;
    case {'Ho','ho'}
        anum = 67;
    case {'Er','er'}
        anum = 68;
    case {'Tm','tm'}
        anum = 69;
    case {'Yb','yb'}
        anum = 70;
    case {'Lu','lu'}
        anum = 71;
    case {'Hf','hf'}
        anum = 72;
    case {'Ta','ta'}
        anum = 73;
    case {'W','w'}
        anum = 74;
    case {'Re','re'}
        anum = 75;
    case {'Os','os'}
        anum = 76;
    case {'Ir','ir'}
        anum = 77;
    case {'Pt','pt'}
        anum = 78;
    case {'Au','au'}
        anum = 79;
    case {'Hg','hg'}
        anum = 80;
    case {'Tl','tl'}
        anum = 81;
    case {'Pb','pb'}
        anum = 82;
    case {'Bi','bi'}
        anum = 83;
    case {'Po','po'}
        anum = 84;
    case {'At','at'}
        anum = 85;
    case {'Rn','rn'}
        anum = 86;
    case {'Fr','fr'}
        anum = 87;
    case {'Ra','ra'}
        anum = 88;
    case {'Ac','ac'}
        anum = 89;
    case {'Th','th'}
        anum = 90;
    case {'Pa','pa'}
        anum = 91;
    case {'U','u'}
        anum = 92;
    case {'Np','np'}
        anum = 93;
    case {'Pu','pu'}
        anum = 94;
    case {'Am','am'}
        anum = 95;
    case {'Cm','cm'}
        anum = 96;
    case {'Bk','bk'}
        anum = 97;
    case {'Cf','cf'}
        anum = 98;
    case {'Es','es'}
        anum = 99;
    case {'Fm','fm'}
        anum = 100;
    case {'Md','md'}
        anum = 101;
    case {'No','no'}
        anum = 102;
    case {'Lr','lr'}
        anum = 103;
    case {'Rf','rf'}
        anum = 104;
    case {'Db','db'}
        anum = 105;
    case {'Sg','sg'}
        anum = 106;
    case {'Bh','bh'}
        anum = 107;
    case {'Hs','hs'}
        anum = 108;
    case {'Mt','mt'}
        anum = 109;
    case {'Ds','ds'}
        anum = 110;
    case {'Rg','rg'}
        anum = 111;
    case {'Cn','cn'}
        anum = 112;
    case {'Uut','uut'}
        anum = 113;
    case {'Fl','fl'}
        anum = 114;
    case {'Uup','uup'}
        anum = 115;
    case {'Lv','lv'}
        anum = 116;
    case {'Uus','uus'}
        anum = 117;
    case {'Uuo','uuo'}
        anum = 118;
    otherwise
        error('The input is not a supported symbol.');
end

end
