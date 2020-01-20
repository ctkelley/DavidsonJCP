function symbol = num2sym(anum)
% NUM2SYM  returns the atom symbol of the given atom number.
%   symbol = NUM2SYM(anum) returns the atom symbol of the given number.
%
%   Example:
%
%       symbol = num2sym(3)
%
%       symbol =
% 
%       Li
%
%
%   See also sym2num.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.
switch anum
    case 1
        symbol = 'H';
    case 2
        symbol = 'He';
    case 3
        symbol = 'Li';
    case 4
        symbol = 'Be';
    case 5
        symbol = 'B';
    case 6
        symbol = 'C';
    case 7
        symbol = 'N';
    case 8
        symbol = 'O';
    case 9
        symbol = 'F';
    case 10
        symbol = 'Ne';
    case 11
        symbol = 'Na';
    case 12
        symbol = 'Mg';
    case 13
        symbol = 'Al';
    case 14
        symbol = 'Si';
    case 15
        symbol = 'P';
    case 16
        symbol = 'S';
    case 17
        symbol = 'Cl';
    case 18
        symbol = 'Ar';
    case 19
        symbol = 'K';
    case 20
        symbol = 'Ca';
    case 21
        symbol = 'Sc';
    case 22
        symbol = 'Ti';
    case 23
        symbol = 'V';
    case 24
        symbol = 'Cr';
    case 25
        symbol = 'Mn';
    case 26
        symbol = 'Fe';
    case 27
        symbol = 'Co';
    case 28
        symbol = 'Ni';
    case 29
        symbol = 'Cu';
    case 30
        symbol = 'Zn';
    case 31
        symbol = 'Ga';
    case 32
        symbol = 'Ge';
    case 33
        symbol = 'As';
    case 34
        symbol = 'Se';
    case 35
        symbol = 'Br';
    case 36
        symbol = 'Kr';
    case 37
        symbol = 'Rb';
    case 38
        symbol = 'Sr';
    case 39
        symbol = 'Y';
    case 40
        symbol = 'Zr';
    case 41
        symbol = 'Nb';
    case 42
        symbol = 'Mo';
    case 43
        symbol = 'Tc';
    case 44
        symbol = 'Ru';
    case 45
        symbol = 'Rh';
    case 46
        symbol = 'Pd';
    case 47
        symbol = 'Ag';
    case 48
        symbol = 'Cd';
    case 49
        symbol = 'In';
    case 50
        symbol = 'Sn';
    case 51
        symbol = 'Sb';
    case 52
        symbol = 'Te';
    case 53
        symbol = 'I';
    case 54
        symbol = 'Xe';
    case 55
        symbol = 'Cs';
    case 56
        symbol = 'Ba';
    case 57
        symbol = 'La';
    case 58
        symbol = 'Ce';
    case 59
        symbol = 'Pr';
    case 60
        symbol = 'Nd';
    case 61
        symbol = 'Pm';
    case 62
        symbol = 'Sm';
    case 63
        symbol = 'Eu';
    case 64
        symbol = 'Gd';
    case 65
        symbol = 'Tb';
    case 66
        symbol = 'Dy';
    case 67
        symbol = 'Ho';
    case 68
        symbol = 'Er';
    case 69
        symbol = 'Tm';
    case 70
        symbol = 'Yb';
    case 71
        symbol = 'Lu';
    case 72
        symbol = 'Hf';
    case 73
        symbol = 'Ta';
    case 74
        symbol = 'W';
    case 75
        symbol = 'Re';
    case 76
        symbol = 'Os';
    case 77
        symbol = 'Ir';
    case 78
        symbol = 'Pt';
    case 79
        symbol = 'Au';
    case 80
        symbol = 'Hg';
    case 81
        symbol = 'Tl';
    case 82
        symbol = 'Pb';
    case 83
        symbol = 'Bi';
    case 84
        symbol = 'Po';
    case 85
        symbol = 'At';
    case 86
        symbol = 'Rn';
    case 87
        symbol = 'Fr';
    case 88
        symbol = 'Ra';
    case 89
        symbol = 'Ac';
    case 90
        symbol = 'Th';
    case 91
        symbol = 'Pa';
    case 92
        symbol = 'U';
    case 93
        symbol = 'Np';
    case 94
        symbol = 'Pu';
    case 95
        symbol = 'Am';
    case 96
        symbol = 'Cm';
    case 97
        symbol = 'Bk';
    case 98
        symbol = 'Cf';
    case 99
        symbol = 'Es';
    case 100
        symbol = 'Fm';
    case 101
        symbol = 'Md';
    case 102
        symbol = 'No';
    case 103
        symbol = 'Lr';
    case 104
        symbol = 'Rf';
    case 105
        symbol = 'Db';
    case 106
        symbol = 'Sg';
    case 107
        symbol = 'Bh';
    case 108
        symbol = 'Hs';
    case 109
        symbol = 'Mt';
    case 110
        symbol = 'Ds';
    case 111
        symbol = 'Rg';
    case 112
        symbol = 'Cn';
    case 113
        symbol = 'Uut';
    case 114
        symbol = 'Fl';
    case 115
        symbol = 'Uup';
    case 116
        symbol = 'Lv';
    case 117
        symbol = 'Uus';
    case 118
        symbol = 'Uuo';
    otherwise
        error('The input is not a supported atom number');
end
