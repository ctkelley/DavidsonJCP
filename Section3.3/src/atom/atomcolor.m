function rgb = atomcolor(anum,scheme)
% ATOMCOLOR returns the color of the atom or atom number.
%   rgb = ATOMCOLOR(anum) returns the atom color of the given number in
%   CPK-Jmol color scheme.
%
%   rgb = ATOMCOLOR(atom) returns the atom color of the given atom in Jmol
%   color scheme.
%
%   rgb = ATOMCOLOR(a,scheme) returns the atom color of the given value in
%   the specified color scheme. Scheme is allowed to be choosen from
%   'CPK'(CPK-Jmol by default), 'CPK-Jmol', 'CPK-Rasmol'.
%
%
%   See also hex2rgb.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

if nargin == 1
    scheme = 'cpk-jmol';
elseif strcmpi(scheme,'cpk')
    scheme = 'cpk-jmol';
end

if isa(anum,'Atom')
    anum = anum.anum;
end

if strcmpi(scheme,'cpk-jmol')
    rgb = hex2rgb(cpk_jmol_hex(anum));
elseif strcmpi(scheme,'cpk-rasmol')
    rgb = hex2rgb(cpk_rasmol_hex(anum));
else
    error(['The scheme, ' scheme ', is not a supported']);
end

    function hex = cpk_jmol_hex(anum)
        switch anum
            case 1
                hex = 'FFFFFF';
            case 2
                hex = 'D9FFFF';
            case 3
                hex = 'CC80FF';
            case 4
                hex = 'C2FF00';
            case 5
                hex = 'FFB5B5';
            case 6
                hex = '909090';
            case 7
                hex = '3050F8';
            case 8
                hex = 'FF0D0D';
            case 9
                hex = '90E050';
            case 10
                hex = 'B3E3F5';
            case 11
                hex = 'AB5CF2';
            case 12
                hex = '8AFF00';
            case 13
                hex = 'BFA6A6';
            case 14
                hex = 'F0C8A0';
            case 15
                hex = 'FF8000';
            case 16
                hex = 'FFFF30';
            case 17
                hex = '1FF01F';
            case 18
                hex = '80D1E3';
            case 19
                hex = '8F40D4';
            case 20
                hex = '3DFF00';
            case 21
                hex = 'E6E6E6';
            case 22
                hex = 'BFC2C7';
            case 23
                hex = 'A6A6AB';
            case 24
                hex = '8A99C7';
            case 25
                hex = '9C7AC7';
            case 26
                hex = 'E06633';
            case 27
                hex = 'F090A0';
            case 28
                hex = '50D050';
            case 29
                hex = 'C88033';
            case 30
                hex = '7D80B0';
            case 31
                hex = 'C28F8F';
            case 32
                hex = '668F8F';
            case 33
                hex = 'BD80E3';
            case 34
                hex = 'FFA100';
            case 35
                hex = 'A62929';
            case 36
                hex = '5CB8D1';
            case 37
                hex = '702EB0';
            case 38
                hex = '00FF00';
            case 39
                hex = '94FFFF';
            case 40
                hex = '94E0E0';
            case 41
                hex = '73C2C9';
            case 42
                hex = '54B5B5';
            case 43
                hex = '3B9E9E';
            case 44
                hex = '248F8F';
            case 45
                hex = '0A7D8C';
            case 46
                hex = '006985';
            case 47
                hex = 'C0C0C0';
            case 48
                hex = 'FFD98F';
            case 49
                hex = 'A67573';
            case 50
                hex = '668080';
            case 51
                hex = '9E63B5';
            case 52
                hex = 'D47A00';
            case 53
                hex = '940094';
            case 54
                hex = '429EB0';
            case 55
                hex = '57178F';
            case 56
                hex = '00C900';
            case 57
                hex = '70D4FF';
            case 58
                hex = 'FFFFC7';
            case 59
                hex = 'D9FFC7';
            case 60
                hex = 'C7FFC7';
            case 61
                hex = 'A3FFC7';
            case 62
                hex = '8FFFC7';
            case 63
                hex = '61FFC7';
            case 64
                hex = '45FFC7';
            case 65
                hex = '30FFC7';
            case 66
                hex = '1FFFC7';
            case 67
                hex = '00FF9C';
            case 68
                hex = '00E675';
            case 69
                hex = '00D452';
            case 70
                hex = '00BF38';
            case 71
                hex = '00AB24';
            case 72
                hex = '4DC2FF';
            case 73
                hex = '4DA6FF';
            case 74
                hex = '2194D6';
            case 75
                hex = '267DAB';
            case 76
                hex = '266696';
            case 77
                hex = '175487';
            case 78
                hex = 'D0D0E0';
            case 79
                hex = 'FFD123';
            case 80
                hex = 'B8B8D0';
            case 81
                hex = 'A6544D';
            case 82
                hex = '575961';
            case 83
                hex = '9E4FB5';
            case 84
                hex = 'AB5C00';
            case 85
                hex = '754F45';
            case 86
                hex = '428296';
            case 87
                hex = '420066';
            case 88
                hex = '007D00';
            case 89
                hex = '70ABFA';
            case 90
                hex = '00BAFF';
            case 91
                hex = '00A1FF';
            case 92
                hex = '008FFF';
            case 93
                hex = '0080FF';
            case 94
                hex = '006BFF';
            case 95
                hex = '545CF2';
            case 96
                hex = '785CE3';
            case 97
                hex = '8A4FE3';
            case 98
                hex = 'A136D4';
            case 99
                hex = 'B31FD4';
            case 100
                hex = 'B31FBA';
            case 101
                hex = 'B30DA6';
            case 102
                hex = 'BD0D87';
            case 103
                hex = 'C70066';
            case 104
                hex = 'CC0059';
            case 105
                hex = 'D1004F';
            case 106
                hex = 'D90045';
            case 107
                hex = 'E00038';
            case 108
                hex = 'E6002E';
            case 109
                hex = 'EB0026';
            case 110
                hex = 'EB0026';
            case 111
                hex = 'EB0026';
            case 112
                hex = 'EB0026';
            case 113
                hex = 'EB0026';
            case 114
                hex = 'EB0026';
            case 115
                hex = 'EB0026';
            case 116
                hex = 'EB0026';
            case 117
                hex = 'EB0026';
            case 118
                hex = 'EB0026';
            otherwise
                error('The input is not a supported atom number');
        end
    end

    function hex = cpk_rasmol_hex(anum)
        switch anum
            case 1
                hex = 'FFFFFF';
            case 2
                hex = 'FFC0CB';
            case 3
                hex = 'B22121';
            case 4
                hex = 'FF1493';
            case 5
                hex = '00FF00';
            case 6
                hex = 'D3D3D3';
            case 7
                hex = '87CEE6';
            case 8
                hex = 'FF0000';
            case 9
                hex = 'DAA520';
            case 10
                hex = 'FF1493';
            case 11
                hex = '0000FF';
            case 12
                hex = '228B22';
            case 13
                hex = '696969';
            case 14
                hex = 'DAA520';
            case 15
                hex = 'FFAA00';
            case 16
                hex = 'FFFF00';
            case 17
                hex = '00FF00';
            case 18
                hex = 'FF1493';
            case 19
                hex = 'FF1493';
            case 20
                hex = '696969';
            case 21
                hex = 'FF1493';
            case 22
                hex = '696969';
            case 23
                hex = 'FF1493';
            case 24
                hex = '696969';
            case 25
                hex = '696969';
            case 26
                hex = 'FFAA00';
            case 27
                hex = 'FF1493';
            case 28
                hex = '802828';
            case 29
                hex = '802828';
            case 30
                hex = '802828';
            case 31
                hex = 'FF1493';
            case 32
                hex = 'FF1493';
            case 33
                hex = 'FF1493';
            case 34
                hex = 'FF1493';
            case 35
                hex = '802828';
            case 36
                hex = 'FF1493';
            case 37
                hex = 'FF1493';
            case 38
                hex = 'FF1493';
            case 39
                hex = 'FF1493';
            case 40
                hex = 'FF1493';
            case 41
                hex = 'FF1493';
            case 42
                hex = 'FF1493';
            case 43
                hex = 'FF1493';
            case 44
                hex = 'FF1493';
            case 45
                hex = 'FF1493';
            case 46
                hex = 'FF1493';
            case 47
                hex = '696969';
            case 48
                hex = 'FF1493';
            case 49
                hex = 'FF1493';
            case 50
                hex = 'FF1493';
            case 51
                hex = 'FF1493';
            case 52
                hex = 'FF1493';
            case 53
                hex = 'A020F0';
            case 54
                hex = 'FF1493';
            case 55
                hex = 'FF1493';
            case 56
                hex = 'FFAA00';
            case 57
                hex = 'FF1493';
            case 58
                hex = 'FF1493';
            case 59
                hex = 'FF1493';
            case 60
                hex = 'FF1493';
            case 61
                hex = 'FF1493';
            case 62
                hex = 'FF1493';
            case 63
                hex = 'FF1493';
            case 64
                hex = 'FF1493';
            case 65
                hex = 'FF1493';
            case 66
                hex = 'FF1493';
            case 67
                hex = 'FF1493';
            case 68
                hex = 'FF1493';
            case 69
                hex = 'FF1493';
            case 70
                hex = 'FF1493';
            case 71
                hex = 'FF1493';
            case 72
                hex = 'FF1493';
            case 73
                hex = 'FF1493';
            case 74
                hex = 'FF1493';
            case 75
                hex = 'FF1493';
            case 76
                hex = 'FF1493';
            case 77
                hex = 'FF1493';
            case 78
                hex = 'FF1493';
            case 79
                hex = 'DAA520';
            case 80
                hex = 'FF1493';
            case 81
                hex = 'FF1493';
            case 82
                hex = 'FF1493';
            case 83
                hex = 'FF1493';
            case 84
                hex = 'FF1493';
            case 85
                hex = 'FF1493';
            case 86
                hex = 'FF1493';
            case 87
                hex = 'FF1493';
            case 88
                hex = 'FF1493';
            case 89
                hex = 'FF1493';
            case 90
                hex = 'FF1493';
            case 91
                hex = 'FF1493';
            case 92
                hex = 'FF1493';
            case 93
                hex = 'FF1493';
            case 94
                hex = 'FF1493';
            case 95
                hex = 'FF1493';
            case 96
                hex = 'FF1493';
            case 97
                hex = 'FF1493';
            case 98
                hex = 'FF1493';
            case 99
                hex = 'FF1493';
            case 100
                hex = 'FF1493';
            case 101
                hex = 'FF1493';
            case 102
                hex = 'FF1493';
            case 103
                hex = 'FF1493';
            case 104
                hex = 'FF1493';
            case 105
                hex = 'FF1493';
            case 106
                hex = 'FF1493';
            case 107
                hex = 'FF1493';
            case 108
                hex = 'FF1493';
            case 109
                hex = 'FF1493';
            case 110
                hex = 'FF1493';
            case 111
                hex = 'FF1493';
            case 112
                hex = 'FF1493';
            case 113
                hex = 'FF1493';
            case 114
                hex = 'FF1493';
            case 115
                hex = 'FF1493';
            case 116
                hex = 'FF1493';
            case 117
                hex = 'FF1493';
            case 118
                hex = 'FF1493';
            otherwise
                error('The input is not a supported atom number');
        end
    end
end
