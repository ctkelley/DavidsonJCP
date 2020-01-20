function rad = atomradius(anum,scheme)
% ATOMRADIUS returns the radius of the atom or atom number.
%   rad = ATOMRADIUS(anum) returns the atomic radius of the given number.
%
%   rad = ATOMRADIUS(atom) returns the atomic radius of the given atom.
%
%   rad = ATOMRADIUS(atom,scheme) returns the atomic radius of the given
%   atom in the given scheme. Currently, we only support 'covalent' and
%   'van_der_waals'.
%
%   See also rad2rgb.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

if nargin == 1
    scheme = 'covalent';
end

if isa(anum,'Atom')
    anum = anum.anum;
end

if strcmpi(scheme,'covalent')
    rad = covalent_rad(anum);
elseif strcmpi(scheme,'van_der_waals')
    rad = van_der_waals_rad(anum);
else
    error(['The scheme, ' scheme ', is not a supported']);
end

    function rad = covalent_rad(anum)
        % atom 1 to 100 is from XCrysDen
        switch anum
            case 1
                rad = 38;
            case 2
                rad = 38;
            case 3
                rad = 123;
            case 4
                rad = 89;
            case 5
                rad = 91;
            case 6
                rad = 77;
            case 7
                rad = 75;
            case 8
                rad = 73;
            case 9
                rad = 71;
            case 10
                rad = 71;
            case 11
                rad = 160;
            case 12
                rad = 140;
            case 13
                rad = 125;
            case 14
                rad = 111;
            case 15
                rad = 100;
            case 16
                rad = 104;
            case 17
                rad = 99;
            case 18
                rad = 98;
            case 19
                rad = 213;
            case 20
                rad = 174;
            case 21
                rad = 160;
            case 22
                rad = 140;
            case 23
                rad = 135;
            case 24
                rad = 140;
            case 25
                rad = 140;
            case 26
                rad = 140;
            case 27
                rad = 135;
            case 28
                rad = 135;
            case 29
                rad = 135;
            case 30
                rad = 135;
            case 31
                rad = 130;
            case 32
                rad = 125;
            case 33
                rad = 115;
            case 34
                rad = 115;
            case 35
                rad = 114;
            case 36
                rad = 112;
            case 37
                rad = 220;
            case 38
                rad = 200;
            case 39
                rad = 185;
            case 40
                rad = 155;
            case 41
                rad = 145;
            case 42
                rad = 145;
            case 43
                rad = 135;
            case 44
                rad = 130;
            case 45
                rad = 135;
            case 46
                rad = 140;
            case 47
                rad = 160;
            case 48
                rad = 155;
            case 49
                rad = 155;
            case 50
                rad = 141;
            case 51
                rad = 145;
            case 52
                rad = 140;
            case 53
                rad = 140;
            case 54
                rad = 131;
            case 55
                rad = 260;
            case 56
                rad = 200;
            case 57
                rad = 175;
            case 58
                rad = 155;
            case 59
                rad = 155;
            case 60
                rad = 155;
            case 61
                rad = 155;
            case 62
                rad = 155;
            case 63
                rad = 155;
            case 64
                rad = 155;
            case 65
                rad = 155;
            case 66
                rad = 155;
            case 67
                rad = 155;
            case 68
                rad = 155;
            case 69
                rad = 155;
            case 70
                rad = 155;
            case 71
                rad = 155;
            case 72
                rad = 155;
            case 73
                rad = 145;
            case 74
                rad = 135;
            case 75
                rad = 135;
            case 76
                rad = 130;
            case 77
                rad = 135;
            case 78
                rad = 135;
            case 79
                rad = 135;
            case 80
                rad = 150;
            case 81
                rad = 190;
            case 82
                rad = 180;
            case 83
                rad = 160;
            case 84
                rad = 155;
            case 85
                rad = 155;
            case 86
                rad = 155;
            case 87
                rad = 280;
            case 88
                rad = 144;
            case 89
                rad = 195;
            case 90
                rad = 155;
            case 91
                rad = 155;
            case 92
                rad = 155;
            case 93
                rad = 155;
            case 94
                rad = 155;
            case 95
                rad = 155;
            case 96
                rad = 155;
            case 97
                rad = 155;
            case 98
                rad = 155;
            case 99
                rad = 155;
            case 100
                rad = 155;
            case 101
                rad = 173;
            case 102
                rad = 176;
            case 103
                rad = 161;
            case 104
                rad = 157;
            case 105
                rad = 149;
            case 106
                rad = 143;
            case 107
                rad = 141;
            case 108
                rad = 134;
            case 109
                rad = 129;
            case 110
                rad = 128;
            case 111
                rad = 121;
            case 112
                rad = 122;
            case 113
                rad = 136;
            case 114
                rad = 143;
            case 115
                rad = 162;
            case 116
                rad = 175;
            case 117
                rad = 165;
            case 118
                rad = 157;
            otherwise
                error('The input is not a supported atom number');
        end
        % pm to bohr
        rad = rad * 0.01889726;
    end

    function rad = van_der_waals_rad(anum)
        % atom 1 to 100 is from XCrysDen
        switch anum
            case 1
                rad = 117;
            case 2
                rad = 140;
            case 3
                rad = 180;
            case 4
                rad = 210;
            case 5
                rad = 210;
            case 6
                rad = 170;
            case 7
                rad = 158;
            case 8
                rad = 152;
            case 9
                rad = 147;
            case 10
                rad = 160;
            case 11
                rad = 230;
            case 12
                rad = 170;
            case 13
                rad = 210;
            case 14
                rad = 170;
            case 15
                rad = 180;
            case 16
                rad = 180;
            case 17
                rad = 178;
            case 18
                rad = 190;
            case 19
                rad = 280;
            case 20
                rad = 210;
            case 21
                rad = 210;
            case 22
                rad = 210;
            case 23
                rad = 210;
            case 24
                rad = 210;
            case 25
                rad = 210;
            case 26
                rad = 210;
            case 27
                rad = 210;
            case 28
                rad = 160;
            case 29
                rad = 140;
            case 30
                rad = 140;
            case 31
                rad = 190;
            case 32
                rad = 210;
            case 33
                rad = 185;
            case 34
                rad = 190;
            case 35
                rad = 185;
            case 36
                rad = 200;
            case 37
                rad = 210;
            case 38
                rad = 210;
            case 39
                rad = 210;
            case 40
                rad = 210;
            case 41
                rad = 210;
            case 42
                rad = 210;
            case 43
                rad = 210;
            case 44
                rad = 210;
            case 45
                rad = 210;
            case 46
                rad = 160;
            case 47
                rad = 170;
            case 48
                rad = 160;
            case 49
                rad = 190;
            case 50
                rad = 100;
            case 51
                rad = 210;
            case 52
                rad = 206;
            case 53
                rad = 196;
            case 54
                rad = 216;
            case 55
                rad = 210;
            case 56
                rad = 210;
            case 57
                rad = 210;
            case 58
                rad = 210;
            case 59
                rad = 210;
            case 60
                rad = 210;
            case 61
                rad = 210;
            case 62
                rad = 210;
            case 63
                rad = 210;
            case 64
                rad = 210;
            case 65
                rad = 210;
            case 66
                rad = 210;
            case 67
                rad = 210;
            case 68
                rad = 210;
            case 69
                rad = 210;
            case 70
                rad = 210;
            case 71
                rad = 210;
            case 72
                rad = 210;
            case 73
                rad = 210;
            case 74
                rad = 210;
            case 75
                rad = 210;
            case 76
                rad = 210;
            case 77
                rad = 210;
            case 78
                rad = 175;
            case 79
                rad = 170;
            case 80
                rad = 150;
            case 81
                rad = 200;
            case 82
                rad = 210;
            case 83
                rad = 210;
            case 84
                rad = 210;
            case 85
                rad = 210;
            case 86
                rad = 210;
            case 87
                rad = 350;
            case 88
                rad = 350;
            case 89
                rad = 350;
            case 90
                rad = 210;
            case 91
                rad = 210;
            case 92
                rad = 210;
            case 93
                rad = 210;
            case 94
                rad = 210;
            case 95
                rad = 210;
            case 96
                rad = 210;
            case 97
                rad = 210;
            case 98
                rad = 210;
            case 99
                rad = 210;
            case 100
                rad = 210;
            case 111
                rad = nan;
            case 112
                rad = nan;
            case 113
                rad = nan;
            case 114
                rad = nan;
            case 115
                rad = nan;
            case 116
                rad = nan;
            case 117
                rad = nan;
            case 118
                rad = nan;
            otherwise
                error('The input is not a supported atom number');
        end
        % pm to bohr
        rad = rad * 0.01889726;
    end

end