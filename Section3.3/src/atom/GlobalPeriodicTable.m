function GlobalPeriodicTable()
% GLOBALPERIODICTABLE  globally define the periodic table. The table
%   currently is of size 110 by 12. Each row is a atom with the
%   corresponding atom number. And the columns represent anum, mass, venum,
%   iloc, occs, occp, occd, iso, ic, isref, ipref, idref in the
%   corresponding order. The explainations of these notations are given in
%   the following tables.
%         Notation    Explaination
%     ------------------------------------------------------------------
%             anum    Atom number
%            amass    Atom mass
%            venum    The number of valence electrons
%             iloc
%             occs
%             occp
%             occd
%              iso
%               ic
%            isref
%            ipref
%            idref
%
%   Remark: the table currently is not full, any contribution to the table
%   is very welcome.
%
%   See also sym2num, num2sym, Atom.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.
global PeriodicTable;
PeriodicTable = [
% anum       amass  venum iloc occs occp occd  iso   ic isref ipref idref
     1      1.0079      1    1    1    0    0    0    0     1     0     1   
     2      4.0026      2    1    2    0    0    0    0     1     1     1
     3      6.9390      1    1    1    0    0    0    0     1     1     1
     4      9.0122      2    1    2    0    0    0    0     1     1     1
     5     10.811       3    1    2    1    0    0    0     1     1     1
     6     12.0112      4    1    2    2    0    0    0     1     1     1
     7     14.0067      5    1    2    3    0    0    0     1     1     1
     8     15.9994      6    2    2    4    0    0    0     1     1     1
     9     18.9984      7    1    2    5    0    0    0     0     1     1
    10     20.183       8    1    2    6    0    0    0     1     1     1
    11     22.9898      1    2    1    0    0    0    0     1     1     1
    12     24.312       2    1    2    0    0    0    0     1     1     1
    13     26.9815      3    1    2    1    0    0    0     1     1     1
    14     28.086       4    1    2    2    0    0    0     1     1     1
    15     30.9738      5    1    2    3    0    0    0     1     1     1
    16     32.064       6    1    2    4    0    0    0     1     1     1
    17     35.453       7    1    2    5    0    0    0     1     1     1
    18     39.948       8    1    2    6    0    0    0     1     1     1
    19     39.0983      1    1    1    0    0    0    1     1     1     1
    20     40.08        2    1    2    0    0    0    0     1     1     1
    21     44.956       3    1    2    0    1    1    0     1     1     1
    22     47.90        4    1    2    0    2    1    0     1     1     1
    23     50.9415      5    1    2    0    3    1    0     1     1     1
    24     51.996       6    1    1    0    5    1    0     1     1     1
    25     54.938       7    1    2    0    5    1    0     1     1     1
    26     55.847       8    1    2    0    6    1    0     1     1     1
    27     58.993       9    1    2    0    7    1    0     1     1     1
    28     58.71       10    1    2    0    8    0    0     1     1     1
    29     63.546      11    1    1    0   10    1    0     1     1     1
    30     65.38       12    1    2    0   10    1    0     0     1     1
    31     69.723       3    1    2    1    0    1    1     1     1     1
    32     72.64        4    2    2    2    0    1    1     1     1     1
    33     74.921       5    3    2    3    0    1    1     1     1     1
    34     78.96        6    3    2    4    0    1    0     1     1     1
    35     79.904       7    3    2    5    0    1    0     1     1     1
    36     83.798       8    3    2    6    0    1    0     1     1     1
    37     85.467     NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    38     87.62      NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    39     88.905     NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    40     91.224     NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    41     92.906     NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    42     95.96      NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    43     97.907     NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    44    101.07      NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    45    102.9       NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    46    106.42      NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    47    107.86      NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    48    112.41       12    1    2    0   10    1    1     1     1     1
    49    114.81        3    1    2    1    0    1    1     1     1     1
    50    118.71      NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    51    121.76      NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    52    127.6         6    1    2    4    0    1    1     1     1     1
    53    126.9       NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    54    131.29      NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    55    132.90      NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    56    137.32      NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    57    NaN         NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    58    NaN          10    1    2    0    1    1    1     1     1     1
    59    NaN         NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    60    NaN         NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    61    NaN         NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    62    NaN         NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    63    NaN         NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    64    NaN         NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    65    NaN         NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    66    NaN         NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    67    NaN         NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    68    NaN         NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    69    NaN         NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    70    NaN         NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    71    NaN         NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    72    178.49      NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    73    180.94      NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    74    183.84      NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    75    186.20      NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    76    190.23      NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    77    192.21      NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    78    195.078      10    1    1    0    9    1    0     1     1     1
    79    196.96      NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    80    200.59      NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    81    204.38      NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    82    207.2       NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    83    208.98        5    2    2    3    0    1    1     1     1     1
    84    208.98      NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    85    208.98      NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    86    222.01      NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    87    NaN         NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    88    NaN         NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    89    NaN         NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    90    NaN         NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    91    NaN         NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    92    NaN         NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    93    NaN         NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    94    NaN         NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    95    NaN         NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    96    NaN         NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    97    NaN         NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    98    NaN         NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
    99    NaN         NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
   100    NaN         NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
   101    NaN         NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
   102    NaN         NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
   103    NaN         NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
   104    NaN         NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
   105    NaN         NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
   106    NaN         NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
   107    NaN         NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
   108    NaN         NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
   109    NaN         NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
   110    NaN         NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
   111    NaN         NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
   112    NaN         NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
   113    NaN         NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
   114    NaN         NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
   115    NaN         NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
   116    NaN         NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
   117    NaN         NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
   118    NaN         NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN   NaN   NaN
];

end