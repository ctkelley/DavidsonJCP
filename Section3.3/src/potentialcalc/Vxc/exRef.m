function exRef
% EXREF reference for exchange correlation.
%    Here attached a exchange correlation name convention and reference
%    list. 
%
%     short name      complete name       Short description
%      "pz"           "sla+pz"            Perdew-Zunger LDA
%      "bp"           "b88+p86"           Becke-Perdew grad.corr.
%      "pw91"         "sla+pw+ggx+ggc"    PW91 (aka GGA)
%      "blyp"         "sla+b88+lyp+blyp"  BLYP
%      "pbe"          "sla+pw+pbx+pbc"    PBE
%      "revpbe"       "sla+pw+rpb+pbc"    revPBE (Zhang-Yang)
%      "pw86pbe"      "sla+pw+pw86+pbc"   PW86 exchange + PBE correlation
%      "b86bpbe"      "sla+pw+b86b+pbc"   B86b exchange + PBE correlation
%      "pbesol"       "sla+pw+psx+psc"    PBEsol
%      "q2d"          "sla+pw+q2dx+q2dc"  PBEQ2D
%      "hcth"         "nox+noc+hcth+hcth" HCTH/120
%      "olyp"         "nox+lyp+optx+blyp" OLYP
%      "wc"           "sla+pw+wcx+pbc"    Wu-Cohen
%      "sogga"        "sla+pw+sox+pbec"   SOGGA
%      "optbk88"      "sla+pw+obk8+p86"   optB88
%      "optb86b"      "sla+pw+ob86+p86"   optB86
%      "ev93"         "sla+pw+evx+nogc"   Engel-Vosko
%      "tpss"         "sla+pw+tpss+tpss"  TPSS Meta-GGA
%      "m06l"         "nox+noc+m6lx+m6lc" M06L Meta-GGA
%      "tb09"         "sla+pw+tb09+tb09"  TB09 Meta-GGA
%      "pbe0"         "pb0x+pw+pb0x+pbc"  PBE0
%      "hse"          "sla+pw+hse+pbc"    Heyd-Scuseria-Ernzerhof (HSE 06)
%      "b3lyp"        "b3lp+b3lp          B3LYP
%                      +b3lp+b3lp"
%      "b3lypv1r"     "b3lp+b3lpv1r       B3LYP-VWN1-RPA
%                      +b3lp+b3lp"
%      "x3lyp"        "x3lp+x3lp          X3LYP
%                      +x3lp+x3lp"
%      "vwn-rpa"      "sla+vwn-rpa"       VWN LDA using vwn1-rpa parametriz
%      "gaupbe"       "sla+pw+gaup+pbc"   Gau-PBE (also "gaup")
%      "vdw-df"       "sla+pw+rpb +vdw1"  vdW-DF1
%      "vdw-df2"      "sla+pw+rw86+vdw2"  vdW-DF2
%      "vdw-df-x"     "sla+pw+????+vdwx"  vdW-DF-x, reserved Thonhauser, not implemented
%      "vdw-df-y"     "sla+pw+????+vdwy"  vdW-DF-y, reserved Thonhauser, not implemented
%      "vdw-df-z"     "sla+pw+????+vdwz"  vdW-DF-z, reserved Thonhauser, not implemented
%      "vdw-df-c09"   "sla+pw+c09x+vdw1"  vdW-DF-C09
%      "vdw-df2-c09"  "sla+pw+c09x+vdw2"  vdW-DF2-C09
%      "vdw-df-cx"    "sla+pw+cx13+vdW1"  vdW-DF-cx
%      "vdw-df-obk8"  "sla+pw+obk8+vdw1"  vdW-DF-obk8 (optB88-vdW)
%      "vdw-df-ob86"  "sla+pw+ob86+vdw1"  vdW-DF-ob86 (optB86b-vdW)
%      "vdw-df2-b86r" "sla+pw+b86r+vdw2"  vdW-DF2-B86R (rev-vdw-df2)
%      "rvv10"        "sla+pw+rw86        rVV10
%                      +pbc+vv10"
%
%    Exchange name list.
%       short name   complete name
%       "sla"        Slater (alpha=2/3)
%       "sl1"        Slater (alpha=1.0)
%       "rxc"        Relativistic Slater
%       "oep"        Optimized Effective Potential
%       "hf"         Hartree-Fock
%       "pb0x"       PBE0 (Slater*0.75+HF*0.25)
%       "b3lp"       B3LYP(Slater*0.80+HF*0.20)
%       "kzk"        Finite-size corrections
%       "x3lp"       X3LYP(Slater*0.782+HF*0.218)
%
%    Correlation name list.
%       short name   complete name
%       "pz"         Perdew-Zunger
%       "vwn"        Vosko-Wilk-Nusair
%       "lyp"        Lee-Yang-Parr
%       "pw"         Perdew-Wang
%       "wig"        Wigner
%       "hl"         Hedin-Lunqvist
%       "obz"        Ortiz-Ballone form for PZ
%       "obw"        Ortiz-Ballone form for PW
%       "gl"         Gunnarson-Lunqvist
%       "kzk"        Finite-size corrections
%       "vwn-rpa"    Vosko-Wilk-Nusair, alt param
%       "b3lp"       B3LYP (0.19*vwn+0.81*lyp)
%       "b3lpv1r"    B3LYP-VWN-1-RPA(0.19*vwn_rpa+0.81*lyp)
%       "x3lp"       X3LYP (0.129*vwn_rpa+0.871*lyp)
%
%    Gradient Correction on Exchange name list.
%       short name   complete name
%       "b88"        Becke88 (beta=0.0042)
%       "ggx"        Perdew-Wang 91
%       "pbx"        Perdew-Burke-Ernzenhof exch
%       "rpb"        revised PBE by Zhang-Yang
%       "hcth"       Cambridge exch, Handy et al
%       "optx"       Handy's exchange functional
%       "pb0x"       PBE0 (PBE exchange*0.75)
%       "b3lp"       B3LYP (Becke88*0.72)
%       "psx"        PBEsol exchange
%       "wcx"        Wu-Cohen
%       "hse"        HSE screened exchange
%       "rw86"       revised PW86
%       "pbe"        same as PBX, back-comp
%       "c09x"       Cooper 09
%       "sox"        sogga
%       "q2dx"       Q2D exchange grad corr
%       "gaup"       Gau-PBE hybrid exchange
%       "pw86"       Perdew-Wang (1986) exchange
%       "b86b"       Becke (1986) exchange
%       "obk8"       optB88  exchange
%       "ob86"       optB86b exchange
%       "evx"        Engel-Vosko exchange
%       "b86r"       revised Becke (b86b)
%       "cx13"       consistent exchange
%       "x3lp"       X3LYP (Becke88*0.542+Perdew-Wang91*0.167)
%
%    Gradient Correction on Correlation name list.
%       short name   complete name
%       "p86"        Perdew86
%       "ggc"        Perdew-Wang 91 corr
%       "blyp"       Lee-Yang-Parr
%       "pbc"        Perdew-Burke-Ernzenhof corr
%       "hcth"       Cambridge corr, Handy et al
%       "b3lp"       B3LYP (Lee-Yang-Parr*0.81)
%       "psc"        PBEsol corr
%       "pbe"        same as PBX, back-comp
%       "q2dc"       Q2D correlation grad corr
%       "x3lp"       X3LYP (Lee-Yang-Parr*0.871)
%
%    Meta-GGA functionals name list.
%       short name   complete name
%       "tpss"       TPSS Meta-GGA
%       "m6lx"       M06L Meta-GGA
%       "tb09"       TB09 Meta-GGA
%
%    Van der Waals functionals (nonlocal term only) name list.
%       short name   complete name
%       "vdw1"       vdW-DF1
%       "vdw2"       vdW-DF2
%       "vv10"       rVV10
%       "vdwx"       vdW-DF-x
%       "vdwy"       vdW-DF-y
%       "vdwz"       vdW-DF-z
%
% References:
%    pz      J.P.Perdew and A.Zunger, PRB 23, 5048 (1981)
%    vwn     S.H.Vosko, L.Wilk, M.Nusair, Can.J.Phys. 58,1200(1980)
%    vwn1-rpa S.H.Vosko, L.Wilk, M.Nusair, Can.J.Phys. 58,1200(1980)
%    wig     E.P.Wigner, Trans. Faraday Soc. 34, 67 (1938)
%    hl      L.Hedin and B.I.Lundqvist, J. Phys. C4, 2064 (1971)
%    gl      O.Gunnarsson and B.I.Lundqvist, PRB 13, 4274 (1976)
%    pw      J.P.Perdew and Y.Wang, PRB 45, 13244 (1992)
%    obpz    G.Ortiz and P.Ballone, PRB 50, 1391 (1994)
%    obpw    as above
%    b88     A.D.Becke, PRA 38, 3098 (1988)
%    p86     J.P.Perdew, PRB 33, 8822 (1986)
%    pw86    J.P.Perdew, PRB 33, 8800 (1986)
%    b86b    A.D.Becke, J.Chem.Phys. 85, 7184 (1986)
%    ob86    Klimes, Bowler, Michaelides, PRB 83, 195131 (2011)
%    b86r    I. Hamada, Phys. Rev. B 89, 121103(R) (2014)
%    pbe     J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996)
%    pw91    J.P.Perdew and Y. Wang, PRB 46, 6671 (1992)
%    blyp    C.Lee, W.Yang, R.G.Parr, PRB 37, 785 (1988)
%    hcth    Handy et al, JCP 109, 6264 (1998)
%    olyp    Handy et al, JCP 116, 5411 (2002)
%    revPBE  Zhang and Yang, PRL 80, 890 (1998)
%    pbesol  J.P. Perdew et al., PRL 100, 136406 (2008)
%    q2d     L. Chiodo et al., PRL 108, 126402 (2012)
%    rw86    E. Amonn D. Murray et al, J. Chem. Theory comp. 5, 2754 (2009)
%    wc      Z. Wu and R. E. Cohen, PRB 73, 235116 (2006)
%    kzk     H.Kwee, S. Zhang, H. Krakauer, PRL 100, 126404 (2008)
%    pbe0    J.P.Perdew, M. Ernzerhof, K.Burke, JCP 105, 9982 (1996)
%    hse     Heyd, Scuseria, Ernzerhof, J. Chem. Phys. 118, 8207 (2003)
%            Heyd, Scuseria, Ernzerhof, J. Chem. Phys. 124, 219906 (2006).
%    b3lyp   P.J. Stephens,F.J. Devlin,C.F. Chabalowski,M.J. Frisch
%            J.Phys.Chem 98, 11623 (1994)
%    x3lyp   X. Xu, W.A Goddard III, PNAS 101, 2673 (2004)
%    vdW-DF  M. Dion et al., PRL 92, 246401 (2004)
%            T. Thonhauser et al., PRL 115, 136402 (2015)
%    vdW-DF2 Lee et al., Phys. Rev. B 82, 081101 (2010)
%    rev-vdW-DF2  I. Hamada, Phys. Rev. B 89, 121103(R) (2014)
%    vdW-DF-cx    K. Berland and P. Hyldgaard, PRB 89, 035412 (2014)
%    vdW-DF-obk8  Klimes et al, J. Phys. Cond. Matter, 22, 022201 (2010)
%    vdW-DF-ob86  Klimes et al, Phys. Rev. B, 83, 195131 (2011)
%    c09x    V. R. Cooper, Phys. Rev. B 81, 161104(R) (2010)
%    tpss    J.Tao, J.P.Perdew, V.N.Staroverov, G.E. Scuseria,
%            PRL 91, 146401 (2003)
%    tb09    F Tran and P Blaha, Phys.Rev.Lett. 102, 226401 (2009)
%    sogga   Y. Zhao and D. G. Truhlar, JCP 128, 184109 (2008)
%    m06l    Y. Zhao and D. G. Truhlar, JCP 125, 194101 (2006)
%    gau-pbe J.-W. Song, K. Yamashita, K. Hirao JCP 135, 071103 (2011)
%    rVV10   R. Sabatini et al. Phys. Rev. B 87, 041108(R) (2013)
%    ev93     Engel-Vosko, Phys. Rev. B 47, 13164 (1993)
%
% NOTE ABOUT HSE: there are two slight deviations with respect to the HSE06
% functional as it is in Gaussian code (that is considered as the reference
% in the chemistry community):
% - The range separation in Gaussian is precisely 0.11 bohr^-1,
%   instead of 0.106 bohr^-1 in this implementation
% - The gradient scaling relation is a bit more complicated
%   [ see: TM Henderson, AF Izmaylov, G Scalmani, and GE Scuseria,
%          J. Chem. Phys. 131, 044108 (2009) ]
% These two modifications accounts only for a 1e-5 Ha difference for a
% single He atom. Info by Fabien Bruneval

% The original source for this list is Quantum ESPRESSO
end