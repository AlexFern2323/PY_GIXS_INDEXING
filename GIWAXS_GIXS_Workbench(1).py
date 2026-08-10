"""GIWAXS/GIXS automatic indexing and manual peak-comparison workbench.

The file intentionally remains self-contained so it can be copied into PyCharm
or run as a normal Python script. Its three data roles are kept separate:

* CIF: crystal structure and calculated reflections.
* NPZ: numerical reciprocal-space data used by the automatic calculator.
* PNG: display-only experimental overlay used for manual clicking.

First run
---------
This is a standalone Python script. Missing pip-installable dependencies are
installed automatically into the interpreter that launches the file. After the
startup installation (if needed), the normal GIWAXS/GIXS GUI opens.
"""

from __future__ import annotations

import base64
import cmath
import contextlib
import ctypes
import gc
import hashlib
import html
import importlib.util
import io
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
import traceback
import zipfile
import zlib
from dataclasses import asdict, dataclass, field, fields, replace
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
from typing import Any, Callable, Iterable, Protocol, Sequence

_REQUIRED_PYTHON_PACKAGES = {
    "numpy": "numpy>=1.24",
    "pandas": "pandas>=2.0",
    "scipy": "scipy>=1.10",
    "matplotlib": "matplotlib>=3.7",
    "PIL": "Pillow>=10.0",
    "gemmi": "gemmi>=0.6.7",
    "cloudpickle": "cloudpickle>=3.0",
    "PyQt6": "PyQt6>=6.6",
}

def _ensure_python_packages() -> None:
    """Install missing third-party packages into the active Python interpreter."""
    missing = [
        requirement
        for module, requirement in _REQUIRED_PYTHON_PACKAGES.items()
        if importlib.util.find_spec(module) is None
    ]
    if not missing:
        return

    print(
        "Installing required packages once:",
        ", ".join(missing),
        flush=True,
    )

    # Use the same interpreter that launched this script.  This is important in
    # IDEs such as PyCharm, where a plain `pip` command may target another Python.
    if importlib.util.find_spec("pip") is None:
        subprocess.check_call([sys.executable, "-m", "ensurepip", "--upgrade"])

    try:
        subprocess.check_call(
            [sys.executable, "-m", "pip", "install", *missing]
        )
    except subprocess.CalledProcessError as exc:
        command = f'"{sys.executable}" -m pip install ' + " ".join(missing)
        raise SystemExit(
            "Automatic dependency installation failed.\n"
            "Check your internet connection and Python package access, then run:\n"
            f"    {command}\n"
            f"Original installation error: {exc}"
        ) from exc


_ensure_python_packages()

import gemmi
import matplotlib.pyplot as plt
from matplotlib import patheffects as matplotlib_patheffects
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont
from scipy.ndimage import (
    binary_closing,
    binary_opening,
    gaussian_filter,
    label,
    map_coordinates,
    maximum_filter,
)
from scipy.optimize import least_squares, linear_sum_assignment
from scipy.spatial import cKDTree
from scipy.spatial.distance import pdist

# =============================== INDEXING ENGINE ===============================
# Self-contained GIWAXS/GIXS indexing backend. The engine compares measured
# reciprocal-space features with CIF-derived reflections to determine orientation,
# calibration, domain structure, and experimental-to-calculated peak assignments.

# region ALL-230 REFLECTION DATABASE - collapse this region/function in the editor
# Compressed affine Wyckoff-orbit data for all 230 crystallographic space groups.
# These symmetry operations support symmetry-consistent structure expansion and
# reflection generation from the selected CIF.
def _all230_database_b85() -> str:
    """Return the compressed all-230 Wyckoff-orbit database."""
    return "".join((
        'c-riJ?QR^olI6SjcM-@)W=f^`x6d&c3@ll)B+GVp-@dauyO?+1vfaLFOG+~0#375U;s7;^MUO)ylOOSM;++5YPTu|7|9kiKj}Nzh{O8BpFSmcZfBMfqzy9Yxe!IPSc>nb7-~RpG`*;7bd-vPde?0v8@b$|-zWran|MtH=|MTh3hu@waZXW*p?|1(u-~LY?zWjOk@Be=IAKA%&?EZav@xKrM&tLyP+59gr|NG&8{^x)GXZLQ-4fEmM?%n<2_wn$%JN#}BzYmAs&Ea=__`ROS!|Q3x<)(4?PyZfr21Cwb_;qEzO*!{1rdNMccJnL4s{iu+FMb66;#*VE|HV7%95ZA2{ugg6t?Ki4{7`9*&wufz(sz9Ri`Ue@sQHfXEB(d4-&8SoQ{``_+xBiI^Tp)f&|>`9H<kV(b3;q@4K3yky-+-amNxKPmsNP}bF8_EZ=c{d_Ttjh3EsUq!tK-Y)`}`dXfycgj=9cmGv@i*8E+7cYxBO+44%INX9j;iCKuFWQuv;5d*EVz1H|b@z5I*ssee%yeqHwP29Un(&(ryB-uv9V-`>~P)vtf*zAkb5=j;XXAA`BvG0;rGc96equPe&2S~@yb4R5W#pL>~dU*_xK;D0&%J|BKhhu^2e@A2@n=b^ryfUii-Y!!K_zfdm>(BrikeC6?eu72qSufbPycmjN0JTsv2d4n|!|I|R>VLS~lr|?+mE9fJSU`@X~G!Q&NU!k}3y)T)5SG`~`()T5j>DSc@2FFO>OST}Y^Xl%}_Uf|UMhjk@CJekKdvy!fBZ_Yr)YMAma(D8sz^iiOeu^iuLu8%7b)|2h_lF3#<y`Kz+|lp*@3(Ij=Yhr8Z+N@=5o5QXHRRTeR2GBTBt5nH;PxGk8|AZ;rK;~{&ZXRESGnmRo?VVKhv;CmxpWkwvkRugHGba5cfcCI<>%)A%l0;j)$B#c-U{2>h%t(+aB7ppYW5;zZ-wol5fJY0+r4k;A_DToaMzW-0U7~G7ZH#rhWj1?W$~l>sOMw1mBr^>-fq;B%y1W~Xf51Y%J-{TnMPwd!wc6@d!E9FrZYHG^*02&fm_O=@!yadh|H2%wEQmm{)Q|?dJsCPU$paf<grHm=KhnF$S!jKJ{)-qx&K~nt^5qW-&$+2%Iki3z3sy|S-Jk?DmIeq_vMxva=e*Y=!KghO-2S?WJuZhJ9dUN=>v3;A!X;F*%{JgWVQ4Y8?v)~B{I10sViJN^4m2ZTYvZQvdx13Jot2f9)ABg{Qfxnen0$vJN$V47XNH7vJEYE`+BKK{w$sg67qh1HHCi^Psj=R(tS0BKZ?h;Le7f5n!@kJgIghA%&(^KTmL;TU;A$u`O<&e&gcGnUY?5g;Z4?TiZ$K#snQgPt%NVX!ejAHLf2Q&oP~!<Qy{`5`|ST3^lc^jzkHj&5BfG=kK0lDzb;R*c52S<^J^OX2=hb#xSN~)aW~h!IzyK7X``ORFS)h-N%7S&a$?S2A#VAw4#u{4Hs*Z;*%((^@B3Nlf_ceSzox(Hsv0}{jC?2io~Sq3e)wKB;3S37A0}i)89Ygq_h6tu$moyv#k*;XI(Lj_+R1v@?3LX)%bVh>+CFuAulG~6ZTVjB%g_H*91oRul5KdX;4v7X0onjjuh#gnI3lWOah?HTV%rOLCW@v@@4FsACWD7!@KBrw4@J`zWuBjH2tf(VmP<z^?}lTUQd|_hu<`MG3T1Fnbk0RFXW97hc>awDNiIe}U%&5<kXCiHK-Vft7>5`<&5mG-pI_Rs4bhkGrdS-wB~8HSCo?sqoj8=oFAZ|+MQnhp!l2Emga(+YT@VKtu}{2haLpLpNFD8MDNuiD>(_|~2HI;ehujexAj6=)mbW4+hz$@c4z9=wKJU_d_bzQPQ7=HEp4SN-JO%|C;EJpuZ^p$-9Up;R+I$)wsm1U(fE1emX@lh1a}Dq)g^WA*E^RQ$Zi_@+l@>;2<fv_;U6d(qsnJs-DfEW=Bdx&+Sm^N??6C+@git4{gB_iFw5nZy-Cx4&<pSuMltn^b;YeQLNM2z{H3ltxh(fUmek;;GF8$Cx?wThz{ayr%Rf}K))Y>ell}RkwPqRW&GWa|Jmpo3pS8F@~m;cVrd%RlXd&u)^hRH~Sw_&6eWB@MP2DJ{r<tqm}FbxjCWgSrKhfcwMry%~-mIQ)yE`i`~qRa}nq~pmCEyR;2PEk-To`f)!AdW2&nWW~N7@{K9xgp1s3pLgOp$#?^mdD<j8lZ<@06}5t;V6WPB8DNFxhO12LH!H);S0AOhQb|Z4DM)u-PgZr3E$twtUQoC|N3wOv{&sVbQKM*WSTx<=UWs2EHX`>u+td|bPk!OPdt%q&Y(`DRLIHEz%+;w#o{n8lqeP#06>XieH)`ZQS2KeWf$j2w@9QEqqCkd$QGM;cbG|yBfUGnWHrv3Qn{c4Bk9DH5~HmFLKZ?M_PBC4QoZ%PBs9MjhFHYjFL&0So{2qy`^phTQCgH0zU|D-ZdtpxLJcRS9wRj9d^bc0t#syl<%2G);e4zk{`admPUp>Q-yYuv!)2-GboKAD<Y(~xoMahz-##PPC%|ie+$yx@oo^wMcO@gfZ?Sn?#ck2e^e-^dzrM(8cQVG|PPTMDe!0Vzafh|ly3%3>%$LLK^$z>xoMiD9I*i5L3G}$Njj`y1S+9NH>9%;!N>Pvpqi4kxLul6cLvgHJQ3wYED`(!v1jI*m@ev^tXvJ#w){MO#lJ;j3WnIt`Ui*Q*ITw$F9<9<Np~-+OZyNs7ut)Msq;q)!T~{ZKO**Rxr(kBKSI@8}RV(5dCUaLyKj^I%HR%+&wWS|)*{*pE<@XnHKXQLP=l)t#cxfNPa$B4esA!KG&(kLI)hN*A{pVMGM{Qyr44sG9VSd%;&;+Kz0e5Q#?shezcG9~F`X*RR+_`T;C$N#qce*BOck6w5?sQEWz^(Vih}&DZ(_5k4>^*L8Bh!-M5YN7e3%6BbKRVWg9QF119(I5aTsrFdV<-z+)8ghhvW454{l3=Bi)$j7Ma~Vs&<h&yByVbA-_+9iEz)mw?U!3*F2b%N7o$)!3N>vNYOFAjiE-yO*^BC%^Jo_EL-Rfqmyta0LJLDi(tMC6@<D2`6<|L~uG}YXBHzvv`~7IXDd-y5rL8xtuwjk@!x9-%p=Ds0t;0H*K>VK4SG7Tme+9CUmJU0G;8juy-0c&UWZ>UUG3O2TnGf4s!|hunVr!@|*kkY-k$sq4!wBp{Yz<opUL&#(1zW?M!Hwr^Uq3I4I3P>K=h5zf#H7|c_qjzJIBa}aOn^Cr&w4M*QY5cE`VSiI^ew5b-9Dz_!{8n%-hXg2xJO(W(f+x+z?yY;++E8$aDcs`OL}-Sg^n)q;Nn~8wXTT)$C7enC6Dg$JA72U6SI9W?&{t!+nwLlGyASSQ_OZGW470U**@lNt@5n1Tv7k)a8v)^!|%TjzyEXi{qMu?zYf2D9)5owe*ZZ9{y6-8Km2|>{JtK3Uk<;|hu_oT_v!F^JpAlx#{cUCUITY>D$w(1fj{fuL59|ub#Ybb|680(eR5SCDHg2;sU!Th%Ku$ym6@kwVsE{PeK@8Ud-?y6JG3inYf}scHoP>vyfOZF|M^?DLDwp04SrkY|2i)wB-J9jNBC-$|9M_qNZ=}a+^TP@{AYK1#;>OEkM0DGUrphU?i7t*P2u<MB#mE9;Wzh2vO2TwTl1we`?b;(h!_7$Pw`7}8NmEkaE270D@}p;FueQ<PuVGWomX?DXfgP03ZE)XffymY{0fiMDNy3|$#@nXDoueHLA?A5_iw=Gxjxd%Mf(W+Mc-k$t2D=?mQ;~BCL2<P>nKinJe6>rqA6cRDl*4p<4=(}{?Z7OW#9JK3O=_`E2s&|N%+V_`hy_!oO588qeGzzd_i{ZC#r5mjfT#^2Ji<D;$9RCn4g-~1`+leNA?;=#Tq{XYaC4|h$-7AV<d*l8x$L?Vp4RXSYQ=@x3}SvlaVH_CBfD>ve!7e*7&|SBk`nZ#hkeYZIr}Qwm&wrk_meoQBsjL7WNtodyOM`ji2}NU2(+JM6M4^<Yc`yTWF2EZ>hZ6noQe7&KXVQG#WcX{QxD!7spNK^*mpm!W-%gG?%~EuX}k{BVA9(n(2KG>rME+(l-$DH;8}l7|tq!l)+v8F>sgOT#Cx#8R0Gq>8Ng$vd<L}SP*ylThn?Xs^dn*fv`C-HsAnD7q}9ha~aC&Y`uwBSh_$^@tDd&X*IUb-6`6-_AvzuXuT_*ilbUJ9y)}6A{*3m09ovz6$I?rUf9?svcae}XkFvS;++PrBwj^>QO!&^(JKfn7}3JUC6NtcYpmJF_r;Of*6Y){OFkAyHCrzy>kRISL(4tiK(mW(i^H-#-#{~i50$<_nTi3Ly=7<6eKm^fN_P>v3RIJ`Y9jmR)2q(<bsDI*ny5FTiu7?h<l8{K71GY`8ui8lsTD%KJ+E=Ms5c%+tq|%h+2>IRku|dsSu+=rr8g)>51<z6O}~$oYfL0;vKq7uuqNsare|Ha4Ny^GkcE1iE2G|E&8>yop%oPdIcxm9k87dccqN#LT$*eV^>&K5LUznv3iYPBI$4-8Pu?VwXCST_)Iz=S;uA%8Uy)tJ^Q~WZ3$9red0)=tsqvgiE3ASZMxX|ACe5%4;#^Xkc%DHvRgGsB#8Maot6-1lE6~%?rFQ9jcTVc3hZgwme4xY{e0NTEs{3lFl_I=Ibc&xHnF&f4)t9FXau7GsIf!;-sw%vizf{tEQ(bts)G09DTJvr4p&BY=E>ugtsfMYT3n3D&8$cxZ*n?WA=X<&mzsIHXiRd@hlp|3A+E0$9h97ZL1gQBMjz{g-kF9s?gT2B_Z5KfhF9%B`Q{sDqhU3zY2T-RP;;UY~Z@NBIuijLu0{7kjh&lhmo(q5I)*k_-i>r2e_mB$*HKi+dde_Ws6;g8GC8ErYiuPRi$F}*#81!I8`NncH)*)6ZkX`#8xWFEvkn=vbmAc-D)rLa!-n?SNl#dX*$V~aNgDv>prMi9U_saB}9sI>RUvl2L&Lh7Y=#&z^>4{_r`P@qgI>pO|*g%aA<*^3YR6N<M4-I4pp@d%53?VY(q!&<~>4VOuLTCE$v8h}aXSb@4A#bbFeTrGwmMZbz{k-)=v#>RRul#=IcDlH0;;`xM^_nj3`_t4pw=K%#2~a>yu~T3*MOK-3ILdYove_W%#Da7-NKMt|wgEOsIzc9#4H5_cw~x~7g>}iNTX5Yp&a@ZSB@fQ@yMKOh`Y26GI}ib-*?A5Yj6#gq)Z#hV4-K#?#d<n-Acc^)bnb2%U{#9s^co|3jWv3W*C)v0(Uynu{sdV(TJmDQ`yZDl;(yeXd8yhZ(nCIwIO{@PBFWgZP(E-WXQ~Y1nO4P9_h+;BiB#Yb^HHl+;1R>w=}-@(9q#I5^OE(_C(?p?W*MLHFArp&^FVCX$@@{A^kzfa=~-=~I!P?V$*7whp*q=Z-K?roolgz0+qzj*6*ZqyU9WYsXOF7ax>+f9liJ<XJ=e{C$y;T$ru81<$5B5!2Zqg>^|OKun|H-CD^Bjmh098)D(-lBSONnDRFTJ&(Rw0&kI~Cx`eD@2>f@BovWC_cj;aYfLrog?qMAuRJi?2teqGuN{LrlzNTe_+mN)Fd3nVh?6v`VW<LH(;05!7FN9_R26-C?c0BocBbDjtvUdpg+uIKx^nAIoQeiLs_lW(mivbcvgu^D8K?@jELH?xgc=af6P+jtR`_w{_*aZJrSo8O(?10xUp&Zc)~Tk=y!&1UDCEsmef&Ni%!m+@SXf-Ss^C2PJHu{cstE>5~};gsp5k3H+6bkdA0=TXwnjlfo!vpz=pxeeHAl5~52(oHm)A<{=E#Dq3y?<>~OH$>_PH{de++9dqWmt<nUbrbzLp2%<0#J&M1_8mP@pOFp7(Y`tXyf5mjKmhL>pzOHdlo~Pxz7euz3O)e$sSPVKk+CBAKsN&8OWy>GqfL}E-;}*|^IxsdeAjB<W6jLB6=+05NCrA~qZMdG%>*DDSRh@L+S-9K6iO&AMVm<4N+wnz+c^DXVjb?6lTRknA(L{J$wdCFkAP)uJVk@07o0c9$)#w7oT8C7Q#9qZCZnw@><C(uVrnz9&e(5nO|2j@;kE215oPGGs8cjYO(fFQI@yUt21QBXM545wNQ62Q8%iV!IdK)wERh(=iNy5?iNv<~_0S|28`C8xo8+pEBjl$MzbCn}3SfV{+r9h!@cZrX`+E3&Is85!eou$rr^D~@@UyFg{n0PFjxP;M6iGc1nXqGjm|oCz-AmRrg0AaHvhMe4=(?(sb-p3eN62Qa3{=|6K8Q7ctF%)q)0g^&NFO0vmoQL;3Hy*v{8~II=3+UkV$kfC5e(u&fs`y3-5`C+uZj#`dPOYj4HAl1CyH3+4KUD0_-dzqE}p_IN_}>8nXy=bDNV>p+~6^|vOme1nen$>`qTjVdqu8xp);o|h&(?vK>ps4#xY&h`P&*lR>L*VTwF7+3&(2D9M3*fT4SPTGm-s<iS0JLjKuD#lpht_lGxRMda<4<D|4h?GIrw<y;aLX#LfXQ)a+b})e8DiA#itFX@%Jkk=7v?Ndoq{9|Cq%P=+sBD2^FiSDFC|<2YSGnO|u4E^;slKP2MwfJBt+5$au=unoc?UV*e9NJM3w%u|I#6x7Ll8inR=AQ2rDnhQcAKJ{7@;l+i8iXyxL7AlJHpzVc<BG0?@(OqOdQc<L;A}%sx9>qkVsts>|h3Ztvbo{=HvQkw|VPwX<gv=3CnVpCowkq>|allR?!l9}RM#zktXGE$hkA=+mnX>6)rLJHDW$Td)Mm1xr-yW9SUf6X(vO)HYlXR^MS531edcn0aeP7Q-HFnvuOkXo!)p(HOwW(@+%H(+U2p;t==tTCPw0n+g@A;*BUrm!|q3&UHg=g{#=X{CpMcvC>QJ2^6{<U8MNU%f!mX=PJ01})v{loN30120R<ShH$`LmY#{qA3=sR>X7u%baw1hAs2V5}#k2pU}3Q_dPhJKmb0D#dy;DP^;hp(;fOD97l{e4@e-u$`JomN0AH@V;88rLIV}mrCdfw>OjiNpc5?1D5f(ph1MXpg}Y>v{0MmBH129s0$lp&p4lMe5|xf_2(@3yZE48s=rIw1N7U2b}4L-J!8psVLVhT@Xn3XKI;*($Dw_O8$YFPCus@V=lfpC$*Z=y`Uco2IeDR8T;48C>KnVtUkep@2xG2iDn{CwOBwysqksAY`=^)M{nITFn6sfHy+A$cK{ME38n4X3Lfyl$nmtB_!wRN(UulJTsndy^>72-k$v67<_c^GC9n$g81f(NT323EX-#EnM^y_bzXqxHQw+-<q{rb{=zrOU)uP^QQ>x-jb|Gw$+BiOlDvHXa2?$s<m&fG5O0_Dg2F8xUEf?mj7&<nS}6wDjh40_r(*?B0iSKo;`@)~@(OYa)!FEx62N9ZrL`A)~^FO}x|OU2P&dRrXC>1t+-m_gg!=j)E-+KcE43ea4*uppe88KYlBr;U6}s>0Sg(OT(jlSVq_O3;CFb&;{6h%`{NX`^hnHeRG!k5E>;E?ykBvWZyVy@~RUVx5q=XO!O?D$t7fJfDf|sZ2z_;;(|iRjy#b)kA(mDlyJO4|(Vz_o9b9!3kkzuQJGY3>Fh6elfvovb|Z;Ka8iaMX9{iujL)5$m%uO-mFif$U1c5hfchrZ@vQM$dKRmzqk|KDVW$V?!@*9Ch|)$(I2yk{ZdTqx9n`$`3uUxfAL?LLkFH~kP=c<H)sx{2nKD18M7-~-5`C+E`t;ouDU@N=MAdz{_@H^D(MZwAc3NWPkQ4^noR7p&cseUP4wK@L=N^Rb~<8WCrl=KPUpzvOEZ~#X>Npk$>U{8J?v6dFqE$-%pRX<!fQQ)D|yCiF?bpYkViMjp7B~2kT<su$eSI2yibUCiBW(pY{o9K!8CeeA)}iY5bH~%r=lBN$unNn1N489mAz#kgX|e&_Ao#^$QnmB7|An!zA7L6SLKnUMR{1sLeiq`N|n0W5?8eU36&OQ>=|DdAz7R?*3irCMBbK%P%b2IYvDp{R5-Mz8ZBLa*H=!M*vbD!l}(>_X)WPUV%mBlXRRl8PIw|`TqkzIbt30mCw4Yd!G+ioV`n{r&4R?Y-ood^F3px9zAeN6iZ;n(Y?2GL?Mt7e@!Lb%+*ZgQ*3^ZxrLE-eA<+0%MW*w%(hjZ3#wpbrJKJa!T~AQBHpynXFT(`6keMJClERDU4rSY8WPSbiFec4~%%r(En&Vq&jIG&d>nmDM<eKm?TdnIg7C-OM>uwpLg-i%iN}Bq8r;N~jP_F!Fj=!98&GEM8_xmQA-xD=uP@mFP68+179>SCz;wy`lm~^+xn2X<_8dbSq?N3YAV1?%Q#9q`C*-U#ek&Irky_hKKtl',
        '(aRdm^@Gtt#D__Tu)Zj@Vw@zIF6~D?&%~irk9^i^?lU53AGOi>n?7$7p_!$i1jHsJs^q@JZdhi0wtN`(;J$MT14<l@qbOxPosqzxh2P_oCjQx_boLi!t*?Z<M7pXsG?BU%ZEK6c&50cEam_cDZTra@$}!<7+V}is*R186EE-ERV(Bt9#wDf|B?q=_IYmq9}57_Fb47;1#k9lUP|nwqa_Bm%ui}mT-j*@)fcTBesO8eV7{KC9)5nFU&jtg{j3h?3Ldw>&N8{UTz!Q9^qxeVBDPXcDgfVPZIX#)Sk1iV5&X|tzaUJB+wEjzOjDSXp!f8<!M*$!j!=m2T%(F(r;|tJ_Hvrl`Vu8afOa@Y7w`GcnPgx%$((og1xVcbDv$5Y+A?QM%a-`nMYccFfsU^W@uLC(YLSQ6=pDZqH|}{vWAza!P&eGlXG5yZJ0EeInnkhM_?a<YY3ja71)S^tzpiCH?LvtZrH$&Xgo|d`#u+<BV+&sQcg$hVZ~Vlf`w@bNE7e@Zm}GLoJlpnk&v1AGrO#52%xxioGfsY4FHpsOsS-HZzQ8iG0a&#L=^VKM6w(Qa>Q3M0&In@TnwGEl8-ood0fd)%nClyb0s+#w)S=`B+c<Af^N^P>sWJQy8jZ@NkCmFaRrS}<@kF`68QldD+rzeuFw~@=I8(ZjW<w;MRuu>cj+a0m#*wxYRF;w&pEVB!4a;|H^%nq&zSr>mQ!+APM*k71xETIax_Ga%848ycqWvUo??*#p`6(wM_!Zd&6@sUUaEY5aktqIEe~-E3J{bxAi4$O#+uk$|78Qt=i-DALqdT?T@O-1R53_r>klvQ)TiPEJA=8N#TX7!>@W<{dFE{n&jt-{7xX1kkshQT9ZDqAm#mG`q2$T}8|<SCUfaj6K8O%5ZIDZ;XS4U}7I-M$+I$4sLZ6p0v3-Jj_Zs^d-k>0cO&O=}&<_@a8@)M8_r!+2C~RU!VRzXz=8)8ZK@&sW&k809OoJN;;H7)&Lm%KSFb!^i24}8cF_E9`m$%&wbp|NMr#SZcH86vE6eu|y<3>4lt_`d5oP+O$x;a;<!wYktPP8!BgR*%Cg}IiL%{wm4ZF9tcTBQlL6gl0(oN3nd5A!l$d_2v<TtmX#9T(;-8*5^3{VAV5+H>8wBc$jcl?vVnefce^n?qy_FtA6cNWfLsW^rn4H<~rAn}rqjrkdH?UchsyQ@voNi5GU7m~2qrS;LT}oT;*u(34L$xRPhQrZ`=SXtRK_&=n10y3m+CB+&zk8do&9l4rbDSH7qTZm0)Fkdp*Au(fD#u*60YVkuFCU~`J|unnaW^z(L+Wa%M6Unb#Ssfe?bDWq{bh1ui`Y<thQ#QSR64ht#J^r?HymF)do;qOJTm@9(m3SM_j6ugWa(~4ZSjflQj_fO@5myPmVF>50ceWC83$^|cI%3uj0@l~Yvf){iG!ODqXAI4Mgx-BlU<;kalitLa<)Lmbwu+~Vnx29H5*aq1%=7|+DAfXEv#~TrYdTJ*20R13XqpCn2xx8Wqwd6upB=k)S30hxN7LJ)v)y%j`K2;p~)G7&DG*HGVmliUQ<?*|_q~io_y@YD*ZZ?!|sb2c9rcCF2lIhKvEF+3B-{R@l@2RULm%75R+mEedw?T+7l*&3`QYHSZ=^v&CQYB$$@8m2=mH4xk`lLz^)bL9%0m6y~L2W;}L0>W}9M#~8>c+wbIjGp8R5cjr)cMw(994~wscPgHR5c$fC0_NxxH`+Xa1*b3vKl||G?MKl^Q6KyIJ@a;bOWr9j4eoyP=O1eCV#_Rm%o8dR={Sg&=E{DT584eMpY}BvkTi<E7@Bwd0Y-XdsAtJ+39&{0zVxwO^iaw!?x4aWzxYlv7L1h;Bzww@b3xoHJ2byC&E9t>6nhsVUqY~I``~tgKyajc7rIL(Q&)TU(Vo-+fulMOfUc9d*CnPJ<Wferi<KZdVPZ4cYA_f-b?Ni^z!j{pP=}66#M14>1f|+=6<JSPf`~*osYP!9cF<##BFTy26^q*?8i;)IN_0f6U*Ok;;Y=Z*}jEuiZef*yT927*IVMBiTToddkug7t02%VdQHwF*7ez2<W*|%o$|c2-|oWwl)DSxz6Y;cxCi-Xc>5kq8C;Mvta<#7-2#~JG&A0V$a^o_Z^`zYYj$X}3G2ksQ_t5=(8RX;te&V1-g<vBv97}@3Ut3Myv$L3>@M1UHHEwG!rfO>xa}_9eKmy-#nY<_@)#+{u#cvY@q$~O6!P>cg}l!)&g-x7k0%o>w|kBE6fQfW+dvP%Ux7mIx(k>?S4_lwaMMsTQ-Yg{D7T=LPlAAl3Te-9XOINv%(cNeIGE;C)N}S(;;^T8JBNbK(<la1C<gQ(pVLFI|EYLVMS-=Z$Xq&BWnu#ero1%F(wj64X_bk`<7F)VaDv6HBdeOcx-ahg&#ZmyKeKk%v<{=@M!c;w12#{~$Znosc5|t!r@0>hO(sC|+@b-Z-y3B3sYi$gdi7IzO<$$_)Fe&D(huOdk{F#cr~I+IyuD{1(`6}D`}nSsX1o-Uey9UvC@F?fvv$IcwraNOm1M=N)4)=;^3X=@#1or+kCfIE!+|*U>`izBdYy#A+E`X~#o(AY5YoOni7IJ$N2RWa5CMygx`J&_<xegsLrfhUNU8ex+yNIfN^#S#xL$2?4b;<;YUf&Y46mtENL?~BnWe|$Y9TY=U9N8=r7wo7MSap;=+kDB&~#VJ=wJFEw`j^^Y5NKl$9)#ru9j%H?hse&xbCMWJz8yzV<|qt*4^2*O7N$s72Ot0<GQ#zTT-Bj?+e!Ei2TbJ%(6t;fe8<-i^U9_^1vw%@)s8QOA}Ht>Tsg%p1}L!;q<wqxfxq)qVN8Ry|X6v=AQ@>90L*GcY7_H*d-M2|GjHKRkj5cShnhI8?^OnWg2{*>TMFkG&<DFcc?Gz^mk*pGrV7%r;Hn2@*C|W{Mb?hOzq)$bojT_;h&{h)}_eqnn;l)f)tQ@`R7|GroziM7voQRZw!9BTK0U+9l365Xtj|ymP$G4@HVAoc*~^NYI?2kg3bsp=(;&)&DME6v32{iteXbc)Ltx+A`>J<CddIyD4(sVy%<g3B}m?3-q+Tv%ea7Me$^xE-)-C=Wd-%bY!4fsN=w<y8<pDk8zINO9d`|Lo?<iA*IH_B+<N(1O)*iEc(|#wL(4So^r^tehky%p$FB{f?$$zzlalz0d(un)<9}a5_O@5}Drrybcc^?f$S12$9Ww3`NAeP{<R!lRYlXk(46^b+zrEYN`+E3&Is85!em<7K-v%#83+Nj<_$Cn1H+=9-09gXRx_Lp0YTEl!(L(l<K=!&X)HS9sG4k5gglrFFUin;UAA@V`L--0gA_dp@snR|M?+Ra2ppfl>cZG+wRoHh$pj*TY3j4NrthCCiogC;ES+Kmp+B#TU2Wx9}RjOm~2i#VgLAn0?!P<H%tgTnm(4>O7i1b&uyZ^)hZZWsBv#RxdWn=D6)BE@JPKYIr1>d?;@cxli)*QGXtoL~7K`y(Hr(4ADyy;)iwO?TE+&QrJ7ndum=$v-5K6LA55lPqccv`6FWx&%Ot>~rUR5$rVNjNP^udi>H;5Y!!Wpb7SaH+9qG5B10r#XUq2Y$o`=C<C8O#LxWW0KciOmUzXoh%`zZOX2(zFEs_@9!)vS6Pwlw_nBU(|i>%@{5I=e298vzZ!?XkdqHlu_uq6d`SASorJ{&k`D_vEdJ#-i35;!hC)@R14vW7-RLd!wVH3fkZ&QCMS)xz+C8}tAZgeX0(B0=RgUace)i5k4{+tlfYO==S2hYtJ1KBwL6Gp{09O_V2@iS3A<sDE8HYUMkY^n7j6<Gr$TRjM&v>;ANMr#><S{9^>5SZgK*>EIkquWR3WD_pB(etZ^g{#jCt{E1AB92&YbpnI#lKqwcG@dub2=$9)T?!1Hc55i$i{5SIbt@2kJ*gXF`FwLvng!MW@KYFF9&8b%d5u0fDX*2e~t6_hVQTA#lUPjvm=p>*%UfvGtx1eGBBHcVm2civ-!(ylGW+N@tL27Uguy8>?0vN^g7SEIN%W1%zTyr>UBOkoHeW;!1nxs7BPq=^!)h@<`gj0@aW`504)mHg%G^$<4x$<esF+t9Z4ajSk<gH42=V%Km~zNH5y^}m3Ha5GJG<3>CE1xGkKTJ>|LtYJ|66yCHJYY_i3P8Y^30SH&nWZN_X=r-9rO;Xdn*_<e`B)G?2T1{p^cE!Q0Y72GhN^9G{D!2IOsnREssO>bz<vKi%6!=^?R<g6JUIa;2KeyVy$)t0q-yDm|=d<y{cAq*OLnw8F~N(!(7!Ri=^dv)32QT6-z|_!-awOvqPjI1bIJb~~aEQ5Q@JLavDw6hN>tp?aS9_|GJaw~lSvZXX7{h6h#13Ls8qte_>?QbCm(9#A2(NOo7q>heFEZ>v_^sq~E+u$>y%Ta~_3%S(5cYU5hCNm!oNc3>yq*jZR1XRls)(UbMCq<mh!e0M!fD!za_H@~(Lf&9nO%o6$v&joVCe`UL_2Jt~(I9vLZ1N*D^GgEBbE|c^Cs~ZQYR)`jo3j`%)Ig|&O7Z;?OTPS3E;2M7@-kxGuV}9QN`+@{l`KGvbXPJ1!r)!Cm2G-HFlvg54GEL{&vER1?0~IKV>X-)htH4$S73g?zyY-z83}4igy|KWH|GN8&|9u16-=2r$$F-gn&&$WAx8rKBAksnav?s-;Ow>f_*mU^xri}&c)H-O!)~U+bOmSr&&)P)gtvsGM4fxN;a~?5K`wi=K&}Lv``PzRtZ%*^!^zBqP1Vi#5^)usGE&Xuz8ngZS?8NcxoJ}pUgLsWc#(64{L?0LIh$PyHB(fSA+t>NcxX$~jtkjJnoRJf!veK4+%f8OrNi>~~f14$T*J0-F#>LoAU7`cbTsc4Yc3|f8jxJc0_%JYY`qzh7`I-Bv{ET?0op|VC#6!J3`=^RD8}U$`{CuT&sC%{PmY$m*m-S^44A~G87!6LX)xw8}z!_sQ8uK1L6eqks?<Oz?;k9|2^6G&?wg+NbZi-W*oO>G`gBr!KPx-Pj1ATRb>Qrqy>C>rB?dj;N5$O#0RcmFb5$QA*=?FX0;pa((9O+0AA98f$GJHs(`;g+{JmwuvsqI6y)79LfgviL~*RuEbuJ9{r`E~H}^LhDeWb$ivwfPp2wT6R}-^$5fCzHQT@-KA%*0YHAHQewTjHJ9m&WkWF%!z!~C+MtCMQ*n`M|~pi(h1t7ZF^p9>gd|U%zt@alBHRD+v~a>UROt67rkzc;wQE(MzYr>!|Q6%>k_utHIu!rmF#uJ)O*&w6aqUlfs8M9@!qmCC%o~kzxV@Le!i~<f8fjE_xbR9I{ZEzevgOW!{K-TZfD?_|9UdOd8TlmyVM!zn8`pdc830V=D(Z(C6U`;phK0y%XR+T9V{j^F9iQ`jh~A5e@@mrh$kowO7VFcN~afV{Hb{Vr$JGTEmCCvoD}@+OO@@%O8XdGW1kEVZH<wA{Jh2wmG&_a_N|wKqTLm^y5hGrzVEIJ_|>!hvAZ_lpimMbQ^c21flKl_6qx~AWi82=O;$IV@>HWzPW!Kpg3FW#u}C@fy`QOsEL$-rjqqip5xPnv5UG4Ze-x7=y!8p#YS#FPii)i6PBkwQtP5|O_V<VeF?UY0k00G@ti{J@D?9-+C{DlPWeFv44D!QUjgW9)8ziL=G<7P0*7)9iAN%)lqF&A3$9lc)z%jTHe-}zOF6?yU!cEgeg+Xi|-xW{OOk_8Jkszoi!>$RQKNP3@^9a9dcF|3xUBs5Ik170cp+wn2O*192ICGh#&q6&1zhXDz;4SVvfnq}<20BT6=ua;9_5|c|-MzGd-tHUh+Lav!_(#slFx2Q9N(Q3702`G+FOF3Ai9IE(AHKjcbd}G|4<M3$gp;5n?7xr#87wF7ax;8@;)Rr6AJCd2mu1spKOWGm2#^n@@H8IRDvxU~jMKy^`#}${y!ub4ED+?|ti%xe@suTkc&qGeD0wK}xMKi`gP?A@!N8oeuExvxi_G?@iQOJ$7$OD~|JZ-R=&paP`)&VN_Yc)k&o|Xk&)3zl<j28TXSK_)h~;raEImam#hIhsQ)a@ngP7aYP99gyhE=~><NJZDW@iImt?|c!t7gNhU#;<7rC23uu#zCQkrojy?1*qF0@QF$0bCUNXA8`oh2kOG!ld4Gk4e4m9?^Lkc*v|Ss67aJ7(vjHf}l3#opJm^0jTAxtj>LLY@{#Mt+D9FMk<ErV?_B_X_c7>)Jc8bc&Tn5f70E+_Y~z7)4Av;kMAj7EU4G8G*cUvN;&B45JPj943KggAm#mV(9hgKFRq~OcF>s^Xi5hCB%K-_W*Nw!vx1<79t161_dJ<*^l>_Pa_E7+%rt|}WoEDPaL~(F+-f=Kk<2I{E=Ga2I1xq@L>o;ISJj|Iy=u^p?uXOve&|Vuo)wu;l<h+?CHk|D82$N;(IetN@e_Ni^1NQ(fuA^=e<QcCA8|J`BD6-B;u~PVEN=fR(??m;L$vRkHN};<AH7QKIg9gxG1T#EKn<Kj6sYWTW<_v;^2OrI7mpm<x#!rnH_+Nij!jND$0nzhW1FWcMevbd8@+lh_UeuNT66g|VdvL`onI4veoZQ8EzTFr%gb8ev7kV)pwMGMam0f8ti|YX<`XD36s$4_DC_|yVEg;gzg${ps-bY7IY2szXQJ_8KYEkdUGJQ1C|;|W9V2AcDrQHmqCCxU$>s}Q9*g&;D%x9?+L^r;a-HRJr`N?Z>?J$=*?GW9nbD=6{(Uy10ngt0cYc9c20C&mSab8i1>lTq#N|LomdM@?5fFw5$V2fs&L;+johR?j=?>%!_(*IY_#-*}?hpJES_~n1@)_JCaAcfV1~;%rDCRnYkfW<0#ltlxUH(MM2F@hB`D88KS`AK|#F)mQKoRd<rCmzwi_j-2z_Kl8m-4zICacfOOQV2ATh=Ze*wUk7OGo~;H1>7WAgDEnfz53Ln>!a#(H+^`qs!<fENZI3wruXpn^>Dc+x{{yE$T988hay~TP&L^_H6E|=%Tb`bJz4<6noBXP3}c8X9bg3hrLI2a@qVZdnBfd3WM9jF#V~^gu$P~!0efEG8?|(9*OBzg~9D%*kZcWH(58AFl%tbAKssERA+eUym?$+{r)HzM<6l0CTDEV<JL#r@QH_)0x=3NYD=nSEbYWR^-wqLxcFDpP0Yn)|2)KDe;$7SIQ;%N{C+?DemnfW9)4dAzt4xCPrUt~rxVi!9uUMGB=(EXX`8<+@TX1xGJ{Qz$EMdRTI|~!`6oH1o4NdDJ`jbFAqu7LO*1NkZqU1ifxd{ZggX9l%1nDgwpU}PQts&Ea|?3m55MiyKTZJk7}!H>Z-k*d(vv>KAoSh(dvS#pGoCC6l-Y3jmXn&Mg{N}oTf{<Z$X>!gTElNU^|wkpl^Ll`>?mkr2Wk^L2%6Yo+C+|lCVG_iQf=U?Lw$vY-7BS1C*<TacpH4Fw8q3ZZX!or6Fca7xsN|n*SImPIU8yY?&GIQYkU->G#e-l?&D8X%`xXb=6jH=CI{zEeyp_0C3(9m&}=*R@<Vr=)6F-jX@uIdXxYYejXxF#{1yelZ1fOyCTx7$y}0xnm_Om31O?wuk-wOfpO$h+oDF1(ltZTcL~1NKqy`6#@nxVfdO~9mQ;MNSw(Eh$I4|yXKN^F#zhCxBejXX_elpx=U5j}tjvF#Cw6X<UFv3?_hMPR&sGy`}H%5aLx13M%i~#>sjVong=dTxX{(52OugNpUb}6>0oJ7o=Px8D=A5VboI<SXWbA&<28<lOZ1;%h4pWm6}o~vqHw3HohQ6@ye<F|B?y{ETBmZ1sSC@=KEoJ80CC?iZ#M_ty_(aB~ZXEQNr-m}{iWKL?nel(Y>=IrX!qvkRduo7JRHQB4m=sVP8uPURT@U#B5o>Tn|VilS+Qj!(wZ(tuI#BT*wPz+vQ*=}fz?WuA%G{oIMDU?#H)Zf6KP$N5<tH2DNX?iy_!h4^@JF7EsyfQvM2dfqUzqinNbJl3I^WI8yjRxfAEK!47*bLDN6mBf!E|Y~U6O1v?%~XQ7)Qls-Xs0oVacqC8v`Y=1nbja?mul<<x@Ng$S7Frh5n4fLmud{!SFmmkZYHaxf{Hlw4iV>V8rVu5*a&4u0&zfXJDU6Sq0&x8{Wb8u@^wCvc52DaBF;_){YiT3ljcrE>ToHSXMyT)i{aT4^(&+MOMyK=O`DIR)>&(f!8M8anlKD1`jD1B)E`Pzb=Xy>IyyXMszj~V^`7@IQirqMb^1`;9_XTVIAjLA%UaQAHgoT4eJqcCgZGtYz`M)e>Vvz;)#32k+ow@`>#6nzws06ATxEOT3TWt#+S~DJZ(#469%EUy_sypj{z%OeXI_4Lo>*g$BL<_fXHRh>9YqT',
        'OSw}MMX@K9t?kKY7i7g0tR?Tf;*W4Cz&23@V+!k&%1KXt$+NF`bOQUO--nY<E1XeQ^b~R&RS2JpM6d^fcjlCKjMVg)Z5$Gs7ZcevQN*4>-rJ}(?S!~`n;pLtbof!L(pKPx$qD+Q**La2>&7e0?2jV3w3)O*qO7?ihc-L%3oA<IPa{NVBSmSxK0$6NC7V;NcD%T*e717M#eYGkIcLu-q0$t{65WXIur@_H@Bo01jw%c9I3G5z;lj#9EOAbC~DbT&a`J1P756~xS?#!u($ghvG)Vuhah)%`9zMG2B9-zme-gO-HF55L|?E(q{<tN0MIc9<ln7H-9ZA!rETXrhZrI^U?^+fk5?6`GiPY35=I6^OwQ@`tDrB!A=zeng}aVmJ-ouFFGk?HGR`CzvVaG}xZ;pkjwIKv)_UvEnb=B9fsnCogROr@%1?;Bvk@@Uok$Q|!<xnq02_+jLWd&(CR1IfU~^_~oFedLR^`C?+Cgui)8<JRv-UZbD9##v|JJ{GTR_z<#<pn}V(wgeW67suSzLIwMDTkcL!XP^yPfWCDGa9selZb8?1pQ;@AKAnLNg9E=X2mZ|oa^SbcZqjSs>^XPYYu@bVH~5h^yBfR?r^owXYm4W1?jHqj>CGe!-iK1W4@TAyI<d(y@RoRzy$<uLkAJ)oyk)tYlg@b(=-KEhc+2b<Q3b?bdrX87^DTv*c~2*x=;I=1hqm}L6LEV9osYb$S5cX`S22U$m)EjzfzEtMdO8_1pV7Jk&d!Mh3~XBP`gI9@If>22OebUJGg?=`Dmu}WflUjJA%pKukkR_c)E!4aQ=U89@`u5;T$pY7rg)EfQ4)1DY26lK^!v^UqrWGu<CIHU>yy+UN0PdyBsGXtGN1{{f;_QQ%OEFiG)WDDJoUhgvLH_+Wo#jUqa{-ob~0rlCsV3|`;ZV7M7n9f7lDqQV-srL0txDq3{k5RYPJ+a_duxmVGwFoC)8}GnXpdM1y2g`sR%2E@Kp&l1-?}_atup95qb}Vn(r&Ea2dkNUE_qGKXWva%hCMtZujo@!|%7l@9W|B<?#FbZs#O<{9`ojnnFKN=qAuxL-~n~8e&sj<Z~Kn+9i^k`J~Dyhnje%1N%^8`JCTJ)9!I-^L&d~Xbss*7)WdQZKwVQ0S;#};BYSV#82exuE(BRO_I2T!H2Q*6>I>D-*)O(cTUH@Q_DG{6ZzdL)51#FspVYJiS9uG7QgM(FD=NU@ykQ(T8UH4tcE;V6H}&tlxMMv&kmAIEN>j7cp}<|;hxw;q)*f}rufA?grJa%5g-ler%L-6Tw`AWM92jxaE%`;?PKt+@a;Pxp9}D=c<7W)w|&}|B%#`3i+D?}F(aNVy2h$l?PCMfHDx2AqO{PeeS8Ojz~@TTH6WI#5OK*ut>Qn=;6taR%k^oa%h2iIl7qFdIla(Ed(bRDs?>9?O3kfK|2C@AJyoZP^~J!h-kl63&Ztgnk05b_^Pvr;DRI8GK#9Yry7%3jq!ALS8v=cZ28;Y_4yDBTGO!f=U@3?f^+bN<C;B7*GM4)6-l@cnCU%IAmB2_@=*XV2rbx%;&jpKqtOD1vU60zAsYfmBdelPA<41+P8oSh2N}i56WA-O`-ld;9rRcB9sNa%OVUwx^>?IUpO|%BZoYNwk6$|vM9xAOeGmI%jk46jzr%@MW-nc3$*HITOD-lvhlDH><h!Qo!+fGGOV&kGI;;JjRRolPmzAAg&Jv!ku>dJA*s?=bt19Qt&2)+26V@K(5?n0Kh5ym>OhFpj5kWR^Xj1Fg?Eb(ItJx(BXeGy$zimb{pdYll}zjzdVTl-hmNvfC!ic2D}sf&5NAxk0he~hTeycuEQqrh6ixcNyml-~E1cB#P=<&j_UI`Ni5Vk#It7iX7hHdSCPA+k#?QYvk-f6~c1tI~<K!T!l+>ny`2ipsl7wOz{QLzNoTipL;-JTg5cdw>dg_p{`v?W68-Y!2PpoVTcu56*U3huyOU_R(aOgv<i$nSqe`PGdED&(VLc9KdYdBA&O!-~g6tD=cJIWcCW@m8O1F!IoSVjMG*5WmLp^s)&K@76$NB*{({gB38Mpa#X|$cU9`KePz2Uk%}0n`-j&AiO+$d;bhnb{F}2z6Qs9J4_=+_(&sCb<Wsf!?M$!V&fNMfh+KEp#R6Nu9a}@XCU#(pVt=Z%Q&FLCr}TIQ^uMz<7TEgj*c#H11~VIJz`(>{*&G_uV4>D<H3>aFA_hx1aca1l1QxBc6070X#=}JRE{)_}D(qb<>|HA4T`KHd8p*r#IWvAgYPd%z3ZG^Tw@ghK?dLm5JWL_@3csIk)X87_t+$m{xTZR5JC%O4l0Ir=`Wo4&q#vb{KKG#Y>g$JW0e|wV8tj^@1ly~=4}(khNL;$k+M>I39@tf}lj(yzOD<g}@WQ>p`H+FqeUSHqo%A^Dq_&MUu`kYvd~Z(lYt!Dl!|c4|RR53Av0ILxG_X2sk*`v>W7o$}dUt}{k8Myj_M+3p=4xbfN!ID#sn6xsl8#vWDqOq9IT6v!<KSwuyjC|S$T!$R_B_eNw)Dfq#%9;sk=g&y9-elX;&#om3I533UYt6U9Y2GMCuhv97xcaze&<CqeLLr=hUU#%yZ3jg3@&kl_os&Tc$#WxzMIpHfH!O``rZ#Ti|whQJ)CB}fsrgH<;scm2I58CpBmbI16kbTs>Pgyg3pmfi%Ki_d~|5b99Ig4-Dy>-Q&WfG?NY&~V(VoaINlRb@Nufkw2)74C;!dP_qUQyZ-IBNoPJL|Ws?%w*;bTI{Q2O^Tt3*A1HK>JO-JHxvQgXkP`u9EB`<a9#j_6Xrc&HZM&>UUiq~1!p`mzjs6WbU{%gO65IpCl(8cNLUlzSL=<GBVf36~KFQN0=boDPQ6Qo&-I&373a%f=(z6*(EtUTtHl%fGax8Mk!FkZIq&fqt%&u`uWuB=SV#8;QB(Uo1FW27Yp4>X(8IzvWFdMxMLnca@kz|VZnK9|qIuzimsuhUarhsdfN*xjoqgY6r6ox<2Y{Sox%b;@D;9$LWmIoBkxNx5rqau?&H4|J{rzVrZ!?Ypm)+o6LRcntDBUScQHF&lovXNwMOy6E~`0(=SxHvFar;>gJ!&nAx#HIQOH3Vw_|{_X@7Yo#8kYE22;s*1I-8n>^8NA`APB8#Rr`N+LdkKk}{4jx&{Yj%BtimXzN<x~_br|*YeXi>e;ADQy2T&Dc%yB$f<^6Tpl8s{_RHYL(in&(V!=nh0q&e0iv3~p26{^su<`UdB3o|4&fM%l$&fCvP3fx#F;U<~Pk^VvIX`j;7uA%z%28smIch@?+N9*a62<wSmb%8*@xUP~Q+a-zFIG=c!e;jHt5Pvpa1rX?ACCF;E26aCo%EN+OS>oIH%vn;V4>|>8G(1&QS*apiI11;iHt*^qTH5k8ms37Ktesy&|m52X!BQZkEK~3!Z&cx0@P2}`WS;{X)L_N`i{)rt>59+xNq;5WOz0Iz*ZltR4%V}!e#A65$z4c#&+!qX6WPcQ?*9*Np=au|bR`!%v_D5O&g?o5(^U4y@t!#Ld$3H5HZB4xc1XVn%!HKe|H1@!U{eikmDPaK*5+-DORd(sj-lZOpo{;T<yY#+z*V$3nPT5^D#rA0nthaqI#(wNnX|XB;Hl3r}_DGBRQLVwb3$mPxr+9SdW@q#BTSPMIxdr|W_Iy*Q>~-^Btz>&_)cC64FKK~!yD(DVs8rrPD_IJt+J~tMqRikuD_Qd1+KK#>DU>YdK02gMig|mPa{VX%zLU?R@1&=`6OeZYWNm^FFa5J|kG>P_5d>+;K-MPc43=Tz9(^aR^_>X4?}V;WApE`)8GR>h_MHg3?}V>XAmqLiVITo%<Z}>N9G%V%q~D^0sGU@fHGr$A|5ulCPv{_O(4vq*)KU>dEuLpQ5G{)584vDd3@@+FPXk2V5r~=?jF-_(W=2yJ**zQ}Y9feQ*dS^lgQ$fLqLyt$p!Lcm0*YA`7wohzT(A{=#u@}QW{))ENd6=m05v9$8`uUpa5arRK786T<CXkLR`!fn@+Wx)Pd}9JZua$FI4!aT#lmstjP;xN>;<~CU5$Pde(%N1ZLy&@0o9!IO8zP<d&(>OqdZ@kT#U4{--OT8UqMNLkC7HSMjEszWQ?>7j5Glw&7(;UYz%TY7-?Z+q&+$ZO(0Nt@}7Y`LN=te7(`Y|U~ZW+7KbMAK=KbHe?27shk=1EjDfysK!>AW)d@}J{(5vu6z-IGUrp&`q3+=buJByR-uibv<wbHReLE`up$SyJQ3hCr&RaHO=9I+Pk|;(l+cEQ&K|@L5XsG<O2C!uj!`rJE9LpF)pZ=|(r3`jAMCYvEm|FCG8u}5ep#~9sXK?$wV3$M{RK9p1!KQt)x|T3l)=)*;baYhy4D6fLd7OjkJg@rqg@ej>eVNz)Y<151jj64^PeVU~HPj$h-x=KgF7``eqmor$4WgV!Igvfe$?=%cy_P!u<V5!;Cy!%QRt=)m3qFxw@X2v`&|ir<@ApK1zbE@cwoy{`H8O3D=5n`(D%Yv`?Gb2BDs{idKpWfNbPpjOd_%-WXoO_AJ~Q@}Irc=-Bix+x%AWHIeHB5r%@ug~d_%-WNV*+#-`W0|H!A#P@c4$*C8Kh(09R(ddTP2qA{&p7v8iG!7=1I|Hi&(aD)g}eR>!OpUR?`lx;j=i9$%l@@UP(NoAI_mfr|LIWXRq^1w-lMXG0m5iR3Xvr?^4S09P<$AjpIf$X+(C+JKtu+aZHPp(b1~kI=zYGxg%bp)fq~(pm)kwE@-Hw?p>lb2ulm!#R<>FFus0jI%V6oh9>MDgt%NCj`7N{mV$!`;NZM26*4mzemP=sFNni#95I@UXk<{6SKXaf&CAt)0~4k=^;)Z2gIo-h!c@IQbi00;AkQr;v@s&)F#A<$mOFr(-6c-$PlMEAWls}oaTzeqY%Vt<wKm}fH<`ZahfX<k3tZqNQO8)M_G5Z>h(ge*RSk)ePq|`FX{FAOLDzFvg`FLyIwEkdi_c-*I%*A^_T2&y$HQtFG4QZU$P7KSL}NI`me;vZ4>-?NcvQhkhh)AnOm=?%d=v1a#-IURIk_63@o37E%zGr{=c8B4XM{_QrZFA)K&QNgxZ_AVlx(My*?@JW`)q7NWmVfghdSeJs|~glvaoUPuG1DlnH_Rjk|_C%v{LB3a#K;5!SLoYq)~n*nZC>OY)wpa@6E~HThQPM@|uSbB>sXIw&JH3aTyC(sgn4XCASDARyQnlS{~+RKu*SWP}DA10!rmN)lAFKV8^Ty)6~HfK@U4LXSNyC@JPZNogR`8&Q2Mi$lr^B%vS(l3?+c1D`wk^Db^4t(g6b*T7$t&A&J+{^EP!FXCjhKey|W+pe#-==*NB==Tpq-J;(=?7tRsJ*UN->5D19EapXSG2gyG-#EBe=Y4^`ogM%6N)pxe_<VlenDbAA@)klQA}8{v#W!mOd5gJ!C8`U0PR&`i=d7{##++q)&N+)4k4UV|XpTI*FL##rS^f21zAOq9zZlm`sKGS(VvHcd*`Wj)Ub<{|z-eQrfu|<YT=O=5qsWTzbUJ)Y;AwjsW!kx1Ai=nigChExyO7E7g>!hpp-wnAy%Emsb1$9lm3s@jx*yNn89sC88a!nE@9Y@E+tW^r(QHw23prsH&nITugraSJ>Unz^=ftO$?N6<;D7l54EEYFhwX5@W)!w#uxw|kpPE$6N-FC#*^j|WM53{%ZujaH%Y7X5m(|^f4zMQ=cSN!UH!!_QH|Nd2oJ3eI(+3_iR$d12E*h6;wB{aUAy^YB2qNB*IPJR)Jj>{x08TMX1^w-66!ue~(w#NBmYg>|+$KrjZ72)mTiM;8iH?rhcx!l3iyWP7_hu`DjXIGJW8kuw37y}rtXA=iRk`)2c9^X#jCr8=a)V`(>!<r#`VP=l-+oC+yS`^KfogsUyMS0#S4?yx6hyE1FpWio70)Tvo=P7*bzXWhs97Cokh86jX$zl^3wWbiJE%F!7u9UM*1GG@X^`h5F5;L=xhw|(u!Uru#@o><#(2Cx1r{{HBX-zcG8|KIbVLM&}xG;bV^<KNIZ1wi*cGK<aRtk?!O)RKz-bRklSS+z{T`Puu&Q7ut3m3%_3-9Y?TOxH){<inXiG?4x&OC0N61Pr^K6MiBf%(*8D6w$-$5crK^QkY2`PA>L#Tk*hCYFtY+&S?a@0?Q2VpLB@=$_7e3Z9OTJ)QLwJe`^B=}2KuXBoFniCd>dPe=Hk&RpHo5wfSVo`9z_lRX_N?CGqd8D*3&nB|LeZARJdy>PX8B>Ecc-7Sgc3txN7c{QVKZ(^X$Bhxs7cDE#UG1A_mYeqS*qNL`g<AL)9vX?*cR*OP`ifG{Z0@=%Nd@B<9L#2J+OW~Qf?=^^jf$XvNz5b!&HDU>5kG1cg@B7i1G8$7xW6Eeu(KM#4SthuA4&innAd57g(0fnlttXO)wXrJ%`Z>(lBG1eb4t~{w{Hp6Eoep^1Tkg~+Kl8~?d@7D8xyT^oNKxpKqT&*19HQLDrL-THQoGL7tCXVGnbNet8Uci}17LsL=FFEfxbChuQ}Wg4#+w{|oxzVQ>FONuDee16FYgY&k7TKsSg)gA{ecrZ%oh{&Fb}ZD`DJR(8|e}FM5aPpM#+@Zg5Th|z&rw<#Z+j^D48OnX%#nK=c}{MyrxeRP4(ob!5O<h=JF$SZvAb8Gj^0_<_L!@%Z4n=h0qk9hKqQ9`KRF``YgTz@_Ajv6pjyCMLjTZxQN>pFH_$cS;kAnXECi}69011Px_-pRG5$vOKecsuRGv}UoclEqLq304`;|kd&ZhR=Ln(r<SrDQ7K9=_Lg%3;BCW-_im#{>9vvIuB-y7k-;Li9g-v=Uj&Orp`r@n!w~6~ecWGgho{1yeAn0B?BIu@bXO=*^eH>c3FK(Z^p_Th$kKFdp*IJ%6w@goy57UP;y6+0G;GYR+?nRbK&Yl~27`c(6xe;P<CD5xY!GO>kxslU?PTtoGi5GH6ywF49g&PV*jX_f=^lWh^j&Q?GzZj~(#|k5`FozzN^brj6T;fm~=6wsSaT9?K7F~mU8Xh1z7Kvez(?aDwc)p8os>PBEH4=NC!FBO~ej@K8MPxUZi6s~641RutACY2~{THjGKaOZ|FmsMd`uUNlR4APF2jvgcIyED>GoAVwOL*TG2UYaP+!^XknGQ$5D(N39?Sh_oEJM91)4>Q>9sO?9c+JA1=T_q#vZ6y)v|fpHEh?ftf?Zm$%dD0Kez|OvM$ebz>aJRZvyda4g&yH7T;Db-rW~+RtB?z$x1%t6=7rG@4HQN}O8AO}Q6HJz%_xi(FN}_&TpLBXiz@KVIJgi;L2kq-gDY>SBA$jSVwU{FKFTTgX@t0x;Kv!vh(iGLL<KN2Cecw&xlbd+r363DIr3oMd-&@Hdgf%%GiOTw8FJ60*)^Ei22)<!?C~wnhV(S^q=B!+%r=;M$Y;CB4J4G6opb==`Ub?Er+9a$H)Wc<xNi;5ORiiMrSPT<j!=kexKaEsbDa-lYpclC^VFmd^`;DtP{?b@x!W|<xB}T@UFU|R6K~4k2!*(Y8z`Quf#P}N@hfo5H(V*Rz3uT~_O}0(Z_8UoK11Wn+1o%ZU7c@F25-k};uOobOxZ(re99iO<1Zl!`OEa)CGPlVm;I(Y5rI%&y$1vm#ztiXf<kO8mTPd~acA(98m7FqMyw6{V@D9UxB3%=cbNm>>40~419&%@gGO`EP#_z?JGD3QGLcRI-dXA!Op^n6H}vKP@b1l03;}p&sc11x4&dEssSDoLgJz6Q-ZndVd9}RDoy^o?-ro$es(Ct-reHj3rfn<e{mmq)q?fb^F{Dar+DJ+2?-p5wJ*7?4#4j}yw~g{r^K}r58|99dcl@aEoo=80sA$%$KVgaN7V24jRB(@7Z^P<=^!lC5g9PHVb;^S{s_J0L^0buHUu8C^HC4XW`<pe|i++em>Qp#es3Fhf3T(0_ebEn*3-61mk_o&xy<7-9L1EjEQv5zwicgo|Ke%c7PO$>-FGbTU(ey&Jf}fX6x-8@nTgWK@1>JAsq{}a7P@X<m-k=Y@%3VpXuc7zdUPDR90CWu{I}Mx?n51hc2|~bRKQEp)(B*Va8)&9)AZdRC33~&rzic4MZJ=$3=w@uBf*XnM5GAYq*iGaCL*ov|j^<#c0#C%jdR~$8VCXUlJZwesd*sJuGoRyS(}UYXzifnk+05NvHbTB^=I$>WAzwDXbY5J>bG&SNNc9wOk!QANGk2cYsL5qsmQVi2Va^H+KOGM8r^62u8h#oF-bCuCn{(h#I`C)xO1tb+_Br(~%-dx?&$>+d#-@NNa6VRA73OZ4=$mEA3l(Q?S3FT&L~ZcCfhP_E4m?lcWB-DQcR>9C!vZ9vLyP>yq',
        '%@~+ciNPf`zi7l&rbOAa;@%=4Kd8~tG1tUv-@zq3q@kRIp3T!FtJlF6Fc>y2skdzJK$jM1ma~lvF}1fkacmscq>LlBCsbM*fZ9IsY^c)vnIt+HF0!2UoOXuusoh`x93LMptoGhhvFHpQhgkWyvAHEudyqT3p$rPNgWa1ksCR`^7P1!cylAn96%6y5+J-R<=<D!jV$!s$U@JJEaZszWzU^EB7VsT<c^4c%Dcbms9y{ck~0rXmVSUqEr-L&3hB}hbgD_3h*xX;p*Vmd-=U2}M$dQgO|@u}+Y|mggX`iUeMz2#Iv~fAef~dV#Vl9+mqR@G`SA1A*1wGSan#nS{Q0U8PQD`ZM2txNMDmtW_~Z$jCsMio+b#Lo5kJypZ)dC!KC?#{*vh<q75A4GXW|II?SrRU`#`gqXRHi9vqu=%lAf%6@N97=j*x@38A{@PD~WS(Hba28Zvo<gu_-aWa0JlKz}O6o&A`|UjLpE<42;eBU~KkOriWf~p@UA3N2f>C5{(_?LX$#WSS~cFzvsl~uyUbd;09-<E*EmjDk*P+<rdPqK2X-`8T`V^UJ5y5yFQEhyzevxqGBC1Z3V*<aHK?7N<<S|!@tvEUY(^YF^E-p#HvJBV(!9~5WXuR9k>!}byp(NU5S<NO7J}iVS5sMSK_jNS7N5S64JaYA!JwL7iR5!mYN+$aD5uC#-QFFt_H)^Q0>PfBe<Sv#UHC>E%<<|3$CZ}D&5o|;Oc^%vA9rI>stj}y}W4ZxEV(UT&Y|u`;<8yfyY8{U68#^)pc-^-lnD-aKs#~Ah<O5A9feqb&B9htAEa*2n!Lxg$=MKg6nw`+!hx+?MZirh-{hJ1`G65#%=6PClIp>DxhjjNxpXhG2`8h6G5`8lp8fcUKGdPn$n7nJ#@kHZl>wAWY=Cs7ISEN71i__*egcCwEZmRd;(d_m7T?u=Chb9JBukNki`^s7ISVsiz)0Z=G=W2bERi7rTHvoDfnH#S<DVVnDi{Bkh7R`_gPFKXEA4X7E_wfV$R)XG0O<r(sP6@XEEpQvzRM6iz(e_F@>DPoV(9ruH-DHbf3i(au)L!E%5z3Sl9;(dodRF!J@}tIB&gkl)=J&+|E&XK}b$C*lOEo)Qqqc2m9d7I(_!E!NWdy*w2-R9p_Eok$K$*4|@wN=Jj~kOEJm+h!$fmTGT~|cOya^JnO7biOx4phk3BM<@rG|KPHwuaUYvmTmRi=Y-7#r5xOHU%Xh;jJrhT`;Z2caLm|h8LXQoF9UBU}mTkQPt-GLny{0ElWz@rMhhVr4!==1txZE!}k2LLxei)wU@I<+uC^2VGPua8UGll{ZSmR=OqDl7R!RJrGdCepkmB(Igd>e0V_hkOa_;{IasjfO$GB5qpw^)a7F>nw>9K`S~PCd(YGqM|{vKz$8w~#B}r{Z6vJ+9;X7s*^hhyKOI@Gpjc(UgC2epR$n%Uz5rntvy-s%X51bA^X7_v2xl&I^8-d%?dRUhXgNz6p$f6PMmOIAOky*WPwLfutV(NeBL{!+^~YyDEs2F2)J-9P4};Z@y)O!PQxW&pGfX{c0<GcBhv_U*sX(hLb*V2{tbH?}dq-T%Fj-)tAi`PnC88^Ie?i&*DUXAtv?%F|pr=mr}N$vQR2oOWRDkD$6M;%-cL_Wfq0^ev~5Qdfq7yK*-gQ=KMy7Y4y05uX5j;g83bQ^z$i>-&XpINu08jQ)%DG?MXRR^)pgI-)B|OeVLi3k(nu&nPJ8v7FGq6sNsB}9%q-C`BW_&U+CfZLJ!9mau60ZW3((|RFHg_R`I#J_pyAAuCI7wcb;+Q1AzQ#BE%{K1+EnXE?vZcHQ`nKAa&^?NaYcW6eTLzh@$5*QFJL*j@wEz;N5o9k+@w_-vH-|smblayA3A-mSu3e@VQ;CMwVpYSL)(d&XgGRyn13j4bZ~o81bDj(b3Bj9lflm$_~&1X(ur6go%D9C`!?1-zO)`0FyE0z@PK!By$RZ97;6HeRP0fjwvE$d5Q?+#POUds?l^#W5|?x`JA~>MyiE;jt@C=Ca99DX(*=*BuRNB$;gqML0|9vK$2{Z(Hj}RgnjrUW0iD@DKZ|~DHxFCjKY)tF%C_H#XW1MgyzLl7}_apj6)yw^`V_|;`Q|_$;FF*9Ad)X55M0IKOS$yA8#z)GGfZ}Ww;({oTuNNz}ew}zQqH5kB0%M84Pe|gr@nue2Tb7IF0|Z(*tw>^b6&xJ9h}txpKYyR_n<3g_~j>d7^tWvAvsf$BLaRsgH~7@nksK`cgM6BIxKHVBoJ|pijCqbzTvAKRy=ni~pS;>JfcwV-|t!(t*W+&%Qa)Ee86usn?vB&a(x^gCZmk{s(Fzldl^_h*kTaD(!;~$Uc~k%0%|1Cbr+Bpz^#pXP`di8rM&mego3DOV1tw7V1rze*e+PT8B!lYka{oK+uV`&A7?#xP};jYf&9V?m^SPK;{iZ*>?8cIL?ajj`l=8h@&%ZbjFR&IJh%z&LJKBsK(C&)i_X%R8%7q2NC!kW1$%5G7$<9=uFd=l^m$XfoeP%GWGN7e|;JdrvY&qsK#?bHGV1&aeCFEPOE+{?CR%2uYNA%>gU3)ezK<=$y1K(DQom6unDs^c4Tktq(TN8+@x<q<!@BFMfQ~G7u->29N9CD<QZ%1u_sKw+`yx*V6hstOQ?2<)qFrI38;2SA7LZ)U(g`33wcY({>BQ{s{L7(r|C*ORcUnS{I_f|BOlMt$Kptup32$B1|VPn2|A$Rq71}iA+PX4UXg|R#UgpWkv!j|Qq|!4D$n(`^n=e}BUW~G6vCGtxX8+Io0VXDyY8H6tR7cdDYmyC2%=bBi)^5Wc<3P_iHf0^EsmIYpoa!0+CUGTCV~Dy51paH;2=&Q5hoBEC#`#*H$c}@PK@Dc#0kW&(j#%gV?7eYp7b(UaKWfh-~i+;h;YG3)(8>0Mu^ZgLR5P(xIh42xS$?I5v;l>54GTeXGxOL30cY+@`CKnj|He_XaT)&Plpz?-<!{yav5krk~`pHU|P=!;aL-2&@;f~v>7w{gBKJ+7Zidk^b9Xx6Jqe5+f#uT#ENi*0WUZOc)>E@1p{8t6TIMTCD#Ty$$%G}I=tZgB-i>PC;2qQSPF`<@JbKDE<LaZc;;d(nwNKG8$^aUAO-_sK!+GaI>aE-AqJ5SF%UMyK-dri;Xn)q^ZqI5FPxv&!0zl+MK<6C7lV83yok64H^3;}zFE9q&Jt);2*xwY<3*b9gQEMsh&n=<Vrxx(hce~on&sV6#nF1@-BJbWALZROSKh_Sx)oAF<@Fo_ok_k_&Fj_pcfk5)>;dch%3?5@G+T6y(Xc<Z|5Pp|WM5#z2!q*TFq#yX4xf-c+#uClW+A!sKkWs!LG#dIdm*(6?kjBqFzx2!e3LfPcu3$W<UlsL2|iX==dF-?MQ#Gk8V9MC9TBp<<R-YQF3uZTuQ|}u+ek<tWP8a?kk4lzc}k(zbwG>uzlG<clk5E91cgX_@ar33MNX01)u*ZH%NQR=Q}ZI1)3D_--cpv{k6gycWsF?LY34GlIBktSfo;c8EanWFIh=yoqGiYZYEZH-Vs?A{YB1z22BS%{Mb{V&Q0(nLWk1o0P|rHTV73^HCdGt(+Wwj#WIFLRDHaKBJH{a@PP8^HqCL^-w21yht4*3My2hxx@iUjl)@jKbG@j_3tPgOa6((6HTBAG#ydhqTmesJIXtNTN>NYid-MI@N@VeC|%@$o_H283KuRCwh=ym5@IKb;xm^9W6tli+Rz#mKwo~SUmJq!$Lx@|9bdNd5{WX>7Kti9lgy>x#nVxM?b7__VtD&e!}?zdM%H-F-o9=iDx$IQ@G&zNE8x`)=fXDsW49lD5|;#FbLygN7St;i|HK5@h@>da9vU#t;0#e#ie>{^JxDNbF~D{_jH6<*JCJxaXYm3Vi3g1WQ`IRc)?XV3!W@Y$o?<se$8RH?K>+@3z;<A!XZ9vOtnE7ff!WP5>aun}nH8(18<7+@SO7{}Ulthvu^<8WFO))S!^B6uSSbY<Gc;j}31TR6}y3Un(0UlX@zUk>{!s;*9iQ8yo~v4X!f{a@hVk7X~HQ^Of`7f0q7Ns_U?n2_z|5-1z+zoi`f&+frU>c4E-Ljpt09>n;&w9^vWKZNxMyz?>#@6^FM@9HxT4_N1bbq-i(Td+<9GJ=A@C>9B@P7E@Ffxs|WCk7coL0}-P6M~FjATSWt2|-3c5Eu&UgdihH!W<71;b=JSuHo38y(SkZgj}FNLzfw;L3v^;Yk>l;L2aa9N>h%GVsi(@W?zaTNVS}YknIJjkd&FJl}Ob=zsc&Ejw)4|Diu;0K<g_`i}acjQl;9+s$QILoB7E49)}3f)29C;tGdr+Ro~8{UN=GX$8Wr&x6<<4bl6$cx@|gO<6Xs8T6z|>YR(01)5+l;bak!@9vUCJK}W#&PzOKm<*VGiyj|hCX`;e4?&fhfkGr|b-CV801re8kP0<Nbf>b1%7^1XM;Yxt_P?`gXxCCs9PA_l^6|VFGZ;mj%;h+i^L|l@hZJW!msV?09j7V;Z*7O|H!29{Oa#rHIK*#Kib<En+Dl5I6GT)GDAG$eff;C0w`%y#cpucY+@qOZ9Rkg%7)SdxC;uQ~_{hy^1-`|4Md75&;>C?M!n$SB?MdWEfGzUa;Ks3(_q8a$xyS)=gza)fsUIP8VSZ@i!2#I$@;2#*`d4X;T{|J>+{^r+y$_KEMA!C7w8bWDky+PNL>!Ie<?0)m$o%ZY(h<}8EI7HjtpKRG6hC$=?K^URZ;<w}eh#Ys`Ye#vHC$`5hRrJ{>@@qWN{e_9`u}}3p_KE!#Pi${tqWkNUP4#cUH4DHs!*O{^%FUYoVOo$U92%7k6Cr&Z8kIw%@@{BU4vos&p;0+BDu;IB&`uoMi9ag3Qd~Q+EeJAQ#Rb%5s;RU_ui#dt<N#G%KwTz!1(&YNR6uFXq2f}lf?HEZsY;)MxReR`6<i;251~(i*rrgpg6o9rLyZdlJ_T_lOgs=|y4*)ST_(C{Pb4e4&~=&MOlF?pl&u0p64|LK5hJ()3tQ-EBPr%%(;^>qX<dMgy$i5|E)@@YBp>ux@1RHWLHBVt2i9e-Iq0{I_ZrZN<3X}~uZTftgSUrAmtU+2#K`@2fs@1cPUg-VADax{_r~$%0p=pAVH3H+TApS}gX6hwfV{pQ@l?Gp<|B$?OI+ZhBtI(zbFM;Q@AG;CLUwgA`n*P;*XZ*aeO{x_YxH@IKCjmMysXC$q#iPs$&W!MR<n|i)#wUT^3^6=BohOvhv;Nt1uFRpldeQ?pdXm}C$O<fzV`^SZwbN(iFZWc9~k0!sZhyx-a1|=n@<3XVv&Hh!&<E|+>4b=(+GQ^pzT1yPdeHz500bJ1A=Hf1=NR#wsVS5quaZ^6G*=#gm_*8{lNIqcIFY5Lfdf;Oj)!Y<JuwKIvfh$W)SgwJT6GejZpZuL&{A@#Djz=2!$_*h&NDZACHvViBv(icY7z0en|-Nyaf7zvEC9j5b^pE2L6E|o)_qr@Q=`dh<~H8ZP^c!qd{^sNRCd8<j8tdtSMni(t2X2)u)Put%;nnp6E&SiCxy3>Xo%7cEWmMr_?8UQOiRO;7J`(l3Noy{WeuhfKTMq)<jRfP3-jdR4@HKu@hSpJM}ivli#efYz;9fLGZ*5)2E7wvWXl5PxK&tVyDTbdTFwW9RN@45PhO2$tGKq`LeIxd+O_2KWgu(uWS9NP39m!8rFO2Yg#{Q@2Rh8{isdmdYq{B3!Bz`g+8&6yMdt&^3~O?XOl96whY)ozWqs7=oCkVJ(7Un+C&HWs`(_AJ+mtI44b?o)FTOC-_!QpP9%KMC?@iFAdpY8k_lf@E4Z>h$x7ZctM6iE@<8yIIY2g-)z!kvB43I1XqzJBQLjGtjr2)svbn76PqLEt%!<6lSk#ocfVOCx;ygq+h6qP{A{^E!JKK~3-6BoKxGuGjbC5M>U5#4E%BJLm;|ID$8qJkrh=D)r!2U4&S<^oZYsya8X${IXKKYo+VrWA#o<D2)hcPOnSeg5WG)mXfD7!wz(6oKa1RX3kjhwzC*e!vCWq?4Ht5Q-f9qQfb8_L{LXK1P(QnM??+$Usvfo-r6Ny88M=X%<4!4uQ27ItFVaj_HAt`>A++Hs*1EN22s%a*FDPvollRIR8!kxMl^Ni*zw%#|pUef8F<TDf&1ms?ri7rjr*RN^LbC2p!#h?~fTxD>weCQG^BGX3xF$Rz8l3)I89A(0)hRwRL6f%SVPS-)<8WIgq-B=ui5?IGy}%$^_=;MGaiPg$t-2Ss_zQIzYc$v-whP5$d;^#XLi_${jYzZQ&iamiaevVduNeYd+t7M)2$AIgF-IxRH?Oc#XQXu%pp7p+Cpj%$IMJZ^4~xIKEkL1ove0Y#r1q9nYhV^EW47I@>9*oNrk);D9&zP|OS$#Xkydp+!3<e3FFj8M0}wW!GxZ^-uG4HlK1b`g7{6FVS6Z#12nJZ9Dxc%x%_f$eW}qE|FRZ!}0cEa?rGSPlCd9W(0-ywMuqdKGH&?G3V~Jd7E6qqnb}6?>!C%}ua^Z#18pd`&N~`HfyP>uZ0bx9_nPe4{~X@|f9i+v{O}qu2BTo8M?vgIWb@@`%LI28)&}6CquM-sp<d<O06Y5s9OdS58FY=+p~HrzVd`9Bmz;(g(7PQRf8Jqb6UGI9hMea%0+QPb2n5SEVM8&>OuXadhh1S&=w;O?#khYVs9{qs?PfIYEXI>U^MD)Z`I~qYV}<SEl`qj##d{A~ksgbKMb%qm!w-h{VwXz0n})a75y0!QSWy=DH&iM^~mMUy(RkZ_sjM+TZ9E%XL?!CSRdk_lm^PscUCN;^+vz(I7SXip0?od!tt<*S#Wfbg~)T$D_7S{TE&G+G=q7UzvxVv1cE8#-MrVu)X9<Vd`bL-27XY|C({`*|*x;{#W*)XAGK$p0Q^iI&4r#%y#SMUs|zKtTAXg3Mp(Vbc$no=;l+mriTvOOXg&39Z88RL!fi4vuHXF33Qi0M|sVR-hA{#GNZ>W3XK~JI?WMylT(Lf1>R(hJ^MwJF=)PsV0+fny#0+YvBSIqk8<*$q>e-GF3Op--$gl#_PYqT$U4m{=3z(RZ5HS>XD%Z9X^xPGy#jA@x}5`fn?b6(+;%dbW|cwUal<W|?;`uj%icu}lagqRxNU-8O@THBEgBhh7L46#aYd~>eziQr6EtYz5mGt%=E2wES5u(sK#L}7oQUm36=tmpzgiy31R69E0<G9y)WC!)@T*1SUQ}7MZz_b{X>mmkJbtwmxfeAC?VHNZX%)$P5hUG;$i1ktXx~)!y%@nan$NGcBKM-kpnX%>_hPQ>vq9J*7qy#hV>N;cl+a@)x2toTi;;G9%$`657YN})+8(!=dWH6%Km;5!b-t0lenKT2S?L5O7(tRa=m>6*E0WzU&*peyb6VBo@Iu)E*qle)Cx*?D#b!fba}-572U;x-DvYO(%>l&bJi<AiHW*v7oI?4ux2?z)H~&@UznIhVSLy#^9v=o%xdPVj3aUg$^7yWx?|TWX-;eLR!FZ!!z1mZ;W0&L|AmkIIKn@t|p)0V0tFOLF8zEM`)RcI(zVc7rpe4r}bpLK=;10U?12GBQK}O~hvqb_GikDpFsN6wb;*#JWW3WwHtzk{)4$2!eg1I?VC{}irqjLux`vlT22_c@BKtC|nTY@k`;vG?#JIHxnpj*N}Lh%3yke)}hE1F`FK!t)eg=u&MSlmHa!<DKjI{3#JY?Dro(cANw;SQRSp2r#EB~V+mc7^eH2`C4q{t2XC5<)yLfqr1Dw*+B?#5*GJ4-E0VK(~Z{gvu#@%QAiR%b$2$REiTlvEBTsVxT>dUFeDK<WKA{d#V>^PizNzV!QYgJ;?URF?pUc%40jR-GQm1Lq3t8+KKKAOl)_2s^^YR?1y$@y8;v48Ba+f+~rXoL%V)x*SD)(k6?&V<!q-#TH68!`1Q^@u*I|av!;KTRMhiETbj%zG@3yH-+Jn<p!8=c0UaURbLFS4Fx|j5xQz>CNh%9lqd5;07qY!%dgDW_gFeLAyAosf^@9YqL5_?wsMO6r=rj4CEAp7VYD7v!4aRz!laZ_uxi+pjnQIE#KpbCqbjJHu92^>!KWZ)YyK>NX>L=E7&0#&$*v=36f>RY;$k|?TG6Zu+!;RF$b`~V+4SvmyU=!L>C#oxAJFDVGGTN219jk_FPJvKNo9K<=1R7<r3gJd-V>?SR>Jt<<j5=Ft+(=byXGPpdYJ_=N9%sB&h1YAgmpD|n(vkE1m1|G=BWj$r;)qtr_dhecUjyp#?2HZ)%VC~;shpiLwprSVZLsLJUi6*Y0Y77}i(KrL@2$ISqPK3uUnBk+@z;pIM*KD6uMvM4<F9JHb<p',
        't>?D!?VSL(`13{l$XtxJ4NRB#kJf8dpqUf>pb>(UFnIl}aYgL><rJ4fpF)<IsfU|zCdmxLHDiDUKFiR!>uZwdSQBB6mO?}%ge)=Bof*moSzTL(o3vJF17mO8C8hbDLr%3HTk!Dmq`eaS$DE>!TDirt%p2aJ$H4n~a`wTj_Z?#}X`SS3Wk*)w`)+)B2Yy|w)Wy%G2#8UkGV&|!K<7QQeybSC68O{p(bg4PweSTzi&<v+1X!h^MEM9-Q$YMP^_*+NaTT0b2$Wx=Mbh~w@5tdV{?=mZLO0tNdpxrjD}2DxqkecuAQy(XGJR+}#-dC?A<A7xV>_=Iw3eowO^56xe%+dx#<`pX=?<I}s{yT`-t;qbF@Dt~GMcRx6l2dDDjRIbLUygf!;W^ZG2zd6L<GU>2r7$auT{I7DGJKz>w27_UY4vQFNomgXj`{El?|3N!MmqF7(`tHX*<_Qff9ADNCXGCve2F3HB8As@Peyp%8#9nYQlv_s_=tC6PS2F|0FoS~oe&!LHEPgxg56E%Xu8Kg#vA{M+IPU5Z5{|om%7{J2aVIJq5Z6$k{7GONbRSN`2u&8h9rt@KU6Cej6W9hhbVZuHS6~^;C+38l=nI&^Mws0?<i>MS#WA_;HfY!PgSxc?p4F;cuq{xxZZHjBlz5fE8-dNR=6ggcyh<2bCd9Dzya{fpn_w_0cE_Z+0cnnSBYi0r*_w^-H^7RVF0I=J7(~jEEo>C`rCf3UXR&slE7pE~w|n<=_~Av{&!cENingO@d+tSB8yWCMm+PC|T)}11VbL%~%%J&S<=OxPZqa2h7{=(Zh-nj9c7V#DWe4cAiRcDR2kE<i@`M#_7jqRy)4gHJtLqIav*2mBa{;P2))rIM>YDB4R@c&2t6E;WLD3J<usKvX=-@T@5*TdR8!CgrZBhP(HU3ocxCi=q59|X^KkqI9SnGKY^aIa6?-h@Gpr7}^KJfJOZawagqi8#dwx(hD2t>htULmIWh#k@NPr}|aCNS1ejXQbOIF42mxPlskLaq-4DPC}ZK8ep@>Y(2b)wrP=*S;uQ$Fg>`gO1D44H(ujVOUOy;i@vot*92?0Sa;{E``~_)>a;Tkr_G63|0eLu~6;Bf2A63sXGAk)gsxRHiP948Y5&61M4hQd(S)I!%&f|4iIZgY2^k`T48g=NCV2WY_1(iDYY&`OPIyBbUI3?S2+M62KsxWo^W+BKz{@DH$Z=bsbw&=45pUB)Y3Xr%g*>lg9<%-fh@X;o^gl@g&u8Loi1k;%$j$>KgRkhdZ0oNT}4kJ;S~sDG^>K?t`-!SZv=TLfcZv?hB2~$`T4D51?Folx(xb&`D%;SdyEFmS6tmj9DRO>E`vdMH6Qbh321mxEE04!jYD*)=vk4%Mp;^mgcn=smLjxilZd3L=o#>IMp@cG7-O)`T6eV=Hn;X_F_?5%bT!Xu3|d}x)#eVkMMm?SYK#tx{^mI|pjEV2H?4V&VlbE|%~o;t>uj1WnziF=sL_|wnXaQoZylpiT-(*pW0i$A>P!c}&~!O9dY>4!{R-wbxNF8}f6AB*?pBz$=r=d1!QE^y7{*9yaM#>pb}x<D;I6vG%+4;i!Tp)?=Z(&oW$$T$%R?MYsPpJjlrj%JW6wTx7=~lRx9343X286C1tVLSHl>`4?7@ibLAI)y^OIOKGBM`Y$X>zKg;cUMH$?0dV^X}x7r_9!5tUW7ui%QEVoW&~*@G+lCGnx<z)!M$NO*eH&?!bx*b&mJuK1FyE`X)&3$Mt<kEu6f6JC8Q$k~G{wg=gYo6b*y-Ta(GC)f(EE)=r)TQOe}pQx(M+gg!}U*n>)CcOGqkbUt7_NFGeM_AOfX7<cRulcecf?h+=YY2LsnxL2Us9148#!8Wp2N>uhbfFpFQeI};3xaJMY=nV6L_Z0duvyjU%ydHMdDJR^GZ+$u%h{?OV4#nX!r3|lXIX16P*^5sHf`nXUqi#Sp&KD)P!M+x93g1Ya4!QcvDn|`i7>!GA7NmylM0#bRLK0x`Vj{D5Cw3*K!hb`P!M+x9HGhL2CGv(q#mNVLo}BfQV-GG3Zl8zqk=`i6Y_NH4H81?33|HG5O^ks)H6M#p4lPwTrZ@a*&+4J4yk8)NIknV$%Nz+!yqN3o~wn_6O<m)+GB*&bR`JN8WIx841<IKc&3Nc+*CJM6k-^pgw%7zka{ME)H6M#p4lPwTq~rWu=JSTpgRDb!1S0#)AbGN**FiT<RMt!o?yLo6+tRz1+4~Xk;Y%G7$Rtq$J$+Pj5R+D2GcYDC4l<eY0<vCb}n_NMbKcVImVhFhP8^CR<b}M%~x5pFE7sPUqj8nyNvu<(?3iK_;~BrwjUxp-HYsSgbYE!G@lM%Bf5S}@v9>0Tpj9>A)0ug&{I9}OP{^KYba|<Te1<+K8~`L&t9^{|3m({W~uuuWOW~;FVolIZ&&x>>1_BK{8`Qf*5z5M_uH<u!?U;2bx&rpY6w@GDMIR3k=&8#>k2mSOsWozrzx2AebM{0Our~gdn=Of1LVB<+C@`EnwzYeWqI8qX+~#=XK$^jtAyBFo6wQc=ssb}TJ<Do{e1$4zXaIJtHZqqR=89S?DX^o|5E;vDX08pQ~bJbUGWQ46RF993Nn?m6~78FmGjj^=q#vo#jlzys3=o8T}>p#?XhhAsh7^gi)Bo(do0_dUlQxMY4I%wKrzPcS#KMZ1OZ9rU+uOBv2~4QCbrvLr|o+&z4p89kGXj<rZ>oLbH8oeW0{#HcEIMvm|pzt9k6waWhS<LkhZVI^x7v5QWJ__+r0iDQ#o5rB+1n))fB%lrg8yQIa6)AAXB+Cq4<^B!fZ8><R=sGjh55`o8RcTeFd-J8y)q}zS82l`IN2L8@(zXef&lJWe_X&Mz4zZFF0P*4L`zUaeIs_^hU?b4zRz`YkCLR-{=*4^c8MI(lI%%;#cO%+_p$o<y5n_sNZP5;+H^iG+Rw1W09=NSwQiNtjZ~na<=T)Y9gsOdPU-B)%zD*p}FO<vcJ(Q5J#u(#TAL8lW+72#L?!!JQ!kqyT@gIqgN!3)<60dFKz1<TY)&*5SRyZ^5EEtU*M}R_BcA7>t-uw^A*47%GtR%56N|lJ&vxd_$Br@Iz2+Va&|6rW6pK2NE{u(H#%au?uf+E5qqOoDA&CradgDq=oQO#M<|XK?2TTbT=$B^(Gh#2!LWA);%EWi=$x~QbF9g3k5Hv-qfst@F{o{LPvdX)W{Qe5zy2Pfn^EKo&C8p+IPJ@uxi`(Lo5^R|R~Ml@O0F(IuaH?@Get#abrI^Q<nnHUX!drhl!%^X$p@Ua-6p13F1y?GmqS2ywwm@R0qT>?`kD{(WqNtdhxrP-Z4!q$s&}^bFt3^QHJ|2HeLE;mb4;(VYHM4Tck6l(=xJWn7ldJXH|$*@PqQinG2fQCc)@;KX2LG>Wof=0R_rt<69?dF2JVN<X|^A}OZ|moxGl}6Ik&vfX%_!;Z&+UYY0fOK{WPzbqrHNU^9&s4OUToX$m0wWnMKso7W?CT37zH@^RyAy{l^Wq!L;U)(x%g5usIX376W=iJ)0K8(QR!P{_`|K1&cvzTFxaQTd#G=3~9xMn^7@xI`awA>cnO~fjAR117tCu<a@<vBD`4pp%FZY6~J9m+Ie+McFV20ppHvvV=~5cavQY4kWFqQ$XhHG$!%QxMK+CfPm^xI9(OSuEl^5(CAVQi>gFq1xEU2Ar}GzaDXLJqAu<aYpT9`#Z?T4KK}dT#B5$!krx=vdUV*n*pi>+%j~bEEer6ZAZJneX<R7j{fC|#R<)hQ&P6&d;r^aSHg|94Fglb-oCP>Qj?N7$N$fihEHGLspCFxrPXy=##tyyT9U~!h1ug>F>lxtxA#w`r;k_vP$O7;mC%+IxLI1pSNk$X`UEeI(wR+pq4+}&co7d7y2o-o3S+=~(1ixGRHBXTcBY%k`JlN+q3a*=~M{#qV1{1+9DX}ijos-YSFm8#-p_*bfmoNkq_(*MQ0ej?BZxdcA&CA5H-nW3wS;fA#}PxumAz{I)>dCNrh3SN=7Ol~dgOSnR(cm=-k<-f|-H^@ba>;mj!i{y)UCGXNp^2NK7H|Zt$;=RI7<%*ml?28EdVvFRPcV%B}m*f>(+4tHdc~f6O6TU+BVEorbc@2^DvMErPO0>@?rRqYaf#*`RMmj06P-{aG&hu|j$nI-X$g*jb*?6`{rjce-+^%FAX}0>x>Y$PSLBT6I3f{+eyLS(VpN)$5(GR~YOT|mEs3iQ+%X9^^#qw0VYUwnEns+doG+T6yu_zU9-k_15Ah$hq7K70wHbxz!X=8K0ImF;H>9A-RBWBS2uX39^;1*p5gJFygi`IM0z5`SSEju80i<u9S&R{T4nnNK}g65A+mtTU4mrEg3f{J&h)BLgN@-dpf0#v+nO)B1*PV<)^qYkjRvAN$IVsM#sSTu|gGid%-xy>DLi!Ot~Fh++(>pf=Q0V;!*9gw@l%m+zl&~%W#=b0W5MN!!#LNd#%OgSPXibX94i#A&uFoW5vShc)%gSq83yg#bdHQN(xfJU(-MY}m<F-0+J6Q5AfTObDY4Vh}10&P|?&@F071q8Wm^aiu+k*ejj8_X@QfsRtOx@LP39iWK}R>kWCSX`5m4i=S`NuxSpH~8FXaqkxc7f0^Yn1RL&iYH%e9HA?_Yn?J;FNh}<YaL;r4^iM!h<SvVLBUHa<`IGxkyMr5$DNUPadr<d&`0RX0$Hbw*b9_6VP*#&gf$eXt{>P21y8}6M`*IRLEQWD>E<w#(D{WlFD;u6DDX%U=1TlR=AsbJ3PERRvewmHyc#}(2*$|xNf0(0eF@A8Gi#L7K(IZ?*JhtvU?$NjQW3aaGzQsZHWk@%Xp~QRgQ+C>UB0j4xSR$~i-Nsh&=5h3putdcj5R+D_UO>rFpkWxQG(N8cC2FHFC%}}^beD|I^Up!k+~rzRm`tuGpRz64(7FDkxt%T;JLwbJ%T>A8%(F!<WhV2>?MH&AMyv?(g%TMAR^gbD)?Kd(P@P2VI3>n4wwkC{jxi1xCx=m$j$av(4IX*$L_&f|5b3jfx0Fa-~=tyo-ZFCHOA=8>Pn2&I6qgiz35_LHIIC&ulEA)nq<Gr?4&iyG0N@}7IN+f?Ze&Zo1JU-&El@;&%Vkz*H@YEp}h5%m<Sz@i8|SYNj%%|Z^imEYW7fK{Ta#D9(7?z-JhYAB>*MmiuO=e?awe{cy$Z>A8DmOp2V|Jv&1Ly92JVp30WiE5Jk-t%9D7`Oe~|k1*|<);yEg0sLeB@N<2q~BD3KjZC|uYPD(FWxew>+63<bq(nwLBn@xMGzAP+V;yI$aZncExR=O%IDDe#G%tlBQR@gF<8sWKoilsuKS>YUl6zE4T?xhP7Y+_zMc|JrEEf`_BD6}5lNVhDtd3h1KNx^t(tI^6n6V<&3R+vrb8LM#tvN=;@mNWH7x+RWu%SgA3bjwJ$xYI4R#6vKoGJ-3E<4E8bQbJh%22YM`#IFHQ;E2_6%H|L@B4S4ZN1*J8Slv(i;(N0r^!$zKAeB$I(6cS<WJ{oDTax)8cD99{YzgdaOQ@J`p=Vpz$(BIRwuF&x=_TC)X%z}z2MN@N6HEt5V68$l-Evt?X`r5L5nC-pT<Mkw%z&*Y!qz)s3s<ebq+6C;y2Z{I{?G(t_(;Eu^vg)UjP%P$zl`+DNWU1<FVz^s1z`*qj4@oc_f%eWR!LlU6vptZq`93jToA@^!5G6u%R_lI#&CtwwSB`ujNyVXhHsB?r6W?;yk7AeCRa8hb;S^=D}+d0D<V=?43U~!``z}(WFk`6#FM%<fk<63MCyt?`nk1F93-4}bRbd(B6T2AdqbpZHa1xM+uj@6`c6_QQ4(WBDyUbfA&T<lfvSJ%a^k8yd%=JBa?6|X{OrqXuP5rawrRhqPM=}B>7=7A`z;}Etz;e@0EYo^7yyU<0S>;`oBnvy-|j$vyp#R~`|+l~-GTmihl;N^{qbhM-GTpj2kYy7$K?&zkVKjI$}pCw?dpY3x86Q<m<Ei+rk5LfX2!a*mTLtyL}liGh3vtF%6yjp#Oi9DR_x4E2(wawWeYVHE-?EZoE=<m_MyXMZ!E&Z+|V;K#txln7#cb{^S`q1LF8J9(4jNu^R+T(sZ{xG6TIOAUp(-|DnLgS-f%(asRg5_R@Bh#|IEI-S30$T;${!2IwOf%V8y5fRt2>H)QGhWb;g{+OTZpR4w@OcnJn`XvWL(er*cgQ-Ej)aewBcpL+*kA1~K;1%+T$8#{wBd2#>Fl%?beFv84L1tnHTmQy#0Gav`7cSnZSx`IJYn#Rxg&m(>pWN<QS5)z0`zKI51FI^x#=Y4~Vb(>ATy-?iL-GC#fby3($ouC%LcSK6C>fBQkR#G+Z+9hnF$vFER_Q*;H*_zLO7lwVpW*Ieu5<GUTuKKa-LdzDcj+iWps@Kv*;dx%G9pHxdp>DwpOCe0RIV?0v(q(R!(vOR$INrg#ljJiUs=tmel799qIVT_<f^UN}*CEx~w%VN-v(P1!H?=kZ(RoOG`(%da(-=;c?mTmgojeYQlEn-I$6IcZ;(vPU`8nHd&*mgx#YO!Y=v`AOd#h!8K>LiTN@{D7baq7_Da}`!K2&up}$Wd6;C@%(<!8{1gh6&5t3z$I`s!mO!;yFtCFA~@W5oI``CEhb2hmF>aixTG?Tw}J(uJ2e@IOo_VfMj-`H^Cis6AU@+?&P#LaAiw%SeetEgxq}4Swo!Y4G>Dqmj&h`XQRBNb?$!LDA>pI6<FP<ZwkIg@&ci67OSM3#Ocq<!$lR+#;c1_br@BLQFRzqhf#GHRfkb^XuayNGrrLvY2$>EHr6CoHV#oi+Nk9i){-_(T8717ra=z=G1ezm21y&~#L82iVKETK7_76_-2<BQXS@CM799qIVT_<f%af`(Ede(eTo!|Vj1GgrdXH^i-P%50(cawZcKo0nqR(Qm&f0Q6gHi>2@dMisec}fyi<URUDlA}BDzFeMv&9dzlUBdc+c&;m1|ZlB^284e7R!ns*h*=lA7SuVbQlbVF@hE?-iPM21l(Y7Sq%CyIt&KuJ!bYKRrXAtWbPKTyOKJK!8&V>_yH<VM6qa=6RIIlgoq!sBT!@=W3WxyK<#W7KQN2`nB|0yCs3rHwZT`n=epH9N1N+58VrUpf)>pg0(-7oYcRMh2K^Wv27~n;GkcQuH`;LiEGY-KMT;xB1-4!G17SZ9_5)#mnh5)y@kBezQGsrep1k%Y*lh(aERh>PUOKG61zMzW=lKd;P*H&}LQ8=QD^p;hHkRtbSQkNy^i-EG9cN32V`(4Y#&o2^L5tuRk&L04b<vd$hZ2T?Fh)Z<9O146Q(b{>kscWPf_Gc;3=0xbrdKR^23n-C(fX2SP{M*RLQC=tOMw{M-8@<cDQ>$3X?2PUgB+~`jXnF&0ro;f>mX+B5#qx<-5~UdSA{{#I-#WIfSkn&ePWCVKW6O_+zZvv?TRH4`os}9#Tt9|-3gyM+b3R;Q;dD$h+WhW6D?iNV&RI?MCcPo=%QYcQ=Du}E+b7eI+GjiOy2d$KwwLhx|l|M1}Aq4(m_Lk-Q^&WN2w0q4oN1gIVlI34C;|Vd_g{glQWg#V!QDfY=lRo4QZki=a4{)hhQ-z&@Ix7RDmuTMH98+%VZTR3<SDGEOm+{$w7<sWGs?`-e85vDyZj!HCFJqrvD2Z{ITria%vkaGYXCR#gX|%l4PteCS-fL1j+_$%~FoF#_qvLR-MxJuV%76L0Ed0ax6W&vtwTR!Af+QV<oaN61m`=gOO-35)DS8!ALY1i3TIlU?ggbkqGIeMA;4~772_*SSKY0@5C61uue)8-U%`iL7kKs+X2W(1cAXowgZ%r2m*s)YzH7C5d`l9*$z-fA_y>wL5)F1BB+xR1qNe`L|7*!#&&=;vS6K*Alm`QqKT2n4euPCl%tb!bW)x|C#CNhXLV_SWjyDj%j={>%6Lx4oY;2-c2YuxL7y{3>$@V3(Ai1paLfz9K8O0)P#+uWW38!=*)f4+*r)8cuurrgJR`_nPsQ2_`*b0Dw+#xdfb7CXb6}*x#-++m93ThTaX}Ly^<3X_#X?7y#)))k_^>^JZIE@B&+NOLs}t$7mpten@(109ltbL;NVbP{cOx!VgzRDXsD)||#ULUM)=IXAHH9LN?+V(phneZlZo9oTg_7Nm4vrcf9KLiL3pMXdTJm}iR1ea1R=D__@T&t6{UcRa8xieMUj2bM-dE8<wdPO<>Uf}z2kN+Q)bT(x4MfvGH1&&Uy04X6ty+d^tdQ4yove9kiVpc(8e+iNfG&Ct`oKJ#2<D+NCaI2Dd>Y_E+5k6=YU$ZkOPeo7c3jSFq*{7AN-3k1QlRj`FV1Y#_Al_#KIV)?vOR?I',
        '8S^9}WDhGvEmV6+f~1Q9?l>24$0(<kJF~v%gze{^emV83IDymE#VDnYQtBwBj#BC<rH)eSD5bVuO5GXXXs13dkHE<`M1$fkkHAS=q0sJab%jDt;G{Nbwn!sz@+t0`-e{1($*Z_)d!zXTPTR{N;2W*A=r9-zV+1YQp49E@2vn9#Tt}!5WPXeegBxSCz0s2TV8|OS$xVg4(F}{hI%~`Q>{PtAUzv82MW==~aATUeTrv)Zoy#TXU;tb$)-eX#qzNR8dL%)x8RU^H8Y~8*N$VH|h@+G5pCECxWJd_RJGn7}7H!_gjJ?r<#L-5JK|e-^LBks@Kpd^LXYwSsZ!tmQXtl*)oi&FD*#H3xbc>ezjY@>9QL#5W7*rx;m4<!rkI_SftjW@veIP0kvUP;PY|&4IJSkp4^D_ax(E{bV1&X6*?2QKfcmd+*8GEAz%ykPAN6*+BEm*EwfH-=_-e>_c(}Kj&Gx$aeLWFG99I~pN1+#8FET3i`dd8l8=w^L4f#lQk`fdR;$kK}!p8&o2UuEnCjzMur<k>*R&?nXxcnB-FZE%G|2Bn2PxVlstS)Vv+KCx{W;+Sy244T(HV=r(F`lL>^F}3RxW88Vw?m^88vYyU}yv3_a8`XB`R^LT!+iA07!r(n<Ir|Yk#X{YKD|Ar{`x6V*n-TUsD9|a!IuN4IALO`O6ZDBA=25S}TO8396v(6YxNm`DP+YR!tcyD9Zf*4EdP3*XpF8?<M}KZ>{kc0Dod(F8<0Efqgbbo-QWg6t7Pfpt-XKe_Lf8TrZ1DkMiw|_nJcO`cfiFotQk+7JDNG3wVM~Yz;B><hBK25zXXq3|7q#?HE7S-B8^_)9*)x080=&hb5RgC~wXnsggsg2KWFa6S3jxJ>P~Hs6oAxMgcE-I3uA(Bq0xi<#9O$YJ0!0Ycad?_a5(Wb}Mq!6B0y*KfRdiYu%&&)rXts!)XoVvlFpH4E0FF`E&bUCXB<ppV(DUa~A`L`~0^Ex;xEBR_qXoJbXKXJDL`(Nl#)QL5hr#gp?BZM=IW{jHhO<QMp)kvt$~(x;f=4pz+)B2G73mRs7<32AXRZzKO7bNN=V}s<h&>dVg=D%`vOPA>CxY7HGvqB)dr98V{1ff#HDF!2$P?phHK_3fE@a}>S}qCp#P~qyJk`~O`UK`N=qwk9{lp0PJ<7n7%c~rhLMQW7lk2%83t2UJA-n8cdp_JM%f-pr+aPw{d#Jt3G)}Q@+*qh3r2Gl@f!$f-i-lmkDJSwhutVD+p57gq8l8>gn&{Q&nT#~U!a7`O`a7}sJmrYR^JL<;et0>7`kPOJMkgxe6Rfbw#Oo5QfFx)YiAw2Y;vi9JU4j*q1Z}&uW!YSgB{jm$-j=v&@uj~Z+Gd$$Yc*nOjN2_+aoeIKdz+GLwcA|zty?U!RJ-k^-@a+nYqi^4`;9s0mTI@X^qcosa;>&^>4v@BQYCij)-9G<s>C+kzG>5IwY^P=BxnY(pqi5kBta8D31LuMoT(&e8nPc?eRI@2_60^r(A3n!UJ^8rsMJ6dY}s%SQ7J?wF4-5*^P&djUA1mnToZ5fYJ2p<+ZKbUdTv+EtB*gaqU}G<ELHCP1J-J{x%L}#UbT;YkZiU_`v7WRRNYb~cInnFmRTzIt8dw+TI>*~@fxIT@Kyd;jn@!mgAd&#<Wi!!d(0Ah>Fp_R77b9)Q{1o*g9;e1MGV?46FkniwVXrhzgy|VEBX3g;5AdtR1&m&0We^cTVkqTK=+^qWoI6SeGjg5%z~Idv9E?7QqT1pI3l^I-aW?elAw3&lt8BEUPQ4vJ=rJGcA+b>k;x~*K9F2asY(x?tJkDXmy)kz0QoH7cTU`>V&cV9MBT~Q$v%+O?%u$6NLIJ`a<`{Cg)3kW19)Uc4gfs%mAW{-kDLgZTod+0Hm8O83s-<W5FwKpx*39U1=)i@aFY5Y01EE=(6=x$`BNUrr(CL?@<=}AQtgyS@+nVksY~)9U#lJROY#|CtDW&n@)2MEI^sV-Mv?;=eR;Qg_xbR9I{ZEzevgOW!{K-T?wjz#HwMd{0*Lgb2_RBV(cWv(q6n{I1q^13<$*|QGEO}ZNo~?>(KW^+0g()pp_c6d0Fe|Xu`%k(2^*XH%^?PtNry$l7%_w9f0f(Z0k`Nf7z|@{ShU___8p)yXxRa|Tg-frbOwWY(t-bAwfCmdpM4>|fPwx5TF(n%gv3`r@Si~X3mDjr1l^Aa^aqgl{<$Ai{+l0!_yW=&K-ZH>fAxiOiPd~l`Om(~qGitC5_g|HKqTylV!W83MevBahe!#vv5kN-cVdnPXpzoWioKA~wL}=9#m~T=3v7)17cCs?1kx`FMyMU=n~}eB_}d%1#|YQh@oq8xMH^Dj_!GuOIFsU1N9V5Ipewr#jnMFjn+NZ-hdRZc5eDKA1qcMx?KTX9f&^2p5eD`U%^m51ZV@*~AEEVja0R3`gC2W<KElFza%~M<kZi1hb-+Ku;2z?)r~8wiqtO4-I+NZLJ1sx$>Xo`Da>{$6C*}Pk48$RNUs`9vdt#^L2_p>bA#%#x6FcKQv2*esgN1U*jBNKr&&N-@21?u$JJmh0v+*v2rXhZNx*rkLt2U^`T)4nCD4uuBI6~L6V~uGMd%;y)wvI5+hbXXzW{znwgMxc}<`J4b^#|m*YgYxk_5<4>;kc_uNI35LDI@k6$DIgHA+90zb_TXV_u({*(EPgJ=Z`x_R^42%VNhVO1euLhxiSjR7<I0>oP=v+gX+2(gDjqy8QFurJwuc0%^OS=&pu`eo`*DS8EIJO5JSx|*8DK3-Fo+Ww_cQZmB1Up#~;pqL@K;W)fLM(<g71B(1(GpiBEe_$o2x;Aj>t$+e>bOH%R27E&EhfFLExj#Hk4V23V1LQV6AP9;3X9;Cm#uB3o1d<>WGcb|%lc&SZNR@<Z!#2pW`wi3xbp5F?v9dQ<`nN9{s3NDXM|`axWk)(1MC96}z?Vw6M3Z8(T(dA3t6)4PyWLsZh5ZkSZIj)2%In?qix3pv9;Hgo%t6gww5M*9Y-(}fJetE_R2Zem`8rh`<k3)!ZyU^ZA^1xfU{w<=x8sKml-R~Q0ymGmM=g}aa`f<Q(YB)<xhPXr0-4blwJ^6seC(Dd$TKJ&4=ko7`oj4tH-)=7=fzCCJnA%l=*bG<~guZzo~;UsO9S~cxU3mHhY=)zdcpkM(VG(uO1s2aS}UJy?_XB{DA5L`pUYUoCY85BG}V;&)BahuRE5Ku@lp&U2+4Jx59dfj`0oLwcTB{z1jMR7QO<`E(W!8NojhH8YMLCbD)j?iLn13ckTGSUcqo8*8n#c)K)NJAxTA3DG#j3^o9?HQGfBCvv*7b4hXX%+w_B^gn8^B&X~H2*7HT*~~h5EMv$#^qNqu};2Ci0whrC)S4?3~3}_1u=;@Q2omHl$3u-kc(RHHCn#a7J95LB_mDQ!T!(eS}%%{8fXz5VqmE)IcH7(Fh@st=5xKE!J0x{ph0km6%?i*Z8GOqv0^TL{UGE>N2*0S!j5z>9|?<e#24w1T#!h1LFjI9B)h>_FblfD5q5(&aJAMxPeul!Oh%);I8v#E#zanNOeDYC7fMc~Ytk2sPYq-<CUP<(1&hAPwtO)`%4}eF9pE_|7|-q^#9``nS{d+Giw)rZMP5*EDdHQJe}n4XDeO1M=htkl-$o<FY{NwOcR++!_-A8A%r$2Cq~e~2RNRn?8&YxSE){2kYaXXb!)VbkM$8}%u4FW6fN)KIlZM)&6?AUj0d};e9)E3_nAc#VO&UWgZb-!qskq}w#Zex>qD>kUH-hO!1p0e*tS6QJ>KjIgJy`TNmH+I!ELxj122uSWsvku4XG~OYOSoHt?1}x}PHd0BW3bRm#&~Ndx;HSj5MWR2hjwE70uw#HPH{3l*6Od+SI~<#Q1SaKK?)v5_aeQ4DtJ4n`2Cd_b2Zq#$Znv5`YR1BUPtS%)Ss1_hmU<d5L%J>hmU<dNWK0_l#{*U;iF#<`u1yUf8|3H4I5|4P7U=1RBqUSI%!m6tl4T%mo~g@j&EkmRMmzJ=W*elmLcFV1U&33=mU-h%%J&S+4o?5sIj><-M(j@5Bf*$L2pDQr7@iurs_&1Sr}V6*R_Hgg92>(344Q&s?_5%m^$cpwGR5wy=YYTf_KnC`IveyI(X38sukctM+9@85Be?8pVdO|gH-G;sn`7y$>XcpLAa1RM3?O*;S0|Woyk3}GufeAsNMtG^Ep#{a1ndYnb{la4#51gNVbP`{bGS)gzRBGsSDNK^A13IG*4G<s`akg5`tq$k7f&f<UeQdqepW&haSx#q&0-JhLF|~g|wjN3a07@UsQvx`%$Pz)11ztOKR{{KMcYoq#jMDq6S^`!)4Ge`Z0vGhLF|}(mHV=EvQG6C!a}IrsRuf7VFWp%4hPmDe3Z=V)eTHJ(^B&O1^%k%c4a;Q`w^_&ibPvFa|FXh(#I<wSDw1wou#VQF;l;u}~X4P;}aLEYyasCBg_Tp*Hqhu#wBCG$a~oLv-ArV_q-R4g|pzz5j=W+92rN<@DM@2QbtIsq3L=2{edO&4c}#v%Mn=WHJCX^t4HRthG24)zH%?qm>5fEv5RDKGNY1X3H(Bt4l>&mY3L)Wk~1?xvCrd;<nlu6Yd7M|5LYw6m8p?#UWCtq113=Lv$73eTbc9ueL3}p;Uh#nyD5ZN<Q7%v!q@W%pk$HP_J@ezhi!Th6^ru?+?(opj?^QGv<C5^;d8p+#ir{L3z4^|L}o~JnB`%PP6#EX2>Z8E$XjWw_lHW7nG+Pq&0HCjD`!$@7E*VMY3PmBW2x)jBsN@%xqcW3cAT`)nBJ9F?-9<g`Jz%IqpzKxG^DSsxEOI;bcgdGLn@u!J8l3gKSB{NLJ#ExDsc9uK%IO+NB<AJ#|v>>SgP@(gkRn*JHgy^jO($G}|R&>wmKKN!!uutfP<1*YD-(<FfU8>*?c)LXS1y@82rLniIALf{*wFf*?HO|F?H_&2bz#^1tYrh-h?8c{o2tp^zn7rdJ_tNP6{uKlS=-ub1qq%tTHP_s$pLebEhwt^$zwDj<#=#3*)rNUcKTg7Cu2@Bj8)r7wM2gL^u~&fEV}#UptM`;=e&{^Qf%8@gSy1pe>r7w(6@UC+z<_k!WS7YzTsVEF&jf?@8B%xxa8lN1RaHG=x1Gsoj~@fn-{@j6UV+uugN<K_>ae;1^`V*77?y#Dus;lCFQ|Gi-Nb_<3r0C?5M>lQb{=|zMcuWkjYcGORc5jzr7{jgu6$bPeX)Jl{mtl6i*ys#l~+Ax?a41t45(iBTPx@WhL!(eVd1P&mn@mLzW5YF)}<$*9a7y^fmY>U#&wrdT6Gfu22?Gc`vZvQlx7d!+`onBD|te1zWFKh^06;$20Z{W~~`r(H_I|5#H2%N=_2s;4!@UIN0w;2Lg#H-s5fit3hT8!9{nEHZ;z}0qr(#As=rf#8V7^Y6L2*lCo5IT0KYLo_*qaerfc&%aT^MP!JLFJ`k>c9vr<1lsntL;-zqcHj$z*PS*b?&Px6W$=<E=)l!5e$OIn+{VKr>gpgsn1y~5e)jH!w*<54^!>1pGE`>I|SB;snZ-@YnZwR0<Swv-Qh^s{fMvwpb!6!7^YtJBEpVWw}MnV>R)Y``qePi5BnvG>^E1plf1B{Wc~#z#NH-lnTTPR2tB9v2206|^X}O7Twz}(LT(TX1CK9kDH*x^pth6@-^mjCSF<+4v6Kv@SRxn%kM<=mw_QqRU!Y<k81zSnAFy5?rrKdYjR+Wa2xRr;h^1u4IhJv4P_^TJox>Nll&pq-WkB^~exwlr?EqMUR6FXY#fTk=seafmQDnY(XDQj`hr^eW6(5I0l||yYSM5pYjh2$pPX_P@WFke|{2BXFvf|Sx{~jpq8IZ=K|D{sgN(1di!tsPU?84}m9KU#LcFbVLZfC{(hBjOpUF@n3J4Lar^R|ma*#^g}Ee-`p+E?z`7Ka*FE2f*n;^NR^bGlt!!mlqc&qz3)XN6rD{gUGsZ_SPw%-HR$nBUOa=6e=9w!%))w%ps(%F4BDX^jaXPKMPQ-$aVG>9w_vts|CuSNq<=zEIh&F5%agmuDm#t>0l6M!)3v#oM*qyV$XHq!sfUTKf!cu^Jq9iefpc%%Y{O%(iEqE~DT@b&;aTqf#*|wPVobV5#gF1|Er`x?NquuP-mpNI2Tc!Y+({$?=P~X2%R>>~>bnZ)jswvi@aWveuVm&wH>e(?6>WhMhyAU|rdft<X~k`xz&GUv-l(^r!U6pOr1%%wZJfr`B2Y%=<1;)Ct--WUWJKsHS7==c1IA4IEy1*7qXSsrl>I$C8!^mRNMp>at7|t1=jN4v8XlWjnn>PaPCwoFH0aksH)CU;h(sg5<0-1Bc;e(KA1^L{TSb=a97yr6FC8F=C5SN;YtK<=LT&RGS7Zwq?v&iC~FE_v{bLbdoBAVds!2Ow&$9k)Aquz&Met#3DENgG4<CL2`DAfx~dK=$S`WqNo$JbI4kU($EjbSdv949~(Hl@+`7Nsz>u~nLak9M6kr7`z^8MwR4rhuyaThqiLt2NKbvpXq<3TVv&1DOlFB41j)D91`fl`qG!HIiK0%>&LL|ZN{@&g&zu*fWNhH@%Ck)tsRqrvW%@Jh5<!Z^tpxq}%mVa#ZpVBUij!o(su9deJ(lFKR6b>aQK1*5ZIPuhB^DM5;ix8(HG+(6R5?7*4cB6=a1+c*4Tt2gRC4+A{|ozMg<h04*03CgoyC-^lmGUocwNQ`re5u&ITR^kcd2G#lE-D7U{><JB!{ImEo;6MdQqC1U0P>lVKL>32B0Ypf~i+XXbwe+gk2U~;^c98F)$c*4oj&W$Pg9jMQL&$%VF4AOgRfmTN{mF>RDBqLyF=SQva-k)QbS!mY!*ipdqAQB<Oaus=DbpLh40R_bYdGYh97yCRh%D7tsyBZYQYLRupc6hLC!Zpxe>j=BDQfsaJPMy&e}*ukMg~b%)feJ*2KTXBm=D4uY1DdOa+pUZnJxr8O(u9EP36dIY?>L+aHXQm^)q+MDXO^mEk9R+a+mwL|LF98$0Lka~58)ay|p^`fQ69PJ!o=a8;af0&2V*u<XS|Kw@(FU&N%-~P7&1ozpXS)b73S)UN3S>K`8ru2d%sGQ*FgInzcJu=ud8O)Ll+OK5)aO*>vXVN5|U<U{?GPrNqnx@731iN&_ZOl5Y#cjg6@dUf3#Rx%<40bn?w9*xeSF*1p^^tHTJ2=rT+4A2MNf4o?)n|uZR}4R57F-|PYbWTD!LG?*mSoU=C7R{;C|znN*d=eK2}qJ|B+;(8Ngj1q3_s$cyFR$rPS9=m!!6V%_ehJH`b1{7(D|pQp>td7Yq7*-En%l<Uo&#V;*#nTojCHV?HA)-{~|6W2s=glnvt+m6bZ(WJRW&5?xE=7@a?cuw2zHGVi0w8D0&?E)%J^VkMkCX%!Zw!eQb2tDT)N+NFI;8n2}5Q@YJPTPkJ>5X|f78!KJuD5vLb!f{Wx>DfGgXsz-Wn!cA}y+yRYVxCt)hO;0GjaAkmymW*%{Tm+QW=!Ki$BDnzyz3`>#SH)1*2rdG6Y4pNPaFO6zp*Nf||9%eF%hl`5*Xo_1IXBX)wN5N|-qs11GojVbu2ofE%vJT}u%h~6F4eeTCdK=>PQ7(8*ISpv%B_pJ-0Het)t+9b61SKuam!(axW!zE8}&E7dkWFO5T4V)onY>sl<U)TYj`u$>t+c}GrixwcazW`pXtEm#zfVlWAtsa_Qam}O9H)|J{x{3oba`BVAsDfckn##BLDUK8L#UP+bMtk^8Vq=n-3449=^Q&_~zT^H(&qx@b2;L$M47d>X%pBP_*q8_ttj4)z7zt=N3h6W_nYE?W$gHH}m!J+vDd?UpI=+Z1T#>%H>4a&0EuyMA#|PtGCK563e`3=3mnb+^gHXqRwM*A8lOcHBVZQ>%4$I7;^ho7<df!tFyY16M6*~C6Wb-1dkd)Wq;=I^<?N2O9YwaAq9%{A>~l=Uz&z0F1e-^sW1&y?Q>B9>*Zmp9rn|RfMJI~7WWMt8c{#|5NJoha61xqKO*b^=)=D<p!zZY3B|DE)vX}aj{4%pu(%zj`eDCBk^N@Inj-qgL6U_cy(+WJBC&jPCid*|x~U%V+c%@mV{jjBT)#Qj6>|M1pnc@F%rNj6>{p|9VDP6@(1%GDC=xts1Ua}GSdY{IaNoeC=tTfL7AWR4tv^Sr97=!!s5z1V1@jFAgZ}971J=vKR6Fda5dp&vfh_JDI5eVu_#x1afZ=u|?0!Vp0nmqkWkB^~{u7E}$E#aGsvY$S4mH5CW2ztaOBC5}W~?c)YL<f}3q^Xhb(uwCSvKtrOy~(d)l9ad&SP*NZCsZvq80C75F7?gkzJV0@|a2@c?hW~Ws',
        '-BsBaag^9yNlT_pyp`bPAbil3)-#7AWFzbSjBylpy2y?OU-F64Nljpg%hNfc5e))eie<M8L2^AdCA34vnZEeh9Q9V7MI#yB`sD0QBKs8BqP0kHyg?yt);n+EJh2Py;MGrut#OM3MbwIxESyW0B+2X^M;m<}n|lNT)gF<?w<}C3q1mkHLL(X<-B8U+{Pejkn>8HH8i@-b)h@HO0Q}lrqVyhN-Ewp{0<IUnTrYv9CLY%=9XH5vetyrIO!YCA3Pd9V?aG`Kn>Qfnd-d9e%)id6;U4{WKzA*ddU`eFKL^)DJ%d+7U3^j)dKh2s;4!@UIN0e$2=Ex=VO<D@e7YKEa^|SawYH!+wb(^Ud1>_{XPjs4qXb_GZrk&GMA0%3#<zBpy&#wo)tf_@`cV4#UVI_e!a>m0BZ6KJ+qh7!?+|=?=`}=@Lbqpq;~%CS-^&lpJrU_=ERwc;$!F3ki$D^2_2IuM$DS-fUmvv!yk=%3#<zBpzf{Hq0yZ_}1|{hhb!q>v}G2nAZrBjpGImqrxILcZ+$9Rida9v~!riY%p8Lv!f@QL{S@uSAK}KkY(R#L-xYfx?LQGoke#em8I3I%3#<zBpSSFr=mztoho2#?=G>(by1P1=O9Qnq8T`h8jGn!w?GEHL{TSb=MceVTL+b+$7`X}_#O_g+}O?Haw@@SsRyk@u*9OfFl#9tFYe-YnM0ysY%v#(7kg>D&S4l?Od;K8EE_NG!nT3KsIizzbh~P(SE?H??#gzB;FT(}j+*geuWA=4`WnR+cfojZ*R)FnDHgX9^v_BO`mf{tvrwFTZXAEO31+1p``i~tUZ;Etoi0spcDhVf6*^v?VB#zi!T~Leh$_MV+>pV@21GE*6!9e)Vb_%-mCtP+?5YxWU2&V6aN8U4i{nPhXEzTCbEthWCK_7j!a#T(1hZ0)oeHn9d<rp7)9Vx7dDGBg4(rwn|7I*z)HKDzcNylL+_~^P!6;K~8ZKLO<w)gon+Ln9gk4wM<|f?sM*QNqk@DHiLz8e>Tb-)-oA_$Ya1+c*J+=~6s4Y2);-TwAqcj?XH4UMBvUrAuXjCCMG7V(s0)WdB#St!di)5FRRTkxQn+Ln9gk4wM<|f?sM*QNqk@DHiLxVZgzPRdd;wzEEO)x9<*h*BPJz;55q6%Hh9i^wnkX^G{FZ`PcGsG*N<`7GdC7OuK+Pd^^k?eA^%A$O3^I%t%u<MH3+=ScSh+iBxQa-zRNC}r`U%Zu|KRhi#|8PPVNz?x_+glNKiu~W2Eb9mx#rdERb;eAHBqfR>k2Zq-i@AmD*u|S9SiI;1?zyNgM`+7@6^&1}ixlnO?6*1lb5qnA)*+IVD2hDV2>LJP7P4a(ZzrTymYy6aIDXj<*Eb_GvAJ6*Qsjk{KSi0Mq<+9W=MYIs6h$751pOCtyzJP;n<iMg=tJ(gq%Oy33-2TQz4%2%m$ef0hLK__0H_Lg$#(N{s)%>mJ)QNi1V>SJ?BY!mEM4>=_gqpJueNL{$e*63AW!#b3UXE=em;vD48aBv+-HMkeL|0C9cmyXSceKM4C~G-`SXiCGJ|x9+ODLTWoB|r*vw<%oF+2kzG3RzyL83-s-QI%4^%<xZd7grZ(jmMhS*7tT}hX;xB>L|?~)0+WM=YN%-rOe%yd9zn9oX0i_tpsHSMb1@OSBoQ_|wDU@@`2g4H^Y%yggw#dam#EHk}kSQSg{(&DCh)M#@D>&?)@<TU)P<Zf}^R<F5lrYmWZ7HfvOyRW31Wu_}wtVK-|1N5)Z`N11HBXiWlO>hyK8#H?1CP=d`)1%|GRS^66P|MH=(!ZL}3pc?<phr91GlE*}+;)OBf{Rd!*XV_tpmrk)JyB$^1w#1Hzt9MB=RT(wZi0(I4_@7i5=_$xA~RaUO>hwk$2KIX5nQCownC5AV#N;5=c$H9kkkED)8Y!jMW9D(v10kt>dh>w9b5&6Rq+QG6-g=t7pV?wa^Ds}WN2wuv^XV?ItA!lH7%}u@uJn86@sIx^ZTcv^Lpae|FKo)<*=gjVlFx_ht-@HbIsWm*lYQ$Yn3Axb2)N3tQxtPbKJhDUb~O$6si_;p=vp-PPLfpRId9~?df$&7>l`tu^d*xSj-iSWsu!|r$Dp<SsKK+bH~R)RxKCpa&i1F#A<~YcUDn@EX9SWWhq@Qj$>R~JwMdV^b{Au4_TG`mCMEPyAZ1tV%(Vs&7LYQL@l%7a&eIH5kEeK=zBeg))3WQ`g0e;jw@;`_QcF;F1lzMsV1T}34JR_Xi8y;*R5{mFuT@adi^l5);s6oPYV;faDKnr&|kmH=(?UaDgFHA{lk|xA09qEe0lrv&9~2QzW(#!-Q(Mj-^KLR<Ew35%r>uQdoZ(g{%k4CEk@dm@FwQ3H>v*m`0erYr>`5e$ESjymz8S{!uG2UOButtTc_!&re|})Cw?^?-o^HpzxbRy{IUDzG+lj8pYQAOYB=G1_f{7TyB|-~oYH?kkUQmoiD}x<`R<S2I!%_t?H%xIdOO_geD`>qNgPzrPOxqEB&~-d`{{O8DW?7`Ch-x637Po+)M)xwvyrs?hu!U-6Vt%;?84uj@qE01-0jNSjvbqP8B+tq4u4+)Uei0(CyKO35ZHuMdl-rC7{_lKwnr<dFYt;^_&EV))*Es$!U4~Ohd$tmc=P8wQ9%BDnVlTq0Z#ynKVN1JDs;dTxZ=-u=9vVCJo8N?7)<2ccj8{(?D)uu_7I)se!TFziBHRa+RUe=Pk#~_@Zupv2D2|CS}L)_$Gl|b(=bn-IQ`<wc49-d5J!(hAhkGp*QR$Zj_!y?GjSff<7yh_bwr=JIIoQlSc#)^c))ZVjX!jN<LDe7Fds*ET%XZ6x{VK*j-!VTOG<I{Y)QNlNAC$i*W&07c-;(K4L#aa;^+?e+*}+D?dl`r=sle1QXJhT5}?G<Ig$WAj_$C}jmFXHymWFI9uW5KV`@i?1R3Y$Q~{|tZ(ldN7)R5)EvIo_8y_$oM|Z^Pka1pukv<wn6VC^efF(!;!8p2&515XlJ8s5M@4o|HGUNRtemogMCF&!C#%~}2U5ukU;L}p>%)r2PDUL=W12T?IDKf~(mq8vV97m@T^#Ya8Kyg38xVJKuhKOqR?zwCAXgCAIgf-ApE-I>-+%ui0*rk_3J$k4uU8_d}+QKI)h$t1xQ!p_qmZ;d}Fok;b5c9iMkKRYruGOP~rC&lUH*uZ#aC4u>r0db6IJ#@~X!6xI6m2qrIT*_k$pb9h+_7U_t4Ajom@}D_eJ>*EdNfhsrV>)vO;ceFbkoTu+ep}7wKzKIdNh<V5fT+M;h_^16Y=Kr6cYvH^At>s3MMKhfW_x2W?(=kDu_rw*>oW4{WmF&=H|!%adcPe(GBpOgm{>VDX~zF81eu{n?z)uaAZKHS=Lg8jZdb`Vhui-GKq(Tj1&}^6B#Kwj-E(>5Jxwyue&Mr=%&TdUF+*c=jh;l-6_<gQ;egV*4N#XdUTWGXnfcY(*U9Lb*E5|PBD&dQeSu1>d~tdN6)+}C5}$5uRG~_bV_k_ll!`pu1BXBN6!K;C5}#^Pql0H=v9lOQww{Ou1BXBM?*8lgmNk0*WI*w^s2?tsr7XyU5{>J91ZCUUoVy4LGXDBC5}#|uN$yBWM6kz>d~pi(M{^>ZdyHhmE!0w^>uTT9D%;>uGFJbi=&&^*WI*w^s2<sP3!BP`?O?VcUS7sDaO%V>+9}HJvy~Gx@mpg^SF1kU*YGcZ>TRnaPuv)@~nt_LJhzL_k?<mXXFd&(z7D=1@#oq=oi#8uwdU$&#DH-=hMvTVtwDt@j{sZ;&ic00CBv~ESSL=qPmPev_v%-erE~b1ozIe^o*ajhnTG4uPg(T4dR()pyF$IW~u2iA#J?o0Z+16+dw6lSzS&d?&fumE{AeYvgY+#@<mpuO__a>1p-Ic%ndr+C&7ND5pWouV+E>vtBA|!kUk8T@wwF?F1PXV^0>^t4CEeepqdQ7+X&6a28~UZaG9`&ii<sB5Cjh5*+vl3R<(zkE~i`srTbyP!{{0=cS)Dg7H=4ryBIH{+DuqVXrv)yDR~&hJG_GsHy)R}6fdU~e;eAn8LtPkHjLW?)ujoq2dhdGZV&HpvM&ud#{t#kF2&0Y)@8<0(lAQcLrj%Z!{w&M%gJ>|7~4QKT}BmU)o?laco}a3jk=GGwv54f>@LR3tAxub<*}2Gms7)Kykm?B==eN#m*VAB!ev9d7~}OwxSUe1D(QGRHC%RV1z;>Ct$AI_V|OWDP6?Mq2NqRjB9EQ4<~0qMyO_uBLcE+BF5|jPWM`B-cJeUVG+dr-aM?=GALb?KWY^H^Q@p{)5yTX4@%=_Ji#@M4#T(rr3QzG4vI!Bm7qX{gelNog3}xaW-6dI?;w??57~YE@arFC*;uP<}n?iMp_dX}$6z{>$j+^3Ld-{sCAXcY%2RDi`9q$*3FVm-Z^AeyBkP@eO3yL(ci`rx%4nv9)La7GIu3_}mCWNx0IOQ-%Sp@Ede7`uqm*HzeWttLg7E2RCDU62hMSGn)x)<g)p`itwX(<jU6eomI7W!HfLWhspr3s;=<Wy(%NfSc3_rr5p6eolZ1V@|@x-U55gwUXX2YEHEO)cOsq&R(yH<^>S5$gL5W%?M4(J92y)#+p6tUlfek|#Bg7+st`CcUa*F}gT?j4>@KVsv%-7-^wb4Wo~qKE`8o<td*weJr#T3=d_KN1MEtrf>HVr;pLjZehrk4Wo;bO~hS1e3++CHW4hQJlTZJb*B(VS0|f@duezlNPY!M=DLfMP58u(kn1i^Hlf{rA=h1<Y{ExT`aHqoCY$iN?xf@B>SPl#scM;dPh<%6$tJv)#$$AGvWejA{;OpD)6<gm>7H9W{nx*FKl&HObj}F&KM3ysie`Pc<K>=lGv)**=R`m_bWBW?{oDVZci9-d5WJt#3+77Be-LEOz3rYw8_D2)Ps<q~NdJnlM$mUnj1hD09Z)(-74`{UBrP6w<V#xI$E?d*$K7y_kBON~In(0(yE<!HJm|PfNsagQ=d7u5CP0zqgNgiSdc1#KY_bstPf8h5(&H|m@ks413K}DZfi^)jy{_qTrU1tyq5gv)g~vlp?gAiBG&u{Vfztn`Y4TY>WKEO1=2v$Ck-MbH6C=cHa+mz-E+F#w8H|9)P1EEiIo2ahp18EMBA->yI(bNLnkJt$-})?I@^9C?1<?;rgXp%_`@Y5equ5_8^?i#uSGB8J>9tiCjrzI=i!GB<uPt|qwaz@PzmYCFn)~Rg-GfR`8Cz++Mwd*fHAw0y(`<y+fui;Ckk@jBevw+2jM6)ycINr2!Nn4|)H~slMX^R$tcyzZYO(65R)duqq56_Tsd~sYLK>shQqSg6zJKac=5(4Z`sx+WM`F*-(P(P&E9e0z@hj;0CwfXLO=ZB|_~I3S6Mer;IzN2e6naXlaDGjYKm6f3aT0>OxX&-G^3o;!?SlNJbbhnv+nb4er<h1|z-`^6r8xY-tK!NLt0s4~il{mW=R~Vem(H)7A!d<Oyi4+zR*2C|TC+p+YBg!i4$&IKq~JpQZX$ntY9jAv7INNZKqhOaDS0l%KQEt>Wua+`pR19ah4VX?N*XTyv_6+EoZoNvyGiI<K|<%vQ=CMdUx38BZ|qzucav=}`2s{fYxSbN8-D%kuivF~UGFzJ{ru(q!<RQ79zH#MdHeCrx6g0B{`2A8<J*tl1@+a>ueSL=TSMhm5849w%|qQx@FwK1H<|wW`0erYr>`5OXO`9EW##%VXM+-f=M#D8V%T;5*n!-$XC6bqt_cBJa8KNr(wM%?y24`ERT)KC8x_P09QFu>5jd=PPe~lT7*-k;1no{z#}Mp~#tr6WH8hC*xKTmFX@C5O9km!%97C{3i(%*cn1X#*cA5>C7N-W#iZpvZAt(|iI@8hQrNywRErum`gAv2Om_41<8BY9SSZqMqd}wJ>wi{jy3(sULZ2)2hl<^V==`<s8kWMEB2JB+kR2~^4OVPA(Vguia9BVN5?t;RdiHVr{P9wjY|FmNk!-``&7I85w5V4^kjU7I}7<Tq)$1H{&;UBjMh!~Iy&J89;1!1$2iIKZt<n?6ukTT8>jZxu|yTEo1G|ppJg~4M83>!D`aeibRov8MKarA6ux)4Wqc!U#X`=Budtoxty=|fZ$g^UWKqNr<)^Ba1!R2+@$v}EJxobR)dM(&~u!)AcmY9n_CQc8{rTBbWN#8E+mA^mFaGZ-=uv%0XD;fY3*ab9UaC$^j<Epry(C>5RsspXw0*e{O1B$sGahlrAKG&-UajH7X#P6`Zw*OI9B-+^iq=l$o@>!v~|Xwd}~LP1_P7e^!Oe5g1Y^1Hb>I;UQ@H2zYpuni(MaX{I=n37TNGT&|!+z`W2T$aLsPHu4@T|H*P8teshvKmc=d8x(G=-^^Jj!vN--HT)K$vC<z^=LLq#YgUviJMi6qf@9y)5|@HvG|PppK=_njm4kG`h2)KbDfSGxy#1U^CAfyN2gSemd5!lu-A=c2(mWLkI$rFVK0%8n#~6k5>lAxP+<)u-b_>ji4GIhVB(#mMr*T>`8Yc1dUPv|s3zm+rqrY97?m0HO=fOZEsm~CH)j$zDaO&QHuReF{#PZA?n*rx@$Q6pxUf=)i=%;lLpqL5p&p&wJSZ&K#G+m%nKFwE_+$zc^->wBNrXjbr06(W9E;zydi1Kr(W&)yr%;b>rLp*IUv~=iXm!0@*ZR7<Qjbn4j>eY@Fbxn&Uw3VFjnDx>_jOBS@h5%VDb=G_C5}eAK>2Fq6~^MTeckgS3EkJ7Qa!q~x<-o&_)H2PM;BMX(g`Ukj?U7|86m8p<LIXKbvLCR-LyEmYkl4D^Ax79JL!6Kmexp+ecer|M>jE!hUOpf9S}-icVz<!-2fr<bvLCRJ?`sHp&q>|aWwRrk}pLHecer|M=x+RL`9oI9GyyEcUS7sDaO%V>+9}HJvy~Gx~YBLDb%Ca{~aCh;prRd3os(>1eL4-_*u{@Jq#CIHefDW02?RQvm)qt89ZCqywt1^!S<%l;1N4rsO_6MT`Y!Wju*-V5T}cEeG|va3pV?bf8KAwWoO8c0mvY6M_78szu+M@mOp|{m&wi04XlBRufZM?#|v36O_M2UgLF9+Ou<J1VxJ77=)k~swj8=3H5Vc;LpHmW+LYO;As}!JZl#O3oM0arCcPa`(I-*Nhjbao<rH*}duvbhWn(e-eC8kQiD<f<pj@Yh%gJ>}az|Kt76@3cRlsF&wg5WmaA;pa<gnO7>@>QF%c*(9D^sUs?P0imY{pXZy4O@(PC8y5eR5agWxDN&jVnOlSaAELIN2yJcO_oN+q*;av2ibOm2f$Qcp09EMz}q^V6!g`xD0`&$&{s}&6h<E#o~~n%BjgU<Mt>>Z?K1i+XGoIO_NbCut~Vwm3UdhWx9+?x*`P9Cmk=RX1<!1w?UPe$YUq3dridU6!O?ziI-EuWqc<Dm7Njt*eS%ztAfkgPJph&%Zrq$YGNKcg?L#rUluuJydD~t@hTh>&<T0$CdA84!(}}8MPz4`Ja#Ja@~Yr+1NJNYG%rC@`?a+#&uQNBoWj0#WGOSxIz3=N&0C&ltWPHPyVkWgWyf<H?09bOW$lZwti6?XG6Z+WP2TWq{Nl6uN{|!eS)ZJsYjgQoS>ml-aRaV+-<SmR#fTwfw79!1e-30&fu>&CgyJ&Pi5P0;gAwSYZsaoFlS5-MA~F=9wK40)D;_U2o2?THPK(ci*~(PGIIWg{-$Z6=Xgphm<FwV|J!$td4$0lm(l=+C@8zr^@*r1|5+rxY&pgj>#dgnU;F>8r-=$b$WTkRy4^0TrB~0!mMA_%hWw1M9uqj>8Fu4g{(nS)P$zAQxyD0JTrL^FUyW&gjiYXhP8YoM)bp*;?2*+}grdhp<D^rP)NB4X#8Eexg4<Wgi@)jJd%^V7&zj>}LFSK2AlQ3C|DG-Cm2-L-m03%FRAnH<Dc<D0O>Bi))<*ZxDOY)>vQp4m{@-tnjkc$#0r-sR`6iiTYa+8ezSxZ{aqAq^(Z@P1uXC+Q<8YV-z0nfEG5t9*Lj`m)h*F1G7jILMHFu6;)>x&R4cQ*p0gvq^xD9ll4276j$+)B+4PfN`|oX}D;{V%iR#^655u`2VC`(#QV7!^o5Kl4S8K_BCE3QBM*&44AnVttc7j~PMXpPUyXo+2{Yi2oUXijoE*<2LAj;WZjsD9Hbe*Xa4`I(j(Nbd(G^wt2p2P7L7&e``1flVK5>$io^xV2zYXqq}Aokv%&>injy?;>C!kn3c?+?G!}~TqV1hZjD|mR?Mf-^Hn4w7VrjcHhRf+iYcK>2gX*4ZV7ArRz;(m;xBf$rvtf8g30t=jC6`3x}@NBL*v=Rh>f#;t5>5JyQh<GyNC=$7VI-`X}jV',
        'tre+tL=5=?)Et$m8Dfx@5SEHNT(@D`TLc>2z7bBlyYNwdqQ<OAt97m_%FRo@!C+!hn`y&PU{%H#GbUz!oc%GIbxjb~6rcqrj%4dncJ}s`h-@0(MkS?5`*2S;!*UbDW_i;AMOYNbHOFzz%{GHAE)~00@S36^-x$bi|GuV+HD5WsvYWfu!Y-R?P8yPYv$xN+1-lZg|H5jD=BW|;BG2b!LLPkrc+)_)*;6*sVv-^gjMn)Hh%9V(Q@M;aUE+sj;wYW>p!dbA?8g{1U?O)7hp&usG7k2|$HzSkA;Ki(SlEf5S@LW1FbuhP}`Q6i?xvjIHc>#xlp}j}wWKH~47dB%IpZ<&uNYtidVf65UI>f)w_kj;HBz{HT#{v$&sedm%<`n%Z<^0e=VA=%ay1tE0J%#g2L4Ug>FWN*bbcsOSUNuSnx^RBu?@iiX;gddzn-}<7t@?`(fUN3o6Y*#j@S{ogZ<X#hCHzh6`E^0x7I-i$A(=-1MU$5fOcDE0aS>BeaWCY%A$2jE3^yb+yUXNtiM0?wJHI^e{JfCR4-9#{kkHNN=hk0MUfk!Ghe3WfL#hV(ixdJV5Z%q^=XQxT{ngG-3n@<45KlY5F5yR)^jDGRZr9(VT;fHvSL5^BT*S9eUBm=-oqAhu&gNE%v%7v{u8yR+_2KNUV@s)Nu`+Hi2x;9==`JE_^@23lnA|lh`xl>!&*5EFq1oH;*T4Sy-CVBgWiyzczr26=^5(<Cr-v_ZKfd|)`OVjVKD>K;`|-QMeD(3w-K&3m|L>RY|2=;H@8?(BM(J&&_x9F$Yj)gfH(PMNS>M;2t^NA=?eX)cuN&o$^}*k>wI?}me<qg*(-da}Bd1-lB#Tcd?(lynobQC5aK5j1R<FY(X^JHt3HgAg8w5uXn25oM2SzbSk{XZL@UR`>ewM{3kt|Roc+?2eznTem=oCu?gW$11k+oyI1V98C37A<iT9Gh<L4S1k0qf;qsvY*zh=5^-KvrK292!wS{19kIz;HVfc0VHQ0O-TNGNAe~KhlVRb^t6vsvY&yV#JQbR6p#OD6-!?oAfqYQ<nJH8E}aqO>ssr5}rl^viOAJ4*z$;`A+Bw=esA{FiDzXiAO>{pvk?#5d@|^Gva|!43ebABR17>N4TG5F-jy06bT+Rg7mLu!W}xr62TyNEKp?a7%u@3K}G^*R*Y68j9}0o9e%)id6;U4{WKzA*ddVB7Xybz)DJ%d+7U3^j)dKh2s;4!@UIN0e$0<FBA^`rOOR?u{j?adBQezv`z4CZH*YOO|CmqjQ3i;MndY!l6jNPg7A+-fEuo@L($*68B1PM@?^>o=O79ssWb*Y|+EwGxpBau@X`tOmI7-%G7e>G2_{CeZV+J#JJ1gclwBa7$FGjIprzpD7Wfm>zaLrj$r@Cv&0Fk0`4M`=6l|lythfFR&3+a5L=W9`C+)4xOM#7N}54$k>CC4w`njJHkvD;ZOzoCs$sQrsFDA(76W$6;DG8lFai9nlXzDQ4P05NukmssQ)u;vfn!cCCuVKZ<TZWcX#;U$VXK|6=6btpAoIoiV(B^)+zc;%k(MT)Qa>(|@YN(4(Rx@)?Ys`KKmi<dbhg6qXxb6)IK@j8cLWHE&znz7=%xNG7D4#Ul&r_XG0mz)=OQM^L%O2u49xp}b{#0wOCJ!Xr$*1WjO;U$6;i(5$jsD#viofMN;Y=jSy2|LAv;z4}OQfLI>+c_z(n$XEoaKwrdtkX-4!I^A8!)c;S@ju7s5W)e`iz>l{N9u-KQx4QTyKN}BkXO<?yKNv;FWkl<+?wK6l(1VHyRpGee&{8j!^y%<F~N8+74wA?e9HZuaJ)Rh+41gL2aI5yUQ-zy%Lb2n62(5yy@_x@^rA`-AN3q=O*v5W?6zS@@J5<vw+)2qh1)oUTT|SM5_W52GB#UC{qxh1dOf+S)LO|S+yt{yvL!hzRan_iR_H}(@Q8hqiSZ1zDo?lx=AGQ!Ky!FtWlKwwhMQnk$~+{8hL9TnWQAT>@^`B{A!s;5E!P}wf~gloX%6XZ_=h=kW9|F+28eJIT!?PHntclw!TFAlG}8->V47}(e26n|PgS<0(63x7o0{gE$H*(06RZqBV)tl-o8bQ7eEKg<=!K8qY*Rp*=|!6Cpn4;#a)O!+N=Q;8sO+nwuY|Y!f5P=o^X(1FlST1CVAv_zpZJEIqDU}~<nhRhaX&E>AHIa0qWvjb*eQwx<47KlyqJ;ce0aumuGffkL@%=zNVo%&10o&v*)vKnj&Orlw;Sws;#Rj4w-y$*rJ3IF&D{IBqF;ZhWwS5IbT9GrE>(3eEp;!EbuZoQW_s==kIK?M_uhE-lD$6*j_;nT`2Epqj@@e;{cA|QYjxdgj@@e;-D^nowYoc-xPRa(;@8AwKL&X-q$w>@P3iouu`7OO(dMd4fTJ$UtqT3u_Pe&sw%>8v?tXp|+wKN6vwez00qbT{`j6aprwtSyzPN?UZ1CI31uP-J73KrC-5olmpYD9e4tR&h-?23zujLQ!2Q6{GQ+(#r+CIr~Puyw$_)j}`Hfh`2`|7q!@XwE?e|5Ld#J*u1e%RgaIqeago&DYZW$qfrv4Gs|%G-_|n|#?o0(K4K`0)27VBandKWtx-_6W|-exOLN58HOP-+MD4h_?*fYD)jn+wPPp{Rr;!=gW`U6#L2!-FA22lz#RB9d^DmxC#b&5JY8g2LgQxcSi5i8i)B$ya?fU6TIR-%>$`3;Gt1}3OtD;AR`1A1TsPZ`4N;RL5>8WNv{vvb|(R9YJfVhz#}IjnU^g$*>*QD&nb=|k`PsNB2s&(cHp+V@c|Pr-@x-_UaDaSyu)Kx0B{(i_I*>;_epjQQ}=PAYp~w|zl8qP3{d$Gn_#|h_ILZ2$-5Ez(AWoP)fGT(8lX-jfXDo<IT5q7A1Ko61Gn8tx#`q$)2YFJ60m=(ZFdPlAdCC_`6>zd)CNNwy6rAO1Y{&m*mfu7<5NmacNhTmu$KUR3U@}4ewVWGUCbnR0I&E@^FS&d>GR?-b4-LU<1w{QF_nb<F6E}Xl;NM{4t8EhBJ%?8m4REMZj^grcEG4#=Df?|C`KIh-`w+WwxN0%eY`VE%pdP0qe-9c#Ip0Js|?&4o~|%EU~;^|o9*H0@^7{$#|tS6IO?Innu1Pn8+H2K3$q9sY;EA@<|je0&u-1@ZrhN!>VXQ@VB^3@8>A_Ov_XnONgHHs7-@ry4JI6-rko|+8}*%L4tP||&1J-5f|;IBU<#J(u9yc6NN2}8iD>xa<wd*MwXaL|Aq9JrIc6~@!NSh@V!8sl0UcHr&54Emwp8rH88i<=6NVjqrd#_iZu`)2?F+>r5p!kC`O<!IgLQ?HH^|~J@-{FyJF^V%HS*k>hLiGIa}gglovaUXhusLP`K8{;`XG14`s9nTMwQYDB5&*HUnw$pDhRc5#c*#6`ueogEy`-0M~jiic;TBMcA_6~+qKSZ7LX+dl-8*yc+%eKp9)Wm>Cd8(WPZg-BGalE>XUVH!lYM{f53-4;dBi;=%aJ;zNl=X(yfy>t3dH`ksFLkCpdbSo>UbNf0rH}aJOsSm>_Q!y5<2RTIbPX<YBjMUdc$eyjt8f4_Gq`k_W6!gJ_ajNKNZfkU?cUN0+piEt{j`Ibu|kerA?cTMr)t2t3y%Cw$GH>mD+CZS*Q_9-R)osDbPY|1^2fF>rd5^~`1pPkyw2_Tl_p8Yw%#yR-ua*iAHWDmt3v(El<M?sNqmz!NQ6TkIn#CBaEeaDc}R{Lgr8>B>bz-EYD_@v<ZO&hw`Tw4wFS;gm2nHB3!Syl1V|qy_>GuO>N~o#ol!tyq`5h(8Bby;pWP{k4&$o_0-vo;m~VzEDb_IhT{JTzktrWjOD@S?>5okNWxT-^VymzPC=DJe%cNx!?1s-$9rHns95QBUN*Ml>{!$xiE9^{wn!f*Zp>>brtPes-NHJK2nkw$=;_~pp^^G9NusWQ!4NG=bp|%I?YC{9R|{Q7|6ImI?W*6BqOffN4LuhZ#0lXu*)e??&B=oZ;D3TYtzn$fpk9*WV|5V%^-con^@}?a<O}LV1sCVIR3?a!(BybUSuwRp!ARPB~4r?b1#K%_p}yu31#-Ryh&DV!a2>~S)$j4`)EUv;=i{`_uD0T=>k1XUIs?vqM6Bybi)@xV#Zb@m2V!Bya>0VN%dK}l(TRre}*$}6+?O4jP%;3)9X+1*E&br3)q99hJ#1#!BBglSL>#h+K$d6wko96y+3Y`NH@!~dymNPoj%|H?|%TrNbFG',
    ))


# endregion

_ALL230_LABEL_RE = re.compile(r"^\s*(?:(\d+)\s*)?([A-Za-z])\s*$")
_ALL230_TWO_PI_I = 2j * cmath.pi


@lru_cache(maxsize=1)
def _all230_database() -> dict:
    compressed = base64.b85decode(_all230_database_b85().encode("ascii"))
    return json.loads(zlib.decompress(compressed).decode("utf-8"))


def _all230_validate_spacegroup(number: int) -> str:
    if not isinstance(number, int) or isinstance(number, bool) or not 1 <= number <= 230:
        raise ValueError("space-group number must be an integer from 1 through 230")
    return str(number)


def _all230_fraction(pair) -> Fraction:
    return Fraction(pair[0], pair[1])


def _all230_resolve_wyckoff(number: int, label: str) -> str:
    match = _ALL230_LABEL_RE.match(label)
    if not match:
        raise ValueError("Wyckoff label must look like '4a' or 'a'")
    multiplicity, letter = match.groups()
    entry = _all230_database()[_all230_validate_spacegroup(number)]
    positions = entry['positions']
    if letter not in positions:
        raise ValueError(f"Space group {number} has no Wyckoff position {label!r}")
    actual = int(positions[letter]['multiplicity'])
    if multiplicity is not None and int(multiplicity) != actual:
        raise ValueError(f"Wyckoff position {letter} has multiplicity {actual}, not {multiplicity}")
    return letter


def _all230_phase_terms(point, h: int, k: int, l: int):
    coefficients = [Fraction(0), Fraction(0), Fraction(0)]
    constant = Fraction(0)
    for index, coordinate in zip((h, k, l), point):
        for variable in range(3):
            coefficients[variable] += index * _all230_fraction(coordinate[variable])
        constant += index * _all230_fraction(coordinate[3])
    return tuple(coefficients), constant


@lru_cache(maxsize=262144)
def _all230_orbit_absent(number: int, letter: str, h: int, k: int, l: int) -> bool:
    orbit = _all230_database()[_all230_validate_spacegroup(number)]['positions'][letter]['points']
    grouped = {}
    for point in orbit:
        variable_key, constant = _all230_phase_terms(point, h, k, l)
        phase = cmath.exp(_ALL230_TWO_PI_I * float(constant % 1))
        grouped[variable_key] = grouped.get(variable_key, 0j) + phase
    return all(abs(amplitude) <= 1e-9 for amplitude in grouped.values())


def all230_reflection_allowed(
        number: int,
        h: int,
        k: int,
        l: int,
        occupied_wyckoff: Sequence[str] | str | None = None,
) -> bool:
    """Exact general or occupied-Wyckoff reflection condition from the all-230 database."""
    if h == k == l == 0:
        return False
    entry = _all230_database()[_all230_validate_spacegroup(number)]
    if occupied_wyckoff is None:
        labels = [entry['order_general_to_special'][0]]
    elif isinstance(occupied_wyckoff, str):
        labels = [_all230_resolve_wyckoff(number, occupied_wyckoff)]
    else:
        labels = [_all230_resolve_wyckoff(number, x) for x in occupied_wyckoff]
        if not labels:
            raise ValueError("occupied_wyckoff cannot be empty")
    return any(not _all230_orbit_absent(number, letter, int(h), int(k), int(l)) for letter in labels)


# ====================== ACCURACY-FOCUSED INDEXING CORE ======================

# ================= ANCHOR-FIRST CROSS-VALIDATED INDEXING CORE =================
# The core is independent of the embedded example data and is GUI-ready:
# configure -> run_gixs_indexing(config) -> DataFrames/dictionaries/plots.

@dataclass
class IndexingConfig:
    cif_path: str
    image_glob: str = "*waxs_stitched*.png"
    search_dirs: tuple[str, ...] = (".", "/content", "/mnt/data")
    input_zip: str | None = None
    colab_upload_if_missing: bool = True

    # Optional dataset manifest. When supplied, this replaces filename parsing and
    # hard-coded scan-block grouping. Relative file paths are resolved from the
    # manifest directory first, then from search_dirs.
    manifest_path: str | None = None
    auto_discover_manifest: bool = True
    manifest_filenames: tuple[str, ...] = (
        "gixs_dataset_manifest.csv",
        "dataset_manifest.csv",
        "measurement_manifest_input.csv",
    )
    manifest_strict: bool = True
    write_manifest_template_if_missing: bool = True
    manifest_template_filename: str = "gixs_dataset_manifest_template.csv"

    # Preview and reporting controls. Preview mode validates measurement inputs,
    # detects reciprocal-space features, and writes diagnostic thumbnails without
    # performing the orientation search.
    preview_only: bool = False
    write_input_preview_images: bool = True
    write_html_report: bool = True
    html_report_filename: str = "indexing_report.html"

    output_dir: str = "gixs_anchor_indexing_results"
    isolate_series_processes: bool = True
    series_worker_timeout_s: int = 600
    # Independent series are launched concurrently in fresh subprocesses.
    # 0 chooses a conservative automatic value based on available CPU cores.
    max_parallel_series_workers: int = 0

    # Numerical reciprocal-space arrays are preferred automatically when available.
    prefer_numerical: bool = True
    numerical_globs: tuple[str, ...] = ("*.npz", "*.npy", "*.h5", "*.hdf5", "*.tif", "*.tiff", "*.csv", "*.mat")
    intensity_keys: tuple[str, ...] = ("intensity", "image", "data", "I", "waxs", "giwaxs")
    qr_keys: tuple[str, ...] = ("qr", "q_r", "qxy", "qx", "x_axis")
    qz_keys: tuple[str, ...] = ("qz", "q_z", "z_axis", "y_axis")

    # PNG q-grid estimate. Numerical axes override these values.
    qr_range: tuple[float, float] = (-1.0, 2.2)
    qz_range: tuple[float, float] = (-0.10, 2.72)
    crop_xyxy: tuple[int, int, int, int] | None = None
    colormap: str = "jet"

    # Explicit acquisition blocks prevent repeat labels from being inferred only
    # from per-angle file order. Unlisted scans retain the deterministic fallback.
    explicit_series_scan_blocks: tuple[tuple[str, int, int], ...] = (
        ("s1:A", 2371795, 2371799),
        ("s1:B", 2371807, 2371811),
        ("s2:A", 2372427, 2372431),
        ("s3:A", 2372437, 2372441),
    )
    use_explicit_series_scan_blocks: bool = True

    # PNGs are cropped independently and registered within each angle series.
    # The transform maps each scan's q coordinates into the median-angle reference.
    reuse_first_png_crop: bool = False
    enable_per_image_registration: bool = True
    registration_png_only: bool = True
    registration_max_features: int = 45
    registration_pair_tolerance_q: float = 0.090
    registration_max_shift_q: float = 0.085
    registration_max_rotation_deg: float = 0.80
    registration_max_scale_change: float = 0.008
    registration_min_matches: int = 4
    registration_min_improvement_fraction: float = 0.08
    registration_iterations: int = 4
    registration_min_series_scatter_improvement_fraction: float = 0.03
    registration_require_no_consensus_support_loss: bool = True

    # Registered visual products. The main overlay uses a high-percentile
    # composite of every angle; separate overlays are also written per scan.
    registered_composite_percentile: float = 100.0
    write_registered_per_angle_overlays: bool = True

    # Fixed-orientation local pixel search. This examines actual registered PNG
    # intensity near unused CIF predictions and requires recurrence across angles.
    enable_local_pixel_completion: bool = True
    local_pixel_completion_png_only: bool = True
    local_pixel_prediction_f2_percentile: float = 0.0
    local_pixel_max_predictions_per_domain: int = 1818
    local_pixel_window_q: float = 0.032
    local_pixel_background_window_q: float = 0.070
    local_pixel_max_centroid_offset_q: float = 0.026
    local_pixel_min_snr: float = 3.2
    local_pixel_min_angle_support: int = 2
    local_pixel_min_three_angle_snr: float = 2.7
    local_pixel_max_prediction_ambiguity: int = 2
    local_pixel_min_prediction_margin_q: float = 0.006
    local_pixel_existing_feature_separation_q: float = 0.020
    local_pixel_candidate_separation_q: float = 0.018

    # X-ray and optional single-layer DWBA optical model. The DWBA implementation
    # is a coherent four-term Fresnel field envelope for intensities; it does not
    # claim to be a full dynamical GIWAXS solver and does not alter Bragg positions.
    xray_wavelength_A: float = 1.0
    enable_dwba: bool = False
    incidence_angle_deg: float = 0.12
    film_delta: float | None = None
    film_beta: float | None = None
    substrate_delta: float | None = None
    substrate_beta: float | None = None
    film_thickness_A: float | None = None
    dwba_scatter_depth_fraction: float = 0.5
    dwba_intensity_floor: float = 0.02

    # GIWAXS physics is selected automatically from measurement metadata.  The
    # program never guesses missing optical constants: when the required data
    # are unavailable it falls back to the kinematic model without stopping the
    # indexing run.
    giwaxs_physics_mode: str = "auto"
    giwaxs_physics_status: str = "auto_pending"
    giwaxs_physics_parameter_source: str = ""
    giwaxs_physics_fallback_reason: str = ""
    giwaxs_incidence_angles_deg: tuple[float, ...] = ()

    # Optional direct-channel refraction correction for predicted qz positions.
    # This maps the internal film reciprocal-vector component to the external
    # detector qz using the supplied wavelength, incidence angle, film delta,
    # and film beta. It is disabled until calibrated optical constants are set.
    enable_refraction_position_correction: bool = False
    refraction_position_channel: str = "direct_direct"

    # Kinematic powder-pattern diagnostic for radial phase/CIF compatibility.
    # Agreement supports structural compatibility but does not uniquely prove that
    # the selected CIF is the correct crystal structure.
    simulate_powder: bool = True
    powder_q_min: float = 0.05
    powder_q_max: float | None = None
    powder_q_step: float = 0.002
    powder_fwhm_q: float = 0.018
    powder_peak_merge_q: float = 0.003
    powder_min_relative_peak: float = 0.01
    powder_apply_lorentz_polarization: bool = False
    alternative_cif_paths: tuple[str, ...] = ()
    cif_radial_match_tolerance_q: float = 0.045
    cif_preference_min_gap: float = 0.08

    # Feature extraction and artifact masking.
    analysis_qz_min: float = 0.15
    exclude_specular_abs_qr: float = 0.07
    exclude_specular_qz_max: float = 1.0
    feature_sigma_px: float = 1.1
    background_sigma_px: float = 13.0
    feature_threshold_mad: float = 3.2
    feature_quantile: float = 0.982
    ridge_threshold_mad: float = 2.5
    min_feature_spacing_px: int = 6
    subpixel_radius_px: int = 3
    component_min_pixels: int = 7
    component_max_pixels: int = 3500
    arc_aspect_ratio: float = 2.4
    feature_merge_tolerance_q: float = 0.020
    max_features_per_image: int = 80

    # Multi-angle consensus.
    consensus_tolerance_q: float = 0.040
    min_angle_support: int = 2
    max_consensus_features: int = 70

    # Reciprocal lattice and dual systematic-absence validation.
    q_max: float = 2.8
    orientation_hkl_max: int = 5
    max_orientation_candidates: int = 70
    all230_compare: bool = True
    all230_policy: str = "gemmi"  # gemmi | agreement | all230
    occupied_wyckoff: tuple[str, ...] | None = None
    structure_factor_zero_percent: float = 1e-8

    # Uncertainty-aware assignment.
    uncertainty_floor_q: float = 0.007
    uncertainty_ceiling_q: float = 0.065
    match_sigma_limit: float = 3.2

    # Anchor/validation separation. Validation membership is selected by a
    # deterministic data-only holdout before any CIF radial prior is evaluated.
    validation_holdout_fraction: float = 0.35
    validation_holdout_min_features: int = 6
    anchor_min_support: int = 3
    anchor_min_qz: float = 0.24
    anchor_min_abs_qr: float = 0.22
    anchor_max_sigma_q: float = 0.040
    anchor_max_major_width_q: float = 0.12
    anchor_strength_quantile: float = 0.45
    anchor_radial_tolerance_q: float = 0.045
    max_anchor_features: int = 8
    max_validation_features: int = 28
    ignored_feature_types: tuple[str, ...] = ("radial_streak",)

    # Conservative post-fit recovery of consensus features that were not used as
    # anchors or holdout validation.  These assignments never influence the
    # selected orientation and are exported as provisional evidence.
    index_ignored_features: bool = True
    ignored_index_f2_percentile: float = 20.0
    ignored_index_max_predictions: int = 450
    ignored_index_tolerance_q: float = 0.050
    ignored_index_sigma_limit: float = 2.6
    ignored_index_min_support: int = 2
    ignored_index_max_sigma_q: float = 0.060
    ignored_index_max_major_width_q: float = 0.20
    ignored_index_max_ambiguity: int = 2
    ignored_index_min_margin_sigma: float = 0.20
    ignored_index_min_support_score: float = 0.015
    ignored_index_include_radial_streaks: bool = True
    ignored_streak_min_support: int = 3
    ignored_streak_sigma_limit: float = 2.0
    ignored_streak_max_ambiguity: int = 2

    # Independent evidence checks for ignored-feature assignments.  These do not
    # change the orientation; they grade each post-fit assignment according to
    # angle-specific corroboration, local calibration stability, and specificity
    # relative to deliberately incorrect orientation decoys.
    ignored_evidence_member_gate_q: float = 0.045
    ignored_evidence_member_sigma_limit: float = 2.4
    ignored_evidence_min_member_fraction: float = 0.60
    ignored_evidence_perturbation_trials: int = 13
    ignored_evidence_tilt_jitter_deg: float = 0.35
    ignored_evidence_scale_jitter: float = 0.004
    ignored_evidence_offset_jitter_q: float = 0.006
    ignored_evidence_min_projection_stability: float = 0.65
    ignored_evidence_decoy_trials: int = 12
    ignored_evidence_decoy_min_angle_deg: float = 12.0
    ignored_evidence_decoy_margin_sigma: float = 0.10
    ignored_evidence_min_decoy_win_fraction: float = 0.70

    # Prediction-guided rescue of angle-specific features that did not survive
    # consensus construction.  A reflection must recur at the same calculated
    # position in multiple independent angles, so isolated single-image noise is
    # not promoted into the indexed table.
    guided_rescue_unclustered_features: bool = True
    guided_rescue_tolerance_q: float = 0.035
    guided_rescue_sigma_limit: float = 1.9
    guided_rescue_min_angle_support: int = 2
    guided_rescue_max_ambiguity: int = 1
    guided_rescue_min_margin_sigma: float = 0.45
    guided_rescue_max_major_width_q: float = 0.16
    # Promote recurring raw detections when they are independently observed in
    # multiple angle images. These points were found from image intensity before
    # CIF assignment, but failed ordinary consensus clustering.
    guided_rescue_supported_min_angle_support: int = 3
    guided_rescue_supported_max_normalized_delta: float = 1.90
    guided_rescue_supported_min_margin_sigma: float = 0.80
    guided_rescue_promote_strong_two_angle: bool = True
    guided_rescue_two_angle_max_normalized_delta: float = 1.75
    guided_rescue_two_angle_max_delta_q: float = 0.014
    guided_rescue_two_angle_min_margin_sigma: float = 1.00
    guided_rescue_two_angle_min_strength: float = 0.001
    guided_rescue_two_angle_excluded_types: tuple[str, ...] = ("radial_streak",)

    # Pair/triplet RANSAC-style hypothesis generation.
    hypothesis_f2_percentile: float = 88.0
    validation_f2_percentile: float = 55.0
    max_hypothesis_predictions: int = 50
    max_validation_predictions: int = 140
    max_pair_hypotheses_per_normal: int = 7
    anchor_match_tolerance_q: float = 0.042
    validation_match_tolerance_q: float = 0.055
    # A second, tighter holdout score prevents a larger loose-tolerance match
    # count from automatically outranking a geometrically cleaner solution.
    strict_validation_tolerance_q: float = 0.042
    min_anchor_matches: int = 3
    min_validation_matches: int = 2
    max_normal_candidates_anchor: int = 45
    refine_hypotheses: int = 10
    expand_top_normals: int = 12
    orientation_family_merge_angle_deg: float = 0.35
    orientation_pattern_merge_q: float = 0.018
    orientation_ambiguity_score_delta: float = 0.05
    orientation_min_unique_separation_deg: float = 3.0
    orientation_min_unique_score_gap: float = 0.06
    orientation_min_bootstrap_iterations: int = 12
    orientation_min_stability_fraction: float = 0.80
    orientation_min_loo_weighted_fraction: float = 0.45
    coarse_tilt_values_deg: tuple[float, ...] = (-6.0, 0.0, 6.0)

    # Tight calibration model: two tilts, separate qr/qz scales, two offsets.
    # Separate detector-axis scales are an empirical correction for anisotropic
    # cell/geometry mismatch; they are not reported as independently refined
    # crystallographic lattice constants.
    enable_anisotropic_q_scale: bool = True
    max_tilt_anchor_deg: float = 7.0
    max_common_scale_change: float = 0.020  # Shared q-scale limit for five-parameter calibration vectors.
    max_qr_scale_change: float = 0.025
    max_qz_scale_change: float = 0.030
    max_anchor_q_offset: float = 0.038
    anchor_regularization: float = 0.025

    # Cross-validated scoring.
    missing_strong_prediction_penalty: float = 0.10
    ambiguity_penalty: float = 0.07
    pair_geometry_penalty: float = 0.12

    # Full independent re-search validation.
    full_leave_one_angle_out: bool = True
    full_bootstrap_iterations: int = 4
    orientation_stability_angle_deg: float = 4.0
    bootstrap_search_normal_limit: int = 22
    loo_search_normal_limit: int = 26

    # Conservative residual second-orientation test.
    test_second_orientation: bool = True
    second_orientation_min_normal_separation_deg: float = 5.0

    # Manual/GUI hooks.
    manual_anchor_feature_ids: tuple[str, ...] = ()
    manual_rejected_feature_ids: tuple[str, ...] = ()
    manual_locked_assignments: tuple[tuple[str, int, int, int], ...] = ()

    # Optional external orientation truth. Each entry is
    # (series_id, h, k, l). Results are called externally tested only when an
    # entry is supplied; synthetic and bootstrap checks remain internal checks.
    external_truth_normals: tuple[tuple[str, int, int, int], ...] = ()
    external_truth_source: str = ""
    external_truth_was_blind: bool = False
    external_truth_tolerance_deg: float = 3.0

    # Robust orientation-selection controls. Usable lower-confidence features can
    # contribute through an explicit outlier-mixture model, while broad or streak-like
    # features receive reduced weight. Features without sufficient positional support
    # may remain unassigned instead of being forced onto a calculated reflection.
    v7_enable_core_solver: bool = True
    v7_fit_f2_percentile: float = 30.0
    v7_max_fit_predictions: int = 320
    v7_all_feature_tolerance_q: float = 0.052
    v7_outlier_cost: float = 0.98
    v7_ignored_feature_weight: float = 0.30
    v7_broad_feature_weight: float = 0.48
    v7_min_feature_quality: float = 0.08
    v7_all_feature_score_weight: float = 0.72
    v7_validation_score_weight: float = 0.28
    v7_angle_consistency_weight: float = 0.10
    v7_max_axis_rotation_deg: float = 2.0
    v7_max_shear: float = 0.035
    v7_affine_regularization: float = 0.070
    v7_refine_hypotheses: int = 14
    v7_refine_cycles: int = 3
    v7_multidomain_max_domains: int = 3
    v7_multidomain_candidate_normals: int = 36
    v7_multidomain_tilt_values_deg: tuple[float, ...] = (-6.0, 0.0, 6.0)
    v7_multidomain_min_gain: float = 0.030
    v7_multidomain_complexity_penalty: float = 0.022
    v7_multidomain_min_matches: int = 3

    # Fixed-orientation reflection completion. Orientation is first selected from
    # an intensity-pruned reflection list for computational efficiency and contrast.
    # Once fixed, residual consensus features are tested against the full visible set
    # of symmetry-allowed reflections so geometrically valid reflections are not
    # excluded solely by the F² intensity cutoff.
    v71_enable_full_reflection_completion: bool = True
    v71_completion_f2_percentile: float = 0.0
    v71_completion_max_predictions_per_domain: int = 1818
    v71_completion_base_tolerance_q: float = 0.052
    v71_completion_max_tolerance_q: float = 0.066
    v71_completion_sigma_limit: float = 3.2
    v71_completion_min_support: int = 2
    v71_completion_min_member_fraction: float = 0.50
    v71_completion_min_member_angles: int = 2
    v71_completion_max_ambiguity: int = 4
    v71_completion_min_margin_sigma: float = 0.10
    v71_completion_max_sigma_q: float = 0.065
    v71_completion_max_major_width_q: float = 0.20
    v71_completion_streak_min_support: int = 3
    v71_completion_streak_sigma_limit: float = 2.4
    v71_completion_streak_max_ambiguity: int = 2

    # Mosaic-reflection completion for weak measurements selected by the configured
    # series prefixes (default: ``s1:``). Residual multi-angle features are compared
    # with small tilt perturbations around each accepted orientation domain. This can
    # recover split or mosaic peaks from one reflection family while limiting
    # unrestricted duplicate assignments.
    v72_enable_s1_mosaic_completion: bool = True
    v72_s1_series_prefixes: tuple[str, ...] = ("s1:",)
    v72_mosaic_only_when_weighted_fraction_below: float = 0.55
    v72_mosaic_tilt_offsets_deg: tuple[tuple[float, float], ...] = (
        (-1.50, 0.00), (-0.75, 0.00), (0.75, 0.00), (1.50, 0.00),
        (0.00, -1.50), (0.00, -0.75), (0.00, 0.75), (0.00, 1.50),
        (-0.75, -0.75), (-0.75, 0.75), (0.75, -0.75), (0.75, 0.75),
    )
    v72_mosaic_f2_percentile: float = 0.0
    v72_mosaic_max_predictions_per_variant: int = 1818
    v72_mosaic_tolerance_q: float = 0.055
    v72_mosaic_sigma_limit: float = 2.8
    v72_mosaic_min_support: int = 3
    v72_mosaic_min_member_angles: int = 2
    v72_mosaic_min_member_fraction: float = 0.50
    v72_mosaic_max_unique_ambiguity: int = 4
    v72_mosaic_min_margin_sigma: float = 0.08
    v72_mosaic_max_sigma_q: float = 0.060
    v72_mosaic_max_major_width_q: float = 0.16
    v72_mosaic_max_assignments_per_reflection: int = 3
    v72_mosaic_variant_separation_q: float = 0.006
    v72_mosaic_tilt_penalty: float = 0.10

    # Sector-aware residual-domain search for weak measurements selected by the
    # configured series prefixes (default: ``s1:A``). Candidate domains are evaluated
    # in local reciprocal-space sectors and admitted to the global one-to-one
    # assignment only when they improve the robust whole-pattern score.
    v73_enable_s1_sector_domains: bool = True
    v73_s1_sector_series_prefixes: tuple[str, ...] = ("s1:A",)
    v73_sector_only_when_weighted_fraction_below: float = 0.50
    v73_sector_max_extra_domains: int = 2
    v73_sector_candidate_normals: int = 30
    v73_sector_tilt_values_deg: tuple[float, ...] = (-6.0, 0.0, 6.0)
    v73_sector_search_f2_percentile: float = 35.0
    v73_sector_search_max_predictions: int = 360
    v73_sector_top_candidates_per_sector: int = 4
    v73_sector_min_features: int = 4
    v73_sector_min_local_matches: int = 4
    v73_sector_min_global_domain_matches: int = 6
    v73_sector_min_new_matches: int = 4
    v73_sector_min_local_global_matches: int = 3
    v73_sector_min_score_gain: float = 0.025
    v73_sector_min_weighted_fraction_gain: float = 0.035
    v73_sector_complexity_penalty: float = 0.012
    v73_sector_max_lost_assignments: int = 1
    v73_preserve_pre_sector_completion: bool = True
    v73_sector_min_normal_separation_deg: float = 8.0
    v73_sector_tolerance_q: float = 0.055
    v73_sector_min_support: int = 2
    v73_sector_max_sigma_q: float = 0.070
    v73_sector_max_major_width_q: float = 0.20
    v73_sector_qr_split: float = 0.12
    v73_sector_qz_edges: tuple[float, float] = (0.90, 1.75)

    # Residual candidate-phase/substrate screening. Additional CIFs may be
    # supplied through alternative_cif_paths or substrate_cif_paths. Automatic
    # discovery looks for other CIFs beside the inputs and in common notebook
    # upload directories. Candidate assignments are diagnostic by default and
    # are never promoted unless the explicit promotion flag is enabled.
    substrate_cif_paths: tuple[str, ...] = ()
    v73_auto_discover_candidate_cifs: bool = True
    v73_candidate_cif_globs: tuple[str, ...] = ("*.cif", "*.CIF")
    v73_residual_phase_min_features: int = 4
    v73_residual_phase_min_support: int = 2
    v73_residual_phase_max_sigma_q: float = 0.080
    v73_residual_phase_max_major_width_q: float = 0.24
    v73_residual_phase_candidate_normals: int = 14
    v73_residual_phase_refine_hypotheses: int = 3
    v73_residual_phase_min_matches: int = 4
    v73_residual_phase_min_score_gap: float = 0.04
    v73_promote_residual_candidate_matches: bool = False

    random_seed: int = 7
    top_candidates_to_print: int = 6
    # Maximum overlay-label count used when labeling every indexed feature is disabled.
    max_labels: int = 28
    # Every indexed consensus feature is labeled on the image.  A numbered key
    # on the right gives the exact experimental reciprocal-space coordinates.
    overlay_label_all_indexed: bool = True
    overlay_show_coordinate_key: bool = True
    overlay_label_include_hkl: bool = True
    overlay_coordinate_decimals: int = 3
    overlay_label_fontsize: float = 5.2
    overlay_coordinate_key_fontsize: float = 5.5
    overlay_coordinate_key_columns: int = 2
    # Presentation-only gap handling. Invalid detector/plot regions remain masked
    # for feature detection, fitting, completion, and validation. The overlay may
    # linearly bridge internal blank runs so the displayed map is continuous.
    overlay_fill_display_gaps: bool = False
    overlay_gap_fill_max_fraction: float = 0.35
    overlay_gap_fill_smoothing_sigma_px: float = 0.0
    overlay_note_display_gap_fill: bool = True
    overlay_panel_spacing: float = 0.08
    # Only evidence-supported post-fit rescues are promoted to indexed in plots.
    overlay_promoted_rescue_tiers: tuple[str, ...] = ("supported", "robust")
    dpi: int = 175


_ENGINE_FILENAME_RE = re.compile(
    r"(?:^|[_-])s(?P<sample>\d+).*?(?:th|angle|ai)(?P<angle>-?\d+(?:\.\d+)?)"
    r".*?(?P<exposure>\d+(?:\.\d+)?)s.*?(?P<scan>\d+)", re.IGNORECASE
)
_SCAN_RE = re.compile(r"(?<!\d)(?P<scan>\d{6,9})(?!\d)")


def _search_roots(config: IndexingConfig) -> list[Path]:
    roots = [Path(p).expanduser() for p in config.search_dirs if Path(p).expanduser().exists()]
    archives = [Path(config.input_zip).expanduser()] if config.input_zip else []
    if not archives:
        for root in roots:
            archives.extend(root.rglob("*gixs*input*.zip"))
    for archive in archives:
        if not archive.is_file():
            continue
        target = Path(config.output_dir) / "_input_cache"
        target.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(archive) as zf:
            base = target.resolve()
            for member in zf.infolist():
                destination = (target / member.filename).resolve()
                if base not in destination.parents and destination != base:
                    raise ValueError(f"Unsafe path in ZIP: {member.filename}")
            zf.extractall(target)
        roots.insert(0, target)
    return list(dict.fromkeys(p.resolve() for p in roots))


def _resolve_file(path: str, roots: Sequence[Path]) -> Path:
    candidate = Path(path).expanduser()
    if candidate.is_file():
        return candidate.resolve()
    for root in roots:
        direct = root / candidate.name
        if direct.is_file():
            return direct.resolve()
        matches = list(root.rglob(candidate.name))
        if matches:
            return matches[0].resolve()
    raise FileNotFoundError(f"File not found: {path}")


def _colab_upload() -> bool:
    try:
        from google.colab import files
    except ImportError:
        return False
    print("Upload the CIF plus q-space arrays/PNGs, or one input ZIP archive.")
    uploaded = files.upload()
    for name, data in uploaded.items():
        Path(name).write_bytes(data)
    return bool(uploaded)


def _parse_measurement_name(path: Path) -> dict | None:
    match = _ENGINE_FILENAME_RE.search(path.name.replace(" ", ""))
    if not match:
        return None
    row = match.groupdict()
    return {
        "sample": f"s{row['sample']}", "angle_deg": float(row["angle"]),
        "exposure_s": float(row["exposure"]), "scan": int(row["scan"]),
    }


def _v739_manifest_template_frame() -> pd.DataFrame:
    """Return a two-row disabled example manifest suitable for editing."""
    columns = [
        "file", "png_file", "numerical_file", "sample", "series", "series_id",
        "angle_deg", "scan", "exposure_s", "qr_min", "qr_max", "qz_min", "qz_max",
        "colormap", "crop_x0", "crop_y0", "crop_x1", "crop_y1", "enabled", "notes",
    ]
    rows = [
        {
            "file": "images/sample1_angle_0p10.png", "sample": "sample1", "series": "A",
            "angle_deg": 0.10, "scan": 1001, "exposure_s": 10.0,
            "qr_min": -1.0, "qr_max": 2.2, "qz_min": -0.10, "qz_max": 2.72,
            "colormap": "jet", "enabled": False,
            "notes": "Disabled example: set enabled=True and replace the path.",
        },
        {
            "numerical_file": "arrays/sample1_angle_0p15.npz",
            "png_file": "images/sample1_angle_0p15.png",
            "sample": "sample1", "series": "A", "angle_deg": 0.15,
            "scan": 1002, "exposure_s": 10.0, "enabled": False,
            "notes": "NPZ axes override qr/qz limits; PNG is optional for display.",
        },
    ]
    return pd.DataFrame(rows, columns=columns)


def write_dataset_manifest_template(path: str | Path) -> Path:
    """Write an editable CSV manifest template and return its absolute path."""
    target = Path(path).expanduser()
    if not target.is_absolute():
        target = Path.cwd() / target
    target.parent.mkdir(parents=True, exist_ok=True)
    _manifest_template_frame().to_csv(target, index=False)
    return target.resolve()


def _manifest_enabled(value) -> bool:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return True
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    return str(value).strip().lower() not in {"0", "false", "no", "off", "disabled"}


def _manifest_text(value) -> str | None:
    if value is None or pd.isna(value):
        return None
    value = str(value).strip()
    return value or None


def _manifest_number(value, default=np.nan):
    if value is None or pd.isna(value) or str(value).strip() == "":
        return default
    return float(value)


def _resolve_manifest_file(value, manifest_dir: Path, roots: Sequence[Path]) -> Path | None:
    text = _manifest_text(value)
    if not text:
        return None
    candidate = Path(text).expanduser()
    attempts = []
    if candidate.is_absolute():
        attempts.append(candidate)
    else:
        attempts.extend([manifest_dir / candidate, *[root / candidate for root in roots]])
    for attempt in attempts:
        if attempt.is_file():
            return attempt.resolve()
    for root in roots:
        matches = list(root.rglob(candidate.name))
        if matches:
            return matches[0].resolve()
    raise FileNotFoundError(f"Manifest file entry not found: {text}")


def _find_dataset_manifest(config: IndexingConfig, roots: Sequence[Path]) -> Path | None:
    if config.manifest_path:
        candidate = Path(config.manifest_path).expanduser()
        attempts = [candidate] if candidate.is_absolute() else [Path.cwd() / candidate, *[r / candidate for r in roots]]
        for attempt in attempts:
            if attempt.is_file():
                return attempt.resolve()
        if config.write_manifest_template_if_missing:
            template_target = candidate if candidate.suffix.lower() == ".csv" else Path(
                config.manifest_template_filename)
            written = write_dataset_manifest_template(template_target)
            raise FileNotFoundError(
                f"Dataset manifest was not found: {config.manifest_path}. "
                f"An editable template was written to {written}"
            )
        raise FileNotFoundError(f"Dataset manifest was not found: {config.manifest_path}")
    if not config.auto_discover_manifest:
        return None
    matches = []
    for root in roots:
        for name in config.manifest_filenames:
            matches.extend(root.rglob(name))
    matches = sorted({p.resolve() for p in matches}, key=lambda p: (len(p.parts), str(p)))
    if len(matches) > 1:
        print(f"Multiple dataset manifests found; using {matches[0]}")
    return matches[0] if matches else None


def _normalize_manifest_columns(frame: pd.DataFrame) -> pd.DataFrame:
    aliases = {
        "filename": "file", "filepath": "file", "data_path": "file", "data_file": "file",
        "image": "png_file", "image_file": "png_file", "png": "png_file", "path_png": "png_file",
        "array": "numerical_file", "array_file": "numerical_file", "numerical_path": "numerical_file",
        "incident_angle": "angle_deg", "incident_angle_deg": "angle_deg", "alpha_i": "angle_deg",
        "repeat": "series", "group": "series", "sample_id": "sample",
        "scan_id": "scan", "exposure": "exposure_s",
    }
    renamed = {}
    for column in frame.columns:
        clean = str(column).strip().lower().replace(" ", "_")
        renamed[column] = aliases.get(clean, clean)
    return frame.rename(columns=renamed)


def _load_dataset_manifest(path: Path, config: IndexingConfig, roots: Sequence[Path]) -> pd.DataFrame:
    raw = _normalize_manifest_columns(pd.read_csv(path))
    if raw.empty:
        raise ValueError(f"Dataset manifest contains no rows: {path}")
    required_any = {"file", "png_file", "numerical_file"}
    if not required_any.intersection(raw.columns):
        raise ValueError("Manifest needs one of: file, png_file, numerical_file")
    if "angle_deg" not in raw.columns:
        raise ValueError("Manifest requires an angle_deg column")
    if "sample" not in raw.columns and "series_id" not in raw.columns:
        raise ValueError("Manifest requires sample or series_id")

    manifest_dir = path.parent
    rows = []
    auto_scan = 1
    numerical_suffixes = {".npz"}
    for source_index, source in raw.iterrows():
        if "enabled" in raw.columns and not _manifest_enabled(source.get("enabled")):
            continue
        generic = _resolve_manifest_file(source.get("file"), manifest_dir, roots) if "file" in raw.columns else None
        png = _resolve_manifest_file(source.get("png_file"), manifest_dir, roots) if "png_file" in raw.columns else None
        numerical = _resolve_manifest_file(source.get("numerical_file"), manifest_dir,
                                           roots) if "numerical_file" in raw.columns else None
        if generic:
            if generic.suffix.lower() == ".png":
                png = png or generic
            elif generic.suffix.lower() in numerical_suffixes:
                numerical = numerical or generic
            else:
                raise ValueError(f"Unsupported manifest file type: {generic}")
        if not png and not numerical:
            raise ValueError(f"Manifest row {source_index + 2} has no usable input file")

        parsed = _parse_measurement_name(png or numerical) or {}
        series_id = _manifest_text(source.get("series_id")) if "series_id" in raw.columns else None
        sample = _manifest_text(source.get("sample")) if "sample" in raw.columns else None
        series = _manifest_text(source.get("series")) if "series" in raw.columns else None
        if series_id:
            if ":" in series_id:
                sample_from_id, series_from_id = series_id.split(":", 1)
                sample = sample or sample_from_id
                series = series or series_from_id
            else:
                sample = sample or series_id
                series = series or "A"
                series_id = f"{sample}:{series}"
        else:
            sample = sample or parsed.get("sample")
            series = series or "A"
            if not sample:
                raise ValueError(f"Manifest row {source_index + 2} needs sample or series_id")
            series_id = f"{sample}:{series}"

        angle = _manifest_number(source.get("angle_deg"))
        if not np.isfinite(angle):
            raise ValueError(f"Manifest row {source_index + 2} has invalid angle_deg")
        scan_value = source.get("scan") if "scan" in raw.columns else None
        if scan_value is None or pd.isna(scan_value) or str(scan_value).strip() == "":
            scan = int(parsed.get("scan", auto_scan))
        else:
            scan = int(float(scan_value))
        auto_scan = max(auto_scan + 1, scan + 1)
        exposure = _manifest_number(source.get("exposure_s") if "exposure_s" in raw.columns else None,
                                    parsed.get("exposure_s", 1.0))
        source_kind = "numerical" if config.prefer_numerical and numerical else "png"
        row = {
            "sample": str(sample), "angle_deg": float(angle), "exposure_s": float(exposure),
            "scan": int(scan), "path": str(png) if png else None,
            "numerical_path": str(numerical) if numerical else None,
            "source_kind": source_kind, "series": str(series), "series_id": str(series_id),
            "series_grouping": "dataset_manifest", "manifest_row": int(source_index + 2),
            "input_manifest_path": str(path),
        }
        for key in ("qr_min", "qr_max", "qz_min", "qz_max"):
            row[key] = _manifest_number(source.get(key) if key in raw.columns else None)
        row["row_colormap"] = _manifest_text(source.get("colormap")) if "colormap" in raw.columns else None
        for key in ("crop_x0", "crop_y0", "crop_x1", "crop_y1"):
            value = source.get(key) if key in raw.columns else None
            row[key] = int(float(value)) if value is not None and not pd.isna(value) and str(value).strip() else np.nan
        # Optional crop-relative pixel anchors for an explicit linear pixel-to-q calibration.
        for key in (
                "qr_pixel_0", "qr_value_0", "qr_pixel_1", "qr_value_1",
                "qz_pixel_0", "qz_value_0", "qz_pixel_1", "qz_value_1",
        ):
            row[key] = _manifest_number(source.get(key) if key in raw.columns else None)
        row["notes"] = _manifest_text(source.get("notes")) if "notes" in raw.columns else None
        rows.append(row)
    if not rows:
        raise ValueError(f"Dataset manifest has no enabled rows: {path}")
    return pd.DataFrame(rows).sort_values(["series_id", "angle_deg", "scan"]).reset_index(drop=True)


def _v738_validate_measurement_manifest(frame: pd.DataFrame, config: IndexingConfig) -> pd.DataFrame:
    """Return a machine-readable preflight report for discovered or manifested inputs."""
    records = []

    def add(severity, check, message, series_id=""):
        records.append({"severity": severity, "check": check, "series_id": series_id, "message": message})

    if frame.empty:
        add("error", "measurements", "No enabled measurements were found.")
        return pd.DataFrame(records)
    duplicates = frame.duplicated(["series_id", "angle_deg"], keep=False)
    for series_id, group in frame.groupby("series_id", sort=False):
        count = len(group)
        if count < 2:
            add("error", "series_size", f"Only {count} measurement; at least two angles are required.", series_id)
        elif count < 3:
            add("warning", "series_size", f"Only {count} measurements; orientation validation will be weak.", series_id)
        else:
            add("ok", "series_size", f"{count} measurements available.", series_id)
        if duplicates.loc[group.index].any():
            repeated = sorted(group.loc[duplicates.loc[group.index], "angle_deg"].unique())
            add("warning", "duplicate_angles", f"Repeated angle values inside one series: {repeated}", series_id)
        if group["source_kind"].eq("png").any():
            add("ok", "png_calibration",
                f"PNG rows use manifest row limits when supplied, otherwise {config.qr_range}/{config.qz_range}.",
                series_id)
        if group["source_kind"].eq("numerical").any():
            add("ok", "numerical_inputs", "Numerical reciprocal-space input is preferred where available.", series_id)
    if frame["scan"].duplicated().any():
        scans = sorted(frame.loc[frame["scan"].duplicated(False), "scan"].unique())
        add("warning", "duplicate_scan_ids", f"Duplicate scan identifiers: {scans}")
    add("ok", "manifest_source", str(frame.get("input_manifest_path", pd.Series(["automatic discovery"])).iloc[0]))
    return pd.DataFrame(records)


def _apply_manifest_coordinate_defaults(config: IndexingConfig, frame: pd.DataFrame) -> IndexingConfig:
    """Use the union of complete manifest q limits as the run-level display/search range."""
    needed = ["qr_min", "qr_max", "qz_min", "qz_max"]
    if not set(needed).issubset(frame.columns):
        return config
    finite = frame[needed].apply(pd.to_numeric, errors="coerce")
    complete = finite.notna().all(axis=1)
    if not complete.any():
        return config
    subset = finite.loc[complete]
    qr_range = (float(subset.qr_min.min()), float(subset.qr_max.max()))
    qz_range = (float(subset.qz_min.min()), float(subset.qz_max.max()))
    return replace(config, qr_range=qr_range, qz_range=qz_range)


def discover_images(config: IndexingConfig, allow_upload: bool = True) -> pd.DataFrame:
    """Discover measurements from an explicit dataset manifest or, when needed, from the supported filename metadata pattern."""
    roots = _search_roots(config)
    manifest_path = _find_dataset_manifest(config, roots)
    if manifest_path:
        frame = _load_dataset_manifest(manifest_path, config, roots)
        print(f"Using dataset manifest: {manifest_path}")
        return frame

    pngs: dict[int, Path] = {}
    metadata: dict[int, dict] = {}
    numerical: dict[int, Path] = {}
    for root in roots:
        for path in root.rglob(config.image_glob):
            row = _parse_measurement_name(path)
            if not row:
                continue
            scan = row["scan"]
            if scan not in pngs or ("(1)" in pngs[scan].name and "(1)" not in path.name):
                pngs[scan], metadata[scan] = path.resolve(), row
        for pattern in config.numerical_globs:
            for path in root.rglob(pattern):
                row = _parse_measurement_name(path)
                scan_match = _SCAN_RE.search(path.name)
                scan = row["scan"] if row else (int(scan_match.group("scan")) if scan_match else None)
                if scan is None:
                    continue
                numerical.setdefault(scan, path.resolve())
                if row:
                    metadata.setdefault(scan, row)
    scans = sorted(set(pngs) | set(numerical))
    if not scans and allow_upload and config.colab_upload_if_missing and _colab_upload():
        return discover_images(config, allow_upload=False)
    columns = [
        "sample", "angle_deg", "exposure_s", "scan", "path", "numerical_path",
        "source_kind", "series", "series_id", "series_grouping",
    ]
    rows = []
    for scan in scans:
        if scan not in metadata:
            continue
        num = numerical.get(scan)
        png = pngs.get(scan)
        source_kind = "numerical" if config.prefer_numerical and num else "png"
        rows.append({
            **metadata[scan], "path": str(png) if png else None,
            "numerical_path": str(num) if num else None, "source_kind": source_kind,
        })
    if not rows:
        return pd.DataFrame(columns=columns)
    frame = pd.DataFrame(rows).sort_values(["sample", "angle_deg", "scan"]).reset_index(drop=True)
    frame["series_id"] = pd.Series(pd.NA, index=frame.index, dtype="string")
    frame["series_grouping"] = "fallback_repeat_order"
    if config.use_explicit_series_scan_blocks:
        for series_id, first_scan, last_scan in config.explicit_series_scan_blocks:
            sample = str(series_id).split(":", 1)[0]
            mask = (frame["sample"].astype(str) == sample) & frame["scan"].between(int(first_scan), int(last_scan))
            frame.loc[mask, "series_id"] = str(series_id)
            frame.loc[mask, "series_grouping"] = "explicit_scan_block"
    remaining = frame["series_id"].isna()
    if remaining.any():
        fallback = frame.loc[remaining].copy()
        fallback["repeat"] = fallback.groupby(["sample", "angle_deg"]).cumcount()
        fallback["series"] = fallback["repeat"].map(lambda i: chr(ord("A") + int(i)))
        frame.loc[remaining, "series_id"] = fallback["sample"] + ":" + fallback["series"]
    frame["series"] = frame["series_id"].astype(str).str.split(":").str[-1]
    return frame.reset_index(drop=True)


def detect_plot_crop(rgb: np.ndarray) -> tuple[int, int, int, int]:
    maximum, minimum = rgb.max(axis=2), rgb.min(axis=2)
    saturation = (maximum - minimum) / np.maximum(maximum, 1e-6)
    components, count = label((saturation > 0.10) & (maximum < 0.995))
    if not count:
        raise ValueError("Could not locate the colored q-space panel.")
    sizes = np.bincount(components.ravel());
    sizes[0] = 0
    y, x = np.where(components == sizes.argmax())
    return int(x.min()), int(y.min()), int(x.max() + 1), int(y.max() + 1)


def _invert_colormap(rgb: np.ndarray, name: str) -> tuple[np.ndarray, np.ndarray]:
    lut = plt.get_cmap(name)(np.linspace(0, 1, 512))[:, :3]
    distance, index = cKDTree(lut).query(rgb.reshape(-1, 3), workers=-1)
    return (index.reshape(rgb.shape[:2]) / (len(lut) - 1)).astype(np.float32), distance.reshape(rgb.shape[:2])


def _axes_from_shape(shape, config, qr=None, qz=None):
    height, width = shape
    qr = np.asarray(qr).squeeze() if qr is not None else np.linspace(*config.qr_range, width)
    qz = np.asarray(qz).squeeze() if qz is not None else np.linspace(config.qz_range[1], config.qz_range[0], height)
    if qr.ndim != 1 or len(qr) != width:
        qr = np.linspace(*config.qr_range, width)
    if qz.ndim != 1 or len(qz) != height:
        qz = np.linspace(config.qz_range[1], config.qz_range[0], height)
    return qr.astype(float), qz.astype(float)



def _first_positive_metadata_value(metadata: dict, keys: Sequence[str]):
    """Return the first finite positive scalar stored under any metadata alias."""
    for key in keys:
        if key not in metadata:
            continue
        try:
            value = float(np.asarray(metadata[key]).reshape(-1)[0])
        except Exception:
            continue
        if np.isfinite(value) and value > 0:
            return value, key
    return None, None


def _automatic_experimental_intensity_normalization(metadata: dict) -> tuple[float, str]:
    """Return a conservative scan-normalization factor from explicit metadata.

    The routine never guesses a beam normalization. An explicitly supplied
    multiplicative normalization factor takes precedence, followed by an explicit
    already-normalized flag. Integrated monitor counts then take precedence
    because they already include exposure duration. If only incident flux is
    available, flux×exposure is used; exposure-only correction is the final
    supported fallback. Absolute units are not claimed.
    """
    explicit_factor, explicit_factor_key = _first_positive_metadata_value(
        metadata,
        ("intensity_normalization_factor", "scan_normalization_factor",
         "multiplicative_intensity_normalization"),
    )
    if explicit_factor is not None:
        return explicit_factor, f"explicit {explicit_factor_key}={explicit_factor:.7g}"

    for flag_key in (
        "intensity_already_normalized", "already_normalized",
        "normalized_by_monitor", "normalized_by_exposure",
    ):
        if flag_key not in metadata:
            continue
        value = metadata[flag_key]
        if isinstance(value, str):
            truthy = value.strip().lower() in {"1", "true", "yes", "y", "on"}
        else:
            try:
                truthy = bool(np.asarray(value).reshape(-1)[0])
            except Exception:
                truthy = False
        if truthy:
            return 1.0, f"no additional scan normalization ({flag_key}=true)"

    monitor, monitor_key = _first_positive_metadata_value(
        metadata,
        (
            "incident_monitor_counts", "beam_monitor_counts", "monitor_counts",
            "i0_counts", "I0_counts", "incident_counts",
        ),
    )
    exposure, exposure_key = _first_positive_metadata_value(
        metadata,
        (
            "exposure_time_s", "exposure_s", "counting_time_s", "dwell_time_s",
            "acquisition_time_s", "integration_time_s",
        ),
    )
    flux, flux_key = _first_positive_metadata_value(
        metadata,
        ("incident_flux", "incident_flux_per_s", "beam_flux", "photons_per_s", "i0_rate"),
    )
    if monitor is not None:
        return 1.0 / monitor, f"normalized by {monitor_key}={monitor:.7g}"
    if flux is not None and exposure is not None:
        denominator = flux * exposure
        if np.isfinite(denominator) and denominator > 0:
            return 1.0 / denominator, (
                f"normalized by {flux_key}×{exposure_key}="
                f"{flux:.7g}×{exposure:.7g}"
            )
    if exposure is not None:
        return 1.0 / exposure, f"normalized by {exposure_key}={exposure:.7g} s"
    return 1.0, "no explicit exposure/incident-beam normalization metadata"


def _robust_local_peak_intensity(image: dict, row: int, column: int,
                                 base_peak_radius: int = 4,
                                 background_radius: int = 9) -> dict:
    """Measure one experimental peak with adaptive aperture/background handling.

    This routine is deliberately independent of the GUI so it can be regression
    tested.  It uses a robust local planar background, estimates an anisotropic
    peak footprint from positive residual moments, excludes pixels assigned more
    naturally to nearby maxima (a conservative local deblend), and propagates a
    background/shot-noise uncertainty estimate.  It never changes the image used
    by the indexing engine.
    """
    raw = np.asarray(image.get("raw_intensity", image.get("intensity")), dtype=float)
    quantitative = np.asarray(image.get("quantitative_intensity", raw), dtype=float)
    valid = np.asarray(image.get("valid", np.isfinite(quantitative)), dtype=bool)
    qr = np.asarray(image.get("qr", np.arange(quantitative.shape[1])), dtype=float)
    qz = np.asarray(image.get("qz", np.arange(quantitative.shape[0])), dtype=float)
    h, w = quantitative.shape
    row = int(np.clip(int(row), 0, h - 1))
    column = int(np.clip(int(column), 0, w - 1))
    base_peak_radius = max(int(base_peak_radius), 1)
    background_radius = max(int(background_radius), base_peak_radius + 3)
    # Give the adaptive aperture room to grow while keeping the fit local.
    patch_radius = max(background_radius, int(math.ceil(base_peak_radius * 3.0)), 8)
    r0, r1 = max(0, row - patch_radius), min(h, row + patch_radius + 1)
    c0, c1 = max(0, column - patch_radius), min(w, column + patch_radius + 1)
    patch = quantitative[r0:r1, c0:c1]
    raw_patch = raw[r0:r1, c0:c1]
    patch_valid = np.isfinite(patch) & valid[r0:r1, c0:c1]
    yy, xx = np.indices(patch.shape, dtype=float)
    gy = yy + r0
    gx = xx + c0
    fit_row, fit_column = float(row), float(column)
    dy = gy - fit_row
    dx = gx - fit_column
    distance = np.hypot(dx, dy)

    annulus_inner = max(float(base_peak_radius + 2), 0.55 * background_radius)
    background_mask = patch_valid & (distance >= annulus_inner) & (distance <= float(background_radius))
    if int(background_mask.sum()) < 12:
        background_mask = patch_valid & (distance > float(base_peak_radius + 1))
    if int(background_mask.sum()) < 6:
        background_mask = patch_valid.copy()

    # Robust local plane: z = a + b*dx + c*dy. Iterative clipping rejects
    # neighboring peaks/ridges from the background estimate.
    bg_y, bg_x = np.where(background_mask)
    bg_values = patch[background_mask]
    design = np.column_stack((np.ones(len(bg_values)), dx[background_mask], dy[background_mask]))
    keep = np.isfinite(bg_values)
    beta = np.array([float(np.nanmedian(bg_values)) if len(bg_values) else 0.0, 0.0, 0.0])
    for _ in range(5):
        if int(keep.sum()) < 6:
            break
        try:
            beta = np.linalg.lstsq(design[keep], bg_values[keep], rcond=None)[0]
        except Exception:
            break
        residual_bg = bg_values - design @ beta
        center = float(np.nanmedian(residual_bg[keep]))
        mad = float(np.nanmedian(np.abs(residual_bg[keep] - center)))
        sigma = max(1.4826 * mad, 1e-12)
        new_keep = np.isfinite(residual_bg) & (np.abs(residual_bg - center) <= 3.5 * sigma)
        if np.array_equal(new_keep, keep):
            break
        keep = new_keep
    plane = beta[0] + beta[1] * dx + beta[2] * dy
    residual = patch - plane
    residual_bg = (bg_values - design @ beta) if len(bg_values) else np.asarray([], float)
    residual_bg_used = residual_bg[keep] if len(residual_bg) and int(keep.sum()) else residual_bg
    if len(residual_bg_used):
        bg_center = float(np.nanmedian(residual_bg_used))
        bg_mad = float(np.nanmedian(np.abs(residual_bg_used - bg_center)))
        background_noise = max(1.4826 * bg_mad, float(np.nanstd(residual_bg_used)), 1e-12)
    else:
        background_noise = max(float(np.nanstd(residual[patch_valid])), 1e-12) if patch_valid.any() else 1e-12

    target_y = row - r0
    target_x = column - c0
    local_search = patch_valid & (distance <= max(float(base_peak_radius * 2.5), 4.0))
    if local_search.any():
        search_values = np.where(local_search, residual, -np.inf)
        local_flat = int(np.nanargmax(search_values))
        target_y, target_x = np.unravel_index(local_flat, search_values.shape)
        row = r0 + int(target_y)
        column = c0 + int(target_x)
        dy = gy - float(row)
        dx = gx - float(column)
        distance = np.hypot(dx, dy)
    peak_height = float(residual[int(target_y), int(target_x)]) if patch_valid[int(target_y), int(target_x)] else np.nan
    if not np.isfinite(peak_height):
        peak_height = float(np.nanmax(residual[patch_valid])) if patch_valid.any() else 0.0
    peak_height = max(peak_height, 0.0)

    # Detect nearby local maxima that could contaminate an aperture.  The target
    # is not counted as a neighbor.  We use a modest threshold so weak shoulders
    # are recognized without allowing background noise to fragment the peak.
    max_radius = max(float(base_peak_radius * 3.0), 5.0)
    maxima = (
        patch_valid
        & (distance <= max_radius)
        & (residual == maximum_filter(np.where(patch_valid, residual, -np.inf), size=3, mode="nearest"))
        & (residual > max(3.0 * background_noise, 0.18 * max(peak_height, background_noise)))
    )
    peak_coords = []
    maxima_coords = list(zip(*np.where(maxima)))
    maxima_coords.sort(key=lambda pair: float(residual[pair[0], pair[1]]), reverse=True)
    for py, px in maxima_coords:
        if math.hypot(float(py - target_y), float(px - target_x)) <= 1.5:
            continue
        if residual[py, px] < max(3.0 * background_noise, 0.22 * max(peak_height, background_noise)):
            continue
        if all(math.hypot(float(py - qy), float(px - qx)) >= 2.0 for qy, qx in peak_coords):
            peak_coords.append((int(py), int(px)))
        if len(peak_coords) >= 4:
            break

    # Estimate the footprint from positive residuals, then integrate the full
    # background-subtracted signal in an ellipse.  The final sum includes
    # negative fluctuations, avoiding the positive bias of summing only >0 data.
    moment_radius = min(max_radius, float(background_radius - 1)) if background_radius > 3 else max_radius
    moment_mask = patch_valid & (distance <= moment_radius)
    positive = np.where(moment_mask, np.maximum(residual, 0.0), 0.0)
    total_positive = float(np.sum(positive))
    if total_positive <= 0:
        positive = np.where(patch_valid & (distance <= float(base_peak_radius)), np.maximum(residual, 0.0), 0.0)
        total_positive = float(np.sum(positive))
    if total_positive > 0:
        cx = float(np.sum(gx * positive) / total_positive)
        cy = float(np.sum(gy * positive) / total_positive)
        ddx = gx - cx
        ddy = gy - cy
        cov_xx = float(np.sum(positive * ddx * ddx) / total_positive)
        cov_yy = float(np.sum(positive * ddy * ddy) / total_positive)
        cov_xy = float(np.sum(positive * ddx * ddy) / total_positive)
    else:
        cx, cy = float(column), float(row)
        cov_xx = cov_yy = max((base_peak_radius / 2.0) ** 2, 0.8 ** 2)
        cov_xy = 0.0
    covariance = np.array([[cov_xx, cov_xy], [cov_xy, cov_yy]], dtype=float)
    try:
        eigval, eigvec = np.linalg.eigh(covariance)
        eigval = np.clip(eigval, 0.8 ** 2, max((max_radius / 2.0) ** 2, 1.0))
        covariance = eigvec @ np.diag(eigval) @ eigvec.T
        inv_cov = np.linalg.inv(covariance)
    except Exception:
        covariance = np.diag([max((base_peak_radius / 2.0) ** 2, 0.8 ** 2)] * 2)
        inv_cov = np.linalg.inv(covariance)
    ddx = gx - cx
    ddy = gy - cy
    mahal = (
        inv_cov[0, 0] * ddx * ddx
        + 2.0 * inv_cov[0, 1] * ddx * ddy
        + inv_cov[1, 1] * ddy * ddy
    )
    geometric_aperture = (mahal <= 9.0) & (np.hypot(ddx, ddy) <= max_radius * 1.35)
    aperture = geometric_aperture & patch_valid

    # Conservative deblending: pixels geometrically closer to another resolved
    # maximum are excluded from this peak.  This prevents obvious double counting
    # without forcing a fragile Gaussian-mixture model on arcs/streaks.
    if peak_coords:
        target_dist2 = (yy - float(target_y)) ** 2 + (xx - float(target_x)) ** 2
        for py, px in peak_coords:
            neighbor_dist2 = (yy - float(py)) ** 2 + (xx - float(px)) ** 2
            aperture &= target_dist2 <= neighbor_dist2
    if int(aperture.sum()) < 3:
        aperture = patch_valid & (distance <= float(base_peak_radius))

    integrated = float(np.nansum(residual[aperture]))
    if not np.isfinite(integrated) or integrated <= 0:
        integrated = float(np.nansum(np.maximum(residual[aperture], 0.0)))
    centroid_weights = np.where(aperture, np.maximum(residual, 0.0), 0.0)
    centroid_sum = float(np.sum(centroid_weights))
    if centroid_sum > 0:
        sub_x = float(np.sum(gx * centroid_weights) / centroid_sum)
        sub_y = float(np.sum(gy * centroid_weights) / centroid_sum)
        effective_pixels = max(
            centroid_sum * centroid_sum / max(float(np.sum(centroid_weights * centroid_weights)), 1e-12),
            1.0,
        )
        variance_x = max(float(np.sum(centroid_weights * (gx - sub_x) ** 2) / centroid_sum), 0.0)
        variance_y = max(float(np.sum(centroid_weights * (gy - sub_y) ** 2) / centroid_sum), 0.0)
        sigma_x_centroid_px = math.sqrt(variance_x / effective_pixels)
        sigma_y_centroid_px = math.sqrt(variance_y / effective_pixels)
    else:
        sub_x, sub_y = float(column), float(row)
        sigma_x_centroid_px = sigma_y_centroid_px = np.nan

    n_aperture = max(int(aperture.sum()), 1)
    n_background = max(int(keep.sum()) if len(bg_values) else 0, 1)
    normalization_factor = float(image.get("intensity_normalization_factor", 1.0) or 1.0)
    native_numerical = not bool(image.get("png_reconstructed", False))
    explicit_sigma_map = image.get("quantitative_intensity_sigma")
    explicit_sigma_integrated = 0.0
    explicit_sigma_used = False
    if explicit_sigma_map is not None:
        try:
            sigma_patch = np.asarray(explicit_sigma_map, dtype=float)[r0:r1, c0:c1]
            sigma_values = sigma_patch[aperture]
            finite_sigma = sigma_values[np.isfinite(sigma_values) & (sigma_values >= 0)]
            if len(finite_sigma) == int(aperture.sum()) and len(finite_sigma):
                explicit_sigma_integrated = math.sqrt(float(np.sum(finite_sigma * finite_sigma)))
                explicit_sigma_used = True
        except Exception:
            explicit_sigma_used = False

    # The fitted background introduces uncertainty even when the detector/processing
    # pipeline provides a per-pixel sigma map.  Avoid double-counting aperture
    # noise: with explicit sigma, add only the uncertainty of the fitted local
    # background level; otherwise use the empirical background noise for both.
    if explicit_sigma_used:
        background_model_sigma = background_noise * n_aperture / math.sqrt(n_background)
        background_sigma_integrated = background_model_sigma
    else:
        background_sigma_integrated = background_noise * math.sqrt(
            n_aperture + (n_aperture * n_aperture) / n_background
        )

    raw_aperture = raw_patch[aperture]
    integer_like = False
    if native_numerical and not explicit_sigma_used and len(raw_aperture):
        finite_raw = raw_aperture[np.isfinite(raw_aperture)]
        if len(finite_raw) and np.nanmin(finite_raw) >= 0:
            integer_like = float(np.nanmedian(np.abs(finite_raw - np.rint(finite_raw)))) < 1e-3
    shot_sigma = 0.0
    if integer_like:
        shot_sigma = math.sqrt(max(float(np.nansum(np.maximum(raw_aperture, 0.0))), 0.0)) * abs(normalization_factor)
    sigma_intensity = math.sqrt(
        background_sigma_integrated ** 2 + explicit_sigma_integrated ** 2 + shot_sigma ** 2
    )
    if explicit_sigma_used:
        uncertainty_source = "NPZ per-pixel uncertainty + local-background fit"
    elif integer_like:
        uncertainty_source = "Poisson counting + local-background estimate"
    else:
        uncertainty_source = "local-background estimate"
    integrated_snr = integrated / max(sigma_intensity, 1e-12)

    qr_step = float(np.nanmedian(np.abs(np.diff(qr)))) if len(qr) > 1 else np.nan
    qz_step = float(np.nanmedian(np.abs(np.diff(qz)))) if len(qz) > 1 else np.nan
    qr_value = float(np.interp(sub_x, np.arange(len(qr), dtype=float), qr))
    qz_value = float(np.interp(sub_y, np.arange(len(qz), dtype=float), qz))
    peak_sigma_x_px = math.sqrt(max(float(covariance[0, 0]), 0.0))
    peak_sigma_y_px = math.sqrt(max(float(covariance[1, 1]), 0.0))

    geometric_count = max(int(geometric_aperture.sum()), 1)
    valid_fraction = float(np.sum(geometric_aperture & patch_valid) / geometric_count)
    edge_truncated = bool(
        r0 == 0 or c0 == 0 or r1 == h or c1 == w or valid_fraction < 0.88
    )
    global_values = quantitative[valid & np.isfinite(quantitative)]
    saturation_fraction = 0.0
    if len(global_values) and aperture.any():
        high = float(np.nanquantile(global_values, 0.9995))
        maximum = float(np.nanmax(global_values))
        threshold = max(high, maximum - max(abs(maximum) * 1e-8, 1e-12))
        saturation_fraction = float(np.mean(patch[aperture] >= threshold))

    snr_score = float(np.clip(np.log1p(max(integrated_snr, 0.0)) / np.log(21.0), 0.0, 1.0))
    background_score = float(np.clip(n_background / 50.0, 0.0, 1.0))
    quality_score = 0.58 * snr_score + 0.17 * background_score + 0.25 * np.clip(valid_fraction, 0.0, 1.0)
    quality_score -= min(0.28, 0.08 * len(peak_coords))
    quality_score -= min(0.25, 1.5 * saturation_fraction)
    if edge_truncated:
        quality_score -= 0.18
    if bool(image.get("png_reconstructed", False)):
        quality_score -= 0.10
    quality_score = float(np.clip(quality_score, 0.0, 1.0))
    quality = "high" if quality_score >= 0.75 else "medium" if quality_score >= 0.45 else "low"
    local_background = float(beta[0] + beta[1] * (sub_x - fit_column) + beta[2] * (sub_y - fit_row))
    relative_uncertainty = sigma_intensity / abs(integrated) if abs(integrated) > 1e-12 else np.nan

    return {
        "QrExp": qr_value,
        "QzExp": qz_value,
        "ExpIntensity": integrated,
        "PixelX": int(column),
        "PixelY": int(row),
        "SubpixelX": sub_x,
        "SubpixelY": sub_y,
        "PeakHeight": peak_height,
        "LocalBackground": local_background,
        "BackgroundNoise": background_noise,
        "BackgroundGradient": float(math.hypot(beta[1], beta[2])),
        "PeakSNR": peak_height / max(background_noise, 1e-12),
        "IntegratedSNR": integrated_snr,
        "SigmaExpIntensity": sigma_intensity,
        "RelativeIntensityUncertainty": relative_uncertainty,
        "PeakAreaPixels": int(aperture.sum()),
        "PeakSigmaQrWidth": peak_sigma_x_px * qr_step if np.isfinite(qr_step) else np.nan,
        "PeakSigmaQzWidth": peak_sigma_y_px * qz_step if np.isfinite(qz_step) else np.nan,
        "SigmaQrExp": sigma_x_centroid_px * qr_step if np.isfinite(qr_step) else np.nan,
        "SigmaQzExp": sigma_y_centroid_px * qz_step if np.isfinite(qz_step) else np.nan,
        "OverlapPeakCount": int(len(peak_coords)),
        "Deblended": bool(len(peak_coords) > 0),
        "SaturationFraction": saturation_fraction,
        "ValidPixelFraction": valid_fraction,
        "EdgeTruncated": edge_truncated,
        "IntensityQualityScore": quality_score,
        "IntensityQuality": quality,
        "IntensityNormalizationFactor": normalization_factor,
        "IntensityNormalization": str(image.get("intensity_normalization_note", "none")),
        "IntensitySource": (
            "PNG-reconstructed NPZ" if bool(image.get("png_reconstructed", False))
            else "native numerical NPZ"
        ),
        "IntensityUncertaintySource": uncertainty_source,
        "PoissonTermUsed": bool(integer_like),
    }


def _v739_load_numerical_qspace(path: str, config: IndexingConfig) -> dict:
    """Load numerical reciprocal-space data with robust NPZ inference.

    Supported NPZ layouts include:
    - intensity + qr + qz
    - image/data/I/waxs/giwaxs/counts + optional axes
    - a single unnamed arr_0 2-D array
    - separable 2-D qr/qz grids
    - RGB/RGBA arrays (converted to a scalar display intensity as a fallback)

    If qr/qz axes are absent, the measurement q ranges supplied by the GUI are
    used only when the explicit-axis requirement is disabled.
    """
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix != ".npz":
        raise ValueError(f"Numerical q-space input must be an NPZ archive, not {path.name}.")
    intensity = qr = qz = mask = intensity_sigma = intensity_variance = None
    mapping = {}
    source_note = ""
    archive_metadata = {}

    def numeric_array(value):
        try:
            array = np.asarray(value)
        except Exception:
            return None
        return array if np.issubdtype(array.dtype, np.number) else None

    def pick_case_insensitive(names, allowed_ndim=None):
        lower = {str(key).lower(): key for key in mapping}
        for name in names:
            key = lower.get(str(name).lower())
            if key is None:
                continue
            array = numeric_array(mapping[key])
            if array is None:
                continue
            if allowed_ndim is None or array.ndim in tuple(allowed_ndim):
                return array
        return None

    def axis_from_grid(grid, wanted, shape):
        if grid is None:
            return None
        array = np.asarray(grid, dtype=float)
        if array.ndim == 1:
            return array
        if array.ndim != 2:
            return None
        if array.shape == shape:
            if wanted == "qr" and np.allclose(array, array[0:1, :], equal_nan=True):
                return array[0, :]
            if wanted == "qz" and np.allclose(array, array[:, 0:1], equal_nan=True):
                return array[:, 0]
        if array.T.shape == shape:
            transposed = array.T
            if wanted == "qr" and np.allclose(transposed, transposed[0:1, :], equal_nan=True):
                return transposed[0, :]
            if wanted == "qz" and np.allclose(transposed, transposed[:, 0:1], equal_nan=True):
                return transposed[:, 0]
        return None

    if suffix == ".npz":
        with np.load(path, allow_pickle=False) as data:
            mapping = {key: data[key] for key in data.files}
        for metadata_key in (
            "source_png", "conversion_note", "colormap", "crop_xyxy",
            "qr_limits", "qz_limits", "color_error",
            "xray_wavelength_A", "wavelength_A", "wavelength_angstrom",
            "xray_energy_keV", "energy_keV", "beam_energy_keV",
            "incidence_angle_deg", "incident_angle_deg", "alpha_i_deg",
            "film_delta", "film_beta", "substrate_delta", "substrate_beta",
            "film_thickness_A", "film_thickness_nm",
            "exposure_time_s", "exposure_s", "counting_time_s", "dwell_time_s",
            "acquisition_time_s", "integration_time_s",
            "incident_monitor_counts", "beam_monitor_counts", "monitor_counts",
            "i0_counts", "I0_counts", "incident_counts",
            "incident_flux", "incident_flux_per_s", "beam_flux", "photons_per_s", "i0_rate",
            "intensity_normalization_factor", "scan_normalization_factor",
            "multiplicative_intensity_normalization", "intensity_already_normalized",
            "already_normalized", "normalized_by_monitor", "normalized_by_exposure",
            "intensity_scale", "display_scale", "normalization", "colorbar_min", "colorbar_max",
            "display_vmin", "display_vmax", "vmin", "vmax", "log_base",
        ):
            if metadata_key not in mapping:
                continue
            try:
                value = np.asarray(mapping[metadata_key])
                archive_metadata[metadata_key] = value.item() if value.ndim == 0 else value.tolist()
            except Exception:
                archive_metadata[metadata_key] = str(mapping[metadata_key])
    # GIWAXS filenames can encode exposure time. Use that value only when the NPZ
    # does not contain explicit exposure metadata, allowing scan-to-scan intensity
    # normalization without requiring an additional manual entry.
        if not any(
            key in archive_metadata for key in (
                "exposure_time_s", "exposure_s", "counting_time_s", "dwell_time_s",
                "acquisition_time_s", "integration_time_s",
            )
        ):
            parsed_name = _parse_measurement_name(path)
            if parsed_name is not None:
                exposure_from_name = float(parsed_name.get("exposure_s", np.nan))
                if np.isfinite(exposure_from_name) and exposure_from_name > 0:
                    archive_metadata["exposure_time_s"] = exposure_from_name
                    archive_metadata["exposure_time_source"] = "measurement filename"

        intensity = pick_case_insensitive(
            ("corrected_intensity", "intensity_corrected", "corrected_counts")
            + tuple(config.intensity_keys) + ("counts", "arr_0"), (2, 3)
        )
        if intensity is not None:
            lower_keys = {str(key).lower() for key in mapping}
            if any(key in lower_keys for key in ("corrected_intensity", "intensity_corrected", "corrected_counts")):
                source_note = "using explicitly corrected numerical intensity array"
        qr_candidate = pick_case_insensitive(tuple(config.qr_keys) + ("q_r_grid", "qr_grid"), (1, 2))
        qz_candidate = pick_case_insensitive(tuple(config.qz_keys) + ("q_z_grid", "qz_grid"), (1, 2))
        mask = pick_case_insensitive(("mask", "valid", "valid_mask"), (2,))
        intensity_sigma = pick_case_insensitive(
            ("intensity_sigma", "sigma_intensity", "intensity_uncertainty",
             "intensity_error", "intensity_errors", "count_sigma"), (2,)
        )
        intensity_variance = pick_case_insensitive(
            ("intensity_variance", "variance_intensity", "count_variance"), (2,)
        )

        if intensity is None:
            candidates = []
            excluded_fallback_keys = {
                "mask", "valid", "valid_mask", "qr", "q_r", "qz", "q_z",
                "qr_grid", "q_r_grid", "qz_grid", "q_z_grid",
                "intensity_sigma", "sigma_intensity", "intensity_uncertainty",
                "intensity_error", "intensity_errors", "count_sigma",
                "intensity_variance", "variance_intensity", "count_variance",
            }
            for key, value in mapping.items():
                if str(key).lower() in excluded_fallback_keys:
                    continue
                array = numeric_array(value)
                if array is not None and array.ndim in (2, 3):
                    candidates.append((key, array))
            if candidates:
                key, intensity = max(candidates, key=lambda item: item[1].size)
                source_note = f"inferred intensity from key {key!r}"

        if intensity is not None and np.asarray(intensity).ndim == 3:
            rgb = np.asarray(intensity, dtype=float)
            if rgb.shape[-1] not in (3, 4):
                raise ValueError(
                    f"The inferred NPZ intensity has unsupported shape {rgb.shape}. "
                    "Expected a 2-D array or RGB/RGBA image."
                )
            rgb = rgb[..., :3]
            finite = np.isfinite(rgb)
            if finite.any() and float(np.nanmax(rgb[finite])) > 1.5:
                rgb = rgb / 255.0
            rgb = np.clip(rgb, 0.0, 1.0)
            try:
                intensity, _ = _invert_colormap(rgb.astype(np.float32), config.colormap)
                source_note = (source_note + "; " if source_note else "") + "RGB colormap inverted"
            except Exception:
                intensity = np.nanmean(rgb, axis=2)
                source_note = (source_note + "; " if source_note else "") + "RGB converted to luminance"

        if intensity is not None:
            shape = np.asarray(intensity).shape
            qr = axis_from_grid(qr_candidate, "qr", shape)
            qz = axis_from_grid(qz_candidate, "qz", shape)
            if qr is None and qr_candidate is not None and np.asarray(qr_candidate).ndim == 1:
                qr = np.asarray(qr_candidate, dtype=float)
            if qz is None and qz_candidate is not None and np.asarray(qz_candidate).ndim == 1:
                qz = np.asarray(qz_candidate, dtype=float)


    if intensity is None:
        available = ", ".join(mapping.keys()) if mapping else "none"
        raise ValueError(
            f"Could not identify a numerical intensity array in {path}. "
            f"Available NPZ keys: {available}."
        )

    intensity = np.asarray(intensity, dtype=float)
    if intensity_sigma is not None:
        intensity_sigma = np.asarray(intensity_sigma, dtype=float)
    elif intensity_variance is not None:
        variance_array = np.asarray(intensity_variance, dtype=float)
        intensity_sigma = np.sqrt(np.clip(variance_array, 0.0, None))

        # PNG-derived NPZ intensity can recover quantitative display units only when
        # the converter stored explicit colorbar normalization. Reconstruct that scale
        # when metadata are available; never infer an absolute intensity scale from
        # image colors alone.
    reconstructed_flag = bool(
        archive_metadata.get("source_png")
        or "rendered png" in str(archive_metadata.get("conversion_note", "")).lower()
        or "reconstructed" in str(archive_metadata.get("conversion_note", "")).lower()
    )
    if reconstructed_flag and np.isfinite(intensity).any():
        finite_intensity = intensity[np.isfinite(intensity)]
        normalized_like = (
            float(np.nanmin(finite_intensity)) >= -1e-6
            and float(np.nanmax(finite_intensity)) <= 1.000001
        )
        def _metadata_number(*names):
            for name in names:
                if name not in archive_metadata:
                    continue
                try:
                    value = float(np.asarray(archive_metadata[name]).reshape(-1)[0])
                except Exception:
                    continue
                if np.isfinite(value):
                    return value
            return None
        display_min = _metadata_number("colorbar_min", "display_vmin", "vmin")
        display_max = _metadata_number("colorbar_max", "display_vmax", "vmax")
        scale_name = str(archive_metadata.get("intensity_scale", archive_metadata.get("display_scale", archive_metadata.get("normalization", "")))).lower()
        if normalized_like and display_min is not None and display_max is not None and display_max > display_min:
            if "log" in scale_name and display_min > 0 and display_max > 0:
                intensity = np.exp(
                    np.log(display_min) + intensity * (np.log(display_max) - np.log(display_min))
                )
                source_note = (source_note + "; " if source_note else "") + "restored explicit logarithmic display intensity scale"
            elif scale_name in ("linear", "lin", "normalize", "normalized", ""):
                intensity = display_min + intensity * (display_max - display_min)
                source_note = (source_note + "; " if source_note else "") + "restored explicit linear display intensity scale"

    # Per-pixel uncertainty is trusted only for native numerical data.  A rendered
    # PNG may have undergone nonlinear display scaling, clipping, interpolation,
    # and color quantization, so an uncertainty map in pre-render units cannot be
    # propagated safely without a complete rendering transform.
    if reconstructed_flag:
        intensity_sigma = None

    if intensity.ndim != 2:
        raise ValueError(f"Numerical q-space intensity must be 2-D after conversion; got {intensity.shape}.")
    if intensity.shape[0] < 4 or intensity.shape[1] < 4:
        raise ValueError(f"Numerical q-space array is too small: {intensity.shape}")

    # Resolve an unambiguous transposed numerical archive *before* axis fallback.
    # Otherwise _axes_from_shape would replace the mismatched explicit axes with
    # GUI-generated axes and silently lose the archive calibration.
    if qr is not None and qz is not None:
        qr_array = np.asarray(qr).squeeze()
        qz_array = np.asarray(qz).squeeze()
        if qr_array.ndim == 1 and qz_array.ndim == 1:
            direct_shape = (len(qz_array), len(qr_array))
            transposed_shape = (len(qr_array), len(qz_array))
            if intensity.shape == transposed_shape and intensity.shape != direct_shape:
                intensity = intensity.T
                if mask is not None:
                    mask = np.asarray(mask).T
                if intensity_sigma is not None:
                    intensity_sigma = np.asarray(intensity_sigma).T

        # Record whether supplied q axes correspond to the actual intensity-array
        # orientation before any GUI-range fallback axes are generated.
    explicit_qr_axis = (
        qr is not None and np.asarray(qr).squeeze().ndim == 1
        and int(np.asarray(qr).size) == int(intensity.shape[1])
    )
    explicit_qz_axis = (
        qz is not None and np.asarray(qz).squeeze().ndim == 1
        and int(np.asarray(qz).size) == int(intensity.shape[0])
    )

    # Axes may be absent. _axes_from_shape then uses the per-row GUI q ranges.
    qr, qz = _axes_from_shape(intensity.shape, config, qr, qz)
    qr = np.asarray(qr, dtype=float).reshape(-1)
    qz = np.asarray(qz, dtype=float).reshape(-1)

    # Accept transposed arrays when the stored axis lengths reveal the layout.
    if intensity.shape == (len(qr), len(qz)) and intensity.shape != (len(qz), len(qr)):
        intensity = intensity.T
        if mask is not None:
            mask = np.asarray(mask).T
        if intensity_sigma is not None:
            intensity_sigma = np.asarray(intensity_sigma).T
    if intensity.shape != (len(qz), len(qr)):
        raise ValueError(
            f"Intensity shape {intensity.shape} does not match qz×qr {(len(qz), len(qr))}."
        )

    if mask is not None:
        mask = np.asarray(mask, dtype=bool)
        if mask.shape == (len(qr), len(qz)) and mask.shape != intensity.shape:
            mask = mask.T
        if mask.shape != intensity.shape:
            raise ValueError(f"Mask shape {mask.shape} does not match intensity {intensity.shape}.")
    if intensity_sigma is not None:
        intensity_sigma = np.asarray(intensity_sigma, dtype=float)
        if intensity_sigma.shape == (len(qr), len(qz)) and intensity_sigma.shape != intensity.shape:
            intensity_sigma = intensity_sigma.T
        if intensity_sigma.shape != intensity.shape:
            # An uncertainty array is optional; a mismatched one must never make a
            # valid intensity archive unusable. Ignore it and estimate uncertainty
            # locally instead.
            intensity_sigma = None

    if qr[0] > qr[-1]:
        qr = qr[::-1]
        intensity = intensity[:, ::-1]
        if mask is not None:
            mask = mask[:, ::-1]
        if intensity_sigma is not None:
            intensity_sigma = intensity_sigma[:, ::-1]
    if qz[0] < qz[-1]:
        qz = qz[::-1]
        intensity = intensity[::-1, :]
        if mask is not None:
            mask = mask[::-1, :]
        if intensity_sigma is not None:
            intensity_sigma = intensity_sigma[::-1, :]

    finite = np.isfinite(intensity)
    valid = finite if mask is None else finite & mask
    values = intensity[valid]
    if not len(values):
        raise ValueError(f"No finite valid intensity pixels were found in {path}.")
    low, high = np.quantile(values, [0.01, 0.995])
    normalized = np.clip((intensity - low) / max(float(high - low), 1e-12), 0, 1).astype(np.float32)
    intensity_normalization_factor, intensity_normalization_note = (
        _automatic_experimental_intensity_normalization(archive_metadata)
    )
    quantitative_intensity = intensity * float(intensity_normalization_factor)
    quantitative_intensity_sigma = (
        np.abs(float(intensity_normalization_factor)) * intensity_sigma
        if intensity_sigma is not None else None
    )
    qr_grid, qz_grid = np.meshgrid(qr, qz)
    valid &= qz_grid >= config.analysis_qz_min
    valid &= ~((np.abs(qr_grid) < config.exclude_specular_abs_qr)
               & (qz_grid < config.exclude_specular_qz_max))
    return {
        "path": str(path), "source_kind": "numerical", "crop": None,
        "rgb": None, "intensity": normalized, "display_intensity": normalized.copy(),
        "raw_intensity": intensity, "quantitative_intensity": quantitative_intensity,
        "quantitative_intensity_sigma": quantitative_intensity_sigma,
        "intensity_normalization_factor": float(intensity_normalization_factor),
        "intensity_normalization_note": intensity_normalization_note,
        "valid": valid, "qr": qr, "qz": qz,
        "qr_grid": qr_grid, "qz_grid": qz_grid,
        "axis_calibration_source": (
            "explicit_npz_q_axes"
            if explicit_qr_axis and explicit_qz_axis
            else "gui_q_range_fallback"
        ),
        "explicit_q_axes": bool(explicit_qr_axis and explicit_qz_axis),
        "numerical_loader_note": source_note or "recognized numerical intensity",
        "archive_metadata": archive_metadata,
        "png_reconstructed": reconstructed_flag,
        "input_provenance": (
            "PNG-reconstructed NPZ: reciprocal-space axes are explicit, but experimental intensities "
            "were estimated from rendered colors."
            if reconstructed_flag
            else "Numerical NPZ with explicit reciprocal-space axes."
            if explicit_qr_axis and explicit_qz_axis
            else "Numerical array using GUI-generated reciprocal-space axes."
        ),
    }


def _row_field(row, name, default=None):
    if isinstance(row, pd.Series):
        value = row.get(name, default)
    else:
        value = getattr(row, name, default)
    return default if value is None or (isinstance(value, float) and np.isnan(value)) else value


def _measurement_specific_config(row, config: IndexingConfig) -> IndexingConfig:
    qr_min, qr_max = _row_field(row, "qr_min"), _row_field(row, "qr_max")
    qz_min, qz_max = _row_field(row, "qz_min"), _row_field(row, "qz_max")
    qr_range = config.qr_range if qr_min is None or qr_max is None else (float(qr_min), float(qr_max))
    qz_range = config.qz_range if qz_min is None or qz_max is None else (float(qz_min), float(qz_max))
    row_colormap = _row_field(row, "row_colormap", config.colormap)
    crop_values = tuple(_row_field(row, key) for key in ("crop_x0", "crop_y0", "crop_x1", "crop_y1"))
    row_crop = config.crop_xyxy
    if all(value is not None for value in crop_values):
        row_crop = tuple(int(value) for value in crop_values)
    row_angle = _row_field(row, "angle_deg", None)
    try:
        row_angle = float(row_angle)
    except Exception:
        row_angle = float(config.incidence_angle_deg)
    if not np.isfinite(row_angle) or row_angle <= 0:
        row_angle = float(config.incidence_angle_deg)
    return replace(
        config, qr_range=qr_range, qz_range=qz_range,
        colormap=str(row_colormap or config.colormap), crop_xyxy=row_crop,
        incidence_angle_deg=row_angle, giwaxs_incidence_angles_deg=(row_angle,),
    )


def _covariance_feature(signal, region, image, threshold, feature_type, source, pixel_x, pixel_y):
    weights = np.clip(signal[region] - threshold, 0.0, None)
    yy, xx = np.where(region)
    total = float(weights.sum())
    if total <= 0:
        return None
    xbar, ybar = float(np.sum(weights * xx) / total), float(np.sum(weights * yy) / total)
    qr_values = np.interp(xx, np.arange(len(image["qr"])), image["qr"])
    qz_values = np.interp(yy, np.arange(len(image["qz"])), image["qz"])
    qr = float(np.sum(weights * qr_values) / total)
    qz = float(np.sum(weights * qz_values) / total)
    dq = np.column_stack((qr_values - qr, qz_values - qz))
    cov = (dq * weights[:, None]).T @ dq / total
    dqr = abs(float(image["qr"][1] - image["qr"][0]))
    dqz = abs(float(image["qz"][1] - image["qz"][0]))
    cov += np.diag([(0.5 * dqr) ** 2, (0.5 * dqz) ** 2])
    eigenvalues, eigenvectors = np.linalg.eigh(cov)
    order = np.argsort(eigenvalues)[::-1]
    major, minor = np.sqrt(np.maximum(eigenvalues[order], 1e-12))
    direction = eigenvectors[:, order[0]]
    orientation = math.degrees(math.atan2(direction[1], direction[0]))
    radial = np.array([qr, qz]);
    radial /= max(np.linalg.norm(radial), 1e-12)
    radial_alignment = abs(float(np.dot(direction, radial)))
    aspect = float(major / max(minor, 1e-12))
    if feature_type == "auto":
        if aspect < 2.0:
            feature_type = "spot"
        elif radial_alignment < 0.55:
            feature_type = "arc"
        elif radial_alignment > 0.82:
            feature_type = "radial_streak"
        else:
            feature_type = "elongated"
    return {
        "qr": qr, "qz": qz, "cov_rr": float(cov[0, 0]), "cov_rz": float(cov[0, 1]),
        "cov_zz": float(cov[1, 1]), "sigma_qr": math.sqrt(max(cov[0, 0], 1e-12)),
        "sigma_qz": math.sqrt(max(cov[1, 1], 1e-12)), "major_width_q": float(major),
        "minor_width_q": float(minor), "aspect_ratio": aspect, "orientation_deg": orientation,
        "feature_type": feature_type, "integrated_signal": total, "source_detector": source,
        "pixel_x": float(pixel_x), "pixel_y": float(pixel_y),
    }


def _v739_detect_features(image: dict, config: IndexingConfig) -> pd.DataFrame:
    intensity, valid = image["intensity"], image["valid"]
    small = gaussian_filter(intensity, config.feature_sigma_px)
    background = gaussian_filter(intensity, config.background_sigma_px)
    signal = small - background
    values = signal[valid]
    columns = [
        "qr", "qz", "cov_rr", "cov_rz", "cov_zz", "sigma_qr", "sigma_qz",
        "major_width_q", "minor_width_q", "aspect_ratio", "orientation_deg",
        "feature_type", "strength", "snr", "integrated_signal", "source_detector",
        "pixel_x", "pixel_y", "experimental_integrated_intensity",
        "experimental_intensity_sigma", "experimental_integrated_snr",
        "experimental_intensity_quality_score", "experimental_intensity_quality",
        "experimental_intensity_source",
    ]
    if not len(values):
        return pd.DataFrame(columns=columns)
    median = float(np.median(values))
    mad = max(float(1.4826 * np.median(np.abs(values - median))), 1e-12)
    peak_threshold = max(median + config.feature_threshold_mad * mad,
                         float(np.quantile(values, config.feature_quantile)))
    ridge_threshold = median + config.ridge_threshold_mad * mad
    candidates = []

    # Local maxima retain distinct spot-like features.
    size = 2 * config.min_feature_spacing_px + 1
    peaks = valid & (signal == maximum_filter(signal, size=size)) & (signal > peak_threshold)
    for y0, x0 in zip(*np.where(peaks)):
        radius = config.subpixel_radius_px
        y1, y2 = max(0, y0 - radius), min(signal.shape[0], y0 + radius + 1)
        x1, x2 = max(0, x0 - radius), min(signal.shape[1], x0 + radius + 1)
        region = np.zeros_like(valid, dtype=bool)
        region[y1:y2, x1:x2] = valid[y1:y2, x1:x2]
        feature = _covariance_feature(signal, region, image, peak_threshold, "spot", "local_maximum", x0, y0)
        if feature:
            feature["snr"] = float((signal[y0, x0] - median) / mad)
            candidates.append(feature)

    # Connected residual components preserve extended arcs or ridges and retain
    # anisotropic shape information for uncertainty-aware feature characterization.
    ridge_mask = valid & (signal > ridge_threshold)
    ridge_mask = binary_opening(ridge_mask, iterations=1)
    ridge_mask = binary_closing(ridge_mask, iterations=1)
    components, count = label(ridge_mask)
    for component in range(1, count + 1):
        region = components == component
        pixels = int(region.sum())
        if pixels < config.component_min_pixels or pixels > config.component_max_pixels:
            continue
        yy, xx = np.where(region)
        feature = _covariance_feature(
            signal, region, image, ridge_threshold, "auto", "connected_component",
            float(xx.mean()), float(yy.mean()),
        )
        if feature and (feature["aspect_ratio"] >= config.arc_aspect_ratio or pixels >= 12):
            feature["snr"] = float((signal[region].max() - median) / mad)
            candidates.append(feature)

    if not candidates:
        return pd.DataFrame(columns=columns)
    frame = pd.DataFrame(candidates)
    maximum = max(float(frame["integrated_signal"].max()), 1e-12)
    frame["strength"] = np.clip(frame["integrated_signal"] / maximum, 0, 1)
    frame["rank"] = frame["strength"] * np.log1p(np.maximum(frame["snr"], 0))
    # Greedy de-duplication favors well-localized high-rank candidates.
    kept = []
    for row in frame.sort_values("rank", ascending=False).itertuples(index=False):
        if all(math.hypot(row.qr - prev.qr, row.qz - prev.qz) > config.feature_merge_tolerance_q for prev in kept):
            kept.append(row)
        if len(kept) >= config.max_features_per_image:
            break
    result = pd.DataFrame([row._asdict() for row in kept]).drop(columns="rank", errors="ignore")
    # Re-measure retained features on the quantitative/native intensity array.
    # These values are output intensity diagnostics; positional indexing weights
    # remain based on the feature-detection measurements.
    measured_rows = []
    for row in result.itertuples(index=False):
        try:
            measurement = _robust_local_peak_intensity(
                image, int(round(float(row.pixel_y))), int(round(float(row.pixel_x))),
                base_peak_radius=max(int(config.subpixel_radius_px), 2),
                background_radius=max(int(round(config.background_sigma_px)), int(config.subpixel_radius_px) + 4),
            )
        except Exception:
            measurement = {}
        measured_rows.append(measurement)
    result["experimental_integrated_intensity"] = [m.get("ExpIntensity", np.nan) for m in measured_rows]
    result["experimental_intensity_sigma"] = [m.get("SigmaExpIntensity", np.nan) for m in measured_rows]
    result["experimental_integrated_snr"] = [m.get("IntegratedSNR", np.nan) for m in measured_rows]
    result["experimental_intensity_quality_score"] = [m.get("IntensityQualityScore", np.nan) for m in measured_rows]
    result["experimental_intensity_quality"] = [m.get("IntensityQuality", "unassessed") for m in measured_rows]
    result["experimental_intensity_source"] = [m.get("IntensitySource", "unreported") for m in measured_rows]
    return result[columns].reset_index(drop=True)


def _registration_feature_subset(features: pd.DataFrame, config: IndexingConfig) -> pd.DataFrame:
    if features is None or features.empty:
        return pd.DataFrame(columns=["qr", "qz", "strength", "snr", "feature_type"])
    frame = features.copy()
    frame = frame[frame["feature_type"].astype(str) != "radial_streak"]
    frame["_registration_rank"] = (
            np.maximum(frame["strength"].to_numpy(float), 1e-9)
            * np.log1p(np.maximum(frame["snr"].to_numpy(float), 0.0))
            / np.sqrt(np.maximum(frame["cov_rr"].to_numpy(float) + frame["cov_zz"].to_numpy(float), 1e-8))
    )
    return frame.nlargest(int(config.registration_max_features), "_registration_rank").reset_index(drop=True)


def _registration_unique_pairs(source_xy, reference_xy, matrix, offset, tolerance):
    if len(source_xy) == 0 or len(reference_xy) == 0:
        return np.empty((0, 2), int), np.empty(0, float)
    transformed = source_xy @ matrix.T + offset
    distance, index = cKDTree(reference_xy).query(transformed, k=1)
    order = np.argsort(distance, kind="mergesort")
    used_ref, pairs, distances = set(), [], []
    for src_index in order:
        ref_index = int(index[src_index])
        if float(distance[src_index]) > tolerance or ref_index in used_ref:
            continue
        used_ref.add(ref_index)
        pairs.append((int(src_index), ref_index))
        distances.append(float(distance[src_index]))
    return np.asarray(pairs, int).reshape(-1, 2), np.asarray(distances, float)


def _fit_similarity_transform(source_xy, reference_xy, weights, config):
    weights = np.maximum(np.asarray(weights, float), 1e-12)
    weights /= weights.sum()
    source_mean = np.sum(source_xy * weights[:, None], axis=0)
    reference_mean = np.sum(reference_xy * weights[:, None], axis=0)
    x = source_xy - source_mean
    y = reference_xy - reference_mean
    covariance = (x * weights[:, None]).T @ y
    u, singular, vt = np.linalg.svd(covariance)
    rotation = vt.T @ u.T
    if np.linalg.det(rotation) < 0:
        vt[-1] *= -1
        rotation = vt.T @ u.T
    angle = math.degrees(math.atan2(rotation[1, 0], rotation[0, 0]))
    angle = float(np.clip(angle, -config.registration_max_rotation_deg, config.registration_max_rotation_deg))
    radians = math.radians(angle)
    rotation = np.array([[math.cos(radians), -math.sin(radians)],
                         [math.sin(radians), math.cos(radians)]], float)
    rotated = x @ rotation.T
    denominator = float(np.sum(weights * np.sum(rotated * rotated, axis=1)))
    scale = float(np.sum(weights * np.sum(rotated * y, axis=1)) / max(denominator, 1e-12))
    scale = float(np.clip(scale, 1.0 - config.registration_max_scale_change,
                          1.0 + config.registration_max_scale_change))
    matrix = scale * rotation
    offset = reference_mean - source_mean @ matrix.T
    norm = float(np.linalg.norm(offset))
    if norm > config.registration_max_shift_q:
        offset *= config.registration_max_shift_q / max(norm, 1e-12)
    return matrix, offset


def _estimate_feature_registration(source_features, reference_features, config):
    source = _registration_feature_subset(source_features, config)
    reference = _registration_feature_subset(reference_features, config)
    identity = np.eye(2, dtype=float)
    zero = np.zeros(2, dtype=float)
    result = {"matrix": identity, "offset": zero, "matches": 0,
              "median_before": np.nan, "median_after": np.nan,
              "accepted": False, "reason": "insufficient_features"}
    if len(source) < config.registration_min_matches or len(reference) < config.registration_min_matches:
        return result
    source_xy = source[["qr", "qz"]].to_numpy(float)
    reference_xy = reference[["qr", "qz"]].to_numpy(float)

    # Translation seed from the densest cluster of plausible pair differences.
    differences = reference_xy[:, None, :] - source_xy[None, :, :]
    flat = differences.reshape(-1, 2)
    flat = flat[np.linalg.norm(flat, axis=1) <= float(config.registration_max_shift_q)]
    offset = zero.copy()
    if len(flat):
        bin_size = max(float(config.registration_pair_tolerance_q) / 4.0, 0.006)
        bins = np.round(flat / bin_size).astype(int)
        _, inverse, counts = np.unique(bins, axis=0, return_inverse=True, return_counts=True)
        best = int(np.argmax(counts))
        offset = np.median(flat[inverse == best], axis=0)
    matrix = identity.copy()

    initial_pairs, initial_distances = _registration_unique_pairs(
        source_xy, reference_xy, matrix, offset, float(config.registration_pair_tolerance_q)
    )
    result["median_before"] = float(np.median(initial_distances)) if len(initial_distances) else np.nan
    for _ in range(max(1, int(config.registration_iterations))):
        pairs, distances = _registration_unique_pairs(
            source_xy, reference_xy, matrix, offset, float(config.registration_pair_tolerance_q)
        )
        if len(pairs) < int(config.registration_min_matches):
            break
        robust_limit = np.median(distances) + 2.5 * max(
            1.4826 * np.median(np.abs(distances - np.median(distances))), 0.003
        )
        pairs = pairs[distances <= robust_limit]
        if len(pairs) < int(config.registration_min_matches):
            break
        src = source_xy[pairs[:, 0]]
        ref = reference_xy[pairs[:, 1]]
        weights = np.sqrt(
            np.maximum(source.iloc[pairs[:, 0]]["strength"].to_numpy(float), 1e-6)
            * np.maximum(reference.iloc[pairs[:, 1]]["strength"].to_numpy(float), 1e-6)
        )
        matrix, offset = _fit_similarity_transform(src, ref, weights, config)

    pairs, distances = _registration_unique_pairs(
        source_xy, reference_xy, matrix, offset, float(config.registration_pair_tolerance_q)
    )
    result["matches"] = int(len(pairs))
    result["median_after"] = float(np.median(distances)) if len(distances) else np.nan
    if len(pairs) < int(config.registration_min_matches):
        result["reason"] = "too_few_matched_features"
        return result
    before = result["median_before"]
    after = result["median_after"]
    improvement = (before - after) / max(before, 1e-9) if np.isfinite(before) else 0.0
    if not np.isfinite(after) or improvement < float(config.registration_min_improvement_fraction):
        result["reason"] = "insufficient_residual_improvement"
        return result
    result.update({"matrix": matrix, "offset": offset, "accepted": True, "reason": "accepted"})
    return result


def _apply_registration_to_features(frame, matrix, offset):
    if frame.empty:
        return frame.copy()
    out = frame.copy()
    out["qr_unregistered"] = out["qr"].to_numpy(float)
    out["qz_unregistered"] = out["qz"].to_numpy(float)
    points = out[["qr", "qz"]].to_numpy(float) @ matrix.T + offset
    out["qr"], out["qz"] = points[:, 0], points[:, 1]
    transformed_covariances = []
    for row in out.itertuples(index=False):
        covariance = np.array([[row.cov_rr, row.cov_rz], [row.cov_rz, row.cov_zz]], float)
        transformed_covariances.append(matrix @ covariance @ matrix.T)
    covariance = np.asarray(transformed_covariances)
    out["cov_rr"] = covariance[:, 0, 0]
    out["cov_rz"] = covariance[:, 0, 1]
    out["cov_zz"] = covariance[:, 1, 1]
    out["sigma_qr"] = np.sqrt(np.maximum(out["cov_rr"], 1e-12))
    out["sigma_qz"] = np.sqrt(np.maximum(out["cov_zz"], 1e-12))
    scale = math.sqrt(abs(float(np.linalg.det(matrix))))
    out["major_width_q"] = out["major_width_q"].to_numpy(float) * scale
    out["minor_width_q"] = out["minor_width_q"].to_numpy(float) * scale
    out["registration_scale"] = scale
    out["registration_rotation_deg"] = math.degrees(math.atan2(matrix[1, 0], matrix[0, 0]))
    out["registration_qr_offset"] = float(offset[0])
    out["registration_qz_offset"] = float(offset[1])
    return out


def _q_to_pixel(image, qr_values, qz_values):
    qr_values = np.asarray(qr_values, float)
    qz_values = np.asarray(qz_values, float)
    x = np.interp(qr_values, image["qr"], np.arange(len(image["qr"]), dtype=float), left=np.nan, right=np.nan)
    qz_axis = np.asarray(image["qz"], float)
    y_indices = np.arange(len(qz_axis), dtype=float)
    if qz_axis[0] > qz_axis[-1]:
        y = np.interp(qz_values, qz_axis[::-1], y_indices[::-1], left=np.nan, right=np.nan)
    else:
        y = np.interp(qz_values, qz_axis, y_indices, left=np.nan, right=np.nan)
    return x, y


def _warp_registered_image(image, matrix, offset, reference_image):
    """Warp scientific and visual rasters separately.

    The scientific raster is masked by ``valid`` for feature analysis. The
    visual raster is warped from ``display_intensity`` without applying that
    analysis mask, so the excluded central/specular region remains visible in
    overlays.
    """
    target_qr = np.asarray(reference_image["qr"], float)
    target_qz = np.asarray(reference_image["qz"], float)
    qr_grid, qz_grid = np.meshgrid(target_qr, target_qz)
    target = np.column_stack((qr_grid.ravel(), qz_grid.ravel()))
    inverse = np.linalg.inv(matrix)
    source = (target - offset) @ inverse.T
    x, y = _q_to_pixel(image, source[:, 0], source[:, 1])
    coordinates = np.vstack((y, x))

    scientific_source = np.asarray(image["intensity"], float)
    display_source = np.asarray(image.get("display_intensity", image["intensity"]), float)

    warped_scientific = map_coordinates(
        scientific_source, coordinates, order=1,
        mode="constant", cval=np.nan,
    ).reshape(qr_grid.shape)
    warped_display = map_coordinates(
        display_source, coordinates, order=1,
        mode="constant", cval=np.nan,
    ).reshape(qr_grid.shape)

    valid = map_coordinates(
        np.asarray(image["valid"], float), coordinates, order=0,
        mode="constant", cval=0.0,
    ).reshape(qr_grid.shape) > 0.5
    valid &= np.isfinite(warped_scientific)

    scientific_intensity = np.where(valid, warped_scientific, np.nan).astype(np.float32)
    display_intensity = np.asarray(warped_display, dtype=np.float32)

    return {
        "path": image.get("path"),
        "source_kind": "registered_" + str(image.get("source_kind", "image")),
        "crop": image.get("crop"),
        "rgb": None,
        "intensity": scientific_intensity,
        "display_intensity": display_intensity,
        "raw_intensity": None,
        "valid": valid,
        "qr": target_qr,
        "qz": target_qz,
    }


def register_measurement_series(all_features, image_records, config):
    """Register each scan to its median-angle reference and update feature q coordinates."""
    features = all_features.copy()
    diagnostics = []
    if not config.enable_per_image_registration:
        for record in image_records:
            record["registered_image"] = record["image"]
        return features, image_records, pd.DataFrame()
    for series_id in sorted({record["series_id"] for record in image_records}):
        records = [record for record in image_records if record["series_id"] == series_id]
        series_indices = features.index[features["series_id"].astype(str) == str(series_id)]
        raw_series_features = features.loc[series_indices].copy()
        reference_record = records[int(np.argmin(np.abs(
            np.asarray([record["angle_deg"] for record in records], float)
            - np.median([record["angle_deg"] for record in records])
        )))]
        reference_scan = int(reference_record["scan"])
        reference_features = features[features["scan"].astype(int) == reference_scan]
        for record in records:
            scan = int(record["scan"])
            source_features = features[features["scan"].astype(int) == scan]
            identity = np.eye(2, dtype=float)
            zero = np.zeros(2, dtype=float)
            if scan == reference_scan:
                registration = {"matrix": identity, "offset": zero, "matches": len(source_features),
                                "median_before": 0.0, "median_after": 0.0,
                                "accepted": True, "reason": "reference_image"}
            elif config.registration_png_only and record["image"].get("source_kind") != "png":
                registration = {"matrix": identity, "offset": zero, "matches": 0,
                                "median_before": np.nan, "median_after": np.nan,
                                "accepted": False, "reason": "numerical_input_not_registered"}
            else:
                registration = _estimate_feature_registration(source_features, reference_features, config)
            matrix = np.asarray(registration["matrix"], float)
            offset = np.asarray(registration["offset"], float)
            if not registration["accepted"]:
                matrix, offset = identity, zero
            transformed = _apply_registration_to_features(source_features, matrix, offset)
            features.loc[transformed.index, transformed.columns] = transformed
            record["registration_matrix"] = matrix
            record["registration_offset"] = offset
            record["registration_reference_scan"] = reference_scan
            record["registration_accepted"] = bool(registration["accepted"])
            record["registration_reason"] = registration["reason"]
            record["registered_image"] = _warp_registered_image(
                record["image"], matrix, offset, reference_record["image"]
            )
        # Downstream scientific analysis uses the registered raster. Both expected
        # dictionary keys therefore reference the same registered image object,
        # avoiding a redundant full-resolution copy in memory.
            record["image"] = record["registered_image"]
            diagnostics.append({
                "series_id": series_id, "scan": scan, "angle_deg": float(record["angle_deg"]),
                "reference_scan": reference_scan, "registration_accepted": bool(registration["accepted"]),
                "registration_reason": registration["reason"], "matched_features": int(registration["matches"]),
                "median_residual_before": registration["median_before"],
                "median_residual_after": registration["median_after"],
                "scale": math.sqrt(abs(float(np.linalg.det(matrix)))),
                "rotation_deg": math.degrees(math.atan2(matrix[1, 0], matrix[0, 0])),
                "qr_offset": float(offset[0]), "qz_offset": float(offset[1]),
            })
        # Registration may improve individual pair residuals yet slightly damage
        # the global consensus. Only use registered feature coordinates for the
        # indexing solver when the full series improves coherently.
        transformed_series = features.loc[series_indices].copy()
        raw_consensus, _ = build_consensus(raw_series_features, config)
        transformed_consensus, _ = build_consensus(transformed_series, config)

        def consensus_metrics(frame):
            if frame.empty:
                return 0, 0, np.inf
            scatter = np.sqrt(np.maximum(frame.cov_rr.to_numpy(float) + frame.cov_zz.to_numpy(float), 1e-12))
            return int(frame.support.sum()), int((frame.support >= 3).sum()), float(np.median(scatter))

        raw_support, raw_support3, raw_scatter = consensus_metrics(raw_consensus)
        reg_support, reg_support3, reg_scatter = consensus_metrics(transformed_consensus)
        scatter_improvement = (raw_scatter - reg_scatter) / max(raw_scatter, 1e-12)
        no_support_loss = (reg_support >= raw_support and reg_support3 >= raw_support3)
        use_for_indexing = (
                np.isfinite(scatter_improvement)
                and scatter_improvement >= float(config.registration_min_series_scatter_improvement_fraction)
                and (no_support_loss or not config.registration_require_no_consensus_support_loss)
        )
        if not use_for_indexing:
            restore_columns = [column for column in (
                "qr", "qz", "cov_rr", "cov_rz", "cov_zz", "sigma_qr", "sigma_qz",
                "major_width_q", "minor_width_q",
            ) if column in raw_series_features]
            features.loc[series_indices, restore_columns] = raw_series_features[restore_columns]
        features.loc[series_indices, "registration_used_for_indexing"] = bool(use_for_indexing)
        for row in diagnostics:
            if row["series_id"] == series_id:
                row.update({
                    "registration_used_for_indexing": bool(use_for_indexing),
                    "raw_consensus_support_sum": raw_support,
                    "registered_consensus_support_sum": reg_support,
                    "raw_consensus_support_ge3": raw_support3,
                    "registered_consensus_support_ge3": reg_support3,
                    "raw_consensus_median_scatter_q": raw_scatter,
                    "registered_consensus_median_scatter_q": reg_scatter,
                    "series_scatter_improvement_fraction": scatter_improvement,
                })
    return features.reset_index(drop=True), image_records, pd.DataFrame(diagnostics)


def registered_composite_image(records, config):
    """Build separate scientific and visual multi-angle composites.

    ``intensity`` remains masked by the scientific validity mask for downstream
    analysis. ``display_intensity`` is composited from the unmasked display
    arrays and is used only for plotting. This prevents the analysis exclusion
    around the specular/central q_r region from appearing as a solid dark-blue
    block in the final overlay.
    """
    images = [record.get("registered_image", record["image"]) for record in records]
    reference = images[int(np.argmin(np.abs(
        np.asarray([record["angle_deg"] for record in records], float)
        - np.median([record["angle_deg"] for record in records])
    )))]
    percentile = float(config.registered_composite_percentile)

    scientific_stack = np.stack([
        np.asarray(image["intensity"], np.float32) for image in images
    ])
    display_stack = np.stack([
        np.asarray(image.get("display_intensity", image["intensity"]), np.float32)
        for image in images
    ])

    def combine(stack):
        with np.errstate(all="ignore"):
            if percentile >= 99.999:
                safe = np.where(np.isfinite(stack), stack, -np.inf)
                combined = np.max(safe, axis=0)
                combined[~np.isfinite(combined)] = np.nan
            else:
                combined = np.nanpercentile(stack, percentile, axis=0)
        return np.asarray(combined, dtype=np.float32)

    scientific_intensity = combine(scientific_stack)
    display_intensity = combine(display_stack)
    valid = np.any(np.stack([
        np.asarray(image["valid"], bool) for image in images
    ]), axis=0)

    # Apply the analysis mask only to the scientific array. The display array is
    # intentionally left unmasked so the central q-space information remains
    # visible beneath the overlay markers.
    scientific_intensity = np.where(valid, scientific_intensity, np.nan).astype(np.float32)

    return {
        "path": "registered_multi_angle_composite",
        "source_kind": "registered_composite",
        "crop": None,
        "rgb": None,
        "intensity": scientific_intensity,
        "display_intensity": display_intensity,
        "raw_intensity": None,
        "valid": valid,
        "qr": reference["qr"],
        "qz": reference["qz"],
    }


def _v739_build_consensus(features: pd.DataFrame, config: IndexingConfig):
    """Return consensus features plus every angle-specific member assignment."""
    columns = [
        "feature_id", "qr", "qz", "cov_rr", "cov_rz", "cov_zz", "sigma_qr", "sigma_qz",
        "strength", "support", "support_fraction", "angles", "feature_type",
        "major_width_q", "minor_width_q", "detection_source", "detection_source_mix",
        "experimental_integrated_intensity", "experimental_intensity_sigma",
        "experimental_integrated_snr", "experimental_intensity_quality_score",
    ]
    if features.empty:
        return pd.DataFrame(columns=columns), features.assign(feature_id=pd.Series(dtype=str))
    coordinates = features[["qr", "qz"]].to_numpy(float)
    if len(coordinates) == 1:
        cluster_ids = np.ones(1, int)
    else:
        parent = np.arange(len(coordinates), dtype=int)
        rank = np.zeros(len(coordinates), dtype=np.int8)
        cluster_members = [{i} for i in range(len(coordinates))]

        def find_root(index):
            while parent[index] != index:
                parent[index] = parent[parent[index]]
                index = parent[index]
            return index

        def union(left, right):
            left, right = find_root(int(left)), find_root(int(right))
            if left == right:
                return
            candidate = cluster_members[left] | cluster_members[right]
            if bool(getattr(config, "consensus_prevent_chaining", True)) and len(candidate) > 2:
                points = coordinates[np.fromiter(candidate, dtype=int)]
                diameter = float(np.max(pdist(points))) if len(points) > 1 else 0.0
                limit = float(config.consensus_tolerance_q) * float(
                    getattr(config, "consensus_cluster_diameter_factor", 1.6)
                )
                if diameter > limit:
                    return
            if rank[left] < rank[right]:
                left, right = right, left
            parent[right] = left
            cluster_members[left] = candidate
            cluster_members[right] = set()
            if rank[left] == rank[right]:
                rank[left] += 1

        pairs = cKDTree(coordinates).query_pairs(config.consensus_tolerance_q, output_type="ndarray")
        for left, right in pairs:
            union(left, right)
        roots = np.array([find_root(i) for i in range(len(coordinates))])
        _, cluster_ids = np.unique(roots, return_inverse=True)
        cluster_ids = cluster_ids + 1
    working = features.copy();
    working["cluster"] = cluster_ids
    rows, members = [], []
    total_angles = max(working["angle_deg"].nunique(), 1)
    for _, group in working.groupby("cluster"):
        group = group.loc[group.groupby("angle_deg")["strength"].idxmax()].copy()
        support = int(group["angle_deg"].nunique())
        if support < config.min_angle_support:
            continue
        weight = np.maximum(group["strength"].to_numpy(float), 1e-6);
        weight /= weight.sum()
        qr, qz = float(np.sum(weight * group["qr"])), float(np.sum(weight * group["qz"]))
        cov = np.zeros((2, 2))
        for w, row in zip(weight, group.itertuples(index=False)):
            local = np.array([[row.cov_rr, row.cov_rz], [row.cov_rz, row.cov_zz]])
            delta = np.array([row.qr - qr, row.qz - qz])
            cov += w * (local + np.outer(delta, delta))
        type_weight = group.groupby("feature_type")["strength"].sum()
        feature_type = str(type_weight.idxmax())
        if "source_detector" in group.columns:
            source_weight = group.groupby(group["source_detector"].astype(str))["strength"].sum()
            detection_source = str(source_weight.idxmax())
            detection_source_mix = ",".join(sorted(set(group["source_detector"].astype(str))))
        else:
            detection_source = "unreported"
            detection_source_mix = "unreported"
        feature_id = f"F{len(rows) + 1:03d}"
        rows.append({
            "feature_id": feature_id, "qr": qr, "qz": qz,
            "cov_rr": float(cov[0, 0]), "cov_rz": float(cov[0, 1]), "cov_zz": float(cov[1, 1]),
            "sigma_qr": math.sqrt(max(cov[0, 0], 1e-12)), "sigma_qz": math.sqrt(max(cov[1, 1], 1e-12)),
            "strength": float(group["strength"].max()), "support": support,
            "support_fraction": support / total_angles,
            "angles": ",".join(f"{x:.3f}" for x in sorted(group["angle_deg"].unique())),
            "feature_type": feature_type,
            "major_width_q": float(np.average(group["major_width_q"], weights=weight)),
            "minor_width_q": float(np.average(group["minor_width_q"], weights=weight)),
            "detection_source": detection_source,
            "detection_source_mix": detection_source_mix,
            "experimental_integrated_intensity": float(np.nanmedian(pd.to_numeric(
                group.get("experimental_integrated_intensity", pd.Series(np.nan, index=group.index)), errors="coerce"
            ))) if np.isfinite(pd.to_numeric(group.get("experimental_integrated_intensity", pd.Series(np.nan, index=group.index)), errors="coerce")).any() else np.nan,
            "experimental_intensity_sigma": float(np.nanmedian(pd.to_numeric(
                group.get("experimental_intensity_sigma", pd.Series(np.nan, index=group.index)), errors="coerce"
            ))) if np.isfinite(pd.to_numeric(group.get("experimental_intensity_sigma", pd.Series(np.nan, index=group.index)), errors="coerce")).any() else np.nan,
            "experimental_integrated_snr": float(np.nanmedian(pd.to_numeric(
                group.get("experimental_integrated_snr", pd.Series(np.nan, index=group.index)), errors="coerce"
            ))) if np.isfinite(pd.to_numeric(group.get("experimental_integrated_snr", pd.Series(np.nan, index=group.index)), errors="coerce")).any() else np.nan,
            "experimental_intensity_quality_score": float(np.nanmedian(pd.to_numeric(
                group.get("experimental_intensity_quality_score", pd.Series(np.nan, index=group.index)), errors="coerce"
            ))) if np.isfinite(pd.to_numeric(group.get("experimental_intensity_quality_score", pd.Series(np.nan, index=group.index)), errors="coerce")).any() else np.nan,
        })
        group["feature_id"] = feature_id;
        members.append(group)
    consensus = pd.DataFrame(rows)
    member_table = pd.concat(members, ignore_index=True) if members else working.iloc[0:0].assign(feature_id="")
    if consensus.empty:
        return pd.DataFrame(columns=columns), member_table
    consensus["rank_weight"] = consensus["strength"] * np.sqrt(consensus["support"])
    keep_ids = set(consensus.nlargest(config.max_consensus_features, "rank_weight")["feature_id"])
    consensus = consensus[consensus["feature_id"].isin(keep_ids)].drop(columns="rank_weight")
    consensus = consensus.sort_values(["support", "strength"], ascending=False).reset_index(drop=True)
    member_table = member_table[member_table["feature_id"].isin(keep_ids)].reset_index(drop=True)
    return consensus[columns], member_table


def reciprocal_basis(cell: gemmi.UnitCell) -> np.ndarray:
    a, b, c = cell.a, cell.b, cell.c
    alpha, beta, gamma = map(math.radians, (cell.alpha, cell.beta, cell.gamma))
    a_vec = np.array([a, 0.0, 0.0])
    b_vec = np.array([b * math.cos(gamma), b * math.sin(gamma), 0.0])
    c_x = c * math.cos(beta)
    c_y = c * (math.cos(alpha) - math.cos(beta) * math.cos(gamma)) / math.sin(gamma)
    c_z = math.sqrt(max(c * c - c_x * c_x - c_y * c_y, 0.0))
    direct = np.column_stack((a_vec, b_vec, [c_x, c_y, c_z]))
    return 2.0 * np.pi * np.linalg.inv(direct).T


def load_reflections(config: IndexingConfig) -> dict:
    path = _resolve_file(config.cif_path, _search_roots(config))
    structure = gemmi.read_small_structure(str(path))
    spacegroup = gemmi.find_spacegroup_by_name(structure.spacegroup_hm)
    if spacegroup is None:
        raise ValueError(f"Could not identify space group {structure.spacegroup_hm!r}")
    basis = reciprocal_basis(structure.cell)
    smallest_reciprocal_singular_value = max(float(np.linalg.svd(basis, compute_uv=False).min()), 1e-12)
    max_index = int(np.ceil(config.q_max / smallest_reciprocal_singular_value)) + 1
    values = np.arange(-max_index, max_index + 1)
    h, k, l = np.meshgrid(values, values, values, indexing="ij")
    hkl = np.column_stack((h.ravel(), k.ravel(), l.ravel()))
    g = hkl @ basis.T;
    q = np.linalg.norm(g, axis=1)
    geometric = (q > 1e-8) & (q <= config.q_max)
    hkl, g, q = hkl[geometric], g[geometric], q[geometric]
    operations = spacegroup.operations()
    gemmi_allowed = np.array([not operations.is_systematically_absent(tuple(map(int, row))) for row in hkl])
    all230_allowed = np.array([
        all230_reflection_allowed(int(spacegroup.number), *map(int, row), occupied_wyckoff=config.occupied_wyckoff)
        for row in hkl
    ]) if config.all230_compare else gemmi_allowed.copy()
    agreement = gemmi_allowed == all230_allowed
    reference_setting = bool(spacegroup.is_reference_setting())
    policy = config.all230_policy.lower()
    if policy not in {"gemmi", "agreement", "all230"}:
        raise ValueError("all230_policy must be 'gemmi', 'agreement', or 'all230'")
    keep = all230_allowed if policy == "all230" else (
        gemmi_allowed & all230_allowed if policy == "agreement" and reference_setting else gemmi_allowed
    )
    hkl_kept, g_kept, q_kept = hkl[keep], g[keep], q[keep]
    calculator = gemmi.StructureFactorCalculatorX(structure.cell)
    f2 = np.array([
        abs(calculator.calculate_sf_from_small_structure(structure, tuple(map(int, row)))) ** 2
        for row in hkl_kept
    ])
    if len(f2) and config.structure_factor_zero_percent > 0:
        f2_keep = f2 >= f2.max() * config.structure_factor_zero_percent / 100.0
        hkl_kept, g_kept, q_kept, f2 = hkl_kept[f2_keep], g_kept[f2_keep], q_kept[f2_keep], f2[f2_keep]
    reflections = pd.DataFrame({
        "h": hkl_kept[:, 0], "k": hkl_kept[:, 1], "l": hkl_kept[:, 2],
        "hkl": [f"({a} {b} {c})" for a, b, c in hkl_kept],
        "gx": g_kept[:, 0], "gy": g_kept[:, 1], "gz": g_kept[:, 2],
        "q": q_kept, "d": 2.0 * np.pi / q_kept, "f2": f2,
    })
    validation = pd.DataFrame({
        "h": hkl[:, 0], "k": hkl[:, 1], "l": hkl[:, 2], "q": q,
        "gemmi_allowed": gemmi_allowed, "all230_allowed": all230_allowed, "agreement": agreement,
    })
    return {
        "path": path, "structure": structure, "spacegroup": spacegroup, "basis": basis,
        "reflections": reflections.sort_values("q").reset_index(drop=True),
        "validation": validation, "validation_disagreements": validation[~agreement].copy(),
        "reference_setting": reference_setting, "all230_policy_used": policy,
    }


def _primitive_hkl(hkl: Iterable[int]) -> tuple[int, int, int]:
    hkl = np.asarray(tuple(hkl), int)
    divisor = math.gcd(math.gcd(abs(int(hkl[0])), abs(int(hkl[1]))), abs(int(hkl[2]))) or 1
    hkl //= divisor
    first = next((x for x in hkl if x), 1)
    if first < 0:
        hkl *= -1
    return tuple(map(int, hkl))


def orientation_candidates(crystal: dict, config: IndexingConfig) -> list[tuple[int, int, int]]:
    subset = crystal["reflections"]
    subset = subset[subset[["h", "k", "l"]].abs().max(axis=1) <= config.orientation_hkl_max].copy()
    subset["primitive"] = [_primitive_hkl(row) for row in subset[["h", "k", "l"]].itertuples(index=False, name=None)]
    subset = subset.drop_duplicates("primitive")
    f2_scale = np.log1p(subset["f2"] / max(float(subset["f2"].median()), 1e-9))
    subset["candidate_rank"] = subset["q"] / (0.35 + f2_scale)
    return subset.nsmallest(config.max_orientation_candidates, "candidate_rank")["primitive"].tolist()


def _align_hkl_to_z(hkl: tuple[int, int, int], basis: np.ndarray) -> np.ndarray:
    """Pure-NumPy Rodrigues rotation that maps the reciprocal hkl vector to +z."""
    vector = basis @ np.asarray(hkl, float)
    vector /= np.linalg.norm(vector)
    target = np.array([0.0, 0.0, 1.0])
    cross = np.cross(vector, target)
    sine = float(np.linalg.norm(cross))
    cosine = float(np.clip(vector @ target, -1.0, 1.0))
    if sine < 1e-12:
        return np.eye(3) if cosine > 0 else np.diag([1.0, -1.0, -1.0])
    kx, ky, kz = cross
    skew = np.array([[0.0, -kz, ky], [kz, 0.0, -kx], [-ky, kx, 0.0]])
    return np.eye(3) + skew + (skew @ skew) * ((1.0 - cosine) / (sine * sine))


def _tilt_matrix(tilt_x_deg: float, tilt_y_deg: float) -> np.ndarray:
    x, y = math.radians(tilt_x_deg), math.radians(tilt_y_deg)
    cx, sx, cy, sy = math.cos(x), math.sin(x), math.cos(y), math.sin(y)
    rotation_x = np.array([[1.0, 0.0, 0.0], [0.0, cx, -sx], [0.0, sx, cx]])
    rotation_y = np.array([[cy, 0.0, sy], [0.0, 1.0, 0.0], [-sy, 0.0, cy]])
    return rotation_y @ rotation_x


def _deduplicate_mirrored_features(frame: pd.DataFrame, tolerance: float = 0.028) -> pd.DataFrame:
    """Use one representative of ±qr duplicate features when constructing anchors.

    Anchor selection must preserve any caller-supplied ``_rank`` values because
    they encode feature priority. Keeping the original DataFrame row indices
    avoids pandas tuple-field renaming and ensures ranking metadata remains
    attached to the correct reciprocal-space feature.
    """
    if frame.empty:
        return frame.copy()

    ranked = frame.copy()
    caller_supplied_rank = "_rank" in ranked.columns
    if not caller_supplied_rank:
        ranked["_rank"] = (
                2 * ranked.support_fraction + ranked.strength
                - 4 * np.sqrt(np.maximum(ranked.cov_rr + ranked.cov_zz, 1e-12))
        )
    ranked = ranked.sort_values("_rank", ascending=False, kind="mergesort")

    kept_indices = []
    kept_positions = []
    for index, row in ranked.iterrows():
        position = (abs(float(row.qr)), float(row.qz))
        if all(math.hypot(position[0] - previous[0], position[1] - previous[1]) > tolerance
               for previous in kept_positions):
            kept_indices.append(index)
            kept_positions.append(position)

    result = ranked.loc[kept_indices].copy()
    if not caller_supplied_rank:
        result = result.drop(columns="_rank", errors="ignore")
    return result.reset_index(drop=True)


def _nearest_radial_residual(features: pd.DataFrame, crystal: dict, percentile: float) -> np.ndarray:
    if features.empty:
        return np.array([])
    reflections = crystal["reflections"]
    strong = reflections[reflections.f2 >= reflections.f2.quantile(percentile / 100.0)]
    q_values = np.sort(strong.q.unique())
    measured = np.hypot(features.qr.to_numpy(float), features.qz.to_numpy(float))
    index = np.searchsorted(q_values, measured)
    left = q_values[np.clip(index - 1, 0, len(q_values) - 1)]
    right = q_values[np.clip(index, 0, len(q_values) - 1)]
    return np.minimum(abs(measured - left), abs(measured - right))


def _stable_holdout_value(feature_id: str, seed: int) -> float:
    digest = hashlib.sha256(f"{seed}:{feature_id}".encode("utf-8")).digest()
    return int.from_bytes(digest[:8], "big") / float(2 ** 64)


def _ignored_reason(row, config: IndexingConfig) -> str:
    reasons = []
    sigma = math.sqrt(max(float(row.cov_rr) + float(row.cov_zz), 1e-12))
    if row.feature_type in config.ignored_feature_types:
        reasons.append("excluded_feature_type")
    if int(row.support) < config.anchor_min_support:
        reasons.append("insufficient_anchor_support")
    if float(row.qz) < config.anchor_min_qz:
        reasons.append("low_qz_for_anchor")
    if abs(float(row.qr)) < config.anchor_min_abs_qr:
        reasons.append("near_specular_for_anchor")
    if sigma > config.anchor_max_sigma_q:
        reasons.append("high_position_uncertainty")
    if float(row.major_width_q) > config.anchor_max_major_width_q:
        reasons.append("broad_or_streak_like")
    radial = getattr(row, "radial_prior_residual", np.nan)
    if np.isfinite(radial) and radial > config.anchor_radial_tolerance_q:
        reasons.append("cif_radial_prior_mismatch")
    if not reasons:
        reasons.append("not_selected_by_role_capacity")
    return ";".join(reasons)


def feature_roles(consensus: pd.DataFrame, crystal: dict, config: IndexingConfig):
    """Create a data-only validation holdout before using any CIF-derived prior."""
    frame = consensus[~consensus.feature_id.isin(config.manual_rejected_feature_ids)].copy()
    if frame.empty:
        return frame, frame, frame
    frame["data_role_score"] = (
            2.0 * frame.support_fraction + 0.9 * frame.strength
            - 7.0 * np.sqrt(np.maximum(frame.cov_rr + frame.cov_zz, 1e-12))
            + 0.35 * (frame.feature_type == "spot").astype(float)
    )
    eligible = frame[
        ~frame.feature_type.isin(config.ignored_feature_types)
        & (frame.support >= config.min_angle_support)
        ].copy()
    eligible["_holdout"] = [
        _stable_holdout_value(feature_id, config.random_seed)
        for feature_id in eligible.feature_id
    ]
    desired = min(
        config.max_validation_features,
        max(config.validation_holdout_min_features,
            int(round(config.validation_holdout_fraction * len(eligible))))
    )
    desired = min(desired, max(0, len(eligible) - config.min_anchor_matches))
    validation = eligible.nsmallest(desired, "_holdout").sort_values(
        "data_role_score", ascending=False
    ).copy()

    # Compute the CIF radial diagnostic only after validation membership has
    # already been chosen, preserving a data-only holdout split.
    frame["radial_prior_residual"] = _nearest_radial_residual(
        frame, crystal, config.hypothesis_f2_percentile
    )
    anchor_source = frame[~frame.feature_id.isin(validation.feature_id)].copy()
    anchor_source["anchor_role_score"] = (
            anchor_source.data_role_score - 5.0 * anchor_source.radial_prior_residual
    )
    anchor_pool = _deduplicate_mirrored_features(anchor_source)
    sigma = np.sqrt(np.maximum(anchor_pool.cov_rr + anchor_pool.cov_zz, 1e-12))
    objective = (
            (anchor_pool.support >= config.anchor_min_support)
            & (anchor_pool.qz >= config.anchor_min_qz)
            & (np.abs(anchor_pool.qr) >= config.anchor_min_abs_qr)
            & (sigma <= config.anchor_max_sigma_q)
            & (anchor_pool.major_width_q <= config.anchor_max_major_width_q)
            & (anchor_pool.radial_prior_residual <= config.anchor_radial_tolerance_q)
            & (~anchor_pool.feature_type.isin(config.ignored_feature_types))
    )
    if objective.any():
        objective &= anchor_pool.strength >= anchor_pool.loc[objective, "strength"].quantile(
            config.anchor_strength_quantile
        )
    if config.manual_anchor_feature_ids:
        objective |= anchor_pool.feature_id.isin(config.manual_anchor_feature_ids)
    anchors = anchor_pool[objective].nlargest(
        config.max_anchor_features, "anchor_role_score"
    ).copy()
    used = set(anchors.feature_id) | set(validation.feature_id)
    ignored = frame[~frame.feature_id.isin(used)].copy()
    if not ignored.empty:
        ignored["ignored_reason"] = [
            _ignored_reason(row, config) for row in ignored.itertuples(index=False)
        ]
    drop_columns = ["data_role_score", "anchor_role_score", "_holdout"]
    return (
        anchors.drop(columns=drop_columns, errors="ignore"),
        validation.drop(columns=drop_columns, errors="ignore"),
        ignored.drop(columns=drop_columns, errors="ignore"),
    )


def normal_vector(hkl, crystal):
    vector = crystal["basis"] @ np.asarray(hkl, float)
    return vector / np.linalg.norm(vector)


def normal_angle(hkl_a, hkl_b, crystal):
    cosine = abs(float(normal_vector(hkl_a, crystal) @ normal_vector(hkl_b, crystal)))
    return math.degrees(math.acos(np.clip(cosine, -1, 1)))


def _unpack(params=None):
    """Return tilt and reciprocal-space calibration parameters.

    Five-parameter vectors are interpreted as a shared q scale for both axes.
    Six-parameter searches permit q_r and q_z scale mismatch to be refined
    independently, which is useful when detector-axis calibration differs.
    """
    if params is None:
        return (0.0, 0.0, 1.0, 1.0, 0.0, 0.0)
    values = tuple(map(float, params))
    if len(values) == 5:
        tx, ty, scale, offset_r, offset_z = values
        return tx, ty, scale, scale, offset_r, offset_z
    if len(values) == 8:
        return values[:6]
    if len(values) != 6:
        raise ValueError(f"Expected 5, 6, or 8 calibration parameters, received {len(values)}")
    return values


def _crystal_arrays(crystal):
    cached = crystal.get("_array_cache")
    if cached is not None:
        return cached
    frame = crystal["reflections"]
    cached = {
        "h": frame.h.to_numpy(np.int32), "k": frame.k.to_numpy(np.int32),
        "l": frame.l.to_numpy(np.int32), "hkl": frame.hkl.to_numpy(object),
        "g": frame[["gx", "gy", "gz"]].to_numpy(float),
        "q": frame.q.to_numpy(float), "d": frame.d.to_numpy(float),
        "f2": frame.f2.to_numpy(float),
    }
    crystal["_array_cache"] = cached
    crystal["_align_cache"] = {}
    return cached


def _base_alignment(crystal, normal_hkl):
    key = tuple(map(int, normal_hkl))
    cache = crystal.setdefault("_align_cache", {})
    if key not in cache:
        cache[key] = _align_hkl_to_z(key, crystal["basis"])
    return cache[key]


def _complex_kz(k0: float, alpha_deg: np.ndarray | float, delta: float, beta: float):
    alpha = np.radians(np.asarray(alpha_deg, float))
    return k0 * np.sqrt(np.sin(alpha) ** 2 - 2.0 * float(delta) + 2.0j * float(beta) + 0j)


def _fresnel_coefficients(kz_a, kz_b):
    denominator = kz_a + kz_b
    denominator = np.where(np.abs(denominator) < 1e-18, 1e-18 + 0j, denominator)
    return (kz_a - kz_b) / denominator, 2.0 * kz_a / denominator


def _layer_field_components(alpha_deg, config: IndexingConfig):
    """Forward/backward film fields for one air/film/substrate layer."""
    required = (
        config.film_delta, config.film_beta,
        config.substrate_delta, config.substrate_beta,
        config.film_thickness_A,
    )
    if any(value is None for value in required):
        raise ValueError(
            "enable_dwba=True requires film_delta, film_beta, substrate_delta, "
            "substrate_beta, and film_thickness_A. Use energy-specific optical constants."
        )
    k0 = 2.0 * np.pi / float(config.xray_wavelength_A)
    alpha = np.asarray(alpha_deg, float)
    kz0 = k0 * np.sin(np.radians(alpha)) + 0j
    kz1 = _complex_kz(k0, alpha, config.film_delta, config.film_beta)
    kz2 = _complex_kz(k0, alpha, config.substrate_delta, config.substrate_beta)
    r01, t01 = _fresnel_coefficients(kz0, kz1)
    r12, _ = _fresnel_coefficients(kz1, kz2)
    thickness = float(config.film_thickness_A)
    z = np.clip(float(config.dwba_scatter_depth_fraction), 0.0, 1.0) * thickness
    phase = np.exp(2.0j * kz1 * thickness)
    denominator = 1.0 + r01 * r12 * phase
    denominator = np.where(np.abs(denominator) < 1e-18, 1e-18 + 0j, denominator)
    forward = t01 * np.exp(1.0j * kz1 * z) / denominator
    backward = t01 * r12 * phase * np.exp(-1.0j * kz1 * z) / denominator
    return forward, backward


def refraction_corrected_qz(qz_internal: np.ndarray, config: IndexingConfig) -> np.ndarray:
    """Approximate external qz for the direct/direct film channel.

    The crystal projection supplies an internal-film reciprocal-vector qz.
    Conservation of the in-plane wavevector and the film refractive index are
    used to map the internal exit kz back into air. This is a single-interface,
    direct/direct position correction, not a complete four-channel DWBA peak
    solver. It is applied only when explicitly enabled with calibrated film
    optical constants.
    """
    qz_internal = np.asarray(qz_internal, float)
    if not config.enable_refraction_position_correction:
        return qz_internal.copy()
    if config.refraction_position_channel != "direct_direct":
        raise ValueError("Only refraction_position_channel='direct_direct' is implemented")
    if config.film_delta is None or config.film_beta is None:
        raise ValueError(
            "enable_refraction_position_correction=True requires film_delta and film_beta "
            "for the measurement wavelength"
        )
    wavelength = float(config.xray_wavelength_A)
    if wavelength <= 0:
        raise ValueError("xray_wavelength_A must be positive")
    k0 = 2.0 * np.pi / wavelength
    angles = tuple(float(value) for value in getattr(config, "giwaxs_incidence_angles_deg", ())
                   if np.isfinite(float(value)) and float(value) > 0)
    if not angles:
        angles = (float(config.incidence_angle_deg),)
    corrected_by_angle = []
    for angle_deg in angles:
        alpha_i = math.radians(angle_deg)
        kz_air_in = k0 * math.sin(alpha_i)
        kz_film_in = _complex_kz(k0, angle_deg, config.film_delta, config.film_beta)
        kz_film_out = kz_film_in + qz_internal
        # kz_film^2 = kz_air^2 - 2*delta*k0^2 + 2i*beta*k0^2.
        kz_air_out = np.sqrt(
            kz_film_out ** 2
            + 2.0 * float(config.film_delta) * k0 ** 2
            - 2.0j * float(config.film_beta) * k0 ** 2
            + 0j
        )
        kz_air_out = np.where(np.real(kz_air_out) < 0, -kz_air_out, kz_air_out)
        corrected_by_angle.append(np.real(kz_air_out) - kz_air_in)
    corrected = np.nanmedian(np.stack(corrected_by_angle, axis=0), axis=0)
    return np.where(np.isfinite(corrected), corrected, qz_internal)


def dwba_intensity_envelope(qz: np.ndarray, config: IndexingConfig) -> np.ndarray:
    """Coherent four-term Fresnel DWBA envelope for projected peak intensities.

    The product (Ei+ + Ei-)(Ef+ + Ef-) expands into four direct/reflected
    channels. This approximation changes intensity ranking only; peak positions
    remain kinematic. A full DWBA/dynamical solver is still needed for quantitative
    line shapes, refraction-shifted peak positions, and multilayer form factors.
    """
    qz = np.asarray(qz, float)
    if not config.enable_dwba:
        return np.ones_like(qz)
    k0 = 2.0 * np.pi / float(config.xray_wavelength_A)
    angles = tuple(float(value) for value in getattr(config, "giwaxs_incidence_angles_deg", ())
                   if np.isfinite(float(value)) and float(value) > 0)
    if not angles:
        angles = (float(config.incidence_angle_deg),)
    intensity_by_angle = []
    for angle_deg in angles:
        sin_exit = qz / k0 - math.sin(math.radians(angle_deg))
        valid = (sin_exit >= 0.0) & (sin_exit <= 1.0)
        alpha_exit = np.degrees(np.arcsin(np.clip(sin_exit, 0.0, 1.0)))
        inc_forward, inc_backward = _layer_field_components(angle_deg, config)
        out_forward, out_backward = _layer_field_components(alpha_exit, config)
        amplitude = (inc_forward + inc_backward) * (out_forward + out_backward)
        intensity = np.abs(amplitude) ** 2
        intensity_by_angle.append(np.where(valid & np.isfinite(intensity), intensity, 0.0))
    intensity = np.mean(np.stack(intensity_by_angle, axis=0), axis=0)
    finite_max = float(np.nanmax(intensity)) if intensity.size else 0.0
    if finite_max > 0:
        intensity /= finite_max
    floor = np.clip(float(config.dwba_intensity_floor), 0.0, 1.0)
    return floor + (1.0 - floor) * intensity


def dwba_model_metadata(config: IndexingConfig) -> dict:
    configured = bool(config.enable_dwba and all(value is not None for value in (
        config.film_delta, config.film_beta, config.substrate_delta,
        config.substrate_beta, config.film_thickness_A,
    )))
    return {
        "mode": str(getattr(config, "giwaxs_physics_mode", "manual")),
        "status": str(getattr(config, "giwaxs_physics_status", "")),
        "parameter_source": str(getattr(config, "giwaxs_physics_parameter_source", "")),
        "fallback_reason": str(getattr(config, "giwaxs_physics_fallback_reason", "")),
        "enabled": bool(config.enable_dwba),
        "configured": configured,
        "model": "single_layer_coherent_four_term_fresnel_envelope",
        "changes_peak_positions": bool(config.enable_refraction_position_correction),
        "position_correction_enabled": bool(config.enable_refraction_position_correction),
        "position_correction_model": "single_interface_direct_direct" if config.enable_refraction_position_correction else "disabled",
        "incidence_angle_deg": float(config.incidence_angle_deg),
        "incidence_angles_deg": list(getattr(config, "giwaxs_incidence_angles_deg", ()) or (config.incidence_angle_deg,)),
        "incidence_angle_policy": (
            "AUTO uses each series' measured angle ensemble; refraction positions use the median "
            "of per-angle corrections and DWBA intensity uses the mean per-angle envelope"
            if len(getattr(config, "giwaxs_incidence_angles_deg", ()) or ()) > 1
            else "single measured incidence angle"
        ),
        "xray_wavelength_A": float(config.xray_wavelength_A),
        "film_delta": config.film_delta,
        "film_beta": config.film_beta,
        "substrate_delta": config.substrate_delta,
        "substrate_beta": config.substrate_beta,
        "film_thickness_A": config.film_thickness_A,
        "quantitative_dwba_claim": False,
        "limitations": [
            "single film on one substrate",
            "intensity envelope uses a coherent four-term Fresnel approximation",
            "position correction, when enabled, uses only the direct/direct channel",
            "no full four-channel qz-dependent crystal form factor",
            "no dynamical Bragg diffraction or multilayer peak-position solver",
            "optical constants must match wavelength/energy",
        ],
    }


def simulate_powder_pattern(crystal: dict, config: IndexingConfig):
    """Kinematic powder pattern in q; summing enumerated hkl supplies multiplicity."""
    reflections = crystal["reflections"].copy()
    q_max = float(config.powder_q_max or config.q_max)
    reflections = reflections[
        (reflections.q >= config.powder_q_min) & (reflections.q <= q_max)
        ].copy()
    if reflections.empty:
        return pd.DataFrame(columns=["q", "intensity"]), pd.DataFrame()
    merge = max(float(config.powder_peak_merge_q), 1e-6)
    reflections["q_bin"] = np.rint(reflections.q / merge).astype(int)
    grouped = reflections.groupby("q_bin", as_index=False).agg(
        q=("q", "mean"), raw_intensity=("f2", "sum"), multiplicity=("h", "size"),
        strongest_f2=("f2", "max"), example_hkl=("hkl", "first"),
    )
    intensity = grouped.raw_intensity.to_numpy(float)
    if config.powder_apply_lorentz_polarization:
        wavelength = float(config.xray_wavelength_A)
        sin_theta = np.clip(grouped.q.to_numpy(float) * wavelength / (4.0 * np.pi), 1e-6, 0.999999)
        theta = np.arcsin(sin_theta)
        lp = (1.0 + np.cos(2.0 * theta) ** 2) / np.maximum(
            np.sin(theta) ** 2 * np.cos(theta), 1e-8
        )
        intensity *= lp
    grouped["relative_intensity"] = intensity / max(float(intensity.max()), 1e-12)
    grouped = grouped[grouped.relative_intensity >= config.powder_min_relative_peak].copy()
    q_grid = np.arange(config.powder_q_min, q_max + config.powder_q_step / 2, config.powder_q_step)
    sigma = max(float(config.powder_fwhm_q) / 2.354820045, config.powder_q_step / 2)
    profile = np.zeros_like(q_grid)
    for row in grouped.itertuples(index=False):
        profile += row.relative_intensity * np.exp(-0.5 * ((q_grid - row.q) / sigma) ** 2)
    if profile.max() > 0:
        profile /= profile.max()
    powder = pd.DataFrame({"q": q_grid, "intensity": profile})
    wavelength = float(config.xray_wavelength_A)
    argument = np.clip(q_grid * wavelength / (4.0 * np.pi), 0.0, 1.0)
    powder["two_theta_deg"] = 2.0 * np.degrees(np.arcsin(argument))
    return powder, grouped.sort_values("q").reset_index(drop=True)


def save_powder_outputs(crystal: dict, config: IndexingConfig, output: Path):
    powder, peaks = simulate_powder_pattern(crystal, config)
    powder.to_csv(output / "powder_pattern_q.csv", index=False)
    peaks.to_csv(output / "powder_peaks.csv", index=False)
    if not powder.empty:
        figure, axis = plt.subplots(figsize=(9.0, 4.2))
        axis.plot(powder.q, powder.intensity, linewidth=1.2)
        if not peaks.empty:
            axis.vlines(peaks.q, 0, peaks.relative_intensity, alpha=0.28, linewidth=0.7)
        axis.set(xlabel=r"$q$ ($\AA^{-1}$)", ylabel="normalized intensity",
                 title=f"Kinematic powder diagnostic: {crystal['path'].name}")
        figure.tight_layout()
        figure.savefig(output / "powder_pattern_q.png", dpi=config.dpi)
        plt.close(figure)
    return powder, peaks


def _powder_peak_compatibility(consensus: pd.DataFrame, peaks: pd.DataFrame, config: IndexingConfig):
    if consensus.empty or peaks.empty:
        return {
            "matched_weight_fraction": 0.0, "matched_features": 0,
            "median_radial_residual_q": np.nan, "score": -1.0,
        }
    observed = np.hypot(consensus.qr.to_numpy(float), consensus.qz.to_numpy(float))
    calculated = peaks.q.to_numpy(float)
    residual = np.abs(observed[:, None] - calculated[None, :])
    rows, cols = linear_sum_assignment(residual)
    accepted = residual[rows, cols] <= config.cif_radial_match_tolerance_q
    rows, cols = rows[accepted], cols[accepted]
    weights = consensus.strength.to_numpy(float) * np.sqrt(consensus.support.to_numpy(float))
    weights /= max(float(weights.sum()), 1e-12)
    weighted = float(weights[rows].sum()) if len(rows) else 0.0
    median = float(np.median(residual[rows, cols])) if len(rows) else np.nan
    score = weighted - (0.25 * median / config.cif_radial_match_tolerance_q if np.isfinite(median) else 0.25)
    return {
        "matched_weight_fraction": weighted,
        "matched_features": int(len(rows)),
        "median_radial_residual_q": median,
        "score": float(score),
    }


def compare_cif_candidates(series_results: dict, crystals: list[dict], config: IndexingConfig):
    rows = []
    peak_cache = {}
    for crystal in crystals:
        _, peaks = simulate_powder_pattern(crystal, config)
        peak_cache[str(crystal["path"])] = peaks
    for series_id, result in series_results.items():
        for candidate_index, crystal in enumerate(crystals):
            metrics = _powder_peak_compatibility(
                result["consensus"], peak_cache[str(crystal["path"])], config
            )
            rows.append({
                "series_id": series_id,
                "candidate_index": candidate_index,
                "cif_path": str(crystal["path"]),
                "spacegroup_number": int(crystal["spacegroup"].number),
                **metrics,
            })
    report = pd.DataFrame(rows)
    conclusions = {}
    for series_id, group in report.groupby("series_id") if not report.empty else []:
        ranked = group.sort_values("score", ascending=False).reset_index(drop=True)
        gap = float(ranked.iloc[0].score - ranked.iloc[1].score) if len(ranked) > 1 else np.nan
        primary_rank = int(ranked.index[ranked.candidate_index == 0][0] + 1) if (
                    ranked.candidate_index == 0).any() else np.nan
        if len(ranked) == 1:
            conclusion = "compatible_only_not_identified_no_alternative_cifs_tested"
        elif primary_rank == 1 and gap >= config.cif_preference_min_gap:
            conclusion = "primary_cif_preferred_among_tested_candidates_not_proven_correct"
        else:
            conclusion = "cif_identity_ambiguous_among_tested_candidates"
        conclusions[series_id] = {
            "cif_conclusion": conclusion,
            "primary_cif_rank": primary_rank,
            "best_vs_second_score_gap": gap,
            "candidate_count": int(len(ranked)),
        }
    return report, conclusions


def _empty_predictions():
    return {
        "h": np.empty(0, np.int32), "k": np.empty(0, np.int32), "l": np.empty(0, np.int32),
        "hkl": np.empty(0, object), "q": np.empty(0), "d": np.empty(0), "f2": np.empty(0),
        "qxy": np.empty(0), "qz": np.empty(0), "qr": np.empty(0),
        "dwba_weight": np.empty(0), "effective_intensity": np.empty(0),
        "prediction_weight": np.empty(0), "domain": np.empty(0, object),
    }


def _predictions_frame(predictions):
    if isinstance(predictions, pd.DataFrame):
        return predictions.copy()
    if not len(predictions["qr"]):
        return pd.DataFrame(
            columns=["h", "k", "l", "hkl", "q", "d", "f2", "qxy", "qz", "qr", "dwba_weight", "effective_intensity",
                     "prediction_weight", "domain"])
    return pd.DataFrame({key: value for key, value in predictions.items()})


def _predictions_dict(predictions):
    if isinstance(predictions, dict):
        return predictions
    if predictions is None or predictions.empty:
        return _empty_predictions()
    result = {key: predictions[key].to_numpy() for key in (
        "h", "k", "l", "hkl", "q", "d", "f2", "qxy", "qz", "qr", "prediction_weight"
    )}
    result["dwba_weight"] = (
        predictions["dwba_weight"].to_numpy(float)
        if "dwba_weight" in predictions else np.ones(len(predictions), float)
    )
    result["effective_intensity"] = (
        predictions["effective_intensity"].to_numpy(float)
        if "effective_intensity" in predictions else result["f2"].astype(float)
    )
    result["domain"] = (
        predictions["domain"].astype(str).to_numpy(object)
        if "domain" in predictions else np.full(len(predictions), "primary", dtype=object)
    )
    return result


def project_reflections(crystal, normal_hkl, config: IndexingConfig, params=None,
                        f2_percentile=0.0, maximum=None) -> pd.DataFrame:
    return _predictions_frame(_prediction_array(
        crystal, normal_hkl, config, params, f2_percentile, maximum
    ))


def _feature_distance_matrices(features, predictions, config):
    experimental = features[["qr", "qz"]].to_numpy(float)
    calculated = np.column_stack((predictions["qr"], predictions["qz"]))
    delta = experimental[:, None, :] - calculated[None, :, :]
    physical = np.linalg.norm(delta, axis=2)
    normalized = np.empty_like(physical)
    floor, ceiling = config.uncertainty_floor_q ** 2, config.uncertainty_ceiling_q ** 2
    for i, row in enumerate(features.itertuples(index=False)):
        covariance = np.array([[row.cov_rr, row.cov_rz], [row.cov_rz, row.cov_zz]], float)
        values, vectors = np.linalg.eigh(covariance)
        covariance = vectors @ np.diag(np.clip(values, floor, ceiling)) @ vectors.T
        inverse = np.linalg.inv(covariance)
        normalized[i] = np.sqrt(np.maximum(np.einsum("ni,ij,nj->n", delta[i], inverse, delta[i]), 0))
    return physical, normalized


def _empty_assignment_metrics(score=-1.0):
    return {
        "score": score, "matches": 0, "weighted_fraction": 0.0, "indexed_fraction": 0.0,
        "median_delta_q": np.nan, "p90_delta_q": np.nan, "ambiguity": np.nan,
        "false_prediction_fraction": 1.0, "pair_geometry_rms": np.nan,
    }


def _greedy_one_to_one(cost, gate):
    """Deterministic low-cost one-to-one assignment used only for hypothesis screening."""
    rows, cols = np.where(gate)
    if not len(rows):
        return np.empty(0, int), np.empty(0, int)
    order = np.argsort(cost[rows, cols], kind="stable")
    used_rows, used_cols, selected_rows, selected_cols = set(), set(), [], []
    for position in order:
        row, col = int(rows[position]), int(cols[position])
        if row in used_rows or col in used_cols:
            continue
        used_rows.add(row);
        used_cols.add(col)
        selected_rows.append(row);
        selected_cols.append(col)
    return np.asarray(selected_rows, int), np.asarray(selected_cols, int)


def _assign_arrays(features, predictions, config: IndexingConfig, tolerance, materialize=True):
    if features.empty or not len(predictions["qr"]):
        return pd.DataFrame(), _empty_assignment_metrics()
    physical, normalized = _feature_distance_matrices(features, predictions, config)
    cost = normalized / config.match_sigma_limit + 0.05 * (1 - predictions["prediction_weight"])[None, :]
    locks = {feature: (h, k, l) for feature, h, k, l in config.manual_locked_assignments}
    for i, feature_id in enumerate(features.feature_id):
        if feature_id in locks:
            hkl = locks[feature_id]
            allowed = (
                    (predictions["h"] == hkl[0]) & (predictions["k"] == hkl[1]) & (predictions["l"] == hkl[2])
            )
            cost[i, ~allowed] = 1e6
    gate = (physical <= tolerance) & (normalized <= config.match_sigma_limit)
    if not materialize:
        row_index, column_index = _greedy_one_to_one(cost, gate)
    else:
        real_cost = np.nan_to_num(cost.copy(), nan=1e4, posinf=1e4, neginf=-1e4)
        real_cost[~gate] = 1e4
        # One private dummy column per feature allows an explicit unassigned result.
        dummy_cost = np.full((len(features), len(features)), 1.08, float)
        augmented = np.hstack((real_cost, dummy_cost))
        row_index, column_index = linear_sum_assignment(augmented)
        accepted = (column_index < len(predictions["qr"]))
        row_index, column_index = row_index[accepted], column_index[accepted]
        accepted = gate[row_index, column_index]
        row_index, column_index = row_index[accepted], column_index[accepted]
    feature_weights = features.strength.to_numpy(float) * np.sqrt(features.support.to_numpy(float))
    feature_weights /= max(float(feature_weights.sum()), 1e-12)
    possible = (physical <= tolerance) & (normalized <= config.match_sigma_limit)
    ambiguity = np.maximum(possible.sum(axis=1), 1)
    selected_norm = normalized[row_index, column_index] if len(row_index) else np.empty(0)
    selected_phys = physical[row_index, column_index] if len(row_index) else np.empty(0)
    weighted = float(feature_weights[row_index].sum()) if len(row_index) else 0.0
    unique_columns = np.unique(column_index)
    false_fraction = (
            1 - float(predictions["prediction_weight"][unique_columns].sum())
            / max(float(predictions["prediction_weight"].sum()), 1e-12)
    ) if len(column_index) else 1.0
    exp_coords = features[["qr", "qz"]].to_numpy(float)
    calc_coords = np.column_stack((predictions["qr"], predictions["qz"]))
    geometry = np.nan
    if len(row_index) >= 3:
        geometry = float(np.sqrt(np.mean((pdist(exp_coords[row_index]) - pdist(calc_coords[column_index])) ** 2)))
    mean_error = float(np.average(selected_norm / config.match_sigma_limit, weights=feature_weights[row_index])) if len(
        row_index) else 1.0
    ambiguity_cost = float(np.average(np.log1p(ambiguity[row_index] - 1), weights=feature_weights[row_index])) if len(
        row_index) else 1.0
    score = weighted - 0.32 * mean_error - config.ambiguity_penalty * ambiguity_cost - config.missing_strong_prediction_penalty * false_fraction
    if np.isfinite(geometry):
        score -= config.pair_geometry_penalty * min(geometry / tolerance, 2.0)
    metrics = {
        "score": score, "matches": len(row_index), "weighted_fraction": weighted,
        "indexed_fraction": len(row_index) / len(features),
        "median_delta_q": float(np.median(selected_phys)) if len(row_index) else np.nan,
        "p90_delta_q": float(np.quantile(selected_phys, 0.9)) if len(row_index) else np.nan,
        "ambiguity": float(np.mean(ambiguity[row_index])) if len(row_index) else np.nan,
        "false_prediction_fraction": false_fraction, "pair_geometry_rms": geometry,
    }
    if not materialize or not len(row_index):
        return pd.DataFrame(), metrics
    rows = []
    for i, j in zip(row_index, column_index):
        exp = features.iloc[i]
        competitors = np.sort(normalized[i][possible[i]])
        margin = float(competitors[1] - competitors[0]) if len(competitors) > 1 else config.match_sigma_limit
        assignment_support_score = (
                math.exp(-0.5 * normalized[i, j] ** 2) * min(1.0, exp.support / 3.0)
                * predictions["prediction_weight"][j] * min(1.0, margin / 1.5) / math.sqrt(ambiguity[i])
        )
        rows.append({
            "feature_id": exp.feature_id, "hkl": predictions["hkl"][j],
            "h": int(predictions["h"][j]), "k": int(predictions["k"][j]), "l": int(predictions["l"][j]),
            "qr_exp": exp.qr, "qz_exp": exp.qz,
            "qr_calc": predictions["qr"][j], "qz_calc": predictions["qz"][j],
            "delta_qr": exp.qr - predictions["qr"][j], "delta_qz": exp.qz - predictions["qz"][j],
            "delta_q": physical[i, j], "normalized_delta": normalized[i, j],
            "strength": exp.strength, "support": int(exp.support), "feature_type": exp.feature_type,
            "f2": predictions["f2"][j], "prediction_weight": predictions["prediction_weight"][j],
            "assignment_ambiguity": int(ambiguity[i]), "assignment_margin_sigma": margin,
            "dwba_weight": predictions.get("dwba_weight", np.ones(len(predictions["qr"])))[j],
            "effective_intensity": predictions.get("effective_intensity", predictions["f2"])[j],
            "orientation_domain": str(
                predictions.get("domain", np.full(len(predictions["qr"]), "primary", dtype=object))[j]),
            "assignment_support_score": assignment_support_score,
            "assignment_support_is_probability": False,
        })
    matches = pd.DataFrame(rows).sort_values("assignment_support_score", ascending=False)
    return matches, metrics


def assign_reflections(features, predictions, config: IndexingConfig, tolerance):
    return _assign_arrays(features, _predictions_dict(predictions), config, tolerance, materialize=True)


def _parameter_bounds(config):
    radial_change = config.max_qr_scale_change if config.enable_anisotropic_q_scale else config.max_common_scale_change
    vertical_change = config.max_qz_scale_change if config.enable_anisotropic_q_scale else config.max_common_scale_change
    return [
        (-config.max_tilt_anchor_deg, config.max_tilt_anchor_deg),
        (-config.max_tilt_anchor_deg, config.max_tilt_anchor_deg),
        (1 - radial_change, 1 + radial_change),
        (1 - vertical_change, 1 + vertical_change),
        (-config.max_anchor_q_offset, config.max_anchor_q_offset),
        (-config.max_anchor_q_offset, config.max_anchor_q_offset),
    ]


def _prediction_pattern_distance(solution_a, solution_b, maximum=80):
    a = solution_a.get("predictions", pd.DataFrame())
    b = solution_b.get("predictions", pd.DataFrame())
    if a.empty or b.empty:
        return np.nan
    a = a.nlargest(min(maximum, len(a)), "effective_intensity" if "effective_intensity" in a else "f2")
    b = b.nlargest(min(maximum, len(b)), "effective_intensity" if "effective_intensity" in b else "f2")
    tree_a = cKDTree(a[["qr", "qz"]].to_numpy(float))
    tree_b = cKDTree(b[["qr", "qz"]].to_numpy(float))
    d_ab = tree_b.query(a[["qr", "qz"]].to_numpy(float))[0]
    d_ba = tree_a.query(b[["qr", "qz"]].to_numpy(float))[0]
    return float(np.median(np.concatenate((d_ab, d_ba))))


def _assignment_jaccard(solution_a, solution_b):
    def keys(solution):
        frames = [solution.get("anchor_matches", pd.DataFrame()), solution.get("validation_matches", pd.DataFrame())]
        frame = pd.concat([item for item in frames if not item.empty], ignore_index=True) if any(
            not item.empty for item in frames) else pd.DataFrame()
        return set(zip(frame.feature_id, frame.h, frame.k, frame.l)) if not frame.empty else set()

    a, b = keys(solution_a), keys(solution_b)
    return len(a & b) / max(len(a | b), 1)


def orientation_family_analysis(solutions, crystal, config):
    """Keep close alternatives and report whether the data can resolve them."""
    unique = {}
    for solution in sorted(solutions, key=lambda item: item["cv_score"], reverse=True):
        key = tuple(solution["hkl"])
        if key not in unique or solution["cv_score"] > unique[key]["cv_score"]:
            unique[key] = solution
    ordered = list(unique.values())
    families = []
    for solution in ordered:
        assigned = False
        for family in families:
            representative = family[0]
            angle = normal_angle(solution["hkl"], representative["hkl"], crystal)
            pattern = _prediction_pattern_distance(solution, representative)
            if (
                    angle <= config.orientation_family_merge_angle_deg
                    and np.isfinite(pattern)
                    and pattern <= config.orientation_pattern_merge_q
            ):
                family.append(solution)
                assigned = True
                break
        if not assigned:
            families.append([solution])
    family_records = [
        (max(family, key=lambda item: item["cv_score"]), family)
        for family in families
    ]
    family_records.sort(key=lambda pair: pair[0]["cv_score"], reverse=True)
    representatives = [pair[0] for pair in family_records]
    best = representatives[0]
    rows = []
    for family_id, (representative, members) in enumerate(family_records, 1):
        delta = float(best["cv_score"] - representative["cv_score"])
        angle = 0.0 if family_id == 1 else normal_angle(best["hkl"], representative["hkl"], crystal)
        pattern = 0.0 if family_id == 1 else _prediction_pattern_distance(best, representative)
        jaccard = 1.0 if family_id == 1 else _assignment_jaccard(best, representative)
        unresolved = bool(
            family_id > 1 and (
                    angle < config.orientation_min_unique_separation_deg
                    or (np.isfinite(pattern) and pattern <= config.orientation_pattern_merge_q)
            )
        )
        rows.append({
            "family_id": family_id,
            "normal_h": representative["hkl"][0],
            "normal_k": representative["hkl"][1],
            "normal_l": representative["hkl"][2],
            "heuristic_cross_validated_score": representative["cv_score"],
            "heuristic_score_delta_from_best": delta,
            "normal_angle_from_best_deg": angle,
            "projected_pattern_distance_from_best_q": pattern,
            "assignment_jaccard_with_best": jaccard,
            "member_count": len(members),
            "in_ambiguity_set": bool(delta <= config.orientation_ambiguity_score_delta or unresolved),
            "unresolved_from_best_at_configured_resolution": unresolved,
            "support_is_statistical_probability": False,
        })
    return pd.DataFrame(rows), representatives


def external_truth_check(series_id: str, best_hkl, crystal, config: IndexingConfig):
    truth = {row[0]: tuple(map(int, row[1:4])) for row in config.external_truth_normals}
    if series_id not in truth:
        return {
            "external_truth_status": "not_provided",
            "external_truth_angle_error_deg": np.nan,
            "external_truth_pass": False,
            "external_truth_source": config.external_truth_source,
            "external_truth_was_blind": bool(config.external_truth_was_blind),
        }
    angle = normal_angle(best_hkl, truth[series_id], crystal)
    return {
        "external_truth_status": "tested",
        "external_truth_normal_hkl": truth[series_id],
        "external_truth_angle_error_deg": angle,
        "external_truth_pass": bool(angle <= config.external_truth_tolerance_deg),
        "external_truth_source": config.external_truth_source,
        "external_truth_was_blind": bool(config.external_truth_was_blind),
    }


def orientation_conclusion(search, leave_one_out, bootstrap_summary, config, repeat_score=np.nan):
    families = search.get("orientation_families", pd.DataFrame())
    ambiguity_count = int(families.in_ambiguity_set.sum()) if not families.empty else 0
    score_gap = search.get("score_margin", np.nan)
    loo = float(leave_one_out.predicted_weighted_fraction.mean()) if not leave_one_out.empty else np.nan
    stability = bootstrap_summary.get("orientation_stability", np.nan)
    completed = int(bootstrap_summary.get("completed_iterations", 0))
    criteria = {
        "one_reported_family_in_ambiguity_set": ambiguity_count == 1,
        "heuristic_score_gap_sufficient": bool(
            np.isfinite(score_gap) and score_gap >= config.orientation_min_unique_score_gap),
        "held_angle_prediction_sufficient": bool(
            np.isfinite(loo) and loo >= config.orientation_min_loo_weighted_fraction),
        "bootstrap_iterations_sufficient": completed >= config.orientation_min_bootstrap_iterations,
        "bootstrap_stability_sufficient": bool(
            np.isfinite(stability) and stability >= config.orientation_min_stability_fraction),
        "repeat_not_contradictory": bool(not np.isfinite(repeat_score) or repeat_score >= 0.40),
    }
    if all(criteria.values()):
        label = "single_orientation_family_supported_not_proven_unique"
    elif ambiguity_count > 1:
        label = "ambiguous_orientation_family_set"
    else:
        label = "insufficient_evidence_for_unique_orientation"
    return label, criteria


def _existing_assignment_keys(search, secondary):
    if secondary and secondary.get("accepted") and secondary.get("joint") is not None:
        tables = [
            secondary["joint"].get("anchor_matches", pd.DataFrame()),
            secondary["joint"].get("validation_matches", pd.DataFrame()),
        ]
        existing = pd.concat([table for table in tables if not table.empty], ignore_index=True) \
            if any(not table.empty for table in tables) else pd.DataFrame()
    else:
        existing = combined_matches(search)
    keys = set()
    for row in existing.itertuples(index=False):
        domain = str(getattr(row, "orientation_domain", "primary"))
        sign = 0 if abs(float(row.qr_calc)) < 1e-10 else int(math.copysign(1, float(row.qr_calc)))
        keys.add((domain, int(row.h), int(row.k), int(row.l), sign))
    return keys


def _drop_used_predictions(predictions, used_keys):
    predictions = _predictions_dict(predictions)
    if not len(predictions["qr"]) or not used_keys:
        return predictions
    keep = []
    for i in range(len(predictions["qr"])):
        sign = 0 if abs(float(predictions["qr"][i])) < 1e-10 else int(math.copysign(1, float(predictions["qr"][i])))
        key = (
            str(predictions["domain"][i]), int(predictions["h"][i]),
            int(predictions["k"][i]), int(predictions["l"][i]), sign,
        )
        keep.append(key not in used_keys)
    keep = np.asarray(keep, dtype=bool)
    return {key: np.asarray(value)[keep] for key, value in predictions.items()}


def _single_prediction_distance(feature_frame, qr_calc, qz_calc, config):
    if feature_frame.empty or not np.isfinite(qr_calc + qz_calc):
        return np.empty(0), np.empty(0)
    physical, normalized = _feature_distance_matrices(
        feature_frame,
        {"qr": np.array([float(qr_calc)]), "qz": np.array([float(qz_calc)])},
        config,
    )
    return physical[:, 0], normalized[:, 0]


def _member_corroboration(assignments, members, config):
    columns = [
        "feature_id", "member_observations", "member_angle_support",
        "member_gate_fraction", "member_median_delta_q", "member_median_normalized_delta",
    ]
    if assignments.empty or members is None or members.empty or "feature_id" not in members:
        return pd.DataFrame(columns=columns)
    rows = []
    for assignment in assignments.itertuples(index=False):
        group = members[members.feature_id == assignment.feature_id].copy()
        physical, normalized = _single_prediction_distance(
            group, assignment.qr_calc, assignment.qz_calc, config
        )
        gate = (
                (physical <= config.ignored_evidence_member_gate_q)
                & (normalized <= config.ignored_evidence_member_sigma_limit)
        ) if len(physical) else np.zeros(0, bool)
        rows.append({
            "feature_id": assignment.feature_id,
            "member_observations": int(len(group)),
            "member_angle_support": int(group.loc[gate, "angle_deg"].nunique()) if len(group) else 0,
            "member_gate_fraction": float(gate.mean()) if len(gate) else np.nan,
            "member_median_delta_q": float(np.median(physical)) if len(physical) else np.nan,
            "member_median_normalized_delta": float(np.median(normalized)) if len(normalized) else np.nan,
        })
    return pd.DataFrame(rows, columns=columns)


def _parameter_perturbations(params, config):
    base = np.asarray(_unpack(params), float)
    steps = np.array([
        config.ignored_evidence_tilt_jitter_deg,
        config.ignored_evidence_tilt_jitter_deg,
        config.ignored_evidence_scale_jitter,
        config.ignored_evidence_scale_jitter,
        config.ignored_evidence_offset_jitter_q,
        config.ignored_evidence_offset_jitter_q,
    ], float)
    trials = [base.copy()]
    for index in range(6):
        for direction in (-1.0, 1.0):
            candidate = base.copy();
            candidate[index] += direction * steps[index]
            trials.append(candidate)
    requested = max(int(config.ignored_evidence_perturbation_trials), 1)
    if requested > len(trials):
        rng = np.random.default_rng(config.random_seed + 1907)
        while len(trials) < requested:
            trials.append(base + rng.normal(0.0, steps / 1.5))
    return trials[:requested]


def _projection_perturbation_stability(assignments, ignored, crystal, search, secondary, config):
    columns = ["feature_id", "projection_stability_fraction", "projection_stability_trials"]
    if assignments.empty:
        return pd.DataFrame(columns=columns)
    feature_lookup = ignored.set_index("feature_id", drop=False)
    rows = []
    for assignment in assignments.itertuples(index=False):
        solution = _domain_orientation_solution(search, secondary, assignment.orientation_domain)
        if solution is None or assignment.feature_id not in feature_lookup.index:
            rows.append({"feature_id": assignment.feature_id,
                         "projection_stability_fraction": np.nan,
                         "projection_stability_trials": 0})
            continue
        feature = feature_lookup.loc[[assignment.feature_id]]
        stable = []
        sign = 1 if float(assignment.qr_calc) >= 0 else -1
        for params in _parameter_perturbations(solution["params"], config):
            qr_calc, qz_calc = _project_exact_reflection(
                crystal, solution["hkl"], config, params,
                assignment.h, assignment.k, assignment.l, sign,
            )
            physical, normalized = _single_prediction_distance(feature, qr_calc, qz_calc, config)
            stable.append(bool(
                len(physical)
                and physical[0] <= config.ignored_index_tolerance_q
                and normalized[0] <= config.ignored_index_sigma_limit
            ))
        rows.append({
            "feature_id": assignment.feature_id,
            "projection_stability_fraction": float(np.mean(stable)) if stable else np.nan,
            "projection_stability_trials": int(len(stable)),
        })
    return pd.DataFrame(rows, columns=columns)


def _orientation_decoy_specificity(assignments, ignored, crystal, search, secondary, config):
    columns = ["feature_id", "decoy_orientation_win_fraction", "decoy_orientation_trials"]
    requested = int(config.ignored_evidence_decoy_trials)
    if assignments.empty or requested <= 0:
        return pd.DataFrame(columns=columns)
    accepted_normals = [search["best"]["hkl"]]
    if secondary and secondary.get("accepted") and secondary.get("search"):
        accepted_normals.append(secondary["search"]["best"]["hkl"])
    candidates = orientation_candidates(
        crystal, replace(config, max_orientation_candidates=max(50, requested * 5))
    )
    candidates = [
        candidate for candidate in candidates
        if all(normal_angle(candidate, normal, crystal) >= config.ignored_evidence_decoy_min_angle_deg
               for normal in accepted_normals)
    ]
    if not candidates:
        return pd.DataFrame(columns=columns)
    rng = np.random.default_rng(config.random_seed + 2711)
    if len(candidates) > requested:
        indices = rng.choice(len(candidates), size=requested, replace=False)
        candidates = [candidates[int(index)] for index in indices]
    feature_lookup = ignored.set_index("feature_id", drop=False)
    domain_predictions = {}
    for domain in sorted(set(assignments.orientation_domain.astype(str))):
        solution = _domain_orientation_solution(search, secondary, domain)
        if solution is None:
            continue
        domain_predictions[domain] = [
            _prediction_array(
                crystal, candidate, config, solution["params"],
                config.ignored_index_f2_percentile, config.ignored_index_max_predictions,
            )
            for candidate in candidates
        ]
    rows = []
    for assignment in assignments.itertuples(index=False):
        feature_id = assignment.feature_id
        if feature_id not in feature_lookup.index:
            rows.append({"feature_id": feature_id, "decoy_orientation_win_fraction": np.nan,
                         "decoy_orientation_trials": 0})
            continue
        feature = feature_lookup.loc[[feature_id]]
        wins = []
        for predictions in domain_predictions.get(str(assignment.orientation_domain), []):
            if not len(predictions["qr"]):
                continue
            _, normalized = _feature_distance_matrices(feature, predictions, config)
            decoy_best = float(np.min(normalized[0]))
            wins.append(decoy_best >= float(assignment.normalized_delta) + config.ignored_evidence_decoy_margin_sigma)
        rows.append({
            "feature_id": feature_id,
            "decoy_orientation_win_fraction": float(np.mean(wins)) if wins else np.nan,
            "decoy_orientation_trials": int(len(wins)),
        })
    return pd.DataFrame(rows, columns=columns)


def _grade_ignored_assignments(assignments, ignored, members, crystal, search, secondary, config):
    if assignments.empty:
        return assignments
    result = assignments.copy()
    extra = ignored[[
        column for column in ("feature_id", "support_fraction", "major_width_q", "minor_width_q")
        if column in ignored
    ]].drop_duplicates("feature_id")
    result = result.merge(extra, on="feature_id", how="left")
    for evidence in (
            _member_corroboration(result, members, config),
            _projection_perturbation_stability(result, ignored, crystal, search, secondary, config),
            _orientation_decoy_specificity(result, ignored, crystal, search, secondary, config),
    ):
        if not evidence.empty:
            result = result.merge(evidence, on="feature_id", how="left")
    for column in (
            "member_gate_fraction", "projection_stability_fraction", "decoy_orientation_win_fraction"
    ):
        if column not in result:
            result[column] = np.nan
    result["member_evidence_pass"] = (
            result.member_gate_fraction >= config.ignored_evidence_min_member_fraction
    )
    result["projection_stability_pass"] = (
            result.projection_stability_fraction >= config.ignored_evidence_min_projection_stability
    )
    result["decoy_specificity_pass"] = (
            result.decoy_orientation_win_fraction >= config.ignored_evidence_min_decoy_win_fraction
    )
    evidence_values = result[[
        "member_gate_fraction", "projection_stability_fraction", "decoy_orientation_win_fraction"
    ]].to_numpy(float)
    evaluated = np.isfinite(evidence_values)
    passed = np.column_stack((
        result.member_evidence_pass.to_numpy(bool),
        result.projection_stability_pass.to_numpy(bool),
        result.decoy_specificity_pass.to_numpy(bool),
    ))
    result["evidence_checks_evaluated"] = evaluated.sum(axis=1)
    result["evidence_checks_passed"] = (passed & evaluated).sum(axis=1)
    geometric = np.exp(-0.5 * np.square(result.normalized_delta.to_numpy(float)))
    margin = np.clip(result.assignment_margin_sigma.to_numpy(float) / 1.5, 0.0, 1.0)
    support = np.clip(result.support.to_numpy(float) / 4.0, 0.0, 1.0)
    evidence_mean = np.nanmean(np.where(evaluated, evidence_values, np.nan), axis=1)
    evidence_mean = np.where(np.isfinite(evidence_mean), evidence_mean, 0.0)
    result["salvage_evidence_score"] = np.clip(
        0.30 * geometric + 0.15 * margin + 0.15 * support + 0.40 * evidence_mean,
        0.0, 1.0,
    )
    robust = (
            (result.evidence_checks_evaluated >= 3)
            & (result.member_gate_fraction >= 0.75)
            & (result.projection_stability_fraction >= 0.75)
            & (result.decoy_orientation_win_fraction >= 0.75)
            & (result.support >= 3)
            & (result.assignment_ambiguity <= 1)
            & (result.assignment_margin_sigma >= 0.45)
    )
    supported = (
            ~robust
            & (result.evidence_checks_evaluated >= 2)
            & (result.evidence_checks_passed == result.evidence_checks_evaluated)
            & (result.assignment_ambiguity <= config.ignored_index_max_ambiguity)
    )
    result["salvage_evidence_tier"] = np.select(
        [robust, supported], ["robust", "supported"], default="provisional"
    )
    result["salvage_status"] = result.salvage_evidence_tier.map({
        "robust": "robust_index", "supported": "supported_index",
        "provisional": "provisional_index",
    })
    return result.sort_values(
        ["salvage_evidence_tier", "salvage_evidence_score"],
        ascending=[True, False], kind="mergesort",
    ).reset_index(drop=True)


def _assignment_prediction_keys(assignments):
    keys = set()
    if assignments is None or assignments.empty:
        return keys
    for row in assignments.itertuples(index=False):
        sign = 0 if abs(float(row.qr_calc)) < 1e-10 else int(math.copysign(1, float(row.qr_calc)))
        keys.add((str(getattr(row, "orientation_domain", "primary")),
                  int(row.h), int(row.k), int(row.l), sign))
    return keys


def prediction_guided_raw_rescue(crystal, raw_features, consensus_members, search, secondary,
                                 ignored_assignments, config):
    """Recover recurring raw detections omitted from consensus using a fixed orientation."""
    assignment_columns = [
        "feature_id", "hkl", "h", "k", "l", "qr_exp", "qz_exp", "qr_calc", "qz_calc",
        "delta_qr", "delta_qz", "delta_q", "normalized_delta", "strength", "support",
        "support_fraction", "feature_type", "f2", "prediction_weight", "assignment_ambiguity",
        "assignment_margin_sigma", "orientation_domain", "assignment_support_score", "role",
        "salvage_status", "salvage_evidence_tier", "raw_member_ids", "angles",
        "salvage_is_primary_orientation_evidence", "salvage_is_statistical_confidence",
    ]
    if not config.guided_rescue_unclustered_features or raw_features is None or raw_features.empty:
        return pd.DataFrame(columns=assignment_columns), pd.DataFrame()
    raw = raw_features.copy().reset_index(drop=True)
    if "raw_feature_id" not in raw:
        raw["raw_feature_id"] = [f"R{i + 1:05d}" for i in range(len(raw))]
    member_ids = set()
    if consensus_members is not None and not consensus_members.empty and "raw_feature_id" in consensus_members:
        member_ids = set(consensus_members.raw_feature_id.astype(str))
    pool = raw[~raw.raw_feature_id.astype(str).isin(member_ids)].copy()
    if pool.empty:
        return pd.DataFrame(columns=assignment_columns), pd.DataFrame()
    sigma = np.sqrt(np.maximum(pool.cov_rr + pool.cov_zz, 1e-12))
    quality = (
            (pool.qz >= config.analysis_qz_min)
            & (sigma <= config.ignored_index_max_sigma_q)
            & (pool.major_width_q <= config.guided_rescue_max_major_width_q)
    )
    pool = pool[quality].copy()
    if pool.empty:
        return pd.DataFrame(columns=assignment_columns), pd.DataFrame()
    primary = _prediction_array(
        crystal, search["best"]["hkl"], config, search["best"]["params"],
        config.ignored_index_f2_percentile, config.ignored_index_max_predictions,
    )
    predictions = primary
    if secondary and secondary.get("accepted") and secondary.get("search"):
        second = secondary["search"]["best"]
        predictions = _combine_domain_predictions(primary, _prediction_array(
            crystal, second["hkl"], config, second["params"],
            config.ignored_index_f2_percentile, config.ignored_index_max_predictions,
        ))
    used = _existing_assignment_keys(search, secondary) | _assignment_prediction_keys(ignored_assignments)
    predictions = _drop_used_predictions(predictions, used)
    if not len(predictions["qr"]):
        return pd.DataFrame(columns=assignment_columns), pd.DataFrame()
    local = replace(config, match_sigma_limit=config.guided_rescue_sigma_limit)
    physical, normalized = _feature_distance_matrices(pool, predictions, local)
    gate = (
            (physical <= config.guided_rescue_tolerance_q)
            & (normalized <= config.guided_rescue_sigma_limit)
    )
    ambiguity = gate.sum(axis=1)
    sorted_norm = np.sort(normalized, axis=1)
    margin = np.where(
        sorted_norm.shape[1] > 1,
        sorted_norm[:, 1] - sorted_norm[:, 0],
        config.guided_rescue_sigma_limit,
    )
    best = np.argmin(normalized, axis=1)
    valid = (
            (ambiguity >= 1)
            & (ambiguity <= config.guided_rescue_max_ambiguity)
            & (margin >= config.guided_rescue_min_margin_sigma)
            & gate[np.arange(len(pool)), best]
    )
    candidates = pool.loc[valid].copy()
    if candidates.empty:
        diagnostics = pool[["raw_feature_id", "qr", "qz", "angle_deg", "feature_type"]].copy()
        diagnostics["guided_rescue_status"] = "no_unique_prediction_within_gate"
        return pd.DataFrame(columns=assignment_columns), diagnostics
    idx = np.flatnonzero(valid)
    chosen = best[valid]
    candidates["prediction_index"] = chosen
    candidates["h"] = predictions["h"][chosen]
    candidates["k"] = predictions["k"][chosen]
    candidates["l"] = predictions["l"][chosen]
    candidates["hkl"] = predictions["hkl"][chosen]
    candidates["qr_calc"] = predictions["qr"][chosen]
    candidates["qz_calc"] = predictions["qz"][chosen]
    candidates["orientation_domain"] = predictions["domain"][chosen]
    candidates["f2"] = predictions["f2"][chosen]
    candidates["prediction_weight"] = predictions["prediction_weight"][chosen]
    candidates["normalized_delta"] = normalized[idx, chosen]
    candidates["delta_q"] = physical[idx, chosen]
    candidates["assignment_ambiguity"] = ambiguity[valid]
    candidates["assignment_margin_sigma"] = margin[valid]
    candidates["qr_sign"] = np.sign(candidates.qr_calc).astype(int)
    candidates = candidates.sort_values(
        ["normalized_delta", "strength"], ascending=[True, False], kind="mergesort"
    )
    candidates = candidates.drop_duplicates(
        ["orientation_domain", "h", "k", "l", "qr_sign", "angle_deg"], keep="first"
    )
    total_angles = max(int(raw.angle_deg.nunique()), 1)
    rows = []
    group_columns = ["orientation_domain", "h", "k", "l", "qr_sign"]
    for key, group in candidates.groupby(group_columns, sort=False):
        support = int(group.angle_deg.nunique())
        if support < config.guided_rescue_min_angle_support:
            continue
        weight = np.maximum(group.strength.to_numpy(float), 1e-6)
        weight /= weight.sum()
        qr_exp = float(np.sum(weight * group.qr.to_numpy(float)))
        qz_exp = float(np.sum(weight * group.qz.to_numpy(float)))
        qr_calc = float(group.qr_calc.iloc[0]);
        qz_calc = float(group.qz_calc.iloc[0])
        delta_q = float(math.hypot(qr_exp - qr_calc, qz_exp - qz_calc))
        median_norm = float(np.median(group.normalized_delta))
        min_margin = float(group.assignment_margin_sigma.min())
        prediction_weight = float(group.prediction_weight.iloc[0])
        assignment_support_score = (
                math.exp(-0.5 * median_norm ** 2) * min(1.0, support / 3.0)
                * prediction_weight * min(1.0, min_margin / 1.5)
        )
        dominant_type = str(group.groupby("feature_type").strength.sum().idxmax())
        maximum_strength = float(group.strength.max())
        strong_multi_angle = (
                support >= int(config.guided_rescue_supported_min_angle_support)
                and median_norm <= float(config.guided_rescue_supported_max_normalized_delta)
                and min_margin >= float(config.guided_rescue_supported_min_margin_sigma)
        )
        strong_two_angle = (
                bool(config.guided_rescue_promote_strong_two_angle)
                and support == 2
                and dominant_type not in set(config.guided_rescue_two_angle_excluded_types)
                and median_norm <= float(config.guided_rescue_two_angle_max_normalized_delta)
                and delta_q <= float(config.guided_rescue_two_angle_max_delta_q)
                and min_margin >= float(config.guided_rescue_two_angle_min_margin_sigma)
                and maximum_strength >= float(config.guided_rescue_two_angle_min_strength)
        )
        tier = "supported" if (strong_multi_angle or strong_two_angle) else "provisional"
        feature_id = f"G{len(rows) + 1:03d}"
        rows.append({
            "feature_id": feature_id, "hkl": str(group.hkl.iloc[0]),
            "h": int(key[1]), "k": int(key[2]), "l": int(key[3]),
            "qr_exp": qr_exp, "qz_exp": qz_exp, "qr_calc": qr_calc, "qz_calc": qz_calc,
            "delta_qr": qr_exp - qr_calc, "delta_qz": qz_exp - qz_calc,
            "delta_q": delta_q, "normalized_delta": median_norm,
            "strength": maximum_strength, "support": support,
            "support_fraction": support / total_angles,
            "feature_type": dominant_type,
            "f2": float(group.f2.iloc[0]), "prediction_weight": prediction_weight,
            "assignment_ambiguity": int(group.assignment_ambiguity.max()),
            "assignment_margin_sigma": min_margin,
            "orientation_domain": str(key[0]),
            "assignment_support_score": assignment_support_score,
            "role": "prediction_guided_raw_rescue",
            "salvage_status": f"{tier}_guided_rescue",
            "salvage_evidence_tier": tier,
            "raw_member_ids": ",".join(group.raw_feature_id.astype(str)),
            "angles": ",".join(f"{value:.3f}" for value in sorted(group.angle_deg.unique())),
            "salvage_is_primary_orientation_evidence": False,
            "salvage_is_statistical_confidence": False,
        })
    assignments = pd.DataFrame(rows, columns=assignment_columns)
    candidate_ids = set()
    if not assignments.empty:
        for value in assignments.raw_member_ids:
            candidate_ids.update(str(value).split(","))
    diagnostics = pool[[
        "raw_feature_id", "qr", "qz", "angle_deg", "feature_type", "strength", "snr"
    ]].copy()
    diagnostics["guided_rescue_status"] = np.where(
        diagnostics.raw_feature_id.astype(str).isin(candidate_ids),
        "used_in_multi_angle_guided_rescue", "not_recovered",
    )
    return assignments.reset_index(drop=True), diagnostics.reset_index(drop=True)


def index_ignored_features(crystal, search, secondary, config: IndexingConfig, members=None):
    """Provisionally index unused consensus peaks after orientation selection.

    The recovered assignments are deliberately excluded from orientation fitting,
    score calculation, bootstrap stability, and uniqueness claims.  This avoids
    circularly improving the solution with the same weak peaks it is attempting
    to explain.
    """
    ignored = search.get("ignored", pd.DataFrame()).copy()
    diagnostic_columns = [
        "feature_id", "qr", "qz", "feature_type", "support", "strength",
        "ignored_reason", "salvage_status", "nearest_hkl", "nearest_domain",
        "nearest_delta_q", "nearest_normalized_delta", "candidate_ambiguity",
        "candidate_margin_sigma",
    ]
    if not config.index_ignored_features or ignored.empty:
        return pd.DataFrame(), pd.DataFrame(columns=diagnostic_columns)

    ignored = ignored[~ignored.feature_id.isin(config.manual_rejected_feature_ids)].copy()
    if ignored.empty:
        return pd.DataFrame(), pd.DataFrame(columns=diagnostic_columns)

    primary_predictions = _prediction_array(
        crystal, search["best"]["hkl"], config, search["best"]["params"],
        config.ignored_index_f2_percentile, config.ignored_index_max_predictions,
    )
    predictions = primary_predictions
    if secondary and secondary.get("accepted") and secondary.get("search"):
        second_best = secondary["search"]["best"]
        secondary_predictions = _prediction_array(
            crystal, second_best["hkl"], config, second_best["params"],
            config.ignored_index_f2_percentile, config.ignored_index_max_predictions,
        )
        predictions = _combine_domain_predictions(primary_predictions, secondary_predictions)
    predictions = _drop_used_predictions(
        predictions, _existing_assignment_keys(search, secondary)
    )
    if not len(predictions["qr"]):
        diagnostics = ignored.copy()
        diagnostics["salvage_status"] = "no_unused_predictions"
        for column in diagnostic_columns:
            if column not in diagnostics:
                diagnostics[column] = np.nan
        return pd.DataFrame(), diagnostics[diagnostic_columns]

    sigma = np.sqrt(np.maximum(ignored.cov_rr + ignored.cov_zz, 1e-12))
    quality = (
            (ignored.support >= config.ignored_index_min_support)
            & (sigma <= config.ignored_index_max_sigma_q)
            & (ignored.major_width_q <= config.ignored_index_max_major_width_q)
    )
    if not config.ignored_index_include_radial_streaks:
        quality &= ~ignored.feature_type.isin(config.ignored_feature_types)
    else:
        streak = ignored.feature_type.isin(config.ignored_feature_types)
        quality &= (~streak) | (ignored.support >= config.ignored_streak_min_support)
    pool = ignored[quality].copy()

    local = replace(config, match_sigma_limit=config.ignored_index_sigma_limit)
    raw_matches, _ = _assign_arrays(
        pool, predictions, local, config.ignored_index_tolerance_q, materialize=True
    ) if not pool.empty else (pd.DataFrame(), _empty_assignment_metrics())

    accepted = raw_matches.copy()
    if not accepted.empty:
        normal_quality = (
                (accepted.assignment_ambiguity <= config.ignored_index_max_ambiguity)
                & (accepted.assignment_margin_sigma >= config.ignored_index_min_margin_sigma)
                & (accepted.assignment_support_score >= config.ignored_index_min_support_score)
        )
        streak_quality = (
                (accepted.feature_type.isin(config.ignored_feature_types))
                & (accepted.support >= config.ignored_streak_min_support)
                & (accepted.normalized_delta <= config.ignored_streak_sigma_limit)
                & (accepted.assignment_ambiguity <= config.ignored_streak_max_ambiguity)
        )
        accepted = accepted[normal_quality & (
                ~accepted.feature_type.isin(config.ignored_feature_types) | streak_quality
        )].copy()
        reason_lookup = ignored.set_index("feature_id").get(
            "ignored_reason", pd.Series(dtype=object)
        )
        accepted["ignored_reason"] = accepted.feature_id.map(reason_lookup).fillna("unused_consensus_feature")
        accepted["role"] = "ignored_salvage"
        accepted["salvage_status"] = "provisional_index"
        accepted["salvage_is_primary_orientation_evidence"] = False
        accepted["salvage_is_statistical_confidence"] = False
        accepted = _grade_ignored_assignments(
            accepted, ignored, members, crystal, search, secondary, config
        )

    physical, normalized = _feature_distance_matrices(ignored, predictions, local)
    diagnostics = ignored.copy()
    if len(predictions["qr"]):
        best_index = np.argmin(normalized, axis=1)
        sorted_norm = np.sort(normalized, axis=1)
        best_norm = normalized[np.arange(len(ignored)), best_index]
        best_phys = physical[np.arange(len(ignored)), best_index]
        gate = (
                (physical <= config.ignored_index_tolerance_q)
                & (normalized <= config.ignored_index_sigma_limit)
        )
        ambiguity = gate.sum(axis=1)
        margin = np.where(
            sorted_norm.shape[1] > 1,
            sorted_norm[:, 1] - sorted_norm[:, 0],
            config.ignored_index_sigma_limit,
        )
        diagnostics["nearest_hkl"] = [predictions["hkl"][j] for j in best_index]
        diagnostics["nearest_domain"] = [predictions["domain"][j] for j in best_index]
        diagnostics["nearest_delta_q"] = best_phys
        diagnostics["nearest_normalized_delta"] = best_norm
        diagnostics["candidate_ambiguity"] = ambiguity
        diagnostics["candidate_margin_sigma"] = margin
    accepted_ids = set(accepted.feature_id) if not accepted.empty else set()
    raw_ids = set(raw_matches.feature_id) if not raw_matches.empty else set()
    statuses = []
    for row_index, row in diagnostics.iterrows():
        feature_id = row.feature_id
        if feature_id in accepted_ids:
            selected = accepted.loc[accepted.feature_id == feature_id, "salvage_status"]
            status = str(selected.iloc[0]) if len(selected) else "provisional_index"
        elif not bool(quality.loc[row_index]):
            status = "rejected_feature_quality"
        elif int(row.get("candidate_ambiguity", 0)) == 0:
            status = "no_prediction_within_gate"
        elif int(row.get("candidate_ambiguity", 0)) > config.ignored_index_max_ambiguity:
            status = "ambiguous_prediction"
        elif float(row.get("candidate_margin_sigma", 0.0)) < config.ignored_index_min_margin_sigma:
            status = "low_candidate_margin"
        elif feature_id in raw_ids:
            status = "rejected_assignment_support"
        else:
            status = "not_selected_one_to_one"
        statuses.append(status)
    diagnostics["salvage_status"] = statuses
    for column in diagnostic_columns:
        if column not in diagnostics:
            diagnostics[column] = np.nan
    return accepted.reset_index(drop=True), diagnostics[diagnostic_columns].reset_index(drop=True)


def _held_angle_features(members, angle):
    frame = members[members.angle_deg == angle].copy()
    if frame.empty:
        return frame
    frame["feature_id"] = [f"H{int(angle * 1000):03d}_{i + 1:03d}" for i in range(len(frame))]
    frame["support"], frame["support_fraction"], frame["angles"] = 1, 1.0, f"{angle:.3f}"
    return frame


def _v739_full_leave_one_angle_out(crystal, members, config):
    if not config.full_leave_one_angle_out:
        return pd.DataFrame()
    rows = []
    for held_angle in sorted(members.angle_deg.unique()):
        training = members[members.angle_deg != held_angle].copy()
        consensus, _ = build_consensus(training, config)
        if len(consensus) < config.min_anchor_matches:
            continue
        local = replace(
            config,
            max_normal_candidates_anchor=min(config.loo_search_normal_limit, config.max_normal_candidates_anchor),
            coarse_tilt_values_deg=config.coarse_tilt_values_deg,
            max_consensus_features=min(36, config.max_consensus_features),
            max_anchor_features=min(6, config.max_anchor_features),
            max_validation_features=min(18, config.max_validation_features),
            refine_hypotheses=2, expand_top_normals=2, max_pair_hypotheses_per_normal=1,
            max_hypothesis_predictions=min(30, config.max_hypothesis_predictions),
            max_validation_predictions=min(75, config.max_validation_predictions),
            full_leave_one_angle_out=False, full_bootstrap_iterations=0, test_second_orientation=False,
        )
        search = anchor_first_search(crystal, consensus, local, fast=True, validation_mode=True)
        held = _held_angle_features(members, held_angle)
        matches, metrics = assign_reflections(held, search["best"]["predictions"], local,
                                              local.validation_match_tolerance_q)
        rows.append({
            "held_angle_deg": held_angle, "trained_normal_hkl": str(search["best"]["hkl"]),
            "held_features": len(held), "predicted_matches": metrics["matches"],
            "predicted_weighted_fraction": metrics["weighted_fraction"],
            "predicted_median_delta_q": metrics["median_delta_q"], "training_score_margin": search["score_margin"],
        })
    return pd.DataFrame(rows)






def _release_working_memory():
    gc.collect()
    try:
        ctypes.CDLL("libc.so.6").malloc_trim(0)
    except Exception:
        pass


def _combine_domain_predictions(primary_predictions, secondary_predictions):
    first = _predictions_dict(primary_predictions)
    second = _predictions_dict(secondary_predictions)
    if not len(first["qr"]):
        combined = {key: np.asarray(value).copy() for key, value in second.items()}
        combined["domain"] = np.full(len(combined["qr"]), "secondary", dtype=object)
        return combined
    if not len(second["qr"]):
        combined = {key: np.asarray(value).copy() for key, value in first.items()}
        combined["domain"] = np.full(len(combined["qr"]), "primary", dtype=object)
        return combined
    combined = {}
    for key in first:
        left = np.asarray(first[key])
        right = np.asarray(second[key])
        combined[key] = np.concatenate((left, right))
    combined["domain"] = np.concatenate((
        np.full(len(first["qr"]), "primary", dtype=object),
        np.full(len(second["qr"]), "secondary", dtype=object),
    ))
    return combined


def repeat_scan_validation(series_results, crystal, config):
    rows = []
    grouped = {}
    for series_id, result in series_results.items():
        grouped.setdefault(series_id.split(":")[0], []).append((series_id, result))
    for sample, items in grouped.items():
        for i in range(len(items)):
            for j in range(i + 1, len(items)):
                name_a, result_a = items[i]
                name_b, result_b = items[j]
                normal_difference = normal_angle(result_a["search"]["best"]["hkl"], result_b["search"]["best"]["hkl"],
                                                 crystal)
                hkl_a = set(combined_matches(result_a["search"])["hkl"]) if not combined_matches(
                    result_a["search"]).empty else set()
                hkl_b = set(combined_matches(result_b["search"])["hkl"]) if not combined_matches(
                    result_b["search"]).empty else set()
                jaccard = len(hkl_a & hkl_b) / max(len(hkl_a | hkl_b), 1)
                features_a = pd.concat([result_a["search"]["anchors"], result_a["search"]["validation"]],
                                       ignore_index=True)
                features_b = pd.concat([result_b["search"]["anchors"], result_b["search"]["validation"]],
                                       ignore_index=True)
                _, cross_ab = assign_reflections(features_b, result_a["search"]["best"]["predictions"], config,
                                                 config.validation_match_tolerance_q)
                _, cross_ba = assign_reflections(features_a, result_b["search"]["best"]["predictions"], config,
                                                 config.validation_match_tolerance_q)
                score = float(np.mean([
                    math.exp(-normal_difference / 10.0), jaccard,
                    cross_ab["weighted_fraction"], cross_ba["weighted_fraction"],
                ]))
                rows.append({
                    "sample": sample, "series_a": name_a, "series_b": name_b,
                    "normal_angle_difference_deg": normal_difference, "indexed_hkl_jaccard": jaccard,
                    "a_predicts_b_weighted_fraction": cross_ab["weighted_fraction"],
                    "b_predicts_a_weighted_fraction": cross_ba["weighted_fraction"],
                    "repeat_agreement_score": score,
                })
    return pd.DataFrame(rows)


def _fit_at_boundary(params, config):
    return any(abs(value - low) < 0.02 * (high - low) or abs(value - high) < 0.02 * (high - low)
               for value, (low, high) in zip(params, _parameter_bounds(config)))


def solution_reliability(search, leave_one_out, bootstrap_summary, config, repeat_score=np.nan):
    flags, best = [], search["best"]
    if best["anchor_metrics"]["matches"] < config.min_anchor_matches:
        flags.append("too_few_anchor_matches")
    if best["validation_metrics"]["matches"] < config.min_validation_matches:
        flags.append("weak_data_only_holdout_validation")
    strict_metrics = best.get("strict_validation_metrics", {})
    if strict_metrics.get("matches", 0) < 1:
        flags.append("no_strict_tolerance_holdout_match")
    families = search.get("orientation_families", pd.DataFrame())
    if families.empty or int(families.in_ambiguity_set.sum()) != 1:
        flags.append("multiple_or_unresolved_orientation_families")
    if not np.isfinite(search.get("score_margin", np.nan)):
        flags.append("orientation_gap_unassessed")
    elif search["score_margin"] < config.orientation_min_unique_score_gap:
        flags.append("small_heuristic_score_gap")
    if leave_one_out.empty:
        flags.append("leave_one_angle_out_not_run")
    elif leave_one_out.predicted_weighted_fraction.mean() < config.orientation_min_loo_weighted_fraction:
        flags.append("weak_held_angle_prediction")
    stability = bootstrap_summary.get("orientation_stability", np.nan)
    completed = int(bootstrap_summary.get("completed_iterations", 0))
    if completed < config.orientation_min_bootstrap_iterations:
        flags.append("too_few_bootstrap_research_runs")
    elif not np.isfinite(stability) or stability < config.orientation_min_stability_fraction:
        flags.append("full_search_bootstrap_instability")
    if np.isfinite(repeat_score) and repeat_score < 0.40:
        flags.append("repeat_scan_disagreement")
    if _fit_at_boundary(best["params"], config):
        flags.append("parameter_at_boundary")
    if len(flags) >= 4:
        label = "weak_hypothesis"
    elif flags:
        label = "provisional_hypothesis"
    else:
        label = "well_supported_internal_hypothesis"
    return label, flags


def _v71_prediction_key(domain, h, k, l, qr):
    sign = 0 if abs(float(qr)) < 1e-10 else int(math.copysign(1, float(qr)))
    return str(domain), int(h), int(k), int(l), sign


def v71_fixed_orientation_completion(crystal, consensus, members, search, secondary, config):
    """Index residual consensus features using the full fixed-orientation lattice.

    The completion assignments do not participate in orientation ranking,
    calibration refinement, bootstrap stability, or validation counts. Existing
    core assignments are locked, and their reflection/sign keys are removed before
    residual one-to-one assignment. This prevents the completion pass from stealing
    or duplicating a primary assignment while restoring low-F2 reflections omitted
    by the fast orientation-search prediction list.
    """
    assignment_columns = [
        "feature_id", "hkl", "h", "k", "l", "qr_exp", "qz_exp", "qr_calc", "qz_calc",
        "delta_qr", "delta_qz", "delta_q", "normalized_delta", "strength", "support",
        "support_fraction", "feature_type", "feature_quality", "f2", "prediction_weight",
        "assignment_ambiguity", "assignment_margin_sigma", "orientation_domain",
        "assignment_support_score", "member_observations", "member_angle_support",
        "member_gate_fraction", "completion_tolerance_q", "role", "index_source",
        "completion_is_primary_orientation_evidence", "completion_is_validation_evidence",
    ]
    diagnostic_columns = [
        "feature_id", "completion_status", "nearest_delta_q", "nearest_normalized_delta",
        "candidate_count", "completion_tolerance_q", "support", "feature_type",
        "member_angle_support", "member_gate_fraction",
    ]
    if not getattr(config, "v71_enable_full_reflection_completion", True):
        return pd.DataFrame(columns=assignment_columns), pd.DataFrame(columns=diagnostic_columns)
    if consensus is None or consensus.empty:
        return pd.DataFrame(columns=assignment_columns), pd.DataFrame(columns=diagnostic_columns)

    if secondary and secondary.get("accepted") and secondary.get("joint") is not None:
        base_tables = [secondary["joint"].get("anchor_matches", pd.DataFrame()),
                       secondary["joint"].get("validation_matches", pd.DataFrame())]
        base = pd.concat([x for x in base_tables if x is not None and not x.empty],
                         ignore_index=True, sort=False) if any(
            x is not None and not x.empty for x in base_tables
        ) else pd.DataFrame()
        domain_solutions = secondary.get("domain_solutions", [])
    else:
        base = combined_matches(search)
        domain_solutions = [{
            "domain": "primary", "hkl": search["best"]["hkl"],
            "params": search["best"]["params"],
        }]
    if not domain_solutions:
        domain_solutions = [{
            "domain": "primary", "hkl": search["best"]["hkl"],
            "params": search["best"]["params"],
        }]

    used_feature_ids = set(base.feature_id.astype(str)) if not base.empty else set()
    residual = consensus[
        ~consensus.feature_id.astype(str).isin(used_feature_ids)
        & ~consensus.feature_id.astype(str).isin(config.manual_rejected_feature_ids)
        ].copy().reset_index(drop=True)
    if residual.empty:
        return pd.DataFrame(columns=assignment_columns), pd.DataFrame(columns=diagnostic_columns)

    sigma = np.sqrt(np.maximum(
        residual.cov_rr.to_numpy(float) + residual.cov_zz.to_numpy(float), 1e-12
    ))
    basic_quality = (
            (residual.support.to_numpy(int) >= int(config.v71_completion_min_support))
            & (sigma <= float(config.v71_completion_max_sigma_q))
            & (residual.major_width_q.to_numpy(float) <= float(config.v71_completion_max_major_width_q))
            & (residual.qz.to_numpy(float) >= float(config.analysis_qz_min))
    )
    residual = residual[basic_quality].copy().reset_index(drop=True)
    sigma = sigma[basic_quality]
    if residual.empty:
        return pd.DataFrame(columns=assignment_columns), pd.DataFrame(columns=diagnostic_columns)

    prediction_sets = []
    for item in domain_solutions:
        prediction_sets.append(_v7_prediction_array(
            crystal, tuple(item["hkl"]), config, np.asarray(item["params"], float),
            config.v71_completion_f2_percentile,
            config.v71_completion_max_predictions_per_domain,
            str(item.get("domain", "primary")),
        ))
    predictions = _v7_concat_predictions(prediction_sets)
    if not len(predictions["qr"]):
        return pd.DataFrame(columns=assignment_columns), pd.DataFrame(columns=diagnostic_columns)

    used_prediction_keys = set()
    if not base.empty:
        for row in base.itertuples(index=False):
            used_prediction_keys.add(_v71_prediction_key(
                getattr(row, "orientation_domain", "primary"), row.h, row.k, row.l, row.qr_calc
            ))
    keep = np.array([
        _v71_prediction_key(predictions["domain"][j], predictions["h"][j],
                            predictions["k"][j], predictions["l"][j], predictions["qr"][j])
        not in used_prediction_keys
        for j in range(len(predictions["qr"]))
    ], bool)
    predictions = {key: np.asarray(value)[keep] for key, value in predictions.items()}
    if not len(predictions["qr"]):
        return pd.DataFrame(columns=assignment_columns), pd.DataFrame(columns=diagnostic_columns)

    physical, normalized = _feature_distance_matrices(residual, predictions, config)
    quality, _ = _v7_feature_quality(residual, config)
    support = residual.support.to_numpy(float)
    spot = residual.feature_type.astype(str).eq("spot").to_numpy(float)
    adaptive_tolerance = np.clip(
        config.v71_completion_base_tolerance_q
        + 0.30 * np.minimum(sigma, 0.035)
        + 0.0025 * np.clip(support - 2.0, 0.0, 4.0)
        + 0.002 * spot,
        config.v71_completion_base_tolerance_q,
        config.v71_completion_max_tolerance_q,
    )
    gate = (
            (physical <= adaptive_tolerance[:, None])
            & (normalized <= float(config.v71_completion_sigma_limit))
    )
    ambiguity = gate.sum(axis=1)
    cost = (
            normalized / max(float(config.v71_completion_sigma_limit), 1e-9)
            + 0.025 * (1.0 - predictions["prediction_weight"])[None, :]
            + 0.050 * np.log1p(np.maximum(ambiguity - 1, 0))[:, None]
    )
    completion_locks = {
        str(feature_id): (int(h), int(k), int(l))
        for feature_id, h, k, l in config.manual_locked_assignments
    }
    for i, feature_id in enumerate(residual.feature_id.astype(str)):
        if feature_id in completion_locks:
            hkl = completion_locks[feature_id]
            allowed = (
                    (predictions["h"] == hkl[0])
                    & (predictions["k"] == hkl[1])
                    & (predictions["l"] == hkl[2])
            )
            cost[i, ~allowed] = 1e6
    real = np.nan_to_num(cost.copy(), nan=1e4, posinf=1e4, neginf=-1e4)
    real[~gate] = 1e4
    dummy_cost = 1.08 + 0.10 * quality
    dummy = np.full((len(residual), len(residual)), 1e4, float)
    np.fill_diagonal(dummy, dummy_cost)
    rows, cols = linear_sum_assignment(np.hstack((real, dummy)))
    selected = cols < len(predictions["qr"])
    rows, cols = rows[selected], cols[selected]
    selected = gate[rows, cols] & (cost[rows, cols] <= dummy_cost[rows])
    rows, cols = rows[selected], cols[selected]

    candidate_rows = []
    for i, j in zip(rows, cols):
        valid_norm = np.sort(normalized[i][gate[i]])
        margin = float(valid_norm[1] - valid_norm[0]) if len(valid_norm) > 1 else float(
            config.v71_completion_sigma_limit)
        candidate_rows.append({
            "feature_id": str(residual.iloc[i].feature_id),
            "hkl": str(predictions["hkl"][j]),
            "h": int(predictions["h"][j]), "k": int(predictions["k"][j]), "l": int(predictions["l"][j]),
            "qr_exp": float(residual.iloc[i].qr), "qz_exp": float(residual.iloc[i].qz),
            "qr_calc": float(predictions["qr"][j]), "qz_calc": float(predictions["qz"][j]),
            "delta_qr": float(residual.iloc[i].qr - predictions["qr"][j]),
            "delta_qz": float(residual.iloc[i].qz - predictions["qz"][j]),
            "delta_q": float(physical[i, j]), "normalized_delta": float(normalized[i, j]),
            "strength": float(residual.iloc[i].strength), "support": int(residual.iloc[i].support),
            "support_fraction": float(residual.iloc[i].support_fraction),
            "feature_type": str(residual.iloc[i].feature_type), "feature_quality": float(quality[i]),
            "f2": float(predictions["f2"][j]), "prediction_weight": float(predictions["prediction_weight"][j]),
            "assignment_ambiguity": int(ambiguity[i]), "assignment_margin_sigma": margin,
            "orientation_domain": str(predictions["domain"][j]),
            "assignment_support_score": float(
                quality[i] * math.exp(-0.5 * normalized[i, j] ** 2)
                * min(1.0, margin / 1.0) / math.sqrt(max(int(ambiguity[i]), 1))
            ),
            "completion_tolerance_q": float(adaptive_tolerance[i]),
        })
    candidates = pd.DataFrame(candidate_rows)
    if candidates.empty:
        nearest = np.min(physical, axis=1) if physical.size else np.full(len(residual), np.nan)
        nearest_norm = np.min(normalized, axis=1) if normalized.size else np.full(len(residual), np.nan)
        diagnostics = pd.DataFrame({
            "feature_id": residual.feature_id.astype(str), "completion_status": "no_assignment",
            "nearest_delta_q": nearest, "nearest_normalized_delta": nearest_norm,
            "candidate_count": ambiguity, "completion_tolerance_q": adaptive_tolerance,
            "support": residual.support.to_numpy(int), "feature_type": residual.feature_type.astype(str),
            "member_angle_support": 0, "member_gate_fraction": np.nan,
        })
        return pd.DataFrame(columns=assignment_columns), diagnostics

    evidence = _member_corroboration(candidates, members, config)
    candidates = candidates.merge(evidence, on="feature_id", how="left")
    for column, default in (("member_observations", 0), ("member_angle_support", 0),
                            ("member_gate_fraction", np.nan)):
        if column not in candidates:
            candidates[column] = default
    is_streak = candidates.feature_type.astype(str).isin(config.ignored_feature_types)
    general_ok = (
            (candidates.assignment_ambiguity <= int(config.v71_completion_max_ambiguity))
            & ((candidates.assignment_ambiguity <= 1)
               | (candidates.assignment_margin_sigma >= float(config.v71_completion_min_margin_sigma)))
            & (candidates.member_angle_support >= int(config.v71_completion_min_member_angles))
            & (candidates.member_gate_fraction >= float(config.v71_completion_min_member_fraction))
    )
    streak_ok = (
            ~is_streak
            | ((candidates.support >= int(config.v71_completion_streak_min_support))
               & (candidates.normalized_delta <= float(config.v71_completion_streak_sigma_limit))
               & (candidates.assignment_ambiguity <= int(config.v71_completion_streak_max_ambiguity)))
    )
    accepted = candidates[general_ok & streak_ok].copy()
    if not accepted.empty:
        accepted["role"] = "fixed_orientation_completion"
        accepted["index_source"] = "full_reflection_fixed_orientation_completion"
        accepted["completion_is_primary_orientation_evidence"] = False
        accepted["completion_is_validation_evidence"] = False
        accepted = accepted.sort_values("assignment_support_score", ascending=False, kind="mergesort")

    accepted_ids = set(accepted.feature_id.astype(str)) if not accepted.empty else set()
    nearest = np.min(physical, axis=1)
    nearest_norm = np.min(normalized, axis=1)
    diagnostics = pd.DataFrame({
        "feature_id": residual.feature_id.astype(str),
        "completion_status": np.where(residual.feature_id.astype(str).isin(accepted_ids),
                                      "indexed_full_reflection_completion", "not_completed"),
        "nearest_delta_q": nearest, "nearest_normalized_delta": nearest_norm,
        "candidate_count": ambiguity, "completion_tolerance_q": adaptive_tolerance,
        "support": residual.support.to_numpy(int), "feature_type": residual.feature_type.astype(str),
    })
    diagnostics = diagnostics.merge(
        candidates[["feature_id", "member_angle_support", "member_gate_fraction"]],
        on="feature_id", how="left"
    )
    return accepted.reindex(columns=assignment_columns), diagnostics.reindex(columns=diagnostic_columns)


def _v739_v72_s1_mosaic_completion(series_id, crystal, consensus, members, search, secondary,
                             existing_completion, config):
    """Recover split/mosaic s1 peaks with a constrained local tilt bank.

    This post-fit stage cannot change the selected normal, affine calibration,
    validation count, bootstrap result, or orientation score. Each additional
    assignment must recur in multiple incidence-angle images. Multiple peaks may
    use the same base reflection only through distinct local tilt variants, with
    a configurable capacity and penalty.
    """
    assignment_columns = [
        "feature_id", "hkl", "h", "k", "l", "qr_exp", "qz_exp", "qr_calc", "qz_calc",
        "delta_qr", "delta_qz", "delta_q", "normalized_delta", "strength", "support",
        "support_fraction", "feature_type", "feature_quality", "f2", "prediction_weight",
        "assignment_ambiguity", "assignment_margin_sigma", "orientation_domain",
        "base_orientation_domain", "mosaic_tilt_x_offset_deg", "mosaic_tilt_y_offset_deg",
        "mosaic_tilt_magnitude_deg", "assignment_support_score", "member_observations",
        "member_angle_support", "member_gate_fraction", "role", "index_source",
        "completion_is_primary_orientation_evidence", "completion_is_validation_evidence",
    ]
    diagnostic_columns = [
        "feature_id", "mosaic_completion_status", "nearest_delta_q",
        "nearest_normalized_delta", "unique_candidate_count", "support", "feature_type",
        "member_angle_support", "member_gate_fraction",
    ]
    enabled = bool(getattr(config, "v72_enable_s1_mosaic_completion", False))
    targeted = any(str(series_id).startswith(prefix)
                   for prefix in getattr(config, "v72_s1_series_prefixes", ("s1:",)))
    weighted_fraction = float(search.get("best", {}).get("all_feature_metrics", {}).get(
        "weighted_fraction", np.nan
    ))
    weak_enough = (not np.isfinite(weighted_fraction) or weighted_fraction <
                   float(config.v72_mosaic_only_when_weighted_fraction_below))
    if not enabled or not targeted or not weak_enough or consensus is None or consensus.empty:
        return pd.DataFrame(columns=assignment_columns), pd.DataFrame(columns=diagnostic_columns)

    if secondary and secondary.get("accepted") and secondary.get("joint") is not None:
        base_tables = [secondary["joint"].get("anchor_matches", pd.DataFrame()),
                       secondary["joint"].get("validation_matches", pd.DataFrame())]
        available = [x for x in base_tables if x is not None and not x.empty]
        base = pd.concat(available, ignore_index=True, sort=False) if available else pd.DataFrame()
        domain_solutions = secondary.get("domain_solutions", [])
    else:
        base = combined_matches(search)
        domain_solutions = [{"domain": "primary", "hkl": search["best"]["hkl"],
                             "params": search["best"]["params"]}]
    if existing_completion is not None and not existing_completion.empty:
        base = pd.concat([base, existing_completion], ignore_index=True, sort=False)
    if not domain_solutions:
        domain_solutions = [{"domain": "primary", "hkl": search["best"]["hkl"],
                             "params": search["best"]["params"]}]

    used_ids = set(base.feature_id.astype(str)) if not base.empty else set()
    residual = consensus[
        ~consensus.feature_id.astype(str).isin(used_ids)
        & ~consensus.feature_id.astype(str).isin(config.manual_rejected_feature_ids)
        ].copy().reset_index(drop=True)
    if residual.empty:
        return pd.DataFrame(columns=assignment_columns), pd.DataFrame(columns=diagnostic_columns)
    sigma = np.sqrt(np.maximum(residual.cov_rr.to_numpy(float)
                               + residual.cov_zz.to_numpy(float), 1e-12))
    allowed_type = ~residual.feature_type.astype(str).isin(config.ignored_feature_types)
    quality_mask = (
            allowed_type.to_numpy(bool)
            & (residual.support.to_numpy(int) >= int(config.v72_mosaic_min_support))
            & (sigma <= float(config.v72_mosaic_max_sigma_q))
            & (residual.major_width_q.to_numpy(float) <= float(config.v72_mosaic_max_major_width_q))
            & (residual.qz.to_numpy(float) >= float(config.analysis_qz_min))
    )
    residual = residual[quality_mask].copy().reset_index(drop=True)
    if residual.empty:
        return pd.DataFrame(columns=assignment_columns), pd.DataFrame(columns=diagnostic_columns)

    banks = []
    offsets = tuple(getattr(config, "v72_mosaic_tilt_offsets_deg", ()))
    for solution in domain_solutions:
        base_domain = str(solution.get("domain", "primary"))
        normal_hkl = tuple(solution["hkl"])
        base_params = np.asarray(solution["params"], float)
        for variant_index, (dx, dy) in enumerate(offsets):
            params = base_params.copy()
            params[0] = np.clip(params[0] + float(dx), -config.max_tilt_anchor_deg,
                                config.max_tilt_anchor_deg)
            params[1] = np.clip(params[1] + float(dy), -config.max_tilt_anchor_deg,
                                config.max_tilt_anchor_deg)
            actual_dx = float(params[0] - base_params[0])
            actual_dy = float(params[1] - base_params[1])
            if math.hypot(actual_dx, actual_dy) < 1e-8:
                continue
            variant = f"{base_domain}:mosaic:{variant_index + 1:02d}"
            prediction = _v7_prediction_array(
                crystal, normal_hkl, config, params, config.v72_mosaic_f2_percentile,
                config.v72_mosaic_max_predictions_per_variant, variant,
            )
            if not len(prediction["qr"]):
                continue
            n = len(prediction["qr"])
            prediction = {key: np.asarray(value) for key, value in prediction.items()}
            prediction["base_domain"] = np.full(n, base_domain, dtype=object)
            prediction["tilt_dx"] = np.full(n, actual_dx, float)
            prediction["tilt_dy"] = np.full(n, actual_dy, float)
            prediction["tilt_magnitude"] = np.full(n, math.hypot(actual_dx, actual_dy), float)
            banks.append(prediction)
    if not banks:
        return pd.DataFrame(columns=assignment_columns), pd.DataFrame(columns=diagnostic_columns)
    keys = banks[0].keys()
    predictions = {key: np.concatenate([bank[key] for bank in banks]) for key in keys}

    sign = np.sign(predictions["qr"]).astype(int)
    order = np.argsort(predictions["tilt_magnitude"], kind="mergesort")
    kept = []
    positions = {}
    for j in order:
        key = (str(predictions["base_domain"][j]), int(predictions["h"][j]),
               int(predictions["k"][j]), int(predictions["l"][j]), int(sign[j]))
        point = (float(predictions["qr"][j]), float(predictions["qz"][j]))
        previous = positions.setdefault(key, [])
        if all(math.hypot(point[0] - p[0], point[1] - p[1]) >=
               float(config.v72_mosaic_variant_separation_q) for p in previous):
            kept.append(int(j));
            previous.append(point)
    kept = np.asarray(kept, int)
    predictions = {key: np.asarray(value)[kept] for key, value in predictions.items()}
    if not len(predictions["qr"]):
        return pd.DataFrame(columns=assignment_columns), pd.DataFrame(columns=diagnostic_columns)

    physical, normalized = _feature_distance_matrices(residual, predictions, config)
    gate = ((physical <= float(config.v72_mosaic_tolerance_q))
            & (normalized <= float(config.v72_mosaic_sigma_limit)))
    quality, _ = _v7_feature_quality(residual, config)
    sign = np.sign(predictions["qr"]).astype(int)
    base_keys = [
        (str(predictions["base_domain"][j]), int(predictions["h"][j]),
         int(predictions["k"][j]), int(predictions["l"][j]), int(sign[j]))
        for j in range(len(predictions["qr"]))
    ]
    unique_ambiguity = np.zeros(len(residual), int)
    margins = np.zeros(len(residual), float)
    for i in range(len(residual)):
        js = np.flatnonzero(gate[i])
        best_by_key = {}
        for j in js:
            key = base_keys[j]
            best_by_key[key] = min(best_by_key.get(key, np.inf), float(normalized[i, j]))
        values = np.sort(np.asarray(list(best_by_key.values()), float))
        unique_ambiguity[i] = len(values)
        margins[i] = (float(values[1] - values[0]) if len(values) > 1
                      else float(config.v72_mosaic_sigma_limit))
    cost = (normalized / max(float(config.v72_mosaic_sigma_limit), 1e-9)
            + float(config.v72_mosaic_tilt_penalty)
            * predictions["tilt_magnitude"][None, :] / 1.5
            + 0.04 * np.log1p(np.maximum(unique_ambiguity - 1, 0))[:, None])
    real = np.nan_to_num(cost, nan=1e4, posinf=1e4, neginf=-1e4)
    real[~gate] = 1e4
    dummy_cost = 1.02 + 0.08 * quality
    dummy = np.full((len(residual), len(residual)), 1e4, float)
    np.fill_diagonal(dummy, dummy_cost)
    rows, cols = linear_sum_assignment(np.hstack((real, dummy)))
    selected = cols < len(predictions["qr"])
    rows, cols = rows[selected], cols[selected]
    selected = gate[rows, cols] & (cost[rows, cols] <= dummy_cost[rows])
    rows, cols = rows[selected], cols[selected]

    records = []
    for i, j in zip(rows, cols):
        exp = residual.iloc[i]
        records.append({
            "feature_id": str(exp.feature_id), "hkl": str(predictions["hkl"][j]),
            "h": int(predictions["h"][j]), "k": int(predictions["k"][j]),
            "l": int(predictions["l"][j]), "qr_exp": float(exp.qr),
            "qz_exp": float(exp.qz), "qr_calc": float(predictions["qr"][j]),
            "qz_calc": float(predictions["qz"][j]),
            "delta_qr": float(exp.qr - predictions["qr"][j]),
            "delta_qz": float(exp.qz - predictions["qz"][j]),
            "delta_q": float(physical[i, j]), "normalized_delta": float(normalized[i, j]),
            "strength": float(exp.strength), "support": int(exp.support),
            "support_fraction": float(exp.support_fraction),
            "feature_type": str(exp.feature_type), "feature_quality": float(quality[i]),
            "f2": float(predictions["f2"][j]),
            "prediction_weight": float(predictions["prediction_weight"][j]),
            "assignment_ambiguity": int(unique_ambiguity[i]),
            "assignment_margin_sigma": float(margins[i]),
            "orientation_domain": str(predictions["domain"][j]),
            "base_orientation_domain": str(predictions["base_domain"][j]),
            "mosaic_tilt_x_offset_deg": float(predictions["tilt_dx"][j]),
            "mosaic_tilt_y_offset_deg": float(predictions["tilt_dy"][j]),
            "mosaic_tilt_magnitude_deg": float(predictions["tilt_magnitude"][j]),
            "assignment_support_score": float(
                quality[i] * math.exp(-0.5 * normalized[i, j] ** 2)
                * min(1.0, margins[i]) / math.sqrt(max(unique_ambiguity[i], 1))
            ),
        })
    candidates = pd.DataFrame(records)
    if candidates.empty:
        diagnostics = pd.DataFrame({
            "feature_id": residual.feature_id.astype(str),
            "mosaic_completion_status": "no_assignment",
            "nearest_delta_q": np.min(physical, axis=1),
            "nearest_normalized_delta": np.min(normalized, axis=1),
            "unique_candidate_count": unique_ambiguity,
            "support": residual.support.to_numpy(int),
            "feature_type": residual.feature_type.astype(str),
            "member_angle_support": 0, "member_gate_fraction": np.nan,
        })
        return pd.DataFrame(columns=assignment_columns), diagnostics.reindex(columns=diagnostic_columns)

    evidence = _member_corroboration(candidates, members, config)
    candidates = candidates.merge(evidence, on="feature_id", how="left")
    general_ok = (
            (candidates.assignment_ambiguity <= int(config.v72_mosaic_max_unique_ambiguity))
            & ((candidates.assignment_ambiguity <= 1)
               | (candidates.assignment_margin_sigma >= float(config.v72_mosaic_min_margin_sigma)))
            & (candidates.member_angle_support >= int(config.v72_mosaic_min_member_angles))
            & (candidates.member_gate_fraction >= float(config.v72_mosaic_min_member_fraction))
            & (candidates.normalized_delta <= float(config.v72_mosaic_sigma_limit))
    )
    candidates = candidates[general_ok].copy()
    if not candidates.empty:
        candidates["_sign"] = np.sign(candidates.qr_calc).astype(int)
        candidates["_capacity_key"] = list(zip(
            candidates.base_orientation_domain.astype(str), candidates.h.astype(int),
            candidates.k.astype(int), candidates.l.astype(int), candidates["_sign"],
        ))
        candidates = candidates.sort_values(
            ["assignment_support_score", "normalized_delta"], ascending=[False, True],
            kind="mergesort",
        )
        candidates["_capacity_rank"] = candidates.groupby("_capacity_key").cumcount() + 1
        candidates = candidates[
            candidates["_capacity_rank"] <= int(config.v72_mosaic_max_assignments_per_reflection)
            ].copy().drop(columns=["_sign", "_capacity_key", "_capacity_rank"])
        candidates["role"] = "s1_mosaic_split_completion"
        candidates["index_source"] = "s1_mosaic_tilt_bank_completion"
        candidates["completion_is_primary_orientation_evidence"] = False
        candidates["completion_is_validation_evidence"] = False

    accepted_ids = set(candidates.feature_id.astype(str)) if not candidates.empty else set()
    diagnostics = pd.DataFrame({
        "feature_id": residual.feature_id.astype(str),
        "mosaic_completion_status": np.where(
            residual.feature_id.astype(str).isin(accepted_ids),
            "indexed_s1_mosaic_completion", "not_completed"
        ),
        "nearest_delta_q": np.min(physical, axis=1),
        "nearest_normalized_delta": np.min(normalized, axis=1),
        "unique_candidate_count": unique_ambiguity,
        "support": residual.support.to_numpy(int),
        "feature_type": residual.feature_type.astype(str),
    })
    support_frame = (candidates[["feature_id", "member_angle_support", "member_gate_fraction"]]
                     if not candidates.empty else pd.DataFrame(columns=[
        "feature_id", "member_angle_support", "member_gate_fraction"
    ]))
    diagnostics = diagnostics.merge(support_frame, on="feature_id", how="left")
    return candidates.reindex(columns=assignment_columns), diagnostics.reindex(columns=diagnostic_columns)


def _overlay_display_intensity(image, config):
    """Return a display copy with internal blank runs bridged for presentation.

    This function never changes ``image["intensity"]`` or ``image["valid"]``.
    Filled pixels are therefore excluded from feature detection, orientation
    fitting, completion, and validation. The returned mask identifies only the
    pixels synthesized for the visualization.
    """
    intensity = np.asarray(image["intensity"], float)
    display = np.asarray(image.get("display_intensity", intensity), float).copy()
    valid = np.asarray(image.get("valid", np.isfinite(intensity)), bool)
    source = valid & np.isfinite(display)
    filled = np.zeros(display.shape, dtype=bool)
    if not bool(getattr(config, "overlay_fill_display_gaps", False)) or not source.any():
        return display, filled

    fraction = float(np.clip(getattr(config, "overlay_gap_fill_max_fraction", 0.35), 0.0, 1.0))
    height, width = display.shape
    max_row_run = max(1, int(round(width * fraction)))
    max_col_run = max(1, int(round(height * fraction)))

    # Bridge only runs bounded by measured pixels on both sides. This avoids
    # inventing an extrapolated detector edge while closing internal module or
    # stitching gaps.
    for row_index in range(height):
        measured = np.flatnonzero(source[row_index])
        if measured.size < 2:
            continue
        left, right = int(measured[0]), int(measured[-1])
        x = np.arange(left, right + 1)
        row_source = source[row_index, left:right + 1]
        for start in np.flatnonzero(~row_source & np.r_[True, row_source[:-1]]):
            stop_candidates = np.flatnonzero(row_source[start:])
            if stop_candidates.size == 0:
                continue
            stop = start + int(stop_candidates[0])
            run_length = stop - start
            if run_length <= 0 or run_length > max_row_run or start == 0:
                continue
            global_start, global_stop = left + start, left + stop
            x0, x1 = global_start - 1, global_stop
            if not (source[row_index, x0] and source[row_index, x1]):
                continue
            positions = np.arange(global_start, global_stop)
            display[row_index, positions] = np.interp(
                positions, [x0, x1], [display[row_index, x0], display[row_index, x1]]
            )
            filled[row_index, positions] = True

    working = source | filled
    for column_index in range(width):
        measured = np.flatnonzero(working[:, column_index])
        if measured.size < 2:
            continue
        top, bottom = int(measured[0]), int(measured[-1])
        column_source = working[top:bottom + 1, column_index]
        for start in np.flatnonzero(~column_source & np.r_[True, column_source[:-1]]):
            stop_candidates = np.flatnonzero(column_source[start:])
            if stop_candidates.size == 0:
                continue
            stop = start + int(stop_candidates[0])
            run_length = stop - start
            if run_length <= 0 or run_length > max_col_run or start == 0:
                continue
            global_start, global_stop = top + start, top + stop
            y0, y1 = global_start - 1, global_stop
            if not (working[y0, column_index] and working[y1, column_index]):
                continue
            positions = np.arange(global_start, global_stop)
            display[positions, column_index] = np.interp(
                positions, [y0, y1], [display[y0, column_index], display[y1, column_index]]
            )
            filled[positions, column_index] = True

    if filled.any():
        sigma = max(0.0, float(getattr(config, "overlay_gap_fill_smoothing_sigma_px", 7.0)))
        if sigma > 0.0:
            numerator = gaussian_filter(
                np.where(source, np.asarray(image.get("display_intensity", intensity), float), 0.0), sigma=sigma)
            denominator = gaussian_filter(source.astype(float), sigma=sigma)
            smooth = np.divide(
                numerator, denominator, out=np.zeros_like(numerator), where=denominator > 1e-8
            )
            usable = filled & (denominator > 1e-8) & np.isfinite(smooth)
            display[usable] = smooth[usable]
    return display, filled


def _indexed_overlay_label_rows(indexed, config):
    """Return one deterministic label row for every indexed consensus feature."""
    if indexed is None or indexed.empty:
        return pd.DataFrame()
    rows = indexed.copy()
    if "qr_exp" not in rows and "qr" in rows:
        rows["qr_exp"] = rows["qr"]
    if "qz_exp" not in rows and "qz" in rows:
        rows["qz_exp"] = rows["qz"]
    rows = rows.sort_values(
        ["qz_exp", "qr_exp", "feature_id"],
        ascending=[False, True, True], kind="mergesort",
    ).reset_index(drop=True)
    rows["overlay_label_id"] = np.arange(1, len(rows) + 1, dtype=int)

    decimals = max(0, int(config.overlay_coordinate_decimals))

    def hkl_text(row):
        value = str(row.get("hkl", "")).strip()
        if value and value.lower() != "nan":
            return value
        try:
            return f"({int(row['h'])} {int(row['k'])} {int(row['l'])})"
        except Exception:
            return "(h k l unavailable)"

    rows["overlay_hkl_text"] = [hkl_text(row) for _, row in rows.iterrows()]
    rows["overlay_coordinate_text"] = [
        f"q=({float(qr):.{decimals}f}, {float(qz):.{decimals}f})"
        for qr, qz in zip(rows.qr_exp, rows.qz_exp)
    ]
    rows["overlay_key_text"] = [
        f"{int(label_id):02d}  {hkl}  {coordinate}"
        for label_id, hkl, coordinate in zip(
            rows.overlay_label_id, rows.overlay_hkl_text, rows.overlay_coordinate_text
        )
    ]
    return rows


def _draw_indexed_coordinate_key(axis, label_rows, config, display_gap_filled=False):
    """Draw a complete, multi-column HKL and reciprocal-coordinate key."""
    axis.set_axis_off()
    axis.set_xlim(0.0, 1.0)
    axis.set_ylim(0.0, 1.0)
    axis.text(
        0.0, 0.995, "Indexed-point coordinate key",
        transform=axis.transAxes, ha="left", va="top", fontsize=8.5,
        fontweight="bold",
    )
    axis.text(
        0.0, 0.965, "ID   (h k l)   q=(qr, qz) [A^-1]",
        transform=axis.transAxes, ha="left", va="top", fontsize=6.2,
    )
    if display_gap_filled and bool(getattr(config, "overlay_note_display_gap_fill", True)):
        axis.text(
            0.0, 0.943,
            "Blank detector/plot gaps are interpolated for display only; indexing uses measured pixels.",
            transform=axis.transAxes, ha="left", va="top", fontsize=5.2, color="0.30",
        )
    if label_rows.empty:
        axis.text(0.0, 0.92, "No indexed points", transform=axis.transAxes,
                  ha="left", va="top", fontsize=7)
        return

    requested_columns = max(1, int(config.overlay_coordinate_key_columns))
    # Add columns automatically when a series contains too many labels for one
    # readable vertical list.  This remains deterministic for identical inputs.
    columns = max(requested_columns, int(math.ceil(len(label_rows) / 36.0)))
    rows_per_column = int(math.ceil(len(label_rows) / columns))
    usable_top = 0.915 if display_gap_filled and bool(getattr(config, "overlay_note_display_gap_fill", True)) else 0.935
    usable_bottom = 0.015
    line_step = (usable_top - usable_bottom) / max(rows_per_column, 1)
    column_width = 1.0 / columns
    fontsize = float(config.overlay_coordinate_key_fontsize)

    for index, row in label_rows.iterrows():
        column = index // rows_per_column
        line = index % rows_per_column
        x = column * column_width
        y = usable_top - line * line_step
        axis.text(
            x, y, str(row.overlay_key_text), transform=axis.transAxes,
            ha="left", va="top", fontsize=fontsize, family="monospace",
            clip_on=False,
        )


def _overlay_render_background(image, config):
    """Return the display background used by the main indexed overlay.

    PNG measurements use the cropped source RGB pixels directly so the actual
    ``indexed_or_ignored_overlay.png`` matches the verified preview. Numerical
    inputs fall back to the finite scientific display/intensity array. The
    scientific validity mask is not painted onto the displayed image.

    Returns
    -------
    tuple[np.ndarray, bool]
        The background array and whether it is an RGB image.
    """
    rgb = image.get("rgb")
    if rgb is not None:
        array = np.asarray(rgb)
        if array.ndim == 3 and array.shape[-1] >= 3:
            if array.dtype.kind in {"u", "i"}:
                array = array.astype(np.float32) / 255.0
            else:
                array = array.astype(np.float32)
                finite = np.isfinite(array)
                if finite.any() and float(np.nanmax(array[finite])) > 1.0:
                    array = array / 255.0
            array = np.clip(array[..., :3], 0.0, 1.0)
            # Reject only a genuinely blank/all-white crop.
            if not np.all(array > 0.985):
                return array, True

    source = image.get("display_intensity", image.get("intensity"))
    if source is None:
        raise ValueError("Overlay image has no RGB or intensity data.")
    intensity = np.asarray(source, dtype=float).copy()
    finite = np.isfinite(intensity)
    if not finite.any():
        raise ValueError("Overlay intensity array contains no finite values.")
    intensity[~finite] = float(np.nanmin(intensity[finite]))
    return intensity, False


def _plot_binary_overlay_tables(image, indexed, ignored, output, title, config,
                                key_filename=None, dpi_override=None):
    """Render the main two-panel overlay using the preview-matched background."""
    label_rows = _indexed_overlay_label_rows(indexed, config)
    show_key = bool(config.overlay_show_coordinate_key)
    if show_key:
        figure = plt.figure(figsize=(15.2, 8.4))
        grid = figure.add_gridspec(1, 2, width_ratios=(1.05, 1.0), wspace=0.08)
        axis = figure.add_subplot(grid[0, 0])
        key_axis = figure.add_subplot(grid[0, 1])
    else:
        figure, axis = plt.subplots(figsize=(8.4, 7.2))
        key_axis = None

    qr_min = float(np.min(image["qr"]))
    qr_max = float(np.max(image["qr"]))
    qz_min = float(np.min(image["qz"]))
    qz_max = float(np.max(image["qz"]))
    extent = [qr_min, qr_max, qz_min, qz_max]

    background, is_rgb = _overlay_render_background(image, config)
    if is_rgb:
        axis.imshow(
            background, origin="upper", extent=extent, aspect="equal",
            interpolation="nearest", resample=False,
        )
    else:
        axis.imshow(
            np.flipud(background), origin="lower", cmap=config.colormap,
            extent=extent, aspect="equal", interpolation="nearest", resample=False,
        )

    if ignored is not None and not ignored.empty:
        axis.scatter(
            ignored.qr, ignored.qz, s=24, marker="x", c="0.72", alpha=0.72,
            linewidths=0.8, label=f"ignored ({len(ignored)})",
        )

    if indexed is not None and not indexed.empty:
        axis.scatter(
            indexed.qr_exp, indexed.qz_exp, s=52, marker="o", facecolors="none",
            edgecolors="lime", linewidths=1.35, label=f"indexed ({len(indexed)})",
        )
        if bool(config.overlay_label_all_indexed):
            shown = label_rows
        else:
            rank_column = next((column for column in (
                "salvage_evidence_score", "assignment_support_score", "prediction_weight"
            ) if column in label_rows), None)
            shown = label_rows
            if rank_column is not None:
                shown = shown.sort_values(rank_column, ascending=False, kind="mergesort")
            shown = shown.head(max(0, int(config.max_labels)))

        offsets = ((4, 4), (4, -10), (-4, 4), (-4, -10),
                   (9, 0), (-9, 0), (0, 9), (0, -13))
        text_effects = [
            matplotlib_patheffects.Stroke(linewidth=1.6, foreground="black"),
            matplotlib_patheffects.Normal(),
        ]
        for display_index, row in enumerate(shown.itertuples(index=False)):
            offset = offsets[display_index % len(offsets)]
            label_id = int(getattr(row, "overlay_label_id"))
            hkl = str(getattr(row, "overlay_hkl_text"))
            label_text = (
                f"{label_id:02d} {hkl}"
                if bool(config.overlay_label_include_hkl) else f"{label_id:02d}"
            )
            annotation = axis.annotate(
                label_text, (float(row.qr_exp), float(row.qz_exp)), color="white",
                fontsize=float(config.overlay_label_fontsize), xytext=offset,
                textcoords="offset points",
                ha="left" if offset[0] >= 0 else "right",
                va="bottom" if offset[1] >= 0 else "top",
                annotation_clip=True, zorder=5,
            )
            annotation.set_path_effects(text_effects)

    axis.set(
        xlabel=r"$q_r$ ($\AA^{-1}$)", ylabel=r"$q_z$ ($\AA^{-1}$)", title=title,
        xlim=(qr_min, qr_max), ylim=(qz_min, qz_max),
    )
    axis.legend(fontsize=8, loc="upper right")
    if key_axis is not None:
        _draw_indexed_coordinate_key(key_axis, label_rows, config, False)

    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    key_columns = [column for column in (
        "overlay_label_id", "feature_id", "hkl", "h", "k", "l",
        "qr_exp", "qz_exp", "orientation_domain", "index_source",
        "overlay_hkl_text", "overlay_coordinate_text", "overlay_key_text",
    ) if column in label_rows.columns]
    if key_filename is None:
        key_path = output.with_name(
            "indexed_overlay_coordinate_key.csv"
            if output.name == "indexed_or_ignored_overlay.png"
            else f"{output.stem}_coordinate_key.csv"
        )
    else:
        key_path = Path(key_filename)
    label_rows[key_columns].to_csv(key_path, index=False)

    figure.tight_layout()
    figure.savefig(
        output,
        dpi=(config.dpi if dpi_override is None else dpi_override),
        bbox_inches="tight",
    )
    plt.close(figure)
    return indexed, ignored


def _plot_fast_per_angle_overlay(image, indexed, ignored, output, title, config):
    """Render a fully labeled per-angle overlay directly with Pillow."""
    intensity, display_fill_mask = _overlay_display_intensity(image, config)
    finite = np.isfinite(intensity)
    normalized = np.where(finite, np.clip(intensity, 0.0, 1.0), 0.0)
    rgb = (plt.get_cmap(config.colormap)(normalized)[..., :3] * 255).astype(np.uint8)
    base = Image.fromarray(rgb, mode="RGB").convert("RGBA")
    width, height = base.size
    top_margin = 28
    canvas = Image.new("RGBA", (width, height + top_margin), (255, 255, 255, 255))
    canvas.alpha_composite(base, (0, top_margin))
    draw = ImageDraw.Draw(canvas, "RGBA")
    font = ImageFont.load_default()
    small_font = ImageFont.load_default()
    draw.rectangle((0, 0, width, top_margin), fill=(255, 255, 255, 245))
    title_text = str(title)
    if display_fill_mask.any() and bool(getattr(config, "overlay_note_display_gap_fill", True)):
        title_text += " | blank gaps filled for display only"
    draw.text((6, 7), title_text, fill=(0, 0, 0, 255), font=font)

    qr_min, qr_max = float(image["qr"][0]), float(image["qr"][-1])
    qz_min, qz_max = float(np.min(image["qz"])), float(np.max(image["qz"]))

    def pixel(qr, qz):
        x = (float(qr) - qr_min) / max(qr_max - qr_min, 1e-12) * (width - 1)
        y = (qz_max - float(qz)) / max(qz_max - qz_min, 1e-12) * (height - 1) + top_margin
        return int(round(x)), int(round(y))

    if ignored is not None and not ignored.empty:
        finite_ignored = ignored[
            np.isfinite(ignored["qr"].to_numpy(float)) & np.isfinite(ignored["qz"].to_numpy(float))]
        for row in finite_ignored.itertuples(index=False):
            x, y = pixel(row.qr, row.qz)
            draw.line((x - 3, y - 3, x + 3, y + 3), fill=(185, 185, 185, 210), width=1)
            draw.line((x - 3, y + 3, x + 3, y - 3), fill=(185, 185, 185, 210), width=1)

    label_rows = _indexed_overlay_label_rows(indexed, config)
    if not label_rows.empty:
        label_rows = label_rows[
            np.isfinite(label_rows["qr_exp"].to_numpy(float))
            & np.isfinite(label_rows["qz_exp"].to_numpy(float))
            ].reset_index(drop=True)
    decimals = max(0, int(config.overlay_coordinate_decimals))
    offsets = ((5, -13), (5, 4), (-76, -13), (-76, 4))
    for index, row in enumerate(label_rows.itertuples(index=False)):
        x, y = pixel(row.qr_exp, row.qz_exp)
        radius = 4
        draw.ellipse((x - radius, y - radius, x + radius, y + radius),
                     outline=(0, 255, 0, 255), width=2)
        label_text = (
            f"{row.overlay_hkl_text} "
            f"({float(row.qr_exp):.{decimals}f},{float(row.qz_exp):.{decimals}f})"
        )
        dx, dy = offsets[index % len(offsets)]
        tx, ty = x + dx, y + dy
        box = draw.textbbox((tx, ty), label_text, font=small_font)
        box = (box[0] - 1, box[1] - 1, box[2] + 1, box[3] + 1)
        draw.rectangle(box, fill=(0, 0, 0, 145))
        draw.text((tx, ty), label_text, fill=(255, 255, 255, 255), font=small_font)

    output = Path(output);
    output.parent.mkdir(parents=True, exist_ok=True)
    canvas.convert("RGB").save(output, format="PNG", optimize=False)
    label_rows.to_csv(output.with_name(f"{output.stem}_coordinate_key.csv"), index=False)
    return indexed, ignored


def plot_indexed_ignored_overlay(image, search, ignored_salvage, guided_rescue,
                                 output, title, config, secondary=None, completion=None,
                                 extra_indexed=None):
    """Create the one canonical old-style indexed/ignored overlay."""
    indexed, ignored = _overlay_indexed_and_ignored(
        search, ignored_salvage, guided_rescue, config,
        secondary, completion, extra_indexed,
    )
    return _plot_binary_overlay_tables(
        image, indexed, ignored, output, title, config,
    )


    # Route overlay entry points through the same indexed-versus-ignored
    # visualization so classification is presented consistently.


def synthetic_recovery_test(crystal, config):
    """Internal regression test using the same forward model; not external validation."""
    candidates = orientation_candidates(crystal, replace(config, max_orientation_candidates=12))
    true_hkl = (0, 1, 1) if (0, 1, 1) in candidates else candidates[0]
    local = replace(
        config, hypothesis_f2_percentile=72.0, validation_f2_percentile=35.0,
        max_hypothesis_predictions=90, max_validation_predictions=180,
        anchor_min_abs_qr=0.15, anchor_strength_quantile=0.0,
        anchor_match_tolerance_q=0.028, validation_match_tolerance_q=0.038,
        max_anchor_features=8, max_validation_features=10,
        max_normal_candidates_anchor=5, refine_hypotheses=7,
        expand_top_normals=5, full_leave_one_angle_out=False,
        full_bootstrap_iterations=0, test_second_orientation=False,
    )
    predictions = project_reflections(crystal, true_hkl, local, [2.0, -1.5, 1.004, 0.998, 0.006, -0.004], 68, 180)
    pool = predictions[
        (abs(predictions.qr) > local.anchor_min_abs_qr) & (predictions.qz > local.anchor_min_qz)].sort_values("f2",
                                                                                                              ascending=False)
    selected = []
    for row in pool.itertuples(index=False):
        if all(math.hypot(row.qr - previous.qr, row.qz - previous.qz) > 0.11 for previous in selected):
            selected.append(row)
        if len(selected) >= 14:
            break
    if len(selected) < 10:
        return {"passed": False, "reason": "too_few_separated_synthetic_predictions"}
    rng = np.random.default_rng(config.random_seed + 404)
    rows = []
    for i, pred in enumerate(selected, 1):
        qr, qz = np.array([pred.qr, pred.qz]) + rng.normal(0, 0.0035, 2)
        rows.append({
            "feature_id": f"S{i:03d}", "qr": qr, "qz": qz,
            "cov_rr": 0.0045 ** 2, "cov_rz": 0.0, "cov_zz": 0.0045 ** 2,
            "sigma_qr": 0.0045, "sigma_qz": 0.0045,
            "strength": max(0.35, 1 - 0.035 * i),
            "support": 5 if i <= 8 else 2, "support_fraction": 1.0 if i <= 8 else 0.4,
            "angles": "synthetic", "feature_type": "spot",
            "major_width_q": 0.012, "minor_width_q": 0.009,
        })
    for i in range(2):
        rows.append({
            "feature_id": f"O{i + 1:03d}", "qr": rng.uniform(0.35, 1.9), "qz": rng.uniform(0.35, 2.2),
            "cov_rr": 0.010 ** 2, "cov_rz": 0.0, "cov_zz": 0.010 ** 2,
            "sigma_qr": 0.010, "sigma_qz": 0.010, "strength": 0.12,
            "support": 2, "support_fraction": 0.4, "angles": "synthetic", "feature_type": "spot",
            "major_width_q": 0.025, "minor_width_q": 0.015,
        })
    synthetic = pd.DataFrame(rows)
    search = anchor_first_search(crystal, synthetic, local, candidate_override=candidates, fast=True)
    angle_error = normal_angle(search["best"]["hkl"], true_hkl, crystal)
    return {
        "passed": bool(angle_error <= config.orientation_stability_angle_deg),
        "true_normal_hkl": true_hkl, "recovered_normal_hkl": search["best"]["hkl"],
        "normal_angle_error_deg": angle_error,
        "anchor_matches": search["best"]["anchor_metrics"]["matches"],
        "validation_matches": search["best"]["validation_metrics"]["matches"],
        "score_margin": search["score_margin"],
    }


def _v73_sector_groups(residual: pd.DataFrame, config: IndexingConfig):
    """Return overlapping left/right and qz-band residual subsets.

    Overlapping sectors are intentional: a physically coherent domain may be
    most visible in one quadrant but still explain peaks elsewhere after global
    reassignment. Each accepted domain must therefore pass both a local-sector
    and a full-pattern improvement test.
    """
    if residual is None or residual.empty:
        return []
    split = float(config.v73_sector_qr_split)
    low_edge, high_edge = map(float, config.v73_sector_qz_edges)
    qz_max = float(config.qz_range[1]) + 1e-9
    minimum = int(config.v73_sector_min_features)
    sectors = []

    def add(name, frame):
        frame = frame.copy()
        if len(frame) >= minimum:
            sectors.append((name, frame))

    add("left", residual[residual.qr < -split])
    add("right", residual[residual.qr > split])
    add("low_qz", residual[(residual.qz >= config.analysis_qz_min) & (residual.qz < low_edge)])
    add("mid_qz", residual[(residual.qz >= low_edge) & (residual.qz < high_edge)])
    add("high_qz", residual[(residual.qz >= high_edge) & (residual.qz <= qz_max)])
    add("left_low", residual[(residual.qr < -split) & (residual.qz < low_edge)])
    add("right_low", residual[(residual.qr > split) & (residual.qz < low_edge)])
    add("left_mid", residual[(residual.qr < -split) & (residual.qz >= low_edge) & (residual.qz < high_edge)])
    add("right_mid", residual[(residual.qr > split) & (residual.qz >= low_edge) & (residual.qz < high_edge)])
    add("left_high", residual[(residual.qr < -split) & (residual.qz >= high_edge)])
    add("right_high", residual[(residual.qr > split) & (residual.qz >= high_edge)])
    return sectors


def _v73_domain_rows(search, secondary):
    if secondary and secondary.get("accepted") and secondary.get("domain_solutions"):
        rows = secondary.get("domain_solutions", [])
    else:
        rows = [{"domain": "primary", "hkl": search["best"]["hkl"],
                 "params": search["best"]["params"]}]
    normalized = []
    for row in rows:
        normalized.append({
            "domain": str(row.get("domain", "primary")),
            "hkl": tuple(row["hkl"]),
            "params": np.asarray(row["params"], float),
        })
    return normalized


def _v739_v73_s1_sector_domain_search(series_id, crystal, consensus, search, secondary, config):
    """Add globally validated domains discovered from coherent s1 residual sectors."""
    assignment_columns = [
        "feature_id", "hkl", "h", "k", "l", "qr_exp", "qz_exp", "qr_calc", "qz_calc",
        "delta_qr", "delta_qz", "delta_q", "normalized_delta", "orientation_domain",
        "sector_discovery_region", "sector_domain_is_new_feature", "role", "index_source",
        "sector_domain_is_primary_orientation_evidence", "sector_domain_is_validation_evidence",
    ]
    diagnostic_columns = [
        "iteration", "sector", "sector_feature_count", "normal_h", "normal_k", "normal_l",
        "tilt_x_deg", "tilt_y_deg", "local_matches", "local_weighted_fraction",
        "global_domain_matches", "new_matches", "lost_matches", "local_global_matches",
        "global_score_gain_after_penalty", "global_weighted_fraction_gain", "accepted",
        "rejection_reason",
    ]
    enabled = bool(getattr(config, "v73_enable_s1_sector_domains", False))
    targeted = any(str(series_id).startswith(prefix) for prefix in
                   getattr(config, "v73_s1_sector_series_prefixes", ("s1:A",)))
    weighted_fraction = float(search.get("best", {}).get("all_feature_metrics", {}).get(
        "weighted_fraction", 0.0
    ))
    weak_enough = weighted_fraction < float(config.v73_sector_only_when_weighted_fraction_below)
    if not enabled or not targeted or not weak_enough or consensus is None or consensus.empty:
        return secondary, pd.DataFrame(columns=assignment_columns), pd.DataFrame(columns=diagnostic_columns)

    features = consensus[
        ~consensus.feature_id.astype(str).isin(config.manual_rejected_feature_ids)
    ].copy().reset_index(drop=True)
    domain_rows = _v73_domain_rows(search, secondary)
    domains = [
        _v7_domain_solution(row["hkl"], row["params"], crystal, config, row["domain"])
        for row in domain_rows
    ]
    combined = _v7_concat_predictions([item["prediction_array"] for item in domains])
    current_matches, current_metrics = _v7_assign_arrays(
        features, combined, config, config.v7_all_feature_tolerance_q, True
    )
    base_feature_ids = set(current_matches.feature_id.astype(str)) if not current_matches.empty else set()
    primary_only = _v7_domain_solution(
        search["best"]["hkl"], search["best"]["params"], crystal, config, "primary"
    )
    _, primary_metrics = _v7_assign_arrays(
        features, primary_only["prediction_array"], config, config.v7_all_feature_tolerance_q, True
    )
    candidate_normals = orientation_candidates(
        crystal, replace(config, max_orientation_candidates=int(config.v73_sector_candidate_normals))
    )
    primary_params = np.asarray(search["best"]["params"], float)
    diagnostics = []
    accepted_records = []

    for iteration in range(1, int(config.v73_sector_max_extra_domains) + 1):
        matched_ids = set(current_matches.feature_id.astype(str)) if not current_matches.empty else set()
        residual = features[~features.feature_id.astype(str).isin(matched_ids)].copy().reset_index(drop=True)
        if residual.empty:
            break
        sigma = np.sqrt(np.maximum(
            residual.cov_rr.to_numpy(float) + residual.cov_zz.to_numpy(float), 1e-12
        ))
        quality = (
                (residual.support.to_numpy(int) >= int(config.v73_sector_min_support))
                & (sigma <= float(config.v73_sector_max_sigma_q))
                & (residual.major_width_q.to_numpy(float) <= float(config.v73_sector_max_major_width_q))
                & (~residual.feature_type.astype(str).isin(config.ignored_feature_types).to_numpy(bool))
                & (residual.qz.to_numpy(float) >= float(config.analysis_qz_min))
        )
        residual = residual[quality].copy().reset_index(drop=True)
        sectors = _v73_sector_groups(residual, config)
        if not sectors:
            break
        used_normals = [item["hkl"] for item in domains]
        local_finalists = []
        for sector_name, sector_features in sectors:
            local_candidates = []
            for normal_hkl in candidate_normals:
                if any(normal_angle(normal_hkl, used, crystal) <
                       float(config.v73_sector_min_normal_separation_deg) for used in used_normals):
                    continue
                for tx in config.v73_sector_tilt_values_deg:
                    for ty in config.v73_sector_tilt_values_deg:
                        params = primary_params.copy()
                        params[0], params[1] = float(tx), float(ty)
                        predictions = _v7_prediction_array(
                            crystal, normal_hkl, config, params,
                            float(config.v73_sector_search_f2_percentile),
                            int(config.v73_sector_search_max_predictions),
                            "sector_candidate",
                        )
                        local_matches, local_metrics = _v7_assign_arrays(
                            sector_features, predictions, config,
                            float(config.v73_sector_tolerance_q), True,
                        )
                        if len(local_matches) < int(config.v73_sector_min_local_matches):
                            continue
                        local_candidates.append({
                            "sector": sector_name,
                            "sector_features": sector_features,
                            "normal_hkl": tuple(normal_hkl),
                            "params": params,
                            "local_matches": local_matches,
                            "local_metrics": local_metrics,
                            "key": (
                                int(len(local_matches)),
                                float(local_metrics.get("weighted_fraction", 0.0)),
                                float(local_metrics.get("score", -np.inf)),
                            ),
                        })
            local_candidates.sort(key=lambda item: item["key"], reverse=True)
            local_finalists.extend(local_candidates[:int(config.v73_sector_top_candidates_per_sector)])

        unique_finalists = []
        seen = set()
        for item in local_finalists:
            key = (item["normal_hkl"], round(float(item["params"][0]), 4),
                   round(float(item["params"][1]), 4))
            if key not in seen:
                seen.add(key)
                unique_finalists.append(item)
        globally_ranked = []
        domain_name = f"s1_sector_{iteration}"
        for item in unique_finalists:
            trial_domain = _v7_domain_solution(
                item["normal_hkl"], item["params"], crystal, config, domain_name
            )
            trial_domain = _v7_refine_domain_tilts(
                trial_domain, item["local_matches"], crystal, config
            )
            trial_predictions = _v7_concat_predictions([
                combined, trial_domain["prediction_array"]
            ])
            trial_matches, trial_metrics = _v7_assign_arrays(
                features, trial_predictions, config, config.v7_all_feature_tolerance_q, True
            )
            domain_matches = trial_matches[
                trial_matches.orientation_domain.astype(str).eq(domain_name)
            ].copy() if not trial_matches.empty else pd.DataFrame()
            previous_ids = set(current_matches.feature_id.astype(str)) if not current_matches.empty else set()
            trial_ids = set(trial_matches.feature_id.astype(str)) if not trial_matches.empty else set()
            new_ids = trial_ids - previous_ids
            lost_ids = previous_ids - trial_ids
            sector_ids = set(item["sector_features"].feature_id.astype(str))
            local_global = len(
                set(domain_matches.feature_id.astype(str)) & sector_ids) if not domain_matches.empty else 0
            score_gain = (
                    float(trial_metrics.get("score", -np.inf))
                    - float(current_metrics.get("score", -np.inf))
                    - float(config.v73_sector_complexity_penalty)
            )
            weighted_gain = (
                    float(trial_metrics.get("weighted_fraction", 0.0))
                    - float(current_metrics.get("weighted_fraction", 0.0))
            )
            reasons = []
            if len(domain_matches) < int(config.v73_sector_min_global_domain_matches):
                reasons.append("too_few_global_domain_matches")
            if len(new_ids) < int(config.v73_sector_min_new_matches):
                reasons.append("too_few_new_matches")
            if local_global < int(config.v73_sector_min_local_global_matches):
                reasons.append("insufficient_sector_specific_support")
            if score_gain < float(config.v73_sector_min_score_gain):
                reasons.append("insufficient_global_score_gain")
            if weighted_gain < float(config.v73_sector_min_weighted_fraction_gain):
                reasons.append("insufficient_weighted_fraction_gain")
            if len(lost_ids) > int(config.v73_sector_max_lost_assignments):
                reasons.append("too_many_displaced_assignments")
            record = {
                "iteration": iteration,
                "sector": item["sector"],
                "sector_feature_count": int(len(item["sector_features"])),
                "normal_h": int(item["normal_hkl"][0]),
                "normal_k": int(item["normal_hkl"][1]),
                "normal_l": int(item["normal_hkl"][2]),
                "tilt_x_deg": float(trial_domain["params"][0]),
                "tilt_y_deg": float(trial_domain["params"][1]),
                "local_matches": int(len(item["local_matches"])),
                "local_weighted_fraction": float(item["local_metrics"].get("weighted_fraction", np.nan)),
                "global_domain_matches": int(len(domain_matches)),
                "new_matches": int(len(new_ids)),
                "lost_matches": int(len(lost_ids)),
                "local_global_matches": int(local_global),
                "global_score_gain_after_penalty": float(score_gain),
                "global_weighted_fraction_gain": float(weighted_gain),
                "accepted": False,
                "rejection_reason": ";".join(reasons),
            }
            diagnostics.append(record)
            globally_ranked.append({
                "key": (score_gain, len(new_ids), local_global, len(domain_matches), weighted_gain),
                "domain": trial_domain,
                "matches": trial_matches,
                "metrics": trial_metrics,
                "domain_matches": domain_matches,
                "predictions": trial_predictions,
                "new_ids": new_ids,
                "lost_ids": lost_ids,
                "sector": item["sector"],
                "accepted": not reasons,
                "diagnostic_index": len(diagnostics) - 1,
            })
        acceptable = [item for item in globally_ranked if item["accepted"]]
        if not acceptable:
            break
        acceptable.sort(key=lambda item: item["key"], reverse=True)
        winner = acceptable[0]
        diagnostics[winner["diagnostic_index"]]["accepted"] = True
        diagnostics[winner["diagnostic_index"]]["rejection_reason"] = ""
        domains.append(winner["domain"])
        combined = winner["predictions"]
        current_matches, current_metrics = winner["matches"], winner["metrics"]
        accepted_records.append({
            "domain": winner["domain"]["domain"],
            "hkl": winner["domain"]["hkl"],
            "params": winner["domain"]["params"],
            "sector": winner["sector"],
            "new_feature_ids": set(winner["new_ids"]),
            "score_gain": float(winner["key"][0]),
        })

    if not accepted_records:
        return secondary, pd.DataFrame(columns=assignment_columns), pd.DataFrame(diagnostics,
                                                                                 columns=diagnostic_columns)

    validation_ids = set(search.get("validation", pd.DataFrame()).feature_id.astype(str))
    joint_validation = current_matches[
        current_matches.feature_id.astype(str).isin(validation_ids)
    ].copy()
    joint_fit = current_matches[
        ~current_matches.feature_id.astype(str).isin(validation_ids)
    ].copy()
    updated = dict(secondary) if secondary else {}
    updated["accepted"] = True
    updated["reason"] = "accepted_v73_s1_sector_aware_domains"
    updated["domain_count"] = int(len(domains))
    updated["domain_solutions"] = [
        {"domain": item["domain"], "hkl": item["hkl"], "params": item["params"]}
        for item in domains
    ]
    updated["s1_sector_domain_count"] = int(len(accepted_records))
    updated["s1_sector_domain_records"] = accepted_records
    updated["joint"] = {
        "score": float(current_metrics.get("score", np.nan)),
        "score_gain": float(
            current_metrics.get("score", np.nan) - primary_metrics.get("score", np.nan)
            - config.v7_multidomain_complexity_penalty * max(len(domains) - 1, 0)
        ),
        "anchor_matches": joint_fit,
        "validation_matches": joint_validation,
        "anchor_metrics": current_metrics,
        "validation_metrics": current_metrics,
        "strict_validation_metrics": _empty_assignment_metrics(score=0.0),
        "secondary_anchor_matches": int((joint_fit.orientation_domain.astype(str) != "primary").sum()),
        "secondary_validation_matches": int((joint_validation.orientation_domain.astype(str) != "primary").sum()),
        "domain_count": int(len(domains)),
        "domain_solutions": updated["domain_solutions"],
    }

    discovery = {item["domain"]: item for item in accepted_records}
    sector_assignments = current_matches[
        current_matches.orientation_domain.astype(str).isin(discovery)
    ].copy()
    if not sector_assignments.empty:
        sector_assignments["sector_discovery_region"] = sector_assignments.orientation_domain.astype(str).map(
            {key: value["sector"] for key, value in discovery.items()}
        )
        new_by_domain = {key: value["new_feature_ids"] for key, value in discovery.items()}
        sector_assignments["sector_domain_is_new_feature"] = [
            str(fid) in new_by_domain.get(str(domain), set())
            for fid, domain in zip(sector_assignments.feature_id, sector_assignments.orientation_domain)
        ]
        sector_assignments["role"] = "s1_sector_aware_domain"
        sector_assignments["index_source"] = "s1_sector_aware_joint_domain_search"
        sector_assignments["sector_domain_is_primary_orientation_evidence"] = True
        sector_assignments["sector_domain_is_validation_evidence"] = sector_assignments.feature_id.astype(str).isin(
            validation_ids)
    return updated, sector_assignments.reindex(columns=assignment_columns, fill_value=np.nan), pd.DataFrame(diagnostics,
                                                                                                            columns=diagnostic_columns)


def _v73_candidate_cif_paths(config: IndexingConfig, primary_path):
    primary = Path(primary_path).resolve()
    try:
        primary_digest = hashlib.sha256(primary.read_bytes()).hexdigest()
    except Exception:
        primary_digest = None
    paths = []
    for value in tuple(getattr(config, "alternative_cif_paths", ())) + tuple(
            getattr(config, "substrate_cif_paths", ())):
        candidate = Path(value).expanduser()
        if candidate.is_file():
            paths.append(candidate.resolve())
    if getattr(config, "v73_auto_discover_candidate_cifs", True):
        roots = [primary.parent, Path.cwd(), Path("/content"), Path("/mnt/data")]
        roots.extend(Path(item).expanduser() for item in getattr(config, "search_dirs", ()))
        for root in roots:
            if not root.is_dir():
                continue
            for pattern in getattr(config, "v73_candidate_cif_globs", ("*.cif", "*.CIF")):
                for candidate in root.glob(pattern):
                    if candidate.is_file():
                        paths.append(candidate.resolve())
    unique = []
    seen = {primary}
    seen_digests = {primary_digest} if primary_digest else set()
    for path in paths:
        if path in seen or path.name == primary.name:
            continue
        try:
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
        except Exception:
            digest = None
        if digest and digest in seen_digests:
            continue
        seen.add(path)
        if digest:
            seen_digests.add(digest)
        unique.append(path)
    return unique


def _v73_candidate_kind(path: Path, config: IndexingConfig):
    explicit_substrates = {Path(item).expanduser().resolve() for item in getattr(config, "substrate_cif_paths", ()) if
                           Path(item).expanduser().is_file()}
    if path.resolve() in explicit_substrates:
        return "substrate"
    name = path.name.lower()
    tokens = ("substrate", "silicon", "sio2", "quartz", "ito", "glass", "sapphire", "al2o3")
    return "substrate_candidate" if any(token in name for token in tokens) else "alternative_phase"


def v73_residual_candidate_phase_screen(series_id, crystal, consensus, indexed_matches, config):
    """Screen residual features against additional phase or substrate CIF candidates.

    The screen is diagnostic unless ``v73_promote_residual_candidate_matches`` is
    explicitly enabled. It combines radial powder compatibility with a compact
    orientation search on features not explained by the primary structure.
    """
    report_columns = [
        "series_id", "status", "candidate_role", "candidate_cif", "spacegroup_number",
        "residual_features", "powder_matched_features", "powder_matched_weight_fraction",
        "powder_median_radial_residual_q", "powder_score", "orientation_normal_h",
        "orientation_normal_k", "orientation_normal_l", "orientation_matches",
        "orientation_weighted_fraction", "orientation_score", "composite_score",
        "score_gap_over_primary_residual", "preferred_over_primary_residual",
        "assignments_promoted",
    ]
    assignment_columns = [
        "feature_id", "hkl", "h", "k", "l", "qr_exp", "qz_exp", "qr_calc", "qz_calc",
        "delta_q", "normalized_delta", "orientation_domain", "candidate_cif",
        "candidate_role", "candidate_rank", "residual_candidate_preferred",
        "role", "index_source", "completion_is_primary_orientation_evidence",
        "completion_is_validation_evidence",
    ]
    if not str(series_id).startswith("s1:"):
        return pd.DataFrame(columns=report_columns), pd.DataFrame(columns=assignment_columns), pd.DataFrame()
    indexed_ids = set(
        indexed_matches.feature_id.astype(str)) if indexed_matches is not None and not indexed_matches.empty else set()
    residual = consensus[
        ~consensus.feature_id.astype(str).isin(indexed_ids)
        & ~consensus.feature_id.astype(str).isin(config.manual_rejected_feature_ids)
        ].copy().reset_index(drop=True)
    if residual.empty:
        row = {"series_id": series_id, "status": "no_residual_features", "residual_features": 0}
        return pd.DataFrame([row], columns=report_columns), pd.DataFrame(columns=assignment_columns), residual
    sigma = np.sqrt(np.maximum(residual.cov_rr.to_numpy(float) + residual.cov_zz.to_numpy(float), 1e-12))
    quality = (
            (residual.support.to_numpy(int) >= int(config.v73_residual_phase_min_support))
            & (sigma <= float(config.v73_residual_phase_max_sigma_q))
            & (residual.major_width_q.to_numpy(float) <= float(config.v73_residual_phase_max_major_width_q))
            & (residual.qz.to_numpy(float) >= float(config.analysis_qz_min))
    )
    residual = residual[quality].copy().reset_index(drop=True)
    candidate_paths = _v73_candidate_cif_paths(config, crystal["path"])
    if len(residual) < int(config.v73_residual_phase_min_features):
        row = {"series_id": series_id, "status": "too_few_quality_residual_features",
               "residual_features": int(len(residual))}
        return pd.DataFrame([row], columns=report_columns), pd.DataFrame(columns=assignment_columns), residual
    if not candidate_paths:
        row = {"series_id": series_id, "status": "no_alternative_or_substrate_cifs_found",
               "candidate_role": "none", "candidate_cif": "",
               "residual_features": int(len(residual)), "assignments_promoted": False}
        return pd.DataFrame([row], columns=report_columns), pd.DataFrame(columns=assignment_columns), residual

    candidates = [("primary_residual_baseline", crystal)]
    for path in candidate_paths:
        try:
            candidate = load_reflections(replace(
                config, cif_path=str(path), all230_compare=False, all230_policy="gemmi"
            ))
            candidates.append((_v73_candidate_kind(path, config), candidate))
        except Exception as exc:
            candidates.append(("load_error", {"path": path, "load_error": str(exc)}))

    rows = []
    assignment_tables = []
    local = replace(
        config,
        max_normal_candidates_anchor=int(config.v73_residual_phase_candidate_normals),
        v7_refine_hypotheses=int(config.v73_residual_phase_refine_hypotheses),
        coarse_tilt_values_deg=(-6.0, 0.0, 6.0),
        validation_holdout_fraction=0.25,
        validation_holdout_min_features=2,
        max_anchor_features=min(6, config.max_anchor_features),
        max_validation_features=min(8, config.max_validation_features),
        v7_max_fit_predictions=min(220, config.v7_max_fit_predictions),
        full_leave_one_angle_out=False,
        full_bootstrap_iterations=0,
        test_second_orientation=False,
    )
    for role, candidate in candidates:
        path = Path(candidate.get("path", "")) if isinstance(candidate, dict) else Path("")
        if role == "load_error":
            rows.append({
                "series_id": series_id, "status": "candidate_load_error",
                "candidate_role": role, "candidate_cif": str(candidate.get("path", "")),
                "residual_features": int(len(residual)),
            })
            continue
        try:
            _, peaks = simulate_powder_pattern(candidate, local)
            powder = _powder_peak_compatibility(residual, peaks, local)
            phase_search = v7_indexing_search(candidate, residual, local, fast=True)
            best = phase_search["best"]
            orientation_matches = best.get("all_feature_matches", pd.DataFrame()).copy()
            orientation_metrics = best.get("all_feature_metrics", _empty_assignment_metrics())
            composite = float(orientation_metrics.get("score", -1.0) + 0.25 * powder.get("score", -1.0))
            rows.append({
                "series_id": series_id, "status": "screened",
                "candidate_role": role, "candidate_cif": str(candidate["path"]),
                "spacegroup_number": int(candidate["spacegroup"].number),
                "residual_features": int(len(residual)),
                "powder_matched_features": int(powder.get("matched_features", 0)),
                "powder_matched_weight_fraction": float(powder.get("matched_weight_fraction", 0.0)),
                "powder_median_radial_residual_q": float(powder.get("median_radial_residual_q", np.nan)),
                "powder_score": float(powder.get("score", -1.0)),
                "orientation_normal_h": int(best["hkl"][0]),
                "orientation_normal_k": int(best["hkl"][1]),
                "orientation_normal_l": int(best["hkl"][2]),
                "orientation_matches": int(orientation_metrics.get("matches", 0)),
                "orientation_weighted_fraction": float(orientation_metrics.get("weighted_fraction", 0.0)),
                "orientation_score": float(orientation_metrics.get("score", -1.0)),
                "composite_score": composite,
                "assignments_promoted": False,
            })
            if not orientation_matches.empty:
                orientation_matches["candidate_cif"] = str(candidate["path"])
                orientation_matches["candidate_role"] = role
                assignment_tables.append(orientation_matches)
        except Exception as exc:
            rows.append({
                "series_id": series_id, "status": "screen_error",
                "candidate_role": role, "candidate_cif": str(candidate.get("path", path)),
                "residual_features": int(len(residual)), "composite_score": np.nan,
            })
    report = pd.DataFrame(rows, columns=report_columns)
    screened = report[report.status.eq("screened")].copy()
    if screened.empty:
        return report, pd.DataFrame(columns=assignment_columns), residual
    screened = screened.sort_values("composite_score", ascending=False, kind="mergesort").reset_index(drop=True)
    report = report.merge(
        screened[["candidate_cif"]].assign(candidate_rank=np.arange(1, len(screened) + 1)),
        on="candidate_cif", how="left",
    )
    primary_rows = screened[screened.candidate_role.eq("primary_residual_baseline")]
    primary_score = float(primary_rows.iloc[0].composite_score) if not primary_rows.empty else -np.inf
    report["score_gap_over_primary_residual"] = report.composite_score - primary_score
    report["preferred_over_primary_residual"] = (
            ~report.candidate_role.eq("primary_residual_baseline")
            & (report.score_gap_over_primary_residual >= float(config.v73_residual_phase_min_score_gap))
            & (report.orientation_matches >= int(config.v73_residual_phase_min_matches))
    )
    best_preferred = report[report.preferred_over_primary_residual.fillna(False)].sort_values(
        "composite_score", ascending=False
    )
    assignments = pd.concat(assignment_tables, ignore_index=True, sort=False) if assignment_tables else pd.DataFrame()
    if assignments.empty:
        return report, pd.DataFrame(columns=assignment_columns), residual
    rank_map = report.set_index("candidate_cif").candidate_rank.to_dict()
    preferred_map = report.set_index("candidate_cif").preferred_over_primary_residual.to_dict()
    assignments["candidate_rank"] = assignments.candidate_cif.map(rank_map)
    assignments["residual_candidate_preferred"] = assignments.candidate_cif.map(preferred_map).fillna(False)
    promote = bool(config.v73_promote_residual_candidate_matches) and not best_preferred.empty
    if promote:
        winning_path = str(best_preferred.iloc[0].candidate_cif)
        promoted = assignments[assignments.candidate_cif.eq(winning_path)].copy()
        promoted["role"] = "residual_candidate_phase"
        promoted["index_source"] = "residual_candidate_cif_screen"
        promoted["completion_is_primary_orientation_evidence"] = False
        promoted["completion_is_validation_evidence"] = False
        report.loc[report.candidate_cif.eq(winning_path), "assignments_promoted"] = True
    else:
        promoted = pd.DataFrame(columns=assignment_columns)
    assignments["role"] = "residual_candidate_phase_diagnostic"
    assignments["index_source"] = "residual_candidate_cif_screen_diagnostic"
    assignments["completion_is_primary_orientation_evidence"] = False
    assignments["completion_is_validation_evidence"] = False
    assignments = assignments.reindex(columns=assignment_columns, fill_value=np.nan)
    if promote:
        promoted = promoted.reindex(columns=assignment_columns, fill_value=np.nan)
        assignments = pd.concat([assignments, promoted], ignore_index=True, sort=False)
    return report, assignments, residual


def _local_pixel_measurement(image, qr_prediction, qz_prediction, config):
    if config.local_pixel_completion_png_only and "png" not in str(image.get("source_kind", "")):
        return None
    x_float, y_float = _q_to_pixel(image, [qr_prediction], [qz_prediction])
    if not np.isfinite(x_float[0]) or not np.isfinite(y_float[0]):
        return None
    x0, y0 = int(round(x_float[0])), int(round(y_float[0]))
    dqr = max(abs(float(image["qr"][1] - image["qr"][0])), 1e-8)
    dqz = max(abs(float(image["qz"][1] - image["qz"][0])), 1e-8)
    inner_x = max(2, int(math.ceil(config.local_pixel_window_q / dqr)))
    inner_y = max(2, int(math.ceil(config.local_pixel_window_q / dqz)))
    outer_x = max(inner_x + 2, int(math.ceil(config.local_pixel_background_window_q / dqr)))
    outer_y = max(inner_y + 2, int(math.ceil(config.local_pixel_background_window_q / dqz)))
    h, w = image["intensity"].shape
    if x0 - outer_x < 0 or x0 + outer_x >= w or y0 - outer_y < 0 or y0 + outer_y >= h:
        return None
    ys = slice(y0 - outer_y, y0 + outer_y + 1)
    xs = slice(x0 - outer_x, x0 + outer_x + 1)
    patch = np.asarray(image["intensity"][ys, xs], float)
    valid = np.asarray(image["valid"][ys, xs], bool) & np.isfinite(patch)
    yy, xx = np.indices(patch.shape)
    center_x, center_y = outer_x, outer_y
    inner = (np.abs(xx - center_x) <= inner_x) & (np.abs(yy - center_y) <= inner_y) & valid
    outer = valid & ~((np.abs(xx - center_x) <= inner_x + 1) & (np.abs(yy - center_y) <= inner_y + 1))
    if inner.sum() < 5 or outer.sum() < 20:
        return None
    background = patch[outer]
    median = float(np.median(background))
    mad = max(float(1.4826 * np.median(np.abs(background - median))), 1e-5)
    inner_values = np.where(inner, patch, -np.inf)
    peak_flat = int(np.argmax(inner_values))
    py, px = np.unravel_index(peak_flat, patch.shape)
    peak = float(patch[py, px])
    snr = (peak - median) / mad
    if not np.isfinite(snr):
        return None
    radius = 2
    y1, y2 = max(0, py - radius), min(patch.shape[0], py + radius + 1)
    x1, x2 = max(0, px - radius), min(patch.shape[1], px + radius + 1)
    local = patch[y1:y2, x1:x2]
    local_valid = valid[y1:y2, x1:x2]
    weights = np.where(local_valid, np.clip(local - median, 0.0, None), 0.0)
    weights = np.nan_to_num(weights, nan=0.0, posinf=0.0, neginf=0.0)
    if weights.sum() <= 0:
        return None
    local_yy, local_xx = np.indices(local.shape)
    x_pixel = (x0 - outer_x + x1) + float(np.sum(weights * local_xx) / weights.sum())
    y_pixel = (y0 - outer_y + y1) + float(np.sum(weights * local_yy) / weights.sum())
    qr = float(np.interp(x_pixel, np.arange(len(image["qr"])), image["qr"]))
    qz_axis = np.asarray(image["qz"], float)
    qz = float(np.interp(y_pixel, np.arange(len(qz_axis)), qz_axis))
    delta = math.hypot(qr - qr_prediction, qz - qz_prediction)
    return {"qr": qr, "qz": qz, "snr": float(snr), "delta_q": delta,
            "peak_intensity": peak, "background": median, "pixel_x": x_pixel, "pixel_y": y_pixel}


def registered_local_pixel_completion(crystal, records, search, secondary, existing_indexed, config):
    columns = [
        "feature_id", "hkl", "h", "k", "l", "qr_exp", "qz_exp", "qr_calc", "qz_calc",
        "delta_qr", "delta_qz", "delta_q", "normalized_delta", "strength", "support",
        "support_fraction", "feature_type", "f2", "prediction_weight", "assignment_ambiguity",
        "assignment_margin_sigma", "orientation_domain", "assignment_support_score", "role",
        "index_source", "angles", "scans", "member_angle_support", "mean_local_snr",
        "minimum_local_snr", "completion_is_primary_orientation_evidence",
        "completion_is_validation_evidence",
    ]
    diagnostics_columns = [
        "orientation_domain", "hkl", "qr_calc", "qz_calc", "angle_support", "angles", "scans",
        "mean_local_snr", "minimum_local_snr", "centroid_delta_q", "prediction_ambiguity",
        "prediction_margin_q", "status",
    ]
    if not config.enable_local_pixel_completion or not records:
        return pd.DataFrame(columns=columns), pd.DataFrame(columns=diagnostics_columns)
    domain_rows = _v73_domain_rows(search, secondary)
    prediction_sets = [_v7_prediction_array(
        crystal, tuple(row["hkl"]), config, np.asarray(row["params"], float),
        config.local_pixel_prediction_f2_percentile,
        config.local_pixel_max_predictions_per_domain, str(row["domain"]),
    ) for row in domain_rows]
    predictions = _v7_concat_predictions(prediction_sets)
    if not len(predictions.get("qr", [])):
        return pd.DataFrame(columns=columns), pd.DataFrame(columns=diagnostics_columns)
    used_keys = set()
    existing_coords = np.empty((0, 2), float)
    if existing_indexed is not None and not existing_indexed.empty:
        for row in existing_indexed.itertuples(index=False):
            used_keys.add(_v71_prediction_key(getattr(row, "orientation_domain", "primary"),
                                              row.h, row.k, row.l, row.qr_calc))
        existing_coords = existing_indexed[["qr_exp", "qz_exp"]].to_numpy(float)
    candidates, diagnostics = [], []
    prediction_coords = np.column_stack((predictions["qr"], predictions["qz"]))
    for j in range(len(predictions["qr"])):
        key = _v71_prediction_key(predictions["domain"][j], predictions["h"][j],
                                  predictions["k"][j], predictions["l"][j], predictions["qr"][j])
        if key in used_keys:
            continue
        qr_calc, qz_calc = float(predictions["qr"][j]), float(predictions["qz"][j])
        if not (config.qr_range[0] <= qr_calc <= config.qr_range[1]
                and config.qz_range[0] <= qz_calc <= config.qz_range[1]
                and qz_calc >= config.analysis_qz_min):
            continue
        detections = []
        for record in records:
            image = record.get("registered_image", record["image"])
            detected = _local_pixel_measurement(image, qr_calc, qz_calc, config)
            if detected is None:
                continue
            threshold = (config.local_pixel_min_three_angle_snr
                         if len(records) >= 3 else config.local_pixel_min_snr)
            if detected["snr"] < threshold or detected["delta_q"] > config.local_pixel_max_centroid_offset_q:
                continue
            detected.update({"angle_deg": float(record["angle_deg"]), "scan": int(record["scan"])})
            detections.append(detected)
        angle_support = len({d["angle_deg"] for d in detections})
        if angle_support < int(config.local_pixel_min_angle_support):
            continue
        weights = np.asarray([max(d["snr"], 0.1) for d in detections], float)
        qr_exp = float(np.average([d["qr"] for d in detections], weights=weights))
        qz_exp = float(np.average([d["qz"] for d in detections], weights=weights))
        centroid_delta = math.hypot(qr_exp - qr_calc, qz_exp - qz_calc)
        if centroid_delta > config.local_pixel_max_centroid_offset_q:
            continue
        if len(existing_coords) and np.min(np.linalg.norm(existing_coords - [qr_exp, qz_exp],
                                                          axis=1)) < config.local_pixel_existing_feature_separation_q:
            continue
        distances = np.linalg.norm(prediction_coords - [qr_exp, qz_exp], axis=1)
        ordered = np.sort(distances)
        ambiguity = int(np.sum(distances <= max(config.local_pixel_max_centroid_offset_q, ordered[0] + 0.006)))
        margin = float(ordered[1] - ordered[0]) if len(ordered) > 1 else np.inf
        if ambiguity > int(
                config.local_pixel_max_prediction_ambiguity) or margin < config.local_pixel_min_prediction_margin_q:
            status = "ambiguous_prediction"
        else:
            status = "accepted"
        diagnostics.append({
            "orientation_domain": str(predictions["domain"][j]),
            "hkl": f"({int(predictions['h'][j])} {int(predictions['k'][j])} {int(predictions['l'][j])})",
            "qr_calc": qr_calc, "qz_calc": qz_calc, "angle_support": angle_support,
            "angles": ",".join(f"{x:.3f}" for x in sorted({d['angle_deg'] for d in detections})),
            "scans": ",".join(str(x) for x in sorted({d['scan'] for d in detections})),
            "mean_local_snr": float(np.mean([d["snr"] for d in detections])),
            "minimum_local_snr": float(np.min([d["snr"] for d in detections])),
            "centroid_delta_q": centroid_delta, "prediction_ambiguity": ambiguity,
            "prediction_margin_q": margin, "status": status,
        })
        if status != "accepted":
            continue
        candidates.append({
            "hkl": f"({int(predictions['h'][j])} {int(predictions['k'][j])} {int(predictions['l'][j])})",
            "h": int(predictions["h"][j]), "k": int(predictions["k"][j]), "l": int(predictions["l"][j]),
            "qr_exp": qr_exp, "qz_exp": qz_exp, "qr_calc": qr_calc, "qz_calc": qz_calc,
            "delta_qr": qr_exp - qr_calc, "delta_qz": qz_exp - qz_calc, "delta_q": centroid_delta,
            "normalized_delta": centroid_delta / max(config.uncertainty_floor_q, 1e-9),
            "strength": float(np.mean([d["peak_intensity"] - d["background"] for d in detections])),
            "support": angle_support, "support_fraction": angle_support / max(len(records), 1),
            "feature_type": "local_pixel_peak", "f2": float(predictions["f2"][j]),
            "prediction_weight": float(predictions["prediction_weight"][j]), "assignment_ambiguity": ambiguity,
            "assignment_margin_sigma": margin / max(config.uncertainty_floor_q, 1e-9),
            "orientation_domain": str(predictions["domain"][j]),
            "assignment_support_score": float(
                angle_support * np.mean([d["snr"] for d in detections]) / (1 + centroid_delta / 0.01)),
            "role": "registered_local_pixel_completion", "index_source": "registered_local_pixel_completion",
            "angles": ",".join(f"{x:.3f}" for x in sorted({d['angle_deg'] for d in detections})),
            "scans": ",".join(str(x) for x in sorted({d['scan'] for d in detections})),
            "member_angle_support": angle_support, "mean_local_snr": float(np.mean([d["snr"] for d in detections])),
            "minimum_local_snr": float(np.min([d["snr"] for d in detections])),
            "completion_is_primary_orientation_evidence": False,
            "completion_is_validation_evidence": False,
        })
    if not candidates:
        return pd.DataFrame(columns=columns), pd.DataFrame(diagnostics).reindex(columns=diagnostics_columns)
    frame = pd.DataFrame(candidates).sort_values(
        ["member_angle_support", "assignment_support_score", "delta_q"],
        ascending=[False, False, True], kind="mergesort",
    )
    kept_indices = []
    kept_positions = []
    for index, row in frame.iterrows():
        position = (float(row.qr_exp), float(row.qz_exp))
        if all(math.hypot(position[0] - prior[0], position[1] - prior[1])
               > config.local_pixel_candidate_separation_q for prior in kept_positions):
            kept_indices.append(index)
            kept_positions.append(position)
    frame = frame.loc[kept_indices].reset_index(drop=True)
    frame.insert(0, "feature_id", [f"P{i + 1:03d}" for i in range(len(frame))])
    return frame.reindex(columns=columns), pd.DataFrame(diagnostics).reindex(columns=diagnostics_columns)


def _row_supported_by_record(row, members, record):
    feature_id = str(row.get("feature_id", ""))
    if members is not None and not members.empty and feature_id in set(members.feature_id.astype(str)):
        subset = members[members.feature_id.astype(str) == feature_id]
        return bool((subset.scan.astype(int) == int(record["scan"])).any())
    scans = str(row.get("scans", ""))
    if scans and scans.lower() != "nan":
        return str(int(record["scan"])) in {token.strip() for token in scans.split(",")}
    angles = str(row.get("angles", ""))
    if angles and angles.lower() != "nan":
        values = []
        for token in angles.split(","):
            try:
                values.append(float(token))
            except ValueError:
                pass
        return any(abs(value - float(record["angle_deg"])) < 5e-4 for value in values)
    return False


def _v739_postprocess_registered_series(series_id, result, records, crystal, config, output):
    series_dir = Path(output) / series_id.replace(":", "_series_")
    _post_start = time.perf_counter()
    print(f"[{series_id}] post-fit image completion and final overlay...", flush=True)
    composite = registered_composite_image(records, config)
    # Use the multi-angle visual composite for the canonical overlay. Its
    # display_intensity is unmasked, while its scientific intensity remains
    # masked for analysis.
    display_image = composite
    base_indexed, _ = _overlay_indexed_and_ignored(
        result["search"], result.get("ignored_salvage", pd.DataFrame()),
        result.get("guided_rescue", pd.DataFrame()), config, result.get("secondary"),
        result.get("completion_assignments", pd.DataFrame()), None,
    )
    local_assignments, local_diagnostics = registered_local_pixel_completion(
        crystal, records, result["search"], result.get("secondary"), base_indexed, config
    )
    print(f"[{series_id}] local-pixel completion finished in {time.perf_counter() - _post_start:.1f}s", flush=True)
    result["local_pixel_completion"] = local_assignments
    result["local_pixel_completion_diagnostics"] = local_diagnostics
    local_assignments.to_csv(series_dir / "registered_local_pixel_completion_assignments.csv", index=False)
    local_diagnostics.to_csv(series_dir / "registered_local_pixel_completion_diagnostics.csv", index=False)
    indexed, ignored = plot_indexed_ignored_overlay(
        display_image, result["search"], result.get("ignored_salvage", pd.DataFrame()),
        result.get("guided_rescue", pd.DataFrame()), series_dir / "indexed_or_ignored_overlay.png",
        f"{series_id} | indexed or ignored | normal {result['search']['best']['hkl']}",
        config, result.get("secondary"), result.get("completion_assignments", pd.DataFrame()),
        local_assignments,
    )
    if bool(getattr(config, "overlay_write_composite_companion", True)):
        plot_indexed_ignored_overlay(
            composite, result["search"], result.get("ignored_salvage", pd.DataFrame()),
            result.get("guided_rescue", pd.DataFrame()), series_dir / "registered_multi_angle_composite_overlay.png",
            f"{series_id} | registered multi-angle composite | normal {result['search']['best']['hkl']}",
            replace(config, overlay_classic_show_labels=False, overlay_static_style="classic"),
            result.get("secondary"), result.get("completion_assignments", pd.DataFrame()),
            local_assignments,
        )
    indexed.assign(overlay_class="indexed").to_csv(series_dir / "overlay_indexed_features.csv", index=False)
    ignored.assign(overlay_class="ignored").to_csv(series_dir / "overlay_ignored_features.csv", index=False)
    pd.concat([indexed.assign(overlay_class="indexed"), ignored.assign(overlay_class="ignored")],
              ignore_index=True, sort=False).to_csv(series_dir / "overlay_feature_classification.csv", index=False)
    shutil.copyfile(series_dir / "indexed_or_ignored_overlay.png", series_dir / "anchor_indexing_overlay.png")
    shutil.copyfile(series_dir / "indexed_or_ignored_overlay.png", series_dir / "postfit_rescue_evidence_overlay.png")
    if config.write_registered_per_angle_overlays:
        angle_dir = series_dir / "registered_per_angle_overlays";
        angle_dir.mkdir(exist_ok=True)
        members = result.get("members", pd.DataFrame())
        for record in sorted(records, key=lambda item: (item["angle_deg"], item["scan"])):
            indexed_mask = indexed.apply(lambda row: _row_supported_by_record(row, members, record), axis=1)
            ignored_mask = ignored.apply(lambda row: _row_supported_by_record(row, members, record), axis=1)
            angle_text = f"{float(record['angle_deg']):.3f}".replace(".", "p")
            output_path = angle_dir / f"angle_{angle_text}_scan_{int(record['scan'])}_indexed_or_ignored.png"
            if (str(getattr(config, "overlay_static_style", "classic")).strip().lower() == "classic"
                    and bool(getattr(config, "overlay_classic_per_angle", True))):
                _plot_classic_overlay_tables(
                    record.get("registered_image", record["image"]),
                    indexed[indexed_mask].copy(), ignored[ignored_mask].copy(),
                    output_path,
                    f"{series_id} | angle {float(record['angle_deg']):.3f} deg | scan {int(record['scan'])}",
                    config, dpi_override=min(int(config.dpi), 125),
                )
            else:
                _plot_fast_per_angle_overlay(
                    record.get("registered_image", record["image"]), indexed[indexed_mask].copy(),
                    ignored[ignored_mask].copy(), output_path,
                    f"{series_id} | angle {float(record['angle_deg']):.3f} deg | scan {int(record['scan'])}",
                    config,
                )
    print(f"[{series_id}] final static overlays finished in {time.perf_counter() - _post_start:.1f}s", flush=True)
    result["summary"].update({
        "overlay_indexed_features": int(len(indexed)), "overlay_ignored_features": int(len(ignored)),
        "registered_local_pixel_completion_features": int(len(local_assignments)),
        "registered_composite_overlay": True,
        "registered_per_angle_overlays": bool(config.write_registered_per_angle_overlays),
    })
    return result


def _coverage_prediction_key(domain, h, k, l, qr_value):
    sign = 0 if abs(float(qr_value)) < 1e-10 else int(math.copysign(1, float(qr_value)))
    return str(domain), int(h), int(k), int(l), sign


def _coverage_prediction_inventory(crystal, search, secondary, config):
    """Return all visible calculated reflections for every accepted domain."""
    solutions = [{
        "domain": "primary",
        "hkl": tuple(search["best"]["hkl"]),
        "params": np.asarray(search["best"]["params"], float),
    }]
    if secondary and secondary.get("accepted"):
        supplied = secondary.get("domain_solutions", []) or []
        for item in supplied:
            domain = str(item.get("domain", "secondary"))
            hkl = tuple(item.get("hkl", ()))
            params = np.asarray(item.get("params", []), float)
            if len(hkl) == 3 and params.size:
                solutions.append({"domain": domain, "hkl": hkl, "params": params})
        if not supplied and secondary.get("search"):
            best = secondary["search"]["best"]
            solutions.append({
                "domain": "secondary",
                "hkl": tuple(best["hkl"]),
                "params": np.asarray(best["params"], float),
            })

    frames = []
    seen = set()
    maximum = int(getattr(config, "v71_completion_max_predictions_per_domain", 1818))
    for item in solutions:
        identity = (item["domain"], item["hkl"], tuple(np.round(item["params"], 10)))
        if identity in seen:
            continue
        seen.add(identity)
        predictions = _prediction_array(
            crystal, item["hkl"], config, item["params"],
            f2_percentile=0.0, maximum=maximum,
        )
        frame = _predictions_frame(predictions)
        if frame.empty:
            continue
        frame["domain"] = item["domain"]
        frame["orientation_normal_hkl"] = str(tuple(item["hkl"]))
        frames.append(frame)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True, sort=False).drop_duplicates(
        ["domain", "h", "k", "l", "qr"], keep="first"
    ).reset_index(drop=True)


def _write_improved_coverage_outputs(
    series_dir, group, consensus, members, guided_rescue,
    guided_rescue_diagnostics, overlay_indexed, overlay_classes,
    crystal, search, secondary, config,
):
    """Write exhaustive experimental-feature and calculated-reflection classifications.

    The tables make the indexing decision auditable by recording raw detections
    rejected before consensus, indexed and ignored consensus features, and visible
    calculated reflections that were not assigned to an experimental peak.
    """
    series_dir = Path(series_dir)

    # Every consensus feature must be classified exactly once.
    consensus_base = consensus.drop_duplicates("feature_id", keep="first").copy()
    class_map = overlay_classes[["feature_id", "overlay_class"]].drop_duplicates(
        "feature_id", keep="first"
    ) if not overlay_classes.empty else pd.DataFrame(columns=["feature_id", "overlay_class"])
    exhaustive = consensus_base.merge(class_map, on="feature_id", how="left")
    exhaustive["overlay_class"] = exhaustive["overlay_class"].fillna("not_indexed")
    exhaustive["classification_status"] = np.where(
        exhaustive.overlay_class.eq("indexed"), "indexed", "not_indexed"
    )
    if exhaustive.feature_id.duplicated().any() or len(exhaustive) != len(consensus_base):
        raise RuntimeError("Consensus feature classification is not one-to-one.")
    exhaustive.to_csv(series_dir / "consensus_feature_classification_exhaustive.csv", index=False)

    # Raw detections that did not survive ordinary multi-angle consensus.
    raw = group.copy().reset_index(drop=True)
    if "raw_feature_id" not in raw:
        raw["raw_feature_id"] = [f"R{i + 1:05d}" for i in range(len(raw))]
    member_ids = set(
        members.get("raw_feature_id", pd.Series(dtype=str)).astype(str)
    ) if members is not None and not members.empty else set()
    raw_only = raw[~raw.raw_feature_id.astype(str).isin(member_ids)].copy()

    rescue_lookup = {}
    if guided_rescue is not None and not guided_rescue.empty:
        for row in guided_rescue.itertuples(index=False):
            tier = str(getattr(row, "salvage_evidence_tier", "provisional"))
            feature_id = str(getattr(row, "feature_id", ""))
            for raw_id in str(getattr(row, "raw_member_ids", "")).split(","):
                raw_id = raw_id.strip()
                if raw_id:
                    rescue_lookup[raw_id] = (feature_id, tier)
    raw_only["guided_rescue_feature_id"] = raw_only.raw_feature_id.astype(str).map(
        lambda value: rescue_lookup.get(value, ("", ""))[0]
    )
    raw_only["guided_rescue_evidence_tier"] = raw_only.raw_feature_id.astype(str).map(
        lambda value: rescue_lookup.get(value, ("", ""))[1]
    )
    promoted_tiers = set(getattr(config, "overlay_promoted_rescue_tiers", ("supported", "robust")))
    raw_only["coverage_status"] = np.where(
        raw_only.guided_rescue_evidence_tier.isin(promoted_tiers),
        "recovered_and_indexed_by_multi_angle_weak_peak_pass",
        np.where(
            raw_only.guided_rescue_feature_id.ne(""),
            "recovered_provisionally_not_promoted",
            "detected_but_not_consensus",
        ),
    )
    if guided_rescue_diagnostics is not None and not guided_rescue_diagnostics.empty:
        diagnostic_columns = [
            column for column in ("raw_feature_id", "guided_rescue_status")
            if column in guided_rescue_diagnostics
        ]
        if len(diagnostic_columns) == 2:
            raw_only = raw_only.merge(
                guided_rescue_diagnostics[diagnostic_columns].drop_duplicates("raw_feature_id"),
                on="raw_feature_id", how="left",
            )
    raw_only.to_csv(series_dir / "detected_but_not_consensus.csv", index=False)

    # Calculated reflections that are visible in the modeled q range but are not
    # used by any displayed experimental-to-calculated assignment.
    predictions = _coverage_prediction_inventory(crystal, search, secondary, config)
    used_keys = set()
    if overlay_indexed is not None and not overlay_indexed.empty:
        for row in overlay_indexed.itertuples(index=False):
            qr_value = getattr(row, "qr_calc", getattr(row, "qr_exp", 0.0))
            used_keys.add(_coverage_prediction_key(
                getattr(row, "orientation_domain", "primary"),
                row.h, row.k, row.l, qr_value,
            ))
    if predictions.empty:
        unused = predictions
    else:
        keys = [
            _coverage_prediction_key(row.domain, row.h, row.k, row.l, row.qr)
            for row in predictions.itertuples(index=False)
        ]
        unused = predictions[[key not in used_keys for key in keys]].copy()
        unused["calculated_status"] = "calculated_visible_but_not_assigned"
        sort_columns = [column for column in ("prediction_weight", "f2") if column in unused]
        if sort_columns:
            unused = unused.sort_values(sort_columns, ascending=False, kind="mergesort")
    unused.to_csv(series_dir / "unused_calculated_reflections.csv", index=False)

    return {
        "consensus_features_exhaustively_classified": int(len(exhaustive)),
        "detected_but_not_consensus_features": int(len(raw_only)),
        "weak_raw_peaks_promoted_to_indexed": int(
            raw_only.coverage_status.eq("recovered_and_indexed_by_multi_angle_weak_peak_pass").sum()
        ) if not raw_only.empty else 0,
        "unused_calculated_reflections": int(len(unused)),
    }


def _v88_analyze_single_series(series_id, group, representative_image, image_count, crystal, config, output):
    """Run every expensive indexing/validation phase for one series in a fresh process."""
    group = group.reset_index(drop=True).copy()
    # Refraction/DWBA predictions describe the shared series solution.  When
    # automatic GIWAXS physics is active, use the median measured incident angle
    # for that series rather than a hard-coded global angle.  With physics
    # disabled this replacement has no effect on the scientific result.
    if str(getattr(config, "giwaxs_physics_mode", "")).lower() == "auto" and "angle_deg" in group:
        series_angles = pd.to_numeric(group["angle_deg"], errors="coerce").to_numpy(float)
        series_angles = series_angles[np.isfinite(series_angles) & (series_angles > 0)]
        if series_angles.size:
            unique_angles = tuple(float(value) for value in np.unique(np.round(series_angles, 8)))
            config = replace(
                config, incidence_angle_deg=float(np.median(series_angles)),
                giwaxs_incidence_angles_deg=unique_angles,
            )
    if "raw_feature_id" not in group:
        group["raw_feature_id"] = [f"R{i + 1:05d}" for i in range(len(group))]
    consensus, members = build_consensus(group, config)
    search = v7_indexing_search(crystal, consensus, config, raw_features=group)
    # Evaluate residual orientation evidence before post-fit completion. For series
    # enabled by the configured sector-search prefixes, accepted residual domains
    # participate in global one-to-one reassignment rather than only visualization.
    base_secondary = second_orientation_test(crystal, consensus, search, config)
    secondary, sector_assignments, sector_diagnostics = v73_s1_sector_domain_search(
        series_id, crystal, consensus, search, base_secondary, config
    )
    completion_assignments, completion_diagnostics = v71_fixed_orientation_completion(
        crystal, consensus, members, search, secondary, config
    )
    preserved_pre_sector_completion = pd.DataFrame(columns=completion_assignments.columns)
    if (getattr(config, "v73_preserve_pre_sector_completion", True)
            and secondary is not base_secondary
            and int((secondary or {}).get("s1_sector_domain_count", 0)) > 0):
        legacy_completion, legacy_completion_diagnostics = v71_fixed_orientation_completion(
            crystal, consensus, members, search, base_secondary, config
        )
        joint_tables = []
        if secondary and secondary.get("accepted") and secondary.get("joint") is not None:
            joint_tables = [
                secondary["joint"].get("anchor_matches", pd.DataFrame()),
                secondary["joint"].get("validation_matches", pd.DataFrame()),
            ]
        joint_tables = [table for table in joint_tables if table is not None and not table.empty]
        occupied = pd.concat(joint_tables + (
            [completion_assignments] if completion_assignments is not None and not completion_assignments.empty else []),
                             ignore_index=True, sort=False) if joint_tables or (
                    completion_assignments is not None and not completion_assignments.empty) else pd.DataFrame()
        occupied_ids = set(occupied.feature_id.astype(str)) if not occupied.empty else set()
        occupied_keys = set()
        if not occupied.empty:
            for row in occupied.itertuples(index=False):
                if all(hasattr(row, name) for name in ("h", "k", "l", "qr_calc")):
                    occupied_keys.add(_v71_prediction_key(
                        getattr(row, "orientation_domain", "primary"), row.h, row.k, row.l, row.qr_calc
                    ))
        if legacy_completion is not None and not legacy_completion.empty:
            keep = []
            for row in legacy_completion.itertuples(index=False):
                key = _v71_prediction_key(
                    getattr(row, "orientation_domain", "primary"), row.h, row.k, row.l, row.qr_calc
                )
                keep.append(str(row.feature_id) not in occupied_ids and key not in occupied_keys)
                if keep[-1]:
                    occupied_ids.add(str(row.feature_id));
                    occupied_keys.add(key)
            preserved_pre_sector_completion = legacy_completion[np.asarray(keep, bool)].copy()
            if not preserved_pre_sector_completion.empty:
                preserved_pre_sector_completion["role"] = "pre_sector_completion_preserved"
                preserved_pre_sector_completion["index_source"] = "pre_sector_domain_completion_preserved"
                completion_assignments = pd.concat(
                    [completion_assignments, preserved_pre_sector_completion], ignore_index=True, sort=False
                ).drop_duplicates("feature_id", keep="first")
        if legacy_completion_diagnostics is not None and not legacy_completion_diagnostics.empty:
            legacy_completion_diagnostics = legacy_completion_diagnostics.copy()
            legacy_completion_diagnostics["completion_diagnostic_source"] = "pre_sector_domain_model"
            completion_diagnostics = pd.concat(
                [completion_diagnostics, legacy_completion_diagnostics], ignore_index=True, sort=False
            )
    mosaic_assignments, mosaic_diagnostics = v72_s1_mosaic_completion(
        series_id, crystal, consensus, members, search, secondary,
        completion_assignments, config
    )
    if mosaic_assignments is not None and not mosaic_assignments.empty:
        completion_assignments = pd.concat(
            [completion_assignments, mosaic_assignments], ignore_index=True, sort=False
        ).drop_duplicates("feature_id", keep="first")
    ignored_salvage, ignored_salvage_diagnostics = index_ignored_features(
        crystal, search, secondary, config, members=members
    )
    guided_rescue, guided_rescue_diagnostics = prediction_guided_raw_rescue(
        crystal, group, members, search, secondary, ignored_salvage, config
    )
    leave_one_out = full_leave_one_angle_out(crystal, group, config)
    bootstrap_runs, bootstrap_assignments, bootstrap_summary = full_search_bootstrap(
        crystal, group, search, config
    )
    matches = combined_matches(search)
    if secondary and secondary.get("accepted") and secondary.get("joint") is not None:
        joint_tables = [
            secondary["joint"].get("anchor_matches", pd.DataFrame()),
            secondary["joint"].get("validation_matches", pd.DataFrame()),
        ]
        joint_tables = [table for table in joint_tables if table is not None and not table.empty]
        if joint_tables:
            matches = pd.concat(joint_tables, ignore_index=True, sort=False)
            matches["role"] = np.where(
                matches.feature_id.isin(set(search["validation"].feature_id)),
                "validation", "joint_all_feature",
            )
            matches = matches.drop_duplicates("feature_id", keep="first")
    if completion_assignments is not None and not completion_assignments.empty:
        matches = pd.concat([matches, completion_assignments], ignore_index=True, sort=False)
        matches = matches.drop_duplicates("feature_id", keep="first")
    if not matches.empty and not bootstrap_assignments.empty:
        matches = matches.merge(
            bootstrap_assignments, on=["feature_id", "h", "k", "l"], how="left"
        )
    residual_phase_report, residual_phase_assignments, residual_unexplained = (
        v73_residual_candidate_phase_screen(series_id, crystal, consensus, matches, config)
    )
    promoted_phase = residual_phase_assignments[
        residual_phase_assignments.get(
            "role", pd.Series(index=residual_phase_assignments.index, dtype=str)
        ).eq("residual_candidate_phase")
    ].copy() if residual_phase_assignments is not None and not residual_phase_assignments.empty else pd.DataFrame()
    if not promoted_phase.empty:
        matches = pd.concat([matches, promoted_phase], ignore_index=True, sort=False).drop_duplicates(
            "feature_id", keep="first"
        )
        completion_assignments = pd.concat(
            [completion_assignments, promoted_phase], ignore_index=True, sort=False
        ).drop_duplicates("feature_id", keep="first")
    series_dir = output / series_id.replace(":", "_series_")
    series_dir.mkdir(parents=True, exist_ok=True)
    roles = pd.concat([
        search["anchors"].assign(role="anchor"),
        search["validation"].assign(role="validation"),
        search["ignored"].assign(role="ignored"),
    ], ignore_index=True)
    roles.to_csv(series_dir / "feature_roles.csv", index=False)
    search["anchors"].to_csv(series_dir / "anchor_features.csv", index=False)
    search["validation"].to_csv(series_dir / "validation_features.csv", index=False)
    search["ignored"].to_csv(series_dir / "ignored_features.csv", index=False)
    search["ranking"].to_csv(series_dir / "ranked_orientation_hypotheses.csv", index=False)
    search["best"].get("angle_consistency_by_image", pd.DataFrame()).to_csv(
        series_dir / "angle_by_angle_consistency.csv", index=False
    )
    search.get("orientation_families", pd.DataFrame()).to_csv(
        series_dir / "orientation_family_diagnostics.csv", index=False
    )
    search.get("ambiguity_set", pd.DataFrame()).to_csv(
        series_dir / "orientation_ambiguity_set.csv", index=False
    )
    search["best"]["anchor_matches"].to_csv(series_dir / "anchor_assignments.csv", index=False)
    search["best"]["validation_matches"].to_csv(series_dir / "validation_assignments.csv", index=False)
    completion_assignments.to_csv(series_dir / "fixed_orientation_completion_assignments.csv", index=False)
    completion_diagnostics.to_csv(series_dir / "fixed_orientation_completion_diagnostics.csv", index=False)
    mosaic_assignments.to_csv(series_dir / "s1_mosaic_completion_assignments.csv", index=False)
    mosaic_diagnostics.to_csv(series_dir / "s1_mosaic_completion_diagnostics.csv", index=False)
    sector_assignments.to_csv(series_dir / "s1_sector_domain_assignments.csv", index=False)
    sector_diagnostics.to_csv(series_dir / "s1_sector_domain_diagnostics.csv", index=False)
    preserved_pre_sector_completion.to_csv(
        series_dir / "s1_pre_sector_completion_preserved.csv", index=False
    )
    residual_phase_report.to_csv(series_dir / "residual_candidate_phase_screening.csv", index=False)
    residual_phase_assignments.to_csv(series_dir / "residual_candidate_phase_assignments.csv", index=False)
    residual_unexplained.to_csv(series_dir / "residual_unexplained_features.csv", index=False)
    matches.to_csv(series_dir / "indexed_reflections.csv", index=False)
    ignored_salvage.to_csv(
        series_dir / "ignored_feature_assignments_provisional.csv", index=False
    )
    ignored_salvage_diagnostics.to_csv(
        series_dir / "ignored_feature_salvage_diagnostics.csv", index=False
    )
    supported_ignored = ignored_salvage[
        ignored_salvage.get("salvage_evidence_tier", pd.Series(index=ignored_salvage.index, dtype=str)).isin(
            ["robust", "supported"])
    ].copy() if not ignored_salvage.empty else pd.DataFrame()
    robust_ignored = ignored_salvage[
        ignored_salvage.get("salvage_evidence_tier", pd.Series(index=ignored_salvage.index, dtype=str)).eq("robust")
    ].copy() if not ignored_salvage.empty else pd.DataFrame()
    supported_ignored.to_csv(series_dir / "ignored_feature_assignments_supported.csv", index=False)
    robust_ignored.to_csv(series_dir / "ignored_feature_assignments_robust.csv", index=False)
    guided_rescue.to_csv(series_dir / "prediction_guided_raw_rescue_assignments.csv", index=False)
    guided_rescue_diagnostics.to_csv(series_dir / "prediction_guided_raw_rescue_diagnostics.csv", index=False)
    indexed_with_salvage = pd.concat(
        [table for table in (matches, ignored_salvage, guided_rescue) if not table.empty],
        ignore_index=True, sort=False,
    ) if (not matches.empty or not ignored_salvage.empty or not guided_rescue.empty) else pd.DataFrame()
    all_rescue_path = series_dir / "indexed_reflections_with_all_postfit_rescue.csv"
    indexed_with_salvage.to_csv(all_rescue_path, index=False)
    # Write an alternate overlay filename for notebook workflows that expect the
    # same byte-identical indexed/ignored classification image under that name.
    shutil.copyfile(all_rescue_path, series_dir / "indexed_reflections_with_provisional_ignored.csv")
    indexed_supported = pd.concat(
        [table for table in (matches, supported_ignored, guided_rescue[
            guided_rescue.get("salvage_evidence_tier", pd.Series(index=guided_rescue.index, dtype=str)).isin(
                ["supported", "robust"])
        ] if not guided_rescue.empty else pd.DataFrame()) if not table.empty],
        ignore_index=True, sort=False,
    ) if (not matches.empty or not supported_ignored.empty or not guided_rescue.empty) else pd.DataFrame()
    indexed_supported.to_csv(
        series_dir / "indexed_reflections_with_supported_postfit_rescue.csv", index=False
    )
    leave_one_out.to_csv(series_dir / "leave_one_angle_out_full_search.csv", index=False)
    bootstrap_runs.to_csv(series_dir / "bootstrap_full_search_runs.csv", index=False)
    bootstrap_assignments.to_csv(series_dir / "bootstrap_assignment_probabilities.csv", index=False)
    if secondary and secondary.get("search"):
        secondary["search"]["ranking"].to_csv(
            series_dir / "second_orientation_hypotheses.csv", index=False
        )
        if secondary.get("joint") is not None:
            secondary["joint"]["anchor_matches"].to_csv(
                series_dir / "joint_two_domain_anchor_assignments.csv", index=False
            )
            secondary["joint"]["validation_matches"].to_csv(
                series_dir / "joint_two_domain_validation_assignments.csv", index=False
            )
            pd.concat([
                secondary["joint"].get("anchor_matches", pd.DataFrame()),
                secondary["joint"].get("validation_matches", pd.DataFrame()),
            ], ignore_index=True, sort=False).drop_duplicates("feature_id", keep="first").to_csv(
                series_dir / "joint_multidomain_assignments.csv", index=False
            )
            pd.DataFrame(secondary.get("domain_solutions", [])).to_csv(
                series_dir / "fitted_orientation_domains.csv", index=False
            )
    # Fresh subprocess workers perform the scientific search and write tables,
    # but defer expensive overlay rendering to the parent process. The parent
    # has the complete registered image set and will replace these preliminary
    # classifications after local-pixel completion. This avoids concurrent
    # Matplotlib rendering and large memory spikes.
    if os.environ.get("GIXS_SERIES_WORKER") == "1":
        overlay_indexed, overlay_ignored = _overlay_indexed_and_ignored(
            search, ignored_salvage, guided_rescue, config, secondary,
            completion_assignments,
        )
    else:
        overlay_indexed, overlay_ignored = plot_indexed_ignored_overlay(
            representative_image, search, ignored_salvage, guided_rescue,
            series_dir / "indexed_or_ignored_overlay.png",
            f"{series_id} | indexed or ignored | normal {search['best']['hkl']}",
            config, secondary, completion_assignments,
        )
    overlay_indexed.assign(overlay_class="indexed").to_csv(
        series_dir / "overlay_indexed_features.csv", index=False
    )
    overlay_ignored.assign(overlay_class="ignored").to_csv(
        series_dir / "overlay_ignored_features.csv", index=False
    )
    overlay_classes = pd.concat([
        overlay_indexed.assign(overlay_class="indexed"),
        overlay_ignored.assign(overlay_class="ignored"),
    ], ignore_index=True, sort=False)
    overlay_classes.to_csv(series_dir / "overlay_feature_classification.csv", index=False)
    coverage_counts = _write_improved_coverage_outputs(
        series_dir, group, consensus, members, guided_rescue,
        guided_rescue_diagnostics, overlay_indexed, overlay_classes,
        crystal, search, secondary, config,
    )
    # Isolated workers return scientific results to the parent process, which renders
    # registered overlays. In-process execution can create equivalent output aliases
    # immediately because the renderer is already available locally.
    if os.environ.get("GIXS_SERIES_WORKER") != "1":
        overlay_path = series_dir / "indexed_or_ignored_overlay.png"
        shutil.copyfile(overlay_path, series_dir / "anchor_indexing_overlay.png")
        shutil.copyfile(overlay_path, series_dir / "postfit_rescue_evidence_overlay.png")
    best = search["best"]
    summary = {
        "series_id": series_id, "images": image_count,
        "giwaxs_physics_status": str(getattr(config, "giwaxs_physics_status", "")),
        "giwaxs_refraction_enabled": bool(config.enable_refraction_position_correction),
        "giwaxs_dwba_enabled": bool(config.enable_dwba),
        "giwaxs_incidence_angles_deg": list(getattr(config, "giwaxs_incidence_angles_deg", ()) or (config.incidence_angle_deg,)),
        "top_ranked_normal_hkl": best["hkl"],
        "normal_hkl": best["hkl"],
        "orientation_output_is_hypothesis": True,
        "heuristic_score_gap": search["score_margin"],
        "score_margin": search["score_margin"],
        "heuristic_values_are_statistical_confidence": False,
        "orientation_ambiguity_set": search.get("ambiguity_set", pd.DataFrame())[
            ["normal_h", "normal_k", "normal_l"]].values.tolist() if not search.get("ambiguity_set",
                                                                                    pd.DataFrame()).empty else [],
        "fallback_anchors": search["fallback_anchors"],
        "anchors": len(search["anchors"]), "validation_features": len(search["validation"]),
        "ignored_features": len(search["ignored"]),
        "overlay_indexed_features": int(len(overlay_indexed)),
        "overlay_ignored_features": int(len(overlay_ignored)),
        "consensus_features_exhaustively_classified": coverage_counts["consensus_features_exhaustively_classified"],
        "detected_but_not_consensus_features": coverage_counts["detected_but_not_consensus_features"],
        "weak_raw_peaks_promoted_to_indexed": coverage_counts["weak_raw_peaks_promoted_to_indexed"],
        "unused_calculated_reflections": coverage_counts["unused_calculated_reflections"],
        "overlay_promoted_rescue_tiers": list(config.overlay_promoted_rescue_tiers),
        "fixed_orientation_completion_features": int(len(completion_assignments)),
        "s1_mosaic_completion_features": int(len(mosaic_assignments)),
        "s1_mosaic_completion_enabled": bool(config.v72_enable_s1_mosaic_completion),
        "s1_sector_domains_added": int((secondary or {}).get("s1_sector_domain_count", 0)) if secondary else 0,
        "s1_sector_domain_assignments": int(len(sector_assignments)),
        "s1_sector_new_features": int(
            sector_assignments.get(
                "sector_domain_is_new_feature", pd.Series(False, index=sector_assignments.index)
            ).fillna(False).astype(bool).sum()
        ) if sector_assignments is not None and not sector_assignments.empty else 0,
        "s1_pre_sector_completion_preserved": int(len(preserved_pre_sector_completion)),
        "residual_candidate_cifs_tested": int(
            residual_phase_report.get(
                "status", pd.Series(index=residual_phase_report.index, dtype=str)
            ).eq("screened").sum()
        ) if residual_phase_report is not None and not residual_phase_report.empty else 0,
        "residual_candidate_promoted_features": int(len(promoted_phase)),
        "completion_is_primary_orientation_evidence": False,
        "completion_is_validation_evidence": False,
        "provisionally_indexed_ignored_features": len(ignored_salvage),
        "supported_indexed_ignored_features": int(len(supported_ignored)),
        "robust_indexed_ignored_features": int(len(robust_ignored)),
        "prediction_guided_raw_rescue_features": int(len(guided_rescue)),
        "supported_prediction_guided_raw_rescue_features": int(
            (guided_rescue.get("salvage_evidence_tier",
                               pd.Series(index=guided_rescue.index, dtype=str)) == "supported").sum()
        ) if not guided_rescue.empty else 0,
        "ignored_salvage_is_primary_orientation_evidence": False,
        "v7_all_feature_ignored_candidates_are_primary_orientation_evidence": True,
        "anchor_matches": best["anchor_metrics"]["matches"],
        "all_feature_matches": best.get("all_feature_metrics", {}).get("matches", best["anchor_metrics"]["matches"]),
        "all_feature_weighted_fraction": best.get("all_feature_metrics", {}).get("weighted_fraction", np.nan),
        "angle_consistency_score": best.get("angle_consistency_score", np.nan),
        "validation_matches": best["validation_metrics"]["matches"],
        "anchor_median_delta_q": best["anchor_metrics"]["median_delta_q"],
        "validation_median_delta_q": best["validation_metrics"]["median_delta_q"],
        "leave_one_angle_out_weighted_fraction": (
            float(leave_one_out.predicted_weighted_fraction.mean()) if not leave_one_out.empty else np.nan
        ),
        "bootstrap_orientation_stability": bootstrap_summary["orientation_stability"],
        "bootstrap_completed_iterations": bootstrap_summary["completed_iterations"],
        "second_orientation_accepted": bool(secondary and secondary.get("accepted")),
        "fitted_domain_count": int((secondary or {}).get("domain_count", 1)) if secondary else 1,
        "joint_two_domain_score_gain": (
            (secondary.get("joint") or {}).get("score_gain", np.nan)
            if secondary else np.nan
        ),
        "tilt_x_deg": best["params"][0], "tilt_y_deg": best["params"][1],
        "qr_scale": _unpack(best["params"])[2], "qz_scale": _unpack(best["params"])[3],
        "common_q_scale": float(np.mean(_unpack(best["params"])[2:4])),
        "qr_offset": _unpack(best["params"])[4], "qz_offset": _unpack(best["params"])[5],
        "axis_rotation_deg": float(best["params"][6]) if len(best["params"]) >= 8 else 0.0,
        "q_space_shear": float(best["params"][7]) if len(best["params"]) >= 8 else 0.0,
        "refraction_position_correction_enabled": bool(config.enable_refraction_position_correction),
    }
    print(
        f"\n[{series_id}] normal={best['hkl']} | anchors "
        f"{best['anchor_metrics']['matches']}/{len(search['anchors'])} | validation "
        f"{best['validation_metrics']['matches']}/{len(search['validation'])} | "
        f"margin={search['score_margin']:.4f} | anchor median="
        f"{best['anchor_metrics']['median_delta_q']:.4f} | validation median="
        f"{best['validation_metrics']['median_delta_q']:.4f} | ignored salvage "
        f"{len(ignored_salvage)}/{len(search['ignored'])} "
        f"(supported={len(supported_ignored)}, robust={len(robust_ignored)}) | "
        f"guided raw rescue={len(guided_rescue)} | full-reflection completion="
        f"{len(completion_assignments)} (s1 mosaic={len(mosaic_assignments)})", flush=True,
    )
    print(search["ranking"].head(config.top_candidates_to_print).to_string(index=False), flush=True)
    return {
        "consensus": consensus, "members": members, "search": search, "indexed": matches,
        "leave_one_out": leave_one_out, "bootstrap_runs": bootstrap_runs,
        "bootstrap_assignments": bootstrap_assignments,
        "bootstrap_summary": bootstrap_summary, "secondary": secondary,
        "ignored_salvage": ignored_salvage,
        "ignored_salvage_diagnostics": ignored_salvage_diagnostics,
        "guided_rescue": guided_rescue,
        "guided_rescue_diagnostics": guided_rescue_diagnostics,
        "completion_assignments": completion_assignments,
        "completion_diagnostics": completion_diagnostics,
        "mosaic_assignments": mosaic_assignments,
        "mosaic_diagnostics": mosaic_diagnostics,
        "sector_assignments": sector_assignments,
        "sector_diagnostics": sector_diagnostics,
        "preserved_pre_sector_completion": preserved_pre_sector_completion,
        "residual_phase_report": residual_phase_report,
        "residual_phase_assignments": residual_phase_assignments,
        "residual_unexplained": residual_unexplained,
        "summary": summary,
    }


def _resolve_series_worker_module_file(output):
    """Return a physical module path for fresh workers, including notebooks.

    Normal script execution uses ``__file__``.  In Jupyter/IPython, the code may
    have been pasted into a cell, so no file exists.  In that case, search recent
    input history for the complete GIXS source and materialize it as a temporary
    Python module.  The worker environment sets ``GIXS_WORKER_IMPORT=1``, so the
    automatic top-level analysis does not run while importing that module.
    """
    # Notebook kernels may retain a previously materialized worker-module path after
    # source changes. Prefer the physical script when available; otherwise rebuild
    # the worker module from recent notebook input so worker configuration classes
    # and analysis functions match the code currently being executed.
    module_file = globals().get("__file__")
    if module_file:
        candidate = Path(str(module_file)).expanduser()
        if candidate.is_file():
            resolved = str(candidate.resolve())
            globals()["_GIXS_RESOLVED_WORKER_MODULE_FILE"] = resolved
            return resolved

    try:
        ip = get_ipython()  # noqa: F821 - defined only in IPython/Jupyter
    except Exception:
        ip = None
    if ip is None:
        return None

    sources = []
    try:
        history = list(ip.user_ns.get("In", []))
        sources.extend(reversed(history[-50:]))
    except Exception:
        pass
    try:
        raw_history = list(ip.history_manager.input_hist_raw)
        sources.extend(reversed(raw_history[-50:]))
    except Exception:
        pass

    source = None
    for candidate_source in sources:
        if not isinstance(candidate_source, str):
            continue
        if (
                "class IndexingConfig" in candidate_source
                and "def run_gixs_indexing" in candidate_source
                and "def _analyze_single_series" in candidate_source
        ):
            source = candidate_source
            break
    if source is None:
        return None

    # Comment out notebook magics/shell escapes so the materialized source is a
    # valid importable Python file. Ordinary Python lines are preserved exactly.
    cleaned_lines = []
    for line in source.splitlines():
        stripped = line.lstrip()
        if stripped.startswith("%") or stripped.startswith("!"):
            cleaned_lines.append("# notebook-only command removed: " + line)
        else:
            cleaned_lines.append(line)
    cleaned_source = "\n".join(cleaned_lines) + "\n"

    import hashlib
    digest = hashlib.sha256(cleaned_source.encode("utf-8")).hexdigest()[:16]
    work_dir = Path(output) / "_series_worker_cache"
    work_dir.mkdir(parents=True, exist_ok=True)
    materialized = work_dir / f"gixs_notebook_worker_{digest}.py"
    if not materialized.is_file() or materialized.read_text(encoding="utf-8") != cleaned_source:
        temporary = materialized.with_suffix(".tmp")
        temporary.write_text(cleaned_source, encoding="utf-8")
        temporary.replace(materialized)
    resolved = str(materialized.resolve())
    globals()["_GIXS_RESOLVED_WORKER_MODULE_FILE"] = resolved
    return resolved


# ======================== ROBUST ORIENTATION SOLVER ========================
# Score all usable consensus features with an explicit outlier option, refine a
# bounded affine q-space calibration, and add orientation domains only when joint
# reassignment improves the robust whole-pattern objective.

def _v7_unpack(params=None):
    if params is None:
        return (0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0)
    values = tuple(map(float, params))
    if len(values) == 5:
        tx, ty, scale, offset_r, offset_z = values
        return tx, ty, scale, scale, offset_r, offset_z, 0.0, 0.0
    if len(values) == 6:
        return (*values, 0.0, 0.0)
    if len(values) != 8:
        raise ValueError(f"Expected 5, 6, or 8 calibration parameters, received {len(values)}")
    return values


def _v7_affine_coordinates(qr, qz, params):
    _, _, scale_r, scale_z, offset_r, offset_z, axis_rotation_deg, shear = _v7_unpack(params)
    theta = math.radians(axis_rotation_deg)
    rotation = np.array([[math.cos(theta), -math.sin(theta)],
                         [math.sin(theta), math.cos(theta)]], float)
    linear = rotation @ np.array([[scale_r, shear], [0.0, scale_z]], float)
    coordinates = np.column_stack((np.asarray(qr, float), np.asarray(qz, float))) @ linear.T
    coordinates[:, 0] += offset_r
    coordinates[:, 1] += offset_z
    return coordinates


def _v7_prediction_array(crystal, normal_hkl, config: IndexingConfig, params=None,
                         f2_percentile=0.0, maximum=None, domain="primary"):
    tx, ty, _, _, _, _, _, _ = _v7_unpack(params)
    arrays = _crystal_arrays(crystal)
    rotation = _tilt_matrix(tx, ty) @ _base_alignment(crystal, normal_hkl)
    lab = arrays["g"] @ rotation.T
    qxy_all = np.hypot(lab[:, 0], lab[:, 1])
    qz_all = refraction_corrected_qz(lab[:, 2], config)
    affine_margin = 0.12
    visible = (
            (qz_all >= config.analysis_qz_min - config.max_anchor_q_offset - affine_margin)
            & (qz_all <= config.qz_range[1] + config.max_anchor_q_offset + affine_margin)
            & (qxy_all <= max(abs(x) for x in config.qr_range) + config.max_anchor_q_offset + affine_margin)
    )
    indices = np.flatnonzero(visible)
    if not len(indices):
        return _empty_predictions()
    if f2_percentile:
        threshold = float(np.quantile(arrays["f2"][indices], f2_percentile / 100.0))
        indices = indices[arrays["f2"][indices] >= threshold]
    ranking_intensity = arrays["f2"][indices] * dwba_intensity_envelope(qz_all[indices], config)
    if maximum and len(indices) > maximum:
        selected = np.argpartition(ranking_intensity, -maximum)[-maximum:]
        indices = indices[selected]
        ranking_intensity = ranking_intensity[selected]
    if not len(indices):
        return _empty_predictions()
    order = np.argsort(ranking_intensity)[::-1]
    indices = indices[order]
    positive_count = len(indices)
    negative_mask = qxy_all[indices] > 1e-8
    source = np.concatenate((np.arange(positive_count), np.flatnonzero(negative_mask)))
    signs = np.concatenate((np.ones(positive_count), -np.ones(int(negative_mask.sum()))))
    raw_qr = signs * qxy_all[indices][source]
    raw_qz = qz_all[indices][source]
    transformed = _v7_affine_coordinates(raw_qr, raw_qz, params)
    qr, qz = transformed[:, 0], transformed[:, 1]
    inside = (
            (qr >= config.qr_range[0]) & (qr <= config.qr_range[1])
            & (qz >= config.qz_range[0]) & (qz <= config.qz_range[1])
    )
    source, raw_qr, raw_qz, qr, qz = source[inside], raw_qr[inside], raw_qz[inside], qr[inside], qz[inside]
    if not len(source):
        return _empty_predictions()
    ref_indices = indices[source]
    # Keep strongest reflection in each near-identical projected bin.
    keys = np.column_stack((np.rint(qr / 0.0035).astype(np.int32),
                            np.rint(qz / 0.0035).astype(np.int32)))
    _, first = np.unique(keys, axis=0, return_index=True)
    first = np.sort(first)
    ref_indices, raw_qr, raw_qz, qr, qz = (
        ref_indices[first], raw_qr[first], raw_qz[first], qr[first], qz[first]
    )
    f2 = arrays["f2"][ref_indices]
    dwba_weight = dwba_intensity_envelope(qz, config)
    effective_intensity = f2 * dwba_weight
    median = max(float(np.median(effective_intensity)), 1e-12)
    weights = np.log1p(effective_intensity / median)
    weights /= max(float(weights.max()), 1e-12)
    return {
        "h": arrays["h"][ref_indices], "k": arrays["k"][ref_indices], "l": arrays["l"][ref_indices],
        "hkl": arrays["hkl"][ref_indices], "q": arrays["q"][ref_indices], "d": arrays["d"][ref_indices],
        "f2": f2, "qxy": np.abs(raw_qr), "qz": qz, "qr": qr,
        "dwba_weight": dwba_weight, "effective_intensity": effective_intensity,
        "prediction_weight": weights,
        "domain": np.full(len(qr), str(domain), dtype=object),
    }


# Make every downstream diagnostic/salvage projection use the same affine model.
def _prediction_array(crystal, normal_hkl, config: IndexingConfig, params=None,
                      f2_percentile=0.0, maximum=None):
    return _v7_prediction_array(crystal, normal_hkl, config, params, f2_percentile, maximum, "primary")


def _v7_feature_quality(features: pd.DataFrame, config: IndexingConfig):
    if features.empty:
        return np.empty(0), np.empty(0)
    sigma = np.sqrt(np.maximum(features.cov_rr.to_numpy(float) + features.cov_zz.to_numpy(float), 1e-12))
    sigma_factor = np.exp(-np.maximum(sigma - 0.018, 0.0) / 0.045)
    width = features.major_width_q.to_numpy(float)
    width_factor = np.exp(-np.maximum(width - 0.070, 0.0) / 0.090)
    support = features.support.to_numpy(float)
    support_factor = np.sqrt(np.clip(support / max(float(np.nanmax(support)), 1.0), 0.0, 1.0))
    type_factor = np.ones(len(features), float)
    types = features.feature_type.astype(str).to_numpy()
    type_factor[np.isin(types, list(config.ignored_feature_types))] *= config.v7_ignored_feature_weight
    type_factor[types != "spot"] *= config.v7_broad_feature_weight
    quality = np.clip(sigma_factor * width_factor * (0.45 + 0.55 * support_factor) * type_factor,
                      config.v7_min_feature_quality, 1.0)
    raw_weight = np.power(np.clip(features.strength.to_numpy(float), 1e-4, 1.0), 0.70) * np.sqrt(
        np.maximum(support, 1.0)) * quality
    return quality, raw_weight


def _v7_assign_arrays(features, predictions, config: IndexingConfig, tolerance, materialize=True):
    predictions = _predictions_dict(predictions)
    if features.empty or not len(predictions["qr"]):
        return pd.DataFrame(), _empty_assignment_metrics()
    physical, normalized = _feature_distance_matrices(features, predictions, config)
    quality, raw_weight = _v7_feature_quality(features, config)
    feature_weight = raw_weight / max(float(raw_weight.sum()), 1e-12)
    gate = (physical <= tolerance) & (normalized <= config.match_sigma_limit)
    ambiguity = np.maximum(gate.sum(axis=1), 1)
    cost = (
            normalized / max(config.match_sigma_limit, 1e-9)
            + 0.07 * (1.0 - predictions["prediction_weight"])[None, :]
            + 0.035 * np.log1p(ambiguity - 1)[:, None]
    )
    locks = {feature: (h, k, l) for feature, h, k, l in config.manual_locked_assignments}
    for i, feature_id in enumerate(features.feature_id.astype(str)):
        if feature_id in locks:
            hkl = locks[feature_id]
            allowed = ((predictions["h"] == hkl[0]) & (predictions["k"] == hkl[1]) & (predictions["l"] == hkl[2]))
            cost[i, ~allowed] = 1e6
    dummy_cost = config.v7_outlier_cost + 0.20 * quality
    if materialize:
        real = np.nan_to_num(cost.copy(), nan=1e4, posinf=1e4, neginf=-1e4)
        real[~gate] = 1e4
        dummy = np.full((len(features), len(features)), 1e4, float)
        np.fill_diagonal(dummy, dummy_cost)
        rows, cols = linear_sum_assignment(np.hstack((real, dummy)))
        accepted = cols < len(predictions["qr"])
        rows, cols = rows[accepted], cols[accepted]
        accepted = gate[rows, cols] & (cost[rows, cols] <= dummy_cost[rows])
        rows, cols = rows[accepted], cols[accepted]
    else:
        eligible = gate & (cost <= dummy_cost[:, None])
        rows, cols = _greedy_one_to_one(cost, eligible)
    selected_norm = normalized[rows, cols] if len(rows) else np.empty(0)
    selected_phys = physical[rows, cols] if len(rows) else np.empty(0)
    utility = np.zeros(len(features), float)
    if len(rows):
        utility[rows] = (
                np.exp(-0.5 * (selected_norm / 1.65) ** 2)
                * (0.62 + 0.38 * predictions["prediction_weight"][cols])
                / np.sqrt(ambiguity[rows])
        )
    matched_utility = float(np.sum(feature_weight * utility))
    matched_weight = float(feature_weight[rows].sum()) if len(rows) else 0.0
    unmatched_weight = 1.0 - matched_weight
    ambiguity_cost = float(np.sum(feature_weight[rows] * np.log1p(ambiguity[rows] - 1))) if len(rows) else 1.0
    selected_predictions = np.unique(cols)
    false_fraction = (
            1.0 - float(predictions["prediction_weight"][selected_predictions].sum())
            / max(float(predictions["prediction_weight"].sum()), 1e-12)
    ) if len(cols) else 1.0
    score = matched_utility - 0.16 * unmatched_weight - 0.045 * ambiguity_cost - 0.025 * false_fraction
    exp_coords = features[["qr", "qz"]].to_numpy(float)
    calc_coords = np.column_stack((predictions["qr"], predictions["qz"]))
    geometry = np.nan
    if len(rows) >= 3:
        geometry = float(np.sqrt(np.mean((pdist(exp_coords[rows]) - pdist(calc_coords[cols])) ** 2)))
        score -= 0.05 * min(geometry / max(tolerance, 1e-9), 2.0)
    metrics = {
        "score": float(score), "matches": int(len(rows)), "weighted_fraction": matched_weight,
        "indexed_fraction": float(len(rows) / len(features)),
        "median_delta_q": float(np.median(selected_phys)) if len(rows) else np.nan,
        "p90_delta_q": float(np.quantile(selected_phys, 0.9)) if len(rows) else np.nan,
        "ambiguity": float(np.mean(ambiguity[rows])) if len(rows) else np.nan,
        "false_prediction_fraction": float(false_fraction), "pair_geometry_rms": geometry,
    }
    if not materialize or not len(rows):
        return pd.DataFrame(), metrics
    output = []
    for i, j in zip(rows, cols):
        exp = features.iloc[i]
        possible = normalized[i][gate[i]]
        ordered = np.sort(possible)
        margin = float(ordered[1] - ordered[0]) if len(ordered) > 1 else config.match_sigma_limit
        support_score = (
                quality[i] * math.exp(-0.5 * normalized[i, j] ** 2)
                * predictions["prediction_weight"][j]
                * min(1.0, margin / 1.5) / math.sqrt(ambiguity[i])
        )
        output.append({
            "feature_id": exp.feature_id, "hkl": predictions["hkl"][j],
            "h": int(predictions["h"][j]), "k": int(predictions["k"][j]), "l": int(predictions["l"][j]),
            "qr_exp": float(exp.qr), "qz_exp": float(exp.qz),
            "qr_calc": float(predictions["qr"][j]), "qz_calc": float(predictions["qz"][j]),
            "delta_qr": float(exp.qr - predictions["qr"][j]), "delta_qz": float(exp.qz - predictions["qz"][j]),
            "delta_q": float(physical[i, j]), "normalized_delta": float(normalized[i, j]),
            "strength": float(exp.strength), "support": int(exp.support), "feature_type": str(exp.feature_type),
            "experimental_integrated_intensity": float(getattr(exp, "experimental_integrated_intensity", np.nan)),
            "experimental_intensity_sigma": float(getattr(exp, "experimental_intensity_sigma", np.nan)),
            "experimental_integrated_snr": float(getattr(exp, "experimental_integrated_snr", np.nan)),
            "experimental_intensity_quality_score": float(getattr(exp, "experimental_intensity_quality_score", np.nan)),
            "feature_quality": float(quality[i]), "f2": float(predictions["f2"][j]),
            "prediction_weight": float(predictions["prediction_weight"][j]),
            "assignment_ambiguity": int(ambiguity[i]), "assignment_margin_sigma": margin,
            "dwba_weight": float(predictions.get("dwba_weight", np.ones(len(predictions["qr"])))[j]),
            "effective_intensity": float(predictions.get("effective_intensity", predictions["f2"])[j]),
            "orientation_domain": str(
                predictions.get("domain", np.full(len(predictions["qr"]), "primary", dtype=object))[j]),
            "assignment_support_score": float(support_score),
            "assignment_support_is_probability": False,
        })
    return pd.DataFrame(output).sort_values("assignment_support_score", ascending=False), metrics


def _v7_bounds(config):
    return np.array([
        [-config.max_tilt_anchor_deg, config.max_tilt_anchor_deg],
        [-config.max_tilt_anchor_deg, config.max_tilt_anchor_deg],
        [1 - config.max_qr_scale_change, 1 + config.max_qr_scale_change],
        [1 - config.max_qz_scale_change, 1 + config.max_qz_scale_change],
        [-config.max_anchor_q_offset, config.max_anchor_q_offset],
        [-config.max_anchor_q_offset, config.max_anchor_q_offset],
        [-config.v7_max_axis_rotation_deg, config.v7_max_axis_rotation_deg],
        [-config.v7_max_shear, config.v7_max_shear],
    ], float)


def _v7_regularization(params, config):
    tx, ty, sr, sz, or_, oz, axis, shear = _v7_unpack(params)
    normalized = np.array([
        tx / max(config.max_tilt_anchor_deg, 1e-9),
        ty / max(config.max_tilt_anchor_deg, 1e-9),
        (sr - 1) / max(config.max_qr_scale_change, 1e-9),
        (sz - 1) / max(config.max_qz_scale_change, 1e-9),
        or_ / max(config.max_anchor_q_offset, 1e-9),
        oz / max(config.max_anchor_q_offset, 1e-9),
        axis / max(config.v7_max_axis_rotation_deg, 1e-9),
        shear / max(config.v7_max_shear, 1e-9),
    ])
    base = config.anchor_regularization * float(np.sum(normalized[:6] ** 2))
    affine = config.v7_affine_regularization * float(np.sum(normalized[6:] ** 2))
    return base + affine


def _v7_prepare_raw_features(raw_features):
    if raw_features is None or raw_features.empty:
        return pd.DataFrame()
    raw = raw_features.copy()
    if "raw_feature_id" not in raw:
        raw["raw_feature_id"] = [f"RAW{i + 1:05d}" for i in range(len(raw))]
    raw["feature_id"] = raw["raw_feature_id"].astype(str)
    raw["support"] = 1
    raw["support_fraction"] = 1.0
    return raw


def _v7_angle_consistency(raw_features, predictions, config):
    raw = _v7_prepare_raw_features(raw_features)
    if raw.empty or "angle_deg" not in raw:
        return np.nan, pd.DataFrame()
    rows = []
    for angle, group in raw.groupby("angle_deg"):
        # Keep the best localized detections per angle so background ridges do not
        # dominate the consistency score.
        group = group.copy()
        quality, weight = _v7_feature_quality(group, config)
        group["_rank"] = weight
        group = group.nlargest(min(80, len(group)), "_rank").drop(columns="_rank")
        _, metrics = _v7_assign_arrays(group, predictions, config, config.v7_all_feature_tolerance_q, materialize=False)
        rows.append({"angle_deg": float(angle), **metrics})
    table = pd.DataFrame(rows)
    if table.empty:
        return np.nan, table
    consistency = float(0.70 * table.weighted_fraction.mean() + 0.30 * table.weighted_fraction.min())
    return consistency, table


def _v7_evaluate(normal_hkl, params, crystal, fit_features, anchors, validation, config,
                 materialize=True, raw_features=None):
    broad = _v7_prediction_array(
        crystal, normal_hkl, config, params,
        config.v7_fit_f2_percentile, config.v7_max_fit_predictions, "primary",
    )
    all_matches, all_metrics = _v7_assign_arrays(
        fit_features, broad, config, config.v7_all_feature_tolerance_q, materialize=materialize
    )
    strong = _v7_prediction_array(
        crystal, normal_hkl, config, params,
        config.hypothesis_f2_percentile, config.max_hypothesis_predictions, "primary",
    )
    anchor_matches, anchor_metrics = _v7_assign_arrays(
        anchors, strong, config, config.anchor_match_tolerance_q, materialize=materialize
    ) if not anchors.empty else (pd.DataFrame(), _empty_assignment_metrics(score=0.0))
    validation_matches, validation_metrics = _v7_assign_arrays(
        validation, broad, config, config.validation_match_tolerance_q, materialize=materialize
    ) if validation is not None and not validation.empty else (pd.DataFrame(), _empty_assignment_metrics(score=0.0))
    strict = _empty_assignment_metrics(score=0.0)
    if validation is not None and not validation.empty:
        _, strict = _v7_assign_arrays(
            validation, broad, config,
            min(config.strict_validation_tolerance_q, config.validation_match_tolerance_q),
            materialize=False,
        )
    angle_score, angle_table = (np.nan, pd.DataFrame())
    if raw_features is not None and materialize:
        angle_score, angle_table = _v7_angle_consistency(raw_features, broad, config)
    hard_penalty = 0.30 * max(0, max(config.min_anchor_matches, 4) - all_metrics["matches"])
    score = (
            config.v7_all_feature_score_weight * all_metrics["score"]
            + config.v7_validation_score_weight * validation_metrics["score"]
            + (config.v7_angle_consistency_weight * angle_score if np.isfinite(angle_score) else 0.0)
            - _v7_regularization(params, config) - hard_penalty
    )
    return {
        "hkl": tuple(map(int, normal_hkl)), "params": np.asarray(params, float),
        "anchor_predictions": _predictions_frame(strong) if materialize else pd.DataFrame(),
        "predictions": _predictions_frame(broad) if materialize else pd.DataFrame(),
        "anchor_matches": anchor_matches, "validation_matches": validation_matches,
        "all_feature_matches": all_matches,
        "anchor_metrics": anchor_metrics, "validation_metrics": validation_metrics,
        "all_feature_metrics": all_metrics, "strict_validation_metrics": strict,
        "angle_consistency_score": angle_score, "angle_consistency_by_image": angle_table,
        "cv_score": float(score),
    }


def _v7_project_fixed_assignments(normal_hkl, params, matches, crystal, config):
    if matches.empty:
        return np.empty((0, 2), float)
    tx, ty, _, _, _, _, _, _ = _v7_unpack(params)
    rotation = _tilt_matrix(tx, ty) @ _base_alignment(crystal, normal_hkl)
    hkl = matches[["h", "k", "l"]].to_numpy(float)
    lab = (hkl @ crystal["basis"].T) @ rotation.T
    sign = np.where(matches.qr_calc.to_numpy(float) < 0, -1.0, 1.0)
    raw_qr = sign * np.hypot(lab[:, 0], lab[:, 1])
    raw_qz = refraction_corrected_qz(lab[:, 2], config)
    return _v7_affine_coordinates(raw_qr, raw_qz, params)


def _v7_refine(solution, crystal, fit_features, anchors, validation, config, raw_features=None, fast=False):
    bounds = _v7_bounds(config)
    lower, upper = bounds[:, 0], bounds[:, 1]
    params = np.clip(np.asarray(_v7_unpack(solution["params"]), float), lower, upper)
    current = _v7_evaluate(solution["hkl"], params, crystal, fit_features, anchors, validation, config, True, None)
    best = current
    cycles = 2 if fast else max(2, int(config.v7_refine_cycles))
    lookup = fit_features.set_index("feature_id")
    for _ in range(cycles):
        matches = current.get("all_feature_matches", pd.DataFrame())
        if len(matches) < 3:
            break
        exp = matches[["qr_exp", "qz_exp"]].to_numpy(float)
        scales = []
        for feature_id in matches.feature_id:
            row = lookup.loc[feature_id]
            covariance = np.array([[row.cov_rr, row.cov_rz], [row.cov_rz, row.cov_zz]], float)
            eigen = np.linalg.eigvalsh(covariance)
            sigma = math.sqrt(
                float(np.clip(eigen.mean(), config.uncertainty_floor_q ** 2, config.uncertainty_ceiling_q ** 2)))
            quality = float(matches.loc[matches.feature_id == feature_id, "feature_quality"].iloc[
                                0]) if "feature_quality" in matches else 1.0
            scales.append(sigma / math.sqrt(max(quality, config.v7_min_feature_quality)))
        scales = np.asarray(scales, float)[:, None]

        def residual(vector):
            calc = _v7_project_fixed_assignments(solution["hkl"], vector, matches, crystal, config)
            data = ((calc - exp) / scales).ravel()
            tx, ty, sr, sz, or_, oz, axis, shear = _v7_unpack(vector)
            regularization = np.sqrt(np.maximum([
                config.anchor_regularization, config.anchor_regularization,
                config.anchor_regularization, config.anchor_regularization,
                config.anchor_regularization, config.anchor_regularization,
                config.v7_affine_regularization, config.v7_affine_regularization,
            ], 1e-12)) * np.array([
                tx / max(config.max_tilt_anchor_deg, 1e-9),
                ty / max(config.max_tilt_anchor_deg, 1e-9),
                (sr - 1) / max(config.max_qr_scale_change, 1e-9),
                (sz - 1) / max(config.max_qz_scale_change, 1e-9),
                or_ / max(config.max_anchor_q_offset, 1e-9),
                oz / max(config.max_anchor_q_offset, 1e-9),
                axis / max(config.v7_max_axis_rotation_deg, 1e-9),
                shear / max(config.v7_max_shear, 1e-9),
            ])
            return np.concatenate((data, regularization))

        fitted = least_squares(
            residual, params, bounds=(lower, upper), method="trf", loss="soft_l1", f_scale=1.0,
            max_nfev=45 if fast else 85, xtol=2e-5, ftol=2e-5, gtol=2e-5,
        )
        params = fitted.x
        current = _v7_evaluate(solution["hkl"], params, crystal, fit_features, anchors, validation, config, True, None)
        if current["cv_score"] > best["cv_score"]:
            best = current
    # The angle-by-angle term is evaluated only once after refinement.
    final = _v7_evaluate(solution["hkl"], best["params"], crystal, fit_features, anchors, validation, config, True,
                         raw_features)
    return final if final["cv_score"] >= best["cv_score"] - 0.02 else best


def v7_indexing_search(crystal, consensus, config: IndexingConfig, candidate_override=None,
                       fast=False, validation_mode=False, raw_features=None):
    anchors, validation, ignored = feature_roles(consensus, crystal, config)
    rejected = set(config.manual_rejected_feature_ids)
    fit_features = consensus[
        ~consensus.feature_id.isin(set(validation.feature_id) | rejected)
    ].copy()
    if fit_features.empty:
        fit_features = consensus[~consensus.feature_id.isin(rejected)].copy()
    candidates = candidate_override or orientation_candidates(
        crystal, replace(config, max_orientation_candidates=config.max_normal_candidates_anchor)
    )
    if validation_mode:
        candidates = candidates[:min(len(candidates), max(8, config.max_normal_candidates_anchor))]
    seeds = []
    for normal_hkl in candidates:
        best = None
        for tx in config.coarse_tilt_values_deg:
            for ty in config.coarse_tilt_values_deg:
                params = np.array([tx, ty, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0])
                evaluated = _v7_evaluate(
                    normal_hkl, params, crystal, fit_features, anchors, validation,
                    config, materialize=False, raw_features=None,
                )
                key = (evaluated["all_feature_metrics"]["matches"], evaluated["cv_score"])
                if best is None or key > (best["all_feature_metrics"]["matches"], best["cv_score"]):
                    best = evaluated
        if best is not None:
            seeds.append(best)
    seeds.sort(key=lambda item: (item["all_feature_metrics"]["matches"], item["cv_score"]), reverse=True)
    if not seeds:
        raise RuntimeError("V7 could not evaluate any orientation hypothesis")
    refine_count = min(
        3 if validation_mode else (6 if fast else config.v7_refine_hypotheses),
        len(seeds),
    )
    refined = [
        _v7_refine(seed, crystal, fit_features, anchors, validation, config,
                   raw_features=raw_features, fast=(fast or validation_mode))
        for seed in seeds[:refine_count]
    ]
    family_table, ranked = orientation_family_analysis(refined, crystal, config)
    if not ranked:
        ranked = sorted(refined, key=lambda item: item["cv_score"], reverse=True)
        family_table = pd.DataFrame()
    best = ranked[0]
    margin = float(best["cv_score"] - ranked[1]["cv_score"]) if len(ranked) > 1 else np.nan
    ranking_rows = []
    for rank, solution in enumerate(ranked, 1):
        ranking_rows.append({
            "rank": rank, "normal_h": solution["hkl"][0], "normal_k": solution["hkl"][1],
            "normal_l": solution["hkl"][2], "heuristic_cross_validated_score": solution["cv_score"],
            "all_feature_matches": solution["all_feature_metrics"]["matches"],
            "all_feature_weighted_fraction": solution["all_feature_metrics"]["weighted_fraction"],
            "all_feature_median_delta_q": solution["all_feature_metrics"]["median_delta_q"],
            "anchor_matches": solution["anchor_metrics"]["matches"],
            "validation_matches": solution["validation_metrics"]["matches"],
            "validation_weighted_fraction": solution["validation_metrics"]["weighted_fraction"],
            "validation_median_delta_q": solution["validation_metrics"]["median_delta_q"],
            "angle_consistency_score": solution.get("angle_consistency_score", np.nan),
            "axis_rotation_deg": float(solution["params"][6]),
            "q_space_shear": float(solution["params"][7]),
            "score_is_statistical_confidence": False,
        })
    alternatives = [{
        "hkl": item["hkl"], "params": item["params"], "cv_score": item["cv_score"],
        "anchor_metrics": item["anchor_metrics"], "validation_metrics": item["validation_metrics"],
        "all_feature_metrics": item["all_feature_metrics"],
    } for item in ranked]
    return {
        "anchors": anchors, "validation": validation, "ignored": ignored,
        "fit_features": fit_features, "best": best, "alternatives": alternatives,
        "ranking": pd.DataFrame(ranking_rows), "orientation_families": family_table,
        "ambiguity_set": family_table[
            family_table.in_ambiguity_set].copy() if not family_table.empty else pd.DataFrame(),
        "score_margin": margin, "heuristic_score_gap": margin,
        "fallback_anchors": False, "solver": "v7_robust_all_feature_affine",
    }


def anchor_first_search(crystal, consensus, config: IndexingConfig, candidate_override=None,
                        fast=False, validation_mode=False):
    if not getattr(config, "v7_enable_core_solver", True):
        raise RuntimeError("The V7 file requires v7_enable_core_solver=True")
    return v7_indexing_search(
        crystal, consensus, config, candidate_override=candidate_override,
        fast=fast, validation_mode=validation_mode, raw_features=None,
    )


def combined_matches(search):
    tables = []
    all_matches = search.get("best", {}).get("all_feature_matches", pd.DataFrame()).copy()
    if not all_matches.empty:
        anchor_ids = set(search.get("anchors", pd.DataFrame()).feature_id) if not search.get("anchors",
                                                                                             pd.DataFrame()).empty else set()
        ignored_ids = set(search.get("ignored", pd.DataFrame()).feature_id) if not search.get("ignored",
                                                                                              pd.DataFrame()).empty else set()
        all_matches["role"] = np.where(
            all_matches.feature_id.isin(anchor_ids), "anchor",
            np.where(all_matches.feature_id.isin(ignored_ids), "all_feature_recovered", "fit_feature")
        )
        tables.append(all_matches)
    validation = search.get("best", {}).get("validation_matches", pd.DataFrame()).copy()
    if not validation.empty:
        validation["role"] = "validation"
        tables.append(validation)
    if not tables:
        return pd.DataFrame()
    result = pd.concat(tables, ignore_index=True, sort=False)
    return result.drop_duplicates(["feature_id"], keep="first")


def _v7_concat_predictions(prediction_sets):
    dictionaries = [_predictions_dict(item) for item in prediction_sets if len(_predictions_dict(item)["qr"])]
    if not dictionaries:
        return _empty_predictions()
    keys = dictionaries[0].keys()
    return {key: np.concatenate([np.asarray(item[key]) for item in dictionaries]) for key in keys}


def _v7_domain_solution(normal_hkl, params, crystal, config, domain):
    predictions = _v7_prediction_array(
        crystal, normal_hkl, config, params,
        config.v7_fit_f2_percentile, config.v7_max_fit_predictions, domain,
    )
    return {"hkl": tuple(normal_hkl), "params": np.asarray(params, float),
            "predictions": _predictions_frame(predictions), "prediction_array": predictions,
            "domain": domain}


def _v7_refine_domain_tilts(domain_solution, assigned, crystal, config):
    if assigned.empty or len(assigned) < 2:
        return domain_solution
    base_params = np.asarray(domain_solution["params"], float).copy()
    exp = assigned[["qr_exp", "qz_exp"]].to_numpy(float)

    def residual(tilts):
        params = base_params.copy();
        params[:2] = tilts
        calc = _v7_project_fixed_assignments(domain_solution["hkl"], params, assigned, crystal, config)
        return ((calc - exp) / max(config.uncertainty_floor_q, 0.012)).ravel()

    fit = least_squares(
        residual, base_params[:2],
        bounds=(-config.max_tilt_anchor_deg, config.max_tilt_anchor_deg),
        loss="soft_l1", max_nfev=40,
    )
    base_params[:2] = fit.x
    return _v7_domain_solution(domain_solution["hkl"], base_params, crystal, config, domain_solution["domain"])


def second_orientation_test(crystal, consensus, primary, config):
    if not config.test_second_orientation or config.v7_multidomain_max_domains < 2:
        return None
    features = consensus[~consensus.feature_id.isin(config.manual_rejected_feature_ids)].copy()
    if len(features) < config.v7_multidomain_min_matches * 2:
        return {"accepted": False, "reason": "too_few_features_for_multidomain", "domain_count": 1}
    primary_solution = _v7_domain_solution(
        primary["best"]["hkl"], primary["best"]["params"], crystal, config, "primary"
    )
    domains = [primary_solution]
    combined_predictions = primary_solution["prediction_array"]
    current_matches, current_metrics = _v7_assign_arrays(
        features, combined_predictions, config, config.v7_all_feature_tolerance_q, True
    )
    single_score = float(current_metrics["score"])
    candidates = orientation_candidates(
        crystal, replace(config, max_orientation_candidates=config.v7_multidomain_candidate_normals)
    )
    accepted_records = []
    domain_names = ["secondary", "tertiary"]
    for domain_index in range(1, min(config.v7_multidomain_max_domains, 3)):
        domain_name = domain_names[domain_index - 1]
        used_normals = [item["hkl"] for item in domains]
        residual_ids = set(features.feature_id) - set(current_matches.feature_id if not current_matches.empty else [])
        residual_count = len(residual_ids)
        if residual_count < config.v7_multidomain_min_matches:
            break
        best_trial = None
        for normal_hkl in candidates:
            if any(normal_angle(normal_hkl, used, crystal) < config.second_orientation_min_normal_separation_deg for
                   used in used_normals):
                continue
            for tx in config.v7_multidomain_tilt_values_deg:
                for ty in config.v7_multidomain_tilt_values_deg:
                    _, _, sr, sz, or_, oz, axis, shear = _v7_unpack(primary["best"]["params"])
                    params = np.array([tx, ty, sr, sz, or_, oz, axis, shear], float)
                    trial_domain = _v7_domain_solution(normal_hkl, params, crystal, config, domain_name)
                    trial_predictions = _v7_concat_predictions([combined_predictions, trial_domain["prediction_array"]])
                    trial_matches, trial_metrics = _v7_assign_arrays(
                        features, trial_predictions, config, config.v7_all_feature_tolerance_q, True
                    )
                    domain_matches = trial_matches[
                        trial_matches.orientation_domain == domain_name] if not trial_matches.empty else pd.DataFrame()
                    gain = float(
                        trial_metrics["score"] - current_metrics["score"] - config.v7_multidomain_complexity_penalty)
                    key = (gain, len(domain_matches), trial_metrics["weighted_fraction"])
                    if best_trial is None or key > best_trial["key"]:
                        best_trial = {"key": key, "domain": trial_domain, "matches": trial_matches,
                                      "metrics": trial_metrics, "domain_matches": domain_matches,
                                      "predictions": trial_predictions, "gain": gain}
        if best_trial is None:
            break
        refined_domain = _v7_refine_domain_tilts(
            best_trial["domain"], best_trial["domain_matches"], crystal, config
        )
        refined_predictions = _v7_concat_predictions([combined_predictions, refined_domain["prediction_array"]])
        refined_matches, refined_metrics = _v7_assign_arrays(
            features, refined_predictions, config, config.v7_all_feature_tolerance_q, True
        )
        refined_domain_matches = refined_matches[
            refined_matches.orientation_domain == domain_name] if not refined_matches.empty else pd.DataFrame()
        refined_gain = float(
            refined_metrics["score"] - current_metrics["score"] - config.v7_multidomain_complexity_penalty)
        if refined_gain >= best_trial["gain"]:
            best_trial.update({"domain": refined_domain, "matches": refined_matches,
                               "metrics": refined_metrics, "domain_matches": refined_domain_matches,
                               "predictions": refined_predictions, "gain": refined_gain})
        if (best_trial["gain"] < config.v7_multidomain_min_gain
                or len(best_trial["domain_matches"]) < config.v7_multidomain_min_matches):
            break
        domains.append(best_trial["domain"])
        combined_predictions = best_trial["predictions"]
        current_matches, current_metrics = best_trial["matches"], best_trial["metrics"]
        accepted_records.append(best_trial)
    if len(domains) == 1:
        return {
            "accepted": False, "reason": "no_multidomain_score_gain", "domain_count": 1,
            "candidate_score_gain": accepted_records[-1]["gain"] if accepted_records else np.nan,
        }
    validation_ids = set(primary.get("validation", pd.DataFrame()).feature_id)
    joint_validation = current_matches[current_matches.feature_id.isin(validation_ids)].copy()
    joint_fit = current_matches[~current_matches.feature_id.isin(validation_ids)].copy()
    non_primary = current_matches[current_matches.orientation_domain != "primary"]
    first_secondary = domains[1]
    secondary_matches = current_matches[current_matches.orientation_domain == "secondary"].copy()
    secondary_search_best = {
        "hkl": first_secondary["hkl"], "params": first_secondary["params"],
        "predictions": first_secondary["predictions"], "anchor_predictions": first_secondary["predictions"],
        "anchor_matches": secondary_matches[~secondary_matches.feature_id.isin(validation_ids)].copy(),
        "validation_matches": secondary_matches[secondary_matches.feature_id.isin(validation_ids)].copy(),
        "anchor_metrics": {**_empty_assignment_metrics(score=0.0),
                           "matches": int((~secondary_matches.feature_id.isin(validation_ids)).sum())},
        "validation_metrics": {**_empty_assignment_metrics(score=0.0),
                               "matches": int(secondary_matches.feature_id.isin(validation_ids).sum())},
        "strict_validation_metrics": _empty_assignment_metrics(score=0.0),
        "cv_score": float(current_metrics["score"]),
    }
    ranking = pd.DataFrame([{
        "rank": index, "domain": item["domain"], "normal_h": item["hkl"][0],
        "normal_k": item["hkl"][1], "normal_l": item["hkl"][2]
    } for index, item in enumerate(domains[1:], 1)])
    search = {
        "anchors": primary.get("anchors", pd.DataFrame()), "validation": primary.get("validation", pd.DataFrame()),
        "ignored": primary.get("ignored", pd.DataFrame()), "best": secondary_search_best,
        "alternatives": [], "ranking": ranking, "score_margin": np.nan, "fallback_anchors": False,
    }
    joint = {
        "score": float(current_metrics["score"]),
        "score_gain": float(
            current_metrics["score"] - single_score - config.v7_multidomain_complexity_penalty * (len(domains) - 1)),
        "anchor_matches": joint_fit, "validation_matches": joint_validation,
        "anchor_metrics": current_metrics, "validation_metrics": current_metrics,
        "strict_validation_metrics": _empty_assignment_metrics(score=0.0),
        "secondary_anchor_matches": int(len(non_primary[~non_primary.feature_id.isin(validation_ids)])),
        "secondary_validation_matches": int(len(non_primary[non_primary.feature_id.isin(validation_ids)])),
        "domain_count": len(domains),
        "domain_solutions": [{"domain": item["domain"], "hkl": item["hkl"], "params": item["params"]} for item in
                             domains],
    }
    return {"accepted": True, "search": search, "joint": joint, "reason": "accepted_v7_joint_multidomain",
            "domain_count": len(domains), "domain_solutions": joint["domain_solutions"]}


def _domain_orientation_solution(search, secondary, domain):
    if str(domain) == "primary":
        return search["best"]
    if secondary and secondary.get("accepted"):
        for item in secondary.get("domain_solutions", []):
            if item.get("domain") == str(domain):
                return {"hkl": tuple(item["hkl"]), "params": np.asarray(item["params"], float)}
        if str(domain) == "secondary" and secondary.get("search"):
            return secondary["search"]["best"]
    return None


def _project_exact_reflection(crystal, normal_hkl, config, params, h, k, l, qr_sign):
    arrays = _crystal_arrays(crystal)
    mask = (arrays["h"] == int(h)) & (arrays["k"] == int(k)) & (arrays["l"] == int(l))
    indices = np.flatnonzero(mask)
    if not len(indices):
        return np.nan, np.nan
    index = int(indices[0])
    tx, ty, _, _, _, _, _, _ = _v7_unpack(params)
    rotation = _tilt_matrix(tx, ty) @ _base_alignment(crystal, normal_hkl)
    lab = arrays["g"][index] @ rotation.T
    raw_qr = float(math.copysign(np.hypot(lab[0], lab[1]), qr_sign if qr_sign else 1.0))
    raw_qz = float(refraction_corrected_qz(np.array([lab[2]]), config)[0])
    point = _v7_affine_coordinates(np.array([raw_qr]), np.array([raw_qz]), params)[0]
    return float(point[0]), float(point[1])


def _overlay_indexed_and_ignored(search, ignored_salvage, guided_rescue, config, secondary=None, completion=None,
                                 extra_indexed=None):
    """Build the binary overlay table with one row per consensus experimental feature.

    When a multidomain model is supported, its joint one-to-one assignments are
    used so an experimental feature is not displayed under multiple domain labels.
    Indexed plus ignored rows therefore account for the full consensus feature set.
    """
    if secondary and secondary.get("accepted") and secondary.get("joint") is not None:
        tables = [secondary["joint"].get("anchor_matches", pd.DataFrame()),
                  secondary["joint"].get("validation_matches", pd.DataFrame())]
        indexed = pd.concat([x for x in tables if x is not None and not x.empty],
                            ignore_index=True, sort=False) if any(
            x is not None and not x.empty for x in tables) else pd.DataFrame()
        if not indexed.empty:
            indexed["index_source"] = "joint_multidomain_index"
    else:
        indexed = combined_matches(search).copy()
        if not indexed.empty:
            indexed["index_source"] = "robust_all_feature_index"
    promoted_tiers = set(config.overlay_promoted_rescue_tiers)
    extras = []
    for table, source in ((ignored_salvage, "supported_ignored_rescue"),
                          (guided_rescue, "supported_raw_rescue")):
        if table is None or table.empty:
            continue
        tier = table.get("salvage_evidence_tier", pd.Series("provisional", index=table.index, dtype=str)).astype(str)
        promoted = table[tier.isin(promoted_tiers)].copy()
        if not promoted.empty:
            promoted["index_source"] = source
            extras.append(promoted)
    if completion is not None and not completion.empty:
        completed = completion.copy()
        completed["index_source"] = completed.get(
            "index_source", pd.Series("full_reflection_fixed_orientation_completion", index=completed.index)
        )
        extras.append(completed)
    if extra_indexed is not None and not extra_indexed.empty:
        additional = extra_indexed.copy()
        additional["index_source"] = additional.get(
            "index_source", pd.Series("registered_local_pixel_completion", index=additional.index)
        )
        extras.append(additional)
    if extras:
        indexed = pd.concat([indexed, *extras], ignore_index=True, sort=False)
    if not indexed.empty:
        if "qr_exp" not in indexed and "qr" in indexed:
            indexed["qr_exp"] = indexed["qr"]
        if "qz_exp" not in indexed and "qz" in indexed:
            indexed["qz_exp"] = indexed["qz"]
        indexed = indexed.sort_values(
            "assignment_support_score" if "assignment_support_score" in indexed else "feature_id",
            ascending=False if "assignment_support_score" in indexed else True,
        ).drop_duplicates("feature_id", keep="first")
    feature_tables = [search.get("anchors", pd.DataFrame()), search.get("validation", pd.DataFrame()),
                      search.get("ignored", pd.DataFrame())]
    all_consensus = pd.concat([x for x in feature_tables if x is not None and not x.empty], ignore_index=True,
                              sort=False)
    if not all_consensus.empty:
        all_consensus = all_consensus.drop_duplicates("feature_id", keep="first")
    indexed_ids = set(indexed.feature_id.astype(str)) if not indexed.empty else set()
    ignored = all_consensus[
        ~all_consensus.feature_id.astype(str).isin(indexed_ids)].copy() if not all_consensus.empty else pd.DataFrame()
    return indexed.reset_index(drop=True), ignored.reset_index(drop=True)


# ====================== END ROBUST ORIENTATION SOLVER ======================

def accuracy_preset(config: IndexingConfig, mode="balanced"):
    mode = mode.lower()
    if mode == "fast":
        return replace(config, max_normal_candidates_anchor=16, refine_hypotheses=4,
                       max_pair_hypotheses_per_normal=4, max_hypothesis_predictions=40,
                       max_validation_predictions=110, expand_top_normals=6,
                       full_leave_one_angle_out=False,
                       full_bootstrap_iterations=0, test_second_orientation=False)
    if mode == "efficient":
    # Use the exhaustive primary orientation search for the scientific decision
    # while limiting repeated diagnostic resampling to reduce unnecessary runtime.
        exhaustive = accuracy_preset(config, "exhaustive")
        return replace(
            exhaustive,
            max_normal_candidates_anchor=20,
            coarse_tilt_values_deg=(-6.0, 0.0, 6.0),
            full_leave_one_angle_out=False,
            full_bootstrap_iterations=0,
            bootstrap_search_normal_limit=14,
            loo_search_normal_limit=16,
            index_ignored_features=False,
            guided_rescue_unclustered_features=True,
            ignored_index_f2_percentile=0.0,
            ignored_index_max_predictions=1818,
            ignored_evidence_perturbation_trials=0,
            ignored_evidence_decoy_trials=0,
            v7_max_fit_predictions=260,
            v7_refine_hypotheses=4,
            v7_multidomain_candidate_normals=10,
            series_worker_timeout_s=max(config.series_worker_timeout_s, 900),
        )
    if mode == "exhaustive":
        return replace(
            config,
            # Moderately more sensitive feature extraction without relaxing the
            # physical assignment gates.
            feature_threshold_mad=3.0, feature_quantile=0.978,
            ridge_threshold_mad=2.35, max_features_per_image=100,
            max_consensus_features=85,
            anchor_strength_quantile=0.35, max_anchor_features=10,
            max_validation_features=36,
            # Include more calculated reflections because thin-film intensities
            # need not follow kinematic CIF ranking.
            hypothesis_f2_percentile=80.0, validation_f2_percentile=35.0,
            max_hypothesis_predictions=90, max_validation_predictions=220,
            max_pair_hypotheses_per_normal=10,
            max_normal_candidates_anchor=80, refine_hypotheses=20,
            expand_top_normals=16,
            coarse_tilt_values_deg=(-8.0, -4.0, 0.0, 4.0, 8.0),
            full_leave_one_angle_out=True, full_bootstrap_iterations=24,
            bootstrap_search_normal_limit=50, loo_search_normal_limit=55,
            series_worker_timeout_s=max(config.series_worker_timeout_s, 1200),
            test_second_orientation=True,
        )
    if mode != "balanced":
        raise ValueError("mode must be fast, balanced, efficient, or exhaustive")
    return config


def completion_preset(config: IndexingConfig, mode="balanced"):
    """Control only the fixed-orientation full-reflection completion pass.

    This does not change the selected orientation or validation statistics.
    ``aggressive`` is intended for visual-review overlays and accepts somewhat
    more ambiguous/broad residual features; use ``conservative`` when minimizing
    false assignments is more important than coverage.
    """
    mode = str(mode).lower()
    if mode == "balanced":
        return config
    if mode == "conservative":
        return replace(
            config,
            v71_completion_min_member_fraction=0.60,
            v71_completion_max_ambiguity=3,
            v71_completion_min_margin_sigma=0.18,
            v71_completion_streak_min_support=3,
            v71_completion_streak_sigma_limit=2.2,
            v71_completion_streak_max_ambiguity=2,
        )
    if mode == "aggressive":
        return replace(
            config,
            v71_completion_min_member_fraction=0.40,
            v71_completion_max_ambiguity=6,
            v71_completion_min_margin_sigma=0.0,
            v71_completion_streak_min_support=2,
            v71_completion_streak_sigma_limit=3.2,
            v71_completion_streak_max_ambiguity=4,
        )
    raise ValueError("completion mode must be conservative, balanced, or aggressive")


# ======================== USER-FACING RUN INTERFACE ========================

_WORKFLOW_PRESETS = {
    "preview": {"accuracy_mode": "fast", "completion_mode": "conservative", "preview_only": True},
    "recommended": {"accuracy_mode": "efficient", "completion_mode": "balanced", "preview_only": False},
    "maximum_coverage": {"accuracy_mode": "efficient", "completion_mode": "aggressive", "preview_only": False},
}


def _v738_apply_workflow_preset(config: IndexingConfig, name: str):
    """Apply one documented run preset and return (config, run_options)."""
    key = str(name or "recommended").strip().lower().replace("-", "_").replace(" ", "_")
    if key not in _WORKFLOW_PRESETS:
        raise ValueError(
            "workflow preset must be preview, recommended, or maximum_coverage"
        )
    settings = _WORKFLOW_PRESETS[key]
    configured = accuracy_preset(config, settings["accuracy_mode"])
    configured = completion_preset(configured, settings["completion_mode"])
    configured = replace(configured, preview_only=bool(settings["preview_only"]))
    return configured, {
        "workflow_preset": key,
        "preview_only": bool(settings["preview_only"]),
        "run_synthetic_test": False,
    }


def _preview_image_rgb(image, config):
    intensity = np.asarray(image["intensity"], float)
    finite = np.isfinite(intensity)
    normalized = np.where(finite, np.clip(intensity, 0.0, 1.0), 0.0)
    return (plt.get_cmap(config.colormap)(normalized)[..., :3] * 255).astype(np.uint8)


def _write_detected_feature_preview(image, features, output, title, config):
    rgb = _preview_image_rgb(image, config)
    canvas = Image.fromarray(rgb, mode="RGB").convert("RGBA")
    draw = ImageDraw.Draw(canvas, "RGBA")
    font = ImageFont.load_default()
    width, height = canvas.size
    qr_min, qr_max = float(image["qr"][0]), float(image["qr"][-1])
    qz_min, qz_max = float(np.min(image["qz"])), float(np.max(image["qz"]))

    def pixel(qr, qz):
        x = (float(qr) - qr_min) / max(qr_max - qr_min, 1e-12) * (width - 1)
        y = (qz_max - float(qz)) / max(qz_max - qz_min, 1e-12) * (height - 1)
        return int(round(x)), int(round(y))

    for row in features.itertuples(index=False):
        qr = getattr(row, "qr", np.nan);
        qz = getattr(row, "qz", np.nan)
        if not np.isfinite(qr) or not np.isfinite(qz):
            continue
        x, y = pixel(qr, qz)
        draw.ellipse((x - 3, y - 3, x + 3, y + 3), outline=(255, 255, 0, 255), width=1)
    draw.rectangle((0, 0, width, 22), fill=(0, 0, 0, 175))
    draw.text((5, 5), f"{title} | detected features: {len(features)}", fill=(255, 255, 255, 255), font=font)
    output = Path(output);
    output.parent.mkdir(parents=True, exist_ok=True)
    canvas.convert("RGB").save(output, format="PNG", optimize=False)


def _relative_link(path: Path, root: Path) -> str:
    try:
        return path.resolve().relative_to(root.resolve()).as_posix()
    except Exception:
        return path.resolve().as_uri()


def write_input_preview_report(output: Path, manifest: pd.DataFrame, preflight: pd.DataFrame,
                               preview_summary: pd.DataFrame) -> Path:
    report = output / "input_preview_report.html"
    cards = []
    for row in preview_summary.itertuples(index=False):
        image_link = html.escape(str(getattr(row, "preview_file", "")))
        cards.append(
            '<div class="card">'
            f'<img src="{image_link}" alt="input preview">'
            f'<div><b>{html.escape(str(row.series_id))}</b> &nbsp; '
            f'angle={float(row.angle_deg):.3f}﷿﷿, scan={html.escape(str(row.scan))}<br>'
            f'source={html.escape(str(row.source_kind))}, features={int(row.detected_features)}</div>'
            '</div>'
        )
    html_text = f"""<!doctype html>
<html><head><meta charset="utf-8"><title>GIXS input preview</title>
<style>body{{font-family:Arial,sans-serif;margin:24px}}table{{border-collapse:collapse;font-size:12px}}
th,td{{border:1px solid #ccc;padding:4px}}.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(320px,1fr));gap:14px}}
.card{{border:1px solid #bbb;padding:8px}}.card img{{max-width:100%;height:auto}}</style></head>
<body><h1>GIXS input preview</h1>
<p>This report validates grouping, q ranges, crops, and detected image features. It does not run orientation indexing.</p>
<h2>Preflight</h2>{preflight.to_html(index=False, escape=True)}
<h2>Measurement summary</h2>{preview_summary.drop(columns=['preview_file'], errors='ignore').to_html(index=False, escape=True)}
<h2>Preview images</h2><div class="grid">{''.join(cards)}</div>
</body></html>"""
    report.write_text(html_text, encoding="utf-8")
    return report.resolve()


def run_input_preview(config: IndexingConfig):
    """Validate inputs and draw detected-feature previews without indexing."""
    output = Path(config.output_dir);
    output.mkdir(parents=True, exist_ok=True)
    manifest = discover_images(config)
    if manifest.empty:
        raise FileNotFoundError("No q-space inputs were found")
    config = _apply_manifest_coordinate_defaults(config, manifest)
    preflight = validate_measurement_manifest(manifest, config)
    preflight.to_csv(output / "input_preflight_report.csv", index=False)
    errors = preflight[preflight.severity == "error"]
    if not errors.empty and config.manifest_strict:
        raise ValueError("Input preflight failed:\n" + "\n".join(errors.message.astype(str)))

    rows = []
    preview_dir = output / "input_preview_images"
    for record in manifest.itertuples(index=False):
        image = load_qspace_measurement(record, config, config.crop_xyxy)
        detected = detect_features(image, config)
        preview_path = preview_dir / f"{record.series_id.replace(':', '_')}_angle_{record.angle_deg:.3f}_scan_{record.scan}.png"
        if config.write_input_preview_images:
            _write_detected_feature_preview(
                image, detected, preview_path,
                f"{record.series_id} angle {record.angle_deg:.3f} scan {record.scan}", config,
            )
        rows.append({
            "series_id": record.series_id, "sample": record.sample, "series": record.series,
            "angle_deg": record.angle_deg, "scan": record.scan,
            "source_kind": image.get("source_kind", "unknown"),
            "detected_features": int(len(detected)),
            "qr_min": float(image["qr"][0]), "qr_max": float(image["qr"][-1]),
            "qz_min": float(np.min(image["qz"])), "qz_max": float(np.max(image["qz"])),
            "crop": str(image.get("crop")),
            "preview_file": _relative_link(preview_path, output) if config.write_input_preview_images else "",
        })
    summary = pd.DataFrame(rows)
    manifest.to_csv(output / "measurement_manifest.csv", index=False)
    summary.to_csv(output / "input_preview_summary.csv", index=False)
    report = write_input_preview_report(output, manifest, preflight, summary)
    print(f"Input preview saved to: {report}")
    return {
        "status": "preview_only", "config": config, "manifest": manifest,
        "preflight": preflight, "preview_summary": summary,
        "output_dir": output.resolve(), "html_report": report,
    }


def _v739_write_indexing_html_report(results) -> Path | None:
    """Write a single navigable HTML report for a completed indexing run."""
    if not results or results.get("status"):
        return None
    output = Path(results["output_dir"])
    summary = results.get("summary", pd.DataFrame())
    preflight_path = output / "input_preflight_report.csv"
    preflight = pd.read_csv(preflight_path) if preflight_path.is_file() else pd.DataFrame()
    sections = []
    for series_id in summary.get("series_id", pd.Series(dtype=str)).astype(str):
        series_dir = output / series_id.replace(":", "_series_")
        overlay = series_dir / "indexed_or_ignored_overlay.png"
        links = []
        for name in (
                "indexed_reflections.csv", "indexed_overlay_coordinate_key.csv",
                "indexing_summary.json", "registered_local_pixel_completion_assignments.csv",
        ):
            path = series_dir / name
            if path.is_file():
                links.append(f'<a href="{html.escape(_relative_link(path, output))}">{html.escape(name)}</a>')
        image_html = (
            f'<img src="{html.escape(_relative_link(overlay, output))}" alt="{html.escape(series_id)} overlay">'
            if overlay.is_file() else '<p>Overlay not found.</p>'
        )
        sections.append(
            f'<section><h2>{html.escape(series_id)}</h2>{image_html}<p>{" | ".join(links)}</p></section>'
        )
    report_path = output / str(results["config"].html_report_filename)
    report_html = f"""<!doctype html>
<html><head><meta charset="utf-8"><title>GIXS indexing report</title>
<style>body{{font-family:Arial,sans-serif;margin:24px;max-width:1500px}}table{{border-collapse:collapse;font-size:12px;display:block;overflow:auto}}
th,td{{border:1px solid #ccc;padding:4px}}section{{margin:28px 0;border-top:2px solid #ddd;padding-top:18px}}section img{{max-width:100%;height:auto}}</style></head>
<body><h1>GIXS indexing report</h1>
<p><a href="measurement_manifest.csv">measurement manifest</a> | <a href="indexing_summary_all_series.csv">summary CSV</a> | <a href="run_configuration.json">run configuration</a></p>
<h2>Series summary</h2>{summary.to_html(index=False, escape=True)}
<h2>Input preflight</h2>{preflight.to_html(index=False, escape=True) if not preflight.empty else '<p>No preflight table found.</p>'}
{''.join(sections)}</body></html>"""
    report_path.write_text(report_html, encoding="utf-8")
    return report_path.resolve()


def _v739_run_user_friendly(config: IndexingConfig, run_options=None):
    """Single public entry point for preview or full indexing."""
    options = dict(run_options or {})
    preview_only = bool(options.get("preview_only", config.preview_only))
    if preview_only:
        return run_input_preview(replace(config, preview_only=True))
    results = run_gixs_indexing(
        replace(config, preview_only=False),
        run_synthetic_test=bool(options.get("run_synthetic_test", False)),
    )
    if results and not results.get("status") and config.write_html_report:
        report = write_indexing_html_report(results)
        results["html_report"] = report
        print(f"HTML report: {report}")
    return results


# ===================== SCIENTIFIC VALIDATION AND REPORTING =====================
# Add calibration checks, held-out validation, standardized evidence tables, and
# reporting around the orientation solver. These diagnostics quantify how strongly
# the proposed indexing solution is supported by independent or perturbed evidence.

@dataclass
class V8Config(IndexingConfig):
    # Explicit PNG pixel-to-q calibration. When all eight anchor fields are supplied
    # in the manifest, the two horizontal and two vertical anchors override q ranges.
    png_axis_warn_if_uncalibrated: bool = True

    # Optional exact HDF5 dataset paths, avoiding ambiguous key inference.

    # Optional stronger image processing used by the supported indexing workflows.
    enable_multiscale_feature_detection: bool = False
    multiscale_feature_sigmas_px: tuple[float, ...] = (0.8, 1.1, 1.6)
    multiscale_merge_tolerance_q: float = 0.018
    enable_covariance_aware_consensus: bool = False
    consensus_mahalanobis_limit: float = 3.0
    consensus_absolute_floor_q: float = 0.014

    # True image-level holdout: remove one angle before consensus and orientation
    # fitting, then predict the raw detections in the excluded image.
    enable_true_image_holdout: bool = True
    true_holdout_min_training_angles: int = 3
    true_holdout_max_normal_candidates: int = 8
    true_holdout_refine_hypotheses: int = 3

    # Apply weak-series recovery by objective fit quality rather than sample name.
    generalize_weak_series_recovery: bool = True
    weak_series_weighted_fraction_threshold: float = 0.56
    weak_series_min_residual_features: int = 8
    weak_series_max_validation_matches: int = 8

    # Standardized evidence output and interactive reporting.
    write_unified_assignment_table: bool = True

    # Static overlay presentation. ``classic`` produces a single-panel q_r/q_z map.
    # The detailed option includes a coordinate key and is written separately as
    # ``indexed_or_ignored_overlay_detailed.png`` when requested.
    overlay_static_style: str = "detailed"  # two-panel overlay with coordinate key
    overlay_write_detailed_companion: bool = False
    overlay_classic_per_angle: bool = True
    overlay_classic_axis_label_fontsize: float = 18.0
    overlay_classic_tick_fontsize: float = 12.0
    # Use the representative single-angle image for the classic overlay so displayed
    # intensity retains the appearance of the measured beamline-style PNG. A
    # multi-angle composite can still be written as a separate diagnostic.
    overlay_classic_use_representative_image: bool = False
    overlay_classic_show_labels: bool = False
    overlay_write_composite_companion: bool = True

    write_run_provenance: bool = True
    provenance_hash_inputs: bool = True

    # Low-memory execution: cache each registered series to disk and release all
    # other image rasters before the expensive orientation/completion stages.
    stream_series_records_to_disk: bool = True
    retain_images_in_result: bool = False
    series_subprocess_streaming: bool = True


# ------------------------ explicit PNG-axis calibration ------------------------












def load_numerical_qspace(path: str, config: IndexingConfig) -> dict:
    """Load numerical NPZ reciprocal-space data for the indexing workbench."""
    result = _v739_load_numerical_qspace(path, config)
    result.setdefault("axis_calibration_source", "numerical_axes_or_config_range")
    return result


# ------------------------ optional multiscale feature detection ------------------------


def _v88_detect_features(image: dict, config: IndexingConfig) -> pd.DataFrame:
    """Run each requested detector scale once and retain its audit information.

    The prior GUI build recomputed the default detector two extra times per
    image for audit bookkeeping. This preserves the same multiscale merge logic
    while reusing the already-computed default-scale result.
    """
    if not bool(getattr(config, "enable_multiscale_feature_detection", False)):
        result = _v739_detect_features(image, config).copy()
        result["detector_mode"] = "single_scale"
        return result

    requested = tuple(float(value) for value in getattr(
        config, "multiscale_feature_sigmas_px", (config.feature_sigma_px,)
    ))
    unique_sigmas = tuple(dict.fromkeys(requested))
    default_sigma = float(config.feature_sigma_px)
    frames_by_sigma = {}
    frames = []
    for sigma in unique_sigmas:
        local = replace(
            config,
            enable_multiscale_feature_detection=False,
            feature_sigma_px=float(sigma),
        )
        frame = _v739_detect_features(image, local)
        frames_by_sigma[float(sigma)] = frame
        if not frame.empty:
            scaled = frame.copy()
            scaled["detection_scale_px"] = float(sigma)
            frames.append(scaled)

    baseline = next(
        (frame for sigma, frame in frames_by_sigma.items()
         if math.isclose(sigma, default_sigma, rel_tol=0.0, abs_tol=1e-12)),
        None,
    )
    if baseline is None:
        baseline = _v739_detect_features(
            image,
            replace(
                config,
                enable_multiscale_feature_detection=False,
                feature_sigma_px=default_sigma,
            ),
        )

    if not frames:
        result = baseline.copy()
        result["detector_mode"] = "multiscale"
        result["single_scale_feature_count"] = int(len(baseline))
        result["multiscale_feature_count"] = int(len(result))
        result["single_multiscale_overlap_fraction"] = 1.0 if len(result) else 0.0
        return result

    pool = pd.concat(frames, ignore_index=True, sort=False)
    pool["_rank"] = (
        pool.strength.to_numpy(float)
        * np.log1p(np.maximum(pool.snr.to_numpy(float), 0.0))
    )
    kept = []
    tolerance = float(getattr(
        config, "multiscale_merge_tolerance_q", config.feature_merge_tolerance_q
    ))
    for row in pool.sort_values("_rank", ascending=False, kind="mergesort").itertuples(index=False):
        if all(
            math.hypot(float(row.qr) - float(previous.qr),
                       float(row.qz) - float(previous.qz)) > tolerance
            for previous in kept
        ):
            kept.append(row)
        if len(kept) >= int(config.max_features_per_image):
            break

    result = pd.DataFrame([row._asdict() for row in kept]).drop(
        columns=["_rank", "detection_scale_px"], errors="ignore"
    )
    result = result.reindex(columns=list(baseline.columns)).reset_index(drop=True)

    if baseline.empty or result.empty:
        overlap = 0.0
    else:
        tree = cKDTree(baseline[["qr", "qz"]].to_numpy(float))
        distance, _ = tree.query(result[["qr", "qz"]].to_numpy(float))
        overlap = float((distance <= tolerance).mean())

    result["detector_mode"] = "multiscale"
    result["single_scale_feature_count"] = int(len(baseline))
    result["multiscale_feature_count"] = int(len(result))
    result["single_multiscale_overlap_fraction"] = overlap
    return result


# ------------------------ covariance-aware consensus option ------------------------


def build_consensus(features: pd.DataFrame, config: IndexingConfig):
    if not bool(getattr(config, "enable_covariance_aware_consensus", False)) or features.empty:
        return _v739_build_consensus(features, config)
    columns = [
        "feature_id", "qr", "qz", "cov_rr", "cov_rz", "cov_zz", "sigma_qr", "sigma_qz",
        "strength", "support", "support_fraction", "angles", "feature_type",
        "major_width_q", "minor_width_q", "detection_source", "detection_source_mix",
        "experimental_integrated_intensity", "experimental_intensity_sigma",
        "experimental_integrated_snr", "experimental_intensity_quality_score",
    ]
    coordinates = features[["qr", "qz"]].to_numpy(float)
    parent = np.arange(len(coordinates), dtype=int)
    rank = np.zeros(len(coordinates), dtype=np.int8)
    cluster_members = [{i} for i in range(len(coordinates))]

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(a, b):
        a, b = root(int(a)), root(int(b))
        if a == b:
            return
        candidate = cluster_members[a] | cluster_members[b]
        if bool(getattr(config, "consensus_prevent_chaining", True)) and len(candidate) > 2:
            points = coordinates[np.fromiter(candidate, dtype=int)]
            diameter = float(np.max(pdist(points))) if len(points) > 1 else 0.0
            limit = float(config.consensus_tolerance_q) * float(
                getattr(config, "consensus_cluster_diameter_factor", 1.6)
            )
            if diameter > limit:
                return
        if rank[a] < rank[b]:
            a, b = b, a
        parent[b] = a
        cluster_members[a] = candidate
        cluster_members[b] = set()
        if rank[a] == rank[b]:
            rank[a] += 1

    max_radius = max(float(config.consensus_tolerance_q), 2.0 * float(config.uncertainty_ceiling_q))
    floor2 = float(getattr(config, "consensus_absolute_floor_q", 0.014)) ** 2
    for left, right in cKDTree(coordinates).query_pairs(max_radius, output_type="ndarray"):
        delta = coordinates[left] - coordinates[right]
        cov = np.array([
            [features.iloc[left].cov_rr + features.iloc[right].cov_rr + floor2,
             features.iloc[left].cov_rz + features.iloc[right].cov_rz],
            [features.iloc[left].cov_rz + features.iloc[right].cov_rz,
             features.iloc[left].cov_zz + features.iloc[right].cov_zz + floor2],
        ], float)
        try:
            mahal = math.sqrt(max(float(delta @ np.linalg.inv(cov) @ delta), 0.0))
        except np.linalg.LinAlgError:
            mahal = np.inf
        physical = float(np.linalg.norm(delta))
        if physical <= max_radius and mahal <= float(config.consensus_mahalanobis_limit):
            union(left, right)
    roots = np.array([root(i) for i in range(len(coordinates))])
    _, cluster_ids = np.unique(roots, return_inverse=True)
    working = features.copy()
    working["cluster"] = cluster_ids + 1
    rows, members = [], []
    total_angles = max(working.angle_deg.nunique(), 1)
    for _, group in working.groupby("cluster"):
        group = group.loc[group.groupby("angle_deg")["strength"].idxmax()].copy()
        support = int(group.angle_deg.nunique())
        if support < int(config.min_angle_support):
            continue
        weight = np.maximum(group.strength.to_numpy(float), 1e-6)
        weight /= weight.sum()
        qr = float(np.sum(weight * group.qr))
        qz = float(np.sum(weight * group.qz))
        cov = np.zeros((2, 2), float)
        for w, row in zip(weight, group.itertuples(index=False)):
            local = np.array([[row.cov_rr, row.cov_rz], [row.cov_rz, row.cov_zz]], float)
            delta = np.array([row.qr - qr, row.qz - qz], float)
            cov += w * (local + np.outer(delta, delta))
        feature_type = str(group.groupby("feature_type").strength.sum().idxmax())
        if "source_detector" in group.columns:
            source_weight = group.groupby(group["source_detector"].astype(str))["strength"].sum()
            detection_source = str(source_weight.idxmax())
            detection_source_mix = ",".join(sorted(set(group["source_detector"].astype(str))))
        else:
            detection_source = "unreported"
            detection_source_mix = "unreported"
        fid = f"F{len(rows) + 1:03d}"
        rows.append({
            "feature_id": fid, "qr": qr, "qz": qz, "cov_rr": float(cov[0, 0]),
            "cov_rz": float(cov[0, 1]), "cov_zz": float(cov[1, 1]),
            "sigma_qr": math.sqrt(max(cov[0, 0], 1e-12)), "sigma_qz": math.sqrt(max(cov[1, 1], 1e-12)),
            "strength": float(group.strength.max()), "support": support,
            "support_fraction": support / total_angles,
            "angles": ",".join(f"{x:.3f}" for x in sorted(group.angle_deg.unique())),
            "feature_type": feature_type,
            "major_width_q": float(np.average(group.major_width_q, weights=weight)),
            "minor_width_q": float(np.average(group.minor_width_q, weights=weight)),
            "detection_source": detection_source,
            "detection_source_mix": detection_source_mix,
            "experimental_integrated_intensity": float(np.nanmedian(pd.to_numeric(
                group.get("experimental_integrated_intensity", pd.Series(np.nan, index=group.index)), errors="coerce"
            ))) if np.isfinite(pd.to_numeric(group.get("experimental_integrated_intensity", pd.Series(np.nan, index=group.index)), errors="coerce")).any() else np.nan,
            "experimental_intensity_sigma": float(np.nanmedian(pd.to_numeric(
                group.get("experimental_intensity_sigma", pd.Series(np.nan, index=group.index)), errors="coerce"
            ))) if np.isfinite(pd.to_numeric(group.get("experimental_intensity_sigma", pd.Series(np.nan, index=group.index)), errors="coerce")).any() else np.nan,
            "experimental_integrated_snr": float(np.nanmedian(pd.to_numeric(
                group.get("experimental_integrated_snr", pd.Series(np.nan, index=group.index)), errors="coerce"
            ))) if np.isfinite(pd.to_numeric(group.get("experimental_integrated_snr", pd.Series(np.nan, index=group.index)), errors="coerce")).any() else np.nan,
            "experimental_intensity_quality_score": float(np.nanmedian(pd.to_numeric(
                group.get("experimental_intensity_quality_score", pd.Series(np.nan, index=group.index)), errors="coerce"
            ))) if np.isfinite(pd.to_numeric(group.get("experimental_intensity_quality_score", pd.Series(np.nan, index=group.index)), errors="coerce")).any() else np.nan,
        })
        group["feature_id"] = fid
        members.append(group)
    consensus = pd.DataFrame(rows)
    member_table = pd.concat(members, ignore_index=True) if members else working.iloc[0:0].assign(feature_id="")
    if consensus.empty:
        return pd.DataFrame(columns=columns), member_table
    consensus["rank_weight"] = consensus.strength * np.sqrt(consensus.support)
    keep = set(consensus.nlargest(int(config.max_consensus_features), "rank_weight").feature_id)
    consensus = consensus[consensus.feature_id.isin(keep)].drop(columns="rank_weight")
    consensus = consensus.sort_values(["support", "strength"], ascending=False).reset_index(drop=True)
    member_table = member_table[member_table.feature_id.isin(keep)].reset_index(drop=True)
    return consensus[columns], member_table


# ------------------------ true raw-image holdout validation ------------------------


def _v803_full_leave_one_angle_out(crystal, members, config):
    """True angle holdout with a fast training-only orientation re-search.

    The held angle is removed before consensus creation and hypothesis selection.
    The excluded raw detections are then scored only after the training solution is
    fixed. This is independent at the image level while remaining practical enough
    for routine runs.
    """
    if not bool(getattr(config, "enable_true_image_holdout", False)):
        return _v739_full_leave_one_angle_out(crystal, members, config)
    if not config.full_leave_one_angle_out or members is None or members.empty:
        return pd.DataFrame()
    rows = []
    for held_angle in sorted(float(x) for x in members.angle_deg.unique()):
        training = members[members.angle_deg != held_angle].copy()
        if training.angle_deg.nunique() < int(config.true_holdout_min_training_angles):
            continue
        consensus, _ = build_consensus(training, config)
        if len(consensus) < int(config.min_anchor_matches):
            continue
        local = replace(
            config,
            max_normal_candidates_anchor=min(int(config.true_holdout_max_normal_candidates),
                                             int(config.max_normal_candidates_anchor)),
            max_orientation_candidates=min(int(config.true_holdout_max_normal_candidates),
                                           int(config.max_orientation_candidates)),
            max_consensus_features=min(34, int(config.max_consensus_features)),
            max_anchor_features=min(7, int(config.max_anchor_features)),
            max_validation_features=min(16, int(config.max_validation_features)),
            refine_hypotheses=min(int(config.true_holdout_refine_hypotheses), int(config.refine_hypotheses)),
            expand_top_normals=min(int(config.true_holdout_refine_hypotheses), int(config.expand_top_normals)),
            max_pair_hypotheses_per_normal=1,
            coarse_tilt_values_deg=(0.0,),
            max_hypothesis_predictions=min(32, int(config.max_hypothesis_predictions)),
            max_validation_predictions=min(80, int(config.max_validation_predictions)),
            full_leave_one_angle_out=False,
            full_bootstrap_iterations=0,
            test_second_orientation=False,
            v72_enable_s1_mosaic_completion=False,
            v73_enable_s1_sector_domains=False,
            v71_enable_full_reflection_completion=False,
            index_ignored_features=False,
            guided_rescue_unclustered_features=False,
            enable_local_pixel_completion=False,
            write_registered_per_angle_overlays=False,
        )
        try:
            search = anchor_first_search(crystal, consensus, local, fast=True, validation_mode=True)
            held = _v7_prepare_raw_features(members[members.angle_deg == held_angle].copy())
            _, metrics = _v7_assign_arrays(
                held, search["best"]["predictions"], local,
                float(local.v7_all_feature_tolerance_q), materialize=True,
            )
        except Exception as error:
            rows.append({
                "held_angle_deg": held_angle,
                "trained_normal_hkl": "",
                "held_features": int((members.angle_deg == held_angle).sum()),
                "predicted_matches": 0,
                "predicted_weighted_fraction": np.nan,
                "predicted_indexed_fraction": np.nan,
                "predicted_median_delta_q": np.nan,
                "predicted_p90_delta_q": np.nan,
                "training_score_margin": np.nan,
                "validation_method": "remove_angle_before_consensus_fast_training_only_research",
                "held_angle_used_in_training": False,
                "holdout_status": "fit_failed",
                "holdout_error": str(error),
            })
            continue
        rows.append({
            "held_angle_deg": held_angle,
            "trained_normal_hkl": str(search["best"]["hkl"]),
            "held_features": int(len(held)),
            "predicted_matches": int(metrics["matches"]),
            "predicted_weighted_fraction": float(metrics["weighted_fraction"]),
            "predicted_indexed_fraction": float(metrics["indexed_fraction"]),
            "predicted_median_delta_q": metrics["median_delta_q"],
            "predicted_p90_delta_q": metrics["p90_delta_q"],
            "training_score_margin": search.get("score_margin", np.nan),
            "validation_method": "remove_angle_before_consensus_fast_training_only_research",
            "held_angle_used_in_training": False,
            "holdout_status": "ok",
            "holdout_error": "",
        })
    return pd.DataFrame(rows)


# ------------------------ generalized weak-series recovery ------------------------


def _series_is_weak(consensus, search, threshold, minimum_residual, maximum_validation_matches):
    weighted = float(search["best"].get("all_feature_metrics", {}).get("weighted_fraction", 0.0))
    validation_matches = int(search["best"].get("validation_metrics", {}).get("matches", 0))
    core = combined_matches(search)
    assigned = set(core.feature_id.astype(str)) if not core.empty else set()
    residual = int((~consensus.feature_id.astype(str).isin(
        assigned)).sum()) if consensus is not None and not consensus.empty else 0
    return (
            weighted < float(threshold)
            and residual >= int(minimum_residual)
            and validation_matches <= int(maximum_validation_matches)
    )


def v72_s1_mosaic_completion(series_id, crystal, consensus, members, search, secondary, completion_assignments, config):
    if bool(getattr(config, "generalize_weak_series_recovery", False)):
        targeted = _series_is_weak(
            consensus, search, config.weak_series_weighted_fraction_threshold,
            config.weak_series_min_residual_features, config.weak_series_max_validation_matches
        )
        local = replace(config, v72_s1_series_prefixes=((str(series_id),) if targeted else ()))
        return _v739_v72_s1_mosaic_completion(
            series_id, crystal, consensus, members, search, secondary, completion_assignments, local
        )
    return _v739_v72_s1_mosaic_completion(
        series_id, crystal, consensus, members, search, secondary, completion_assignments, config
    )


def v73_s1_sector_domain_search(series_id, crystal, consensus, search, base_secondary, config):
    if bool(getattr(config, "generalize_weak_series_recovery", False)):
        targeted = _series_is_weak(
            consensus, search, config.weak_series_weighted_fraction_threshold,
            config.weak_series_min_residual_features, config.weak_series_max_validation_matches
        )
        local = replace(config, v73_s1_sector_series_prefixes=((str(series_id),) if targeted else ()))
        return _v739_v73_s1_sector_domain_search(series_id, crystal, consensus, search, base_secondary, local)
    return _v739_v73_s1_sector_domain_search(series_id, crystal, consensus, search, base_secondary, config)


# ------------------------ unified assignment evidence table ------------------------
def _unified_evidence_table(result):
    sources = []

    def add(frame, stage, tier, primary=False, validation=False):
        if frame is None or frame.empty:
            return
        table = frame.copy()
        table["evidence_stage"] = stage
        table["evidence_tier"] = tier
        table["is_primary_orientation_evidence"] = bool(primary)
        table["is_independent_validation_evidence"] = bool(validation)
        sources.append(table)

    search = result.get("search", {})
    add(search.get("best", {}).get("all_feature_matches", pd.DataFrame()), "core_robust_all_feature_fit", "core", True,
        False)
    add(search.get("best", {}).get("validation_matches", pd.DataFrame()), "model_selection_subset", "selection", False, False)
    secondary = result.get("secondary") or {}
    if secondary.get("accepted") and secondary.get("joint") is not None:
        joint = secondary["joint"]
        add(joint.get("all_feature_matches", joint.get("anchor_matches", pd.DataFrame())), "joint_multidomain", "core",
            True, False)
        add(joint.get("validation_matches", pd.DataFrame()), "joint_multidomain_selection_subset", "selection", False, False)
    add(result.get("completion_assignments", pd.DataFrame()), "fixed_orientation_completion", "supported")
    add(result.get("mosaic_assignments", pd.DataFrame()), "mosaic_completion", "supported")
    add(result.get("sector_assignments", pd.DataFrame()), "sector_domain", "core", True, False)
    ignored = result.get("ignored_salvage", pd.DataFrame())
    if ignored is not None and not ignored.empty:
        tiers = ignored.get("salvage_evidence_tier", pd.Series("provisional", index=ignored.index))
        for tier, group in ignored.groupby(tiers):
            add(group, "ignored_consensus_rescue", str(tier))
    guided = result.get("guided_rescue", pd.DataFrame())
    if guided is not None and not guided.empty:
        tiers = guided.get("salvage_evidence_tier", pd.Series("supported", index=guided.index))
        for tier, group in guided.groupby(tiers):
            add(group, "raw_detection_rescue", str(tier))
    add(result.get("local_pixel_completion", pd.DataFrame()), "registered_local_pixel", "image_supported")
    if not sources:
        return pd.DataFrame()
    merged = pd.concat(sources, ignore_index=True, sort=False)
    priority = {
        "selection": 0,
        "validation": 0,
        "core": 1,
        "robust": 2,
        "supported": 3,
        "image_supported": 4,
        "provisional": 8,
    }
    merged["_priority"] = merged.evidence_tier.map(priority).fillna(6)
    support_values = merged["assignment_support_score"] if "assignment_support_score" in merged else pd.Series(0.0,
                                                                                                               index=merged.index)
    merged["_support"] = pd.to_numeric(support_values, errors="coerce").fillna(0.0)
    if "feature_id" in merged:
        merged = merged.sort_values(["_priority", "_support"], ascending=[True, False], kind="mergesort")
        merged = merged.drop_duplicates("feature_id", keep="first")
    merged["display_as_indexed"] = merged.evidence_tier.astype(str).isin(
        {"core", "selection", "validation", "supported", "robust", "image_supported"}
    )
    return merged.drop(columns=["_priority", "_support"], errors="ignore").reset_index(drop=True)


# ------------------------ interactive SVG overlay ------------------------


def _write_interactive_overlay(series_id, image, indexed, ignored, series_dir, config):
    """Return the optional output hook used by reporting code; supported GUI workflows leave this analysis disabled."""
    return None


def postprocess_registered_series(series_id, result, records, crystal, config, output):
    _v8_post_start = time.perf_counter()
    result = _v739_postprocess_registered_series(series_id, result, records, crystal, config, output)
    series_dir = Path(output) / series_id.replace(":", "_series_")
    unified = _unified_evidence_table(result)
    result["unified_assignments"] = unified
    if bool(getattr(config, "write_unified_assignment_table", False)):
        unified.to_csv(series_dir / "unified_assignment_evidence.csv", index=False)
        approved = unified[unified.get("display_as_indexed", pd.Series(False, index=unified.index)).astype(bool)].copy()
        approved.to_csv(series_dir / "unified_indexed_assignments.csv", index=False)
        if not approved.empty:
            completion_mask = (
                    ~approved.is_primary_orientation_evidence.astype(bool)
                    & ~approved.is_independent_validation_evidence.astype(bool)
            )
            approved[completion_mask].to_csv(series_dir / "unified_completion_assignments.csv", index=False)
    indexed_path = series_dir / "overlay_indexed_features.csv"
    ignored_path = series_dir / "overlay_ignored_features.csv"
    indexed = pd.read_csv(indexed_path) if indexed_path.is_file() else pd.DataFrame()
    ignored = pd.read_csv(ignored_path) if ignored_path.is_file() else pd.DataFrame()
    composite = registered_composite_image(records, config)
    interactive = _write_interactive_overlay(series_id, composite, indexed, ignored, series_dir, config)
    result["interactive_overlay"] = interactive
    result["summary"]["unified_assignment_evidence"] = int(len(unified))
    result["summary"]["unified_indexed_assignments"] = int(
        unified.get("display_as_indexed", pd.Series(False, index=unified.index)).astype(bool).sum()
    )
    result["summary"]["interactive_overlay"] = str(interactive) if interactive else ""
    print(
        f"[{series_id}] unified evidence and interactive overlay finished in {time.perf_counter() - _v8_post_start:.1f}s",
        flush=True)
    return result


# ------------------------ provenance and enhanced report ------------------------
def _sha256_file(path, block=1024 * 1024):
    import hashlib
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            chunk = handle.read(block)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def write_run_provenance(results):
    if not results or results.get("status") or not bool(getattr(results.get("config"), "write_run_provenance", False)):
        return None
    import platform
    import datetime
    import importlib.metadata

    output = Path(results["output_dir"])
    config = results["config"]
    inputs = []
    manifest = results.get("manifest", pd.DataFrame())
    paths = [config.cif_path]
    if not manifest.empty:
        for column in ("path", "numerical_path"):
            if column in manifest:
                paths.extend(x for x in manifest[column].dropna().astype(str) if x)
    for raw in dict.fromkeys(paths):
        path = Path(str(raw)).expanduser()
        record = {"path": str(path.resolve()) if path.exists() else str(path), "exists": path.is_file()}
        if path.is_file() and bool(getattr(config, "provenance_hash_inputs", True)):
            record["sha256"] = _sha256_file(path)
            record["size_bytes"] = path.stat().st_size
        inputs.append(record)
    packages = {}
    for name in ("numpy", "pandas", "scipy", "matplotlib", "Pillow", "gemmi", "cloudpickle"):
        try:
            packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            packages[name] = None
    module_file = globals().get("__file__")
    code = {"path": str(module_file) if module_file else "notebook_cell", "sha256": None}
    if module_file and Path(module_file).is_file():
        code["sha256"] = _sha256_file(module_file)
    payload = {
        "workflow_version": "8.0",
        "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "packages": packages,
        "random_seed": int(config.random_seed),
        "code": code,
        "inputs": inputs,
        "config": asdict(config),
    }
    target = output / "run_provenance.json"
    target.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
    return target.resolve()


def write_indexing_html_report(results) -> Path | None:
    report = _v739_write_indexing_html_report(results)
    if report is None:
        return None
    output = Path(results["output_dir"])
    text = Path(report).read_text(encoding="utf-8")
    cards = []
    summary = results.get("summary", pd.DataFrame())
    for series_id in summary.get("series_id", pd.Series(dtype=str)).astype(str):
        series_dir = output / series_id.replace(":", "_series_")
        interactive = series_dir / "interactive_indexed_or_ignored_overlay.html"
        unified = series_dir / "unified_indexed_assignments.csv"
        holdout = series_dir / "leave_one_angle_out_full_search.csv"
        links = []
        if interactive.is_file():
            links.append(f'<a href="{html.escape(_relative_link(interactive, output))}">interactive overlay</a>')
        if unified.is_file():
            links.append(f'<a href="{html.escape(_relative_link(unified, output))}">unified assignments</a>')
        if holdout.is_file():
            links.append(f'<a href="{html.escape(_relative_link(holdout, output))}">true image holdout</a>')
        cards.append(f'<li><b>{html.escape(series_id)}</b>: {" | ".join(links)}</li>')
    addition = (
            '<h2>V8 defensibility outputs</h2><ul>' + ''.join(cards)
            + '</ul><p><a href="run_provenance.json">run provenance</a></p>'
    )
    text = text.replace('</body></html>', addition + '</body></html>')
    Path(report).write_text(text, encoding="utf-8")
    return report


def _v803_run_user_friendly(config: IndexingConfig, run_options=None):
    results = _v739_run_user_friendly(config, run_options)
    if results and not results.get("status"):
    # Write provenance after the indexing report. The report references the expected
    # provenance filename, so the link becomes valid as soon as the file is created
    # and no additional HTML rewrite is required.
        results["provenance"] = write_run_provenance(results)
    return results


    # Build the compact configuration used by the supported GUI workflows. Additional
    # scientific controls remain available through the configuration dataclass.


def _ensure_assignment_compatible_column(frame: pd.DataFrame, column: str, values) -> None:
    """Prepare a DataFrame column for dtype-safe assignment.

    Pandas 2.2+ and 3.x reject assigning strings into a column that was first
    created from ``np.nan`` and therefore inferred as ``float64``.  This occurs
    when a fresh series worker returns metadata such as ``input_manifest_path``
    to a parent manifest that did not previously contain that column.

    The helper creates new columns with a nullable dtype compatible with the
    incoming values and promotes an existing numeric column to ``object`` only
    when non-numeric values must be stored.
    """
    incoming = pd.Series(values)
    non_missing = incoming.dropna()

    if column not in frame.columns:
        if pd.api.types.is_bool_dtype(incoming.dtype):
            frame[column] = pd.Series(pd.NA, index=frame.index, dtype="boolean")
        elif pd.api.types.is_integer_dtype(incoming.dtype):
            frame[column] = pd.Series(pd.NA, index=frame.index, dtype="Int64")
        elif pd.api.types.is_float_dtype(incoming.dtype):
            frame[column] = pd.Series(np.nan, index=frame.index, dtype="float64")
        elif pd.api.types.is_datetime64_any_dtype(incoming.dtype):
            frame[column] = pd.Series(pd.NaT, index=frame.index, dtype=incoming.dtype)
        else:
            # Object is intentionally used instead of pandas' StringDtype here:
            # worker metadata can contain strings, Paths, numbers, or mixed values.
            frame[column] = pd.Series([None] * len(frame), index=frame.index, dtype="object")
        return

    target = frame[column]
    incoming_has_text = bool(
        not non_missing.empty
        and non_missing.map(lambda value: isinstance(value, (str, Path))).any()
    )
    if incoming_has_text and pd.api.types.is_numeric_dtype(target.dtype):
        frame[column] = target.astype("object")


def _run_v8_series_manifest_subprocess(module_file, series_id, series_manifest, config, output):
    # Run one manifest-defined series in a fresh interpreter with a small payload.
    import cloudpickle
    safe = str(series_id).replace(":", "_series_")
    work = Path(output) / "_v8_series_subprocess"
    work.mkdir(parents=True, exist_ok=True)
    manifest_path = work / f"{safe}_manifest.csv"
    payload_path = work / f"{safe}_payload.pkl"
    result_path = work / f"{safe}_result.pkl"
    error_path = work / f"{safe}_error.txt"
    exported_manifest = series_manifest.copy().rename(columns={
        "path": "png_file", "numerical_path": "numerical_file",
    })
    exported_manifest.to_csv(manifest_path, index=False)
    child_config = asdict(config)
    child_config.update({
        "manifest_path": str(manifest_path),
        "output_dir": str(output),
        "series_subprocess_streaming": False,
        "isolate_series_processes": False,
        "stream_series_records_to_disk": False,
        "retain_images_in_result": False,
        "write_html_report": False,
        "full_bootstrap_iterations": int(config.full_bootstrap_iterations),
    })
    with open(payload_path, "wb") as handle:
        cloudpickle.dump((str(series_id), child_config), handle)
    worker_code = r'''import importlib.util, os, sys, traceback
from pathlib import Path
import cloudpickle
module_path, payload_path, result_path, error_path = sys.argv[1:5]
os.environ["GIXS_WORKER_IMPORT"] = "1"
try:
    spec = importlib.util.spec_from_file_location("gixs_v8_series_module", module_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    with open(payload_path, "rb") as handle:
        series_id, config_dict = cloudpickle.load(handle)
    config_class = getattr(module, "V8Config", module.IndexingConfig)
    # Notebook kernels may retain/materialize an older source revision. Construct
    # the worker config using only fields accepted by that dataclass, then attach
    # newer wrapper settings as ordinary attributes. The active code consistently
    # reads optional wrapper settings through getattr(), so this is backward-safe.
    import dataclasses
    accepted_fields = (
        {field.name for field in dataclasses.fields(config_class)}
        if dataclasses.is_dataclass(config_class) else set(config_dict)
    )
    constructor_values = {
        key: value for key, value in config_dict.items() if key in accepted_fields
    }
    config = config_class(**constructor_values)
    for key, value in config_dict.items():
        if key not in accepted_fields:
            try:
                setattr(config, key, value)
            except Exception:
                pass
    run = module.run_gixs_indexing(config, run_synthetic_test=False)
    import pandas as pd
    registration_path = Path(config.output_dir) / "per_image_registration_diagnostics.csv"
    bundle = {
        "series_result": run["series_results"][series_id],
        "features": run.get("features"),
        "manifest": run.get("manifest"),
        "registration": pd.read_csv(registration_path) if registration_path.is_file() else pd.DataFrame(),
    }
    with open(result_path, "wb") as handle:
        cloudpickle.dump(bundle, handle)
except Exception:
    Path(error_path).write_text(traceback.format_exc(), encoding="utf-8")
    raise
'''
    env = os.environ.copy()
    env["GIXS_WORKER_IMPORT"] = "1"
    env["PYTHONUNBUFFERED"] = "1"
    env.setdefault("MPLBACKEND", "Agg")
    completed = subprocess.run(
        [sys.executable, "-c", worker_code, str(module_file), str(payload_path), str(result_path), str(error_path)],
        capture_output=True, text=True, env=env, timeout=int(config.series_worker_timeout_s),
    )
    if completed.stdout:
        print(completed.stdout, end="")
    if completed.returncode != 0 or not result_path.is_file():
        detail = error_path.read_text(encoding="utf-8") if error_path.is_file() else completed.stderr
        raise RuntimeError(f"V8 series subprocess {series_id} failed:\n{detail}")
    with open(result_path, "rb") as handle:
        bundle = cloudpickle.load(handle)
    for item in (payload_path, result_path, error_path, manifest_path):
        item.unlink(missing_ok=True)
    return bundle


# ------------------------ true low-memory series-streaming runner ------------------------


def _v88_run_gixs_indexing(config: IndexingConfig, run_synthetic_test=True):
    """Process one measurement series at a time to limit full-resolution raster memory.

    Each series is loaded, analyzed, rendered, and released independently. This
    keeps feature detection and registration deterministic while avoiding the
    memory cost of retaining every series raster simultaneously.
    """
    output = Path(config.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    manifest = discover_images(config)
    if manifest.empty:
        print("No q-space inputs were found.")
        return {"status": "missing_images", "config": config}
    config = _apply_manifest_coordinate_defaults(config, manifest)
    preflight = validate_measurement_manifest(manifest, config)
    preflight.to_csv(output / "input_preflight_report.csv", index=False)
    errors = preflight[preflight.severity == "error"]
    if not errors.empty and config.manifest_strict:
        raise ValueError("Input preflight failed:\n" + "\n".join(errors.message.astype(str)))

    crystal = load_reflections(config)
    print(f"CIF: {crystal['path'].name} | space group {crystal['spacegroup'].number} {crystal['spacegroup'].hm}")
    print(
        f"Allowed reflections: {len(crystal['reflections'])} | Gemmi/all-230 disagreements: {len(crystal['validation_disagreements'])}")
    (output / "dwba_model_metadata.json").write_text(
        json.dumps(dwba_model_metadata(config), indent=2, default=str)
    )
    if config.enable_dwba:
        dwba_intensity_envelope(np.array([max(config.analysis_qz_min, 0.2)]), config)
    if config.simulate_powder:
        save_powder_outputs(crystal, config, output)

    alternative_crystals = []
    for candidate_path in config.alternative_cif_paths:
        alternative_crystals.append(load_reflections(replace(
            config, cif_path=candidate_path, all230_compare=False, all230_policy="gemmi"
        )))
    synthetic = synthetic_recovery_test(crystal, config) if run_synthetic_test else {"passed": np.nan}
    synthetic["validation_type"] = "internal_same_model_regression_only"
    print("Internal synthetic regression test:", synthetic)
    (output / "internal_synthetic_regression_test.json").write_text(
        json.dumps(synthetic, indent=2, default=str)
    )

    series_results = {}
    all_feature_frames = []
    registration_frames = []
    retained_images = []
    use_series_subprocesses = (
            bool(getattr(config, "series_subprocess_streaming", False))
            and manifest["series_id"].astype(str).nunique() > 1
    )
    module_file = None
    if use_series_subprocesses:
        module_file = _resolve_series_worker_module_file(output)
        if not module_file:
            print("Fresh series module unavailable; using in-process series streaming.")
            use_series_subprocesses = False
        else:
            print("V8 multi-series run: sequential fresh-process streaming enabled.")

    for series_id, series_manifest in manifest.groupby("series_id", sort=True):
        print(f"[{series_id}] loading and registering series...", flush=True)
        series_manifest = series_manifest.sort_values(["angle_deg", "scan"]).copy()
        if use_series_subprocesses:
            bundle = _run_v8_series_manifest_subprocess(
                module_file, str(series_id), series_manifest, config, output
            )
            series_results[str(series_id)] = bundle["series_result"]
            if bundle.get("features") is not None and not bundle["features"].empty:
                all_feature_frames.append(bundle["features"])
            if bundle.get("registration") is not None and not bundle["registration"].empty:
                registration_frames.append(bundle["registration"])
            child_manifest = bundle.get("manifest")
            if child_manifest is not None and not child_manifest.empty:
                child_lookup = child_manifest.set_index("scan")
                mask = manifest["scan"].isin(child_lookup.index)
                for column in child_manifest.columns:
                    if column == "scan":
                        continue
                    mapping = child_lookup[column].to_dict()
                    mapped_values = manifest.loc[mask, "scan"].map(mapping)
                    _ensure_assignment_compatible_column(manifest, column, mapped_values)
                    manifest.loc[mask, column] = mapped_values.to_numpy()
            print(f"[{series_id}] series subprocess complete.", flush=True)
            _release_working_memory()
            continue
        records, feature_tables = [], []
        shared_crop = config.crop_xyxy
        series_rows = list(series_manifest.itertuples(index=False))
        for image_number, row in enumerate(series_rows, start=1):
            image_name = Path(str(row.path)).name
            print(
                f"[{series_id}] image {image_number}/{len(series_rows)}: loading {image_name}...",
                flush=True,
            )
            crop_for_image = shared_crop if config.reuse_first_png_crop else config.crop_xyxy
            image = load_qspace_measurement(row, config, crop_for_image)
            if config.reuse_first_png_crop and image["source_kind"] == "png" and shared_crop is None:
                shared_crop = image["crop"]
            print(
                f"[{series_id}] image {image_number}/{len(series_rows)}: detecting features...",
                flush=True,
            )
            features = detect_features(image, config)
            print(
                f"[{series_id}] image {image_number}/{len(series_rows)}: "
                f"detected {len(features)} features.",
                flush=True,
            )
            features["sample"], features["series"], features["series_id"] = row.sample, row.series, row.series_id
            features["angle_deg"], features["scan"] = row.angle_deg, row.scan
            for key in ("qr_grid", "qz_grid", "rgb", "raw_intensity"):
                image.pop(key, None)
            image["intensity"] = np.asarray(image["intensity"], dtype=np.float32)
            image["display_intensity"] = image["intensity"]
            image["valid"] = np.asarray(image["valid"], dtype=bool)
            records.append({**row._asdict(), "image": image})
            feature_tables.append(features)

        group = pd.concat(feature_tables, ignore_index=True) if feature_tables else pd.DataFrame()
        print(f"[{series_id}] registering {len(records)} measurement images...", flush=True)
        group, records, registration = register_measurement_series(group, records, config)
        print(f"[{series_id}] registration complete; building consensus...", flush=True)
        if not registration.empty:
            registration_frames.append(registration)
            lookup = registration.set_index("scan")
            for column in (
                    "reference_scan", "registration_accepted", "registration_reason", "matched_features",
                    "median_residual_before", "median_residual_after", "scale", "rotation_deg",
                    "qr_offset", "qz_offset", "registration_used_for_indexing",
                    "series_scatter_improvement_fraction",
            ):
                if column in lookup:
                    mapping = lookup[column].to_dict()
                    mask = manifest["scan"].isin(mapping)
                    mapped_values = manifest.loc[mask, "scan"].map(mapping)
                    _ensure_assignment_compatible_column(manifest, column, mapped_values)
                    manifest.loc[mask, column] = mapped_values.to_numpy()
        all_feature_frames.append(group.copy())
        representative = registered_composite_image(records, config)

        previous_worker_flag = os.environ.get("GIXS_SERIES_WORKER")
        os.environ["GIXS_SERIES_WORKER"] = "1"
        try:
            result = _analyze_single_series(
                str(series_id), group.reset_index(drop=True), representative,
                len(records), crystal, config, output,
            )
        finally:
            if previous_worker_flag is None:
                os.environ.pop("GIXS_SERIES_WORKER", None)
            else:
                os.environ["GIXS_SERIES_WORKER"] = previous_worker_flag
        result = postprocess_registered_series(
            str(series_id), result, records, crystal, config, output
        )
        result["_registered_postprocessed"] = True
        series_results[str(series_id)] = result
        if bool(getattr(config, "retain_images_in_result", False)):
            retained_images.extend(records)
        del records, feature_tables, group, representative, result
        _release_working_memory()
        print(f"[{series_id}] series complete; memory released.", flush=True)

    all_features = pd.concat(all_feature_frames, ignore_index=True) if all_feature_frames else pd.DataFrame()
    registration_diagnostics = (
        pd.concat(registration_frames, ignore_index=True)
        if registration_frames else pd.DataFrame()
    )
    registration_diagnostics.to_csv(output / "per_image_registration_diagnostics.csv", index=False)
    if bool(getattr(config, "series_subprocess_streaming", False)):
        shutil.rmtree(output / "_v8_series_subprocess", ignore_errors=True)

    repeats = repeat_scan_validation(series_results, crystal, config)
    repeat_scores = {}
    for row in repeats.itertuples(index=False):
        repeat_scores[row.series_a] = row.repeat_agreement_score
        repeat_scores[row.series_b] = row.repeat_agreement_score
    cif_report, cif_conclusions = compare_cif_candidates(
        series_results, [crystal, *alternative_crystals], config
    )
    cif_report.to_csv(output / "cif_candidate_compatibility.csv", index=False)

    summaries = []
    for series_id, result in series_results.items():
        repeat_score = repeat_scores.get(series_id, np.nan)
        status, flags = solution_reliability(
            result["search"], result["leave_one_out"], result["bootstrap_summary"], config,
            repeat_score,
        )
        conclusion, criteria = orientation_conclusion(
            result["search"], result["leave_one_out"], result["bootstrap_summary"], config,
            repeat_score,
        )
        external = external_truth_check(
            series_id, result["search"]["best"]["hkl"], crystal, config
        )
        result["summary"].update({
            "status": status,
            "reliability_flags": flags,
            "orientation_conclusion": conclusion,
            "orientation_uniqueness_criteria": criteria,
            "repeat_agreement_score": repeat_score,
            "validation_scope": (
                "external_orientation_test_included" if external["external_truth_status"] == "tested"
                else "internal_consistency_only"
            ),
            "validated_end_to_end": False,
            **external,
            **cif_conclusions.get(series_id, {}),
        })
        series_dir = output / series_id.replace(":", "_series_")
        cif_report[cif_report.series_id == series_id].to_csv(
            series_dir / "cif_candidate_compatibility.csv", index=False
        )
        (series_dir / "indexing_summary.json").write_text(
            json.dumps(result["summary"], indent=2, default=str)
        )
        summaries.append(result["summary"])

    manifest.to_csv(output / "measurement_manifest.csv", index=False)
    all_features.to_csv(output / "detected_features_all_measurements.csv", index=False)
    crystal["validation_disagreements"].to_csv(output / "gemmi_all230_disagreements.csv", index=False)
    repeats.to_csv(output / "repeat_scan_validation.csv", index=False)
    summary_table = pd.DataFrame(summaries)
    summary_table.to_csv(output / "indexing_summary_all_series.csv", index=False)
    (output / "run_configuration.json").write_text(json.dumps(asdict(config), indent=2, default=str))
    claims = {
        "unique_orientation_claim": "never based solely on rank 1; see orientation_conclusion and ambiguity set",
        "heuristic_scores": "ranking/support quantities are not calibrated probabilities or confidence intervals",
        "cif_claim": "compatibility or preference among supplied candidates only; never proof of correctness",
        "external_validation": "absent unless external_truth_normals are supplied",
        "end_to_end_validation": False,
        "dwba_scope": dwba_model_metadata(config),
        "powder_scope": "kinematic diagnostic for radial/CIF compatibility",
        "true_image_holdout": bool(getattr(config, "enable_true_image_holdout", False)),
        "completion_evidence": "post-fit completion is never counted as independent orientation validation",
    }
    (output / "scientific_claim_boundaries.json").write_text(
        json.dumps(claims, indent=2, default=str)
    )
    print("\nAMBIGUITY-AWARE INDEXING SUMMARY")
    print(summary_table[[
        "series_id", "status", "orientation_conclusion", "top_ranked_normal_hkl",
        "anchor_matches", "validation_matches", "heuristic_score_gap",
        "leave_one_angle_out_weighted_fraction", "bootstrap_orientation_stability",
        "repeat_agreement_score", "validation_scope", "cif_conclusion",
    ]].to_string(index=False))
    if not repeats.empty:
        print("\nREPEAT-SCAN PREDICTION")
        print(repeats.to_string(index=False))
    print(f"\nSaved to: {output.resolve()}")
    return {
        "config": config, "manifest": manifest, "crystal": crystal,
        "images": retained_images, "features": all_features,
        "series_results": series_results, "repeat_scan_validation": repeats,
        "internal_synthetic_regression_test": synthetic,
        "cif_candidate_compatibility": cif_report,
        "summary": summary_table, "output_dir": output.resolve(),
    }


# --------------------- workflow presets, preflight, and templates ---------------------


def _v803_apply_workflow_preset(config: IndexingConfig, name: str):
    configured, options = _v738_apply_workflow_preset(config, name)
    key = options["workflow_preset"]
    if key == "preview":
        configured = replace(
            configured,
            enable_true_image_holdout=False,
            write_run_provenance=True,
        )
    elif key == "recommended":
        configured = replace(
            configured,
            enable_true_image_holdout=True,
            full_leave_one_angle_out=True,
            max_parallel_series_workers=1,
            isolate_series_processes=False,
            stream_series_records_to_disk=True,
            true_holdout_max_normal_candidates=min(8, configured.max_normal_candidates_anchor),
            true_holdout_refine_hypotheses=3,
            enable_multiscale_feature_detection=False,
            enable_covariance_aware_consensus=False,
        )
    elif key == "maximum_coverage":
        configured = replace(
            configured,
            enable_true_image_holdout=True,
            full_leave_one_angle_out=True,
            max_parallel_series_workers=1,
            isolate_series_processes=False,
            stream_series_records_to_disk=True,
            true_holdout_max_normal_candidates=min(8, configured.max_normal_candidates_anchor),
            true_holdout_refine_hypotheses=2,
            enable_multiscale_feature_detection=True,
            enable_covariance_aware_consensus=False,
        )
    return configured, options


def validate_measurement_manifest(frame: pd.DataFrame, config: IndexingConfig) -> pd.DataFrame:
    report = _v738_validate_measurement_manifest(frame, config)
    rows = []
    anchor_columns = (
        "qr_pixel_0", "qr_value_0", "qr_pixel_1", "qr_value_1",
        "qz_pixel_0", "qz_value_0", "qz_pixel_1", "qz_value_1",
    )
    if not frame.empty and bool(getattr(config, "png_axis_warn_if_uncalibrated", True)):
        for series_id, group in frame.groupby("series_id", sort=False):
            png = group[group.source_kind.eq("png")]
            if png.empty:
                continue
            complete = pd.Series(False, index=png.index)
            if set(anchor_columns).issubset(png.columns):
                complete = png[list(anchor_columns)].apply(pd.to_numeric, errors="coerce").notna().all(axis=1)
            count = int(complete.sum())
            if count == len(png):
                rows.append({"severity": "ok", "check": "png_axis_anchors", "series_id": series_id,
                             "message": f"All {count} PNG measurements use explicit pixel-to-q anchors."})
            else:
                rows.append({"severity": "warning", "check": "png_axis_anchors", "series_id": series_id,
                             "message": f"{len(png) - count}/{len(png)} PNG measurements rely on crop limits rather than explicit tick anchors."})
    if rows:
        report = pd.concat([report, pd.DataFrame(rows)], ignore_index=True)
    return report


# ===================== VALIDATION AND PRESENTATION CONTROLS =====================
@dataclass
class V81Config(V8Config):
    """Configuration controls for stronger validation and consistent scientific reporting."""
    strict_numerical_input: bool = False
    true_holdout_exact_primary_search: bool = True
    consensus_prevent_chaining: bool = True
    consensus_cluster_diameter_factor: float = 1.6
    overlay_classic_axis_label_fontsize: float = 28.0
    overlay_classic_tick_fontsize: float = 22.0
    write_output_inventory: bool = True
    write_scientific_parameter_audit: bool = True


def load_qspace_measurement(row, config: IndexingConfig, crop=None) -> dict:
    """Load the NPZ numerical reciprocal-space input required by the supported GUI workflows."""
    measurement_config = _measurement_specific_config(row, config)
    numerical_path = _row_field(row, "numerical_path")
    if not numerical_path:
        raise FileNotFoundError(
            f"No NPZ numerical q-space data for scan {_row_field(row, 'scan', '?')}"
        )
    if Path(str(numerical_path)).suffix.lower() != ".npz":
        raise ValueError(
            f"Automatic calculator input must be NPZ, not {Path(str(numerical_path)).name}."
        )
    return load_numerical_qspace(numerical_path, measurement_config)


def full_leave_one_angle_out(crystal, members, config):
    """Perform image-level holdout validation with the primary orientation solver when exact holdout mode is enabled."""
    if not bool(getattr(config, "true_holdout_exact_primary_search", False)):
        return _v803_full_leave_one_angle_out(crystal, members, config)
    if not bool(getattr(config, "enable_true_image_holdout", False)):
        return _v739_full_leave_one_angle_out(crystal, members, config)
    if not config.full_leave_one_angle_out or members is None or members.empty:
        return pd.DataFrame()
    rows = []
    method = "remove_angle_before_consensus_exact_primary_v7_solver"
    for held_angle in sorted(float(x) for x in members.angle_deg.unique()):
        training = members[members.angle_deg != held_angle].copy()
        if training.angle_deg.nunique() < int(config.true_holdout_min_training_angles):
            continue
        consensus, _ = build_consensus(training, config)
        if len(consensus) < int(config.min_anchor_matches):
            continue
        local = replace(
            config,
            full_leave_one_angle_out=False,
            full_bootstrap_iterations=0,
            test_second_orientation=False,
            v72_enable_s1_mosaic_completion=False,
            v73_enable_s1_sector_domains=False,
            v71_enable_full_reflection_completion=False,
            index_ignored_features=False,
            guided_rescue_unclustered_features=False,
            enable_local_pixel_completion=False,
            write_registered_per_angle_overlays=False,
        )
        try:
            search = v7_indexing_search(
                crystal, consensus, local, fast=False, validation_mode=False, raw_features=training
            )
            held = _v7_prepare_raw_features(members[members.angle_deg == held_angle].copy())
            _, metrics = _v7_assign_arrays(
                held, search["best"]["predictions"], local,
                float(local.v7_all_feature_tolerance_q), materialize=True,
            )
            error_text, status = "", "ok"
        except Exception as error:
            search, held = None, members[members.angle_deg == held_angle]
            metrics = {
                "matches": 0, "weighted_fraction": np.nan, "indexed_fraction": np.nan,
                "median_delta_q": np.nan, "p90_delta_q": np.nan,
            }
            error_text, status = str(error), "fit_failed"
        rows.append({
            "held_angle_deg": held_angle,
            "trained_normal_hkl": "" if search is None else str(search["best"]["hkl"]),
            "held_features": int(len(held)),
            "predicted_matches": int(metrics["matches"]),
            "predicted_weighted_fraction": float(metrics["weighted_fraction"]),
            "predicted_indexed_fraction": float(metrics["indexed_fraction"]),
            "predicted_median_delta_q": metrics["median_delta_q"],
            "predicted_p90_delta_q": metrics["p90_delta_q"],
            "training_score_margin": np.nan if search is None else search.get("score_margin", np.nan),
            "validation_method": method,
            "validation_scope": "primary_orientation_solver_only",
            "held_angle_used_in_training": False,
            "holdout_status": status,
            "holdout_error": error_text,
        })
    return pd.DataFrame(rows)


def _plot_classic_overlay_tables(image, indexed, ignored, output, title, config,
                                 key_filename=None, dpi_override=None):
    """Write the alternate filename for the detailed two-panel indexed/ignored overlay."""
    return _plot_binary_overlay_tables(
        image, indexed, ignored, output, title, config,
        key_filename=key_filename, dpi_override=dpi_override,
    )


def _parameter_audit(config, output):
    rows = []
    tokens = ("tolerance", "threshold", "percentile", "penalty", "weight", "margin", "max_", "min_")
    for field in fields(config):
        if any(token in field.name.lower() for token in tokens):
            rows.append({"parameter": field.name, "value": repr(getattr(config, field.name)),
                         "audit_note": "threshold-like configuration; validate against independent datasets"})
    frame = pd.DataFrame(rows)
    path = Path(output) / "scientific_parameter_audit.csv"
    frame.to_csv(path, index=False)
    return path.resolve(), frame


def _write_output_inventory(results):
    output = Path(results["output_dir"])
    rows = []
    for path in sorted(p for p in output.rglob("*") if p.is_file()):
        rel = path.relative_to(output).as_posix()
        suffix = path.suffix.lower()
        category = "image" if suffix in {".png", ".jpg", ".jpeg"} else "table" if suffix == ".csv" else "report" if suffix in {".html", ".json"} else "other"
        rows.append({"relative_path": rel, "category": category, "size_bytes": path.stat().st_size})
    frame = pd.DataFrame(rows)
    target = output / "output_inventory.csv"
    self_row = pd.DataFrame([{"relative_path": target.name, "category": "table", "size_bytes": 0}])
    frame = pd.concat([frame, self_row], ignore_index=True)
    frame.to_csv(target, index=False)
    frame.loc[frame.relative_path.eq(target.name), "size_bytes"] = target.stat().st_size
    frame.to_csv(target, index=False)
    return target.resolve(), frame


def _append_inventory_to_html(results, inventory):
    report = results.get("html_report")
    if not report or not Path(report).is_file():
        return
    path = Path(report)
    text = path.read_text(encoding="utf-8")
    start, end = "<!-- V81_OUTPUT_INDEX_START -->", "<!-- V81_OUTPUT_INDEX_END -->"
    if start in text and end in text:
        text = text.split(start)[0] + text.split(end, 1)[1]
    links = []
    root = Path(results["output_dir"])
    for row in inventory.itertuples(index=False):
        rel = html.escape(str(row.relative_path))
        links.append(f'<tr><td>{html.escape(str(row.category))}</td><td><a href="{rel}">{rel}</a></td><td>{int(row.size_bytes):,}</td></tr>')
    section = (start + '<h2>Complete output inventory</h2><p>This table links every file generated by the run.</p>'
               '<table><tr><th>type</th><th>artifact</th><th>bytes</th></tr>' + ''.join(links) + '</table>' + end)
    text = text.replace('</body></html>', section + '</body></html>')
    path.write_text(text, encoding="utf-8")


def run_user_friendly(config: IndexingConfig, run_options=None):
    results = _v803_run_user_friendly(config, run_options)
    if not results or results.get("status"):
        return results
    summary = results.get("summary", pd.DataFrame()).copy()
    if not summary.empty:
        summary["model_selection_subset_is_independent_validation"] = False
        summary["true_image_holdout_is_independent_validation"] = bool(getattr(config, "enable_true_image_holdout", False))
        if "reliability_flags" in summary:
            summary["fit_parameter_at_boundary"] = summary.reliability_flags.astype(str).str.contains("parameter_at_boundary")
        results["summary"] = summary
        summary.to_csv(Path(results["output_dir"]) / "indexing_summary_all_series.csv", index=False)
    if bool(getattr(config, "write_scientific_parameter_audit", True)):
        audit_path, audit = _parameter_audit(config, results["output_dir"])
        results["scientific_parameter_audit"] = audit_path
        results["threshold_like_parameter_count"] = int(len(audit))
    if bool(getattr(config, "write_output_inventory", True)):
        inventory_path, inventory = _write_output_inventory(results)
        results["output_inventory"] = inventory_path
        _append_inventory_to_html(results, inventory)
    return results


def _v88_apply_workflow_preset(config: IndexingConfig, name: str):
    requested = str(name).strip().lower()
    if requested not in {"preview", "recommended", "maximum_coverage"}:
        raise ValueError("internal workflow preset must be preview, recommended, or maximum_coverage")
    configured, options = _v803_apply_workflow_preset(config, requested)
    configured = replace(
        configured,
        true_holdout_exact_primary_search=(requested == "recommended"),
        consensus_prevent_chaining=True,
        strict_numerical_input=False,
        overlay_static_style="detailed",
        overlay_write_detailed_companion=False,
        overlay_show_coordinate_key=True,
        overlay_label_all_indexed=True,
        overlay_panel_spacing=0.08,
        overlay_fill_display_gaps=False,
        overlay_classic_use_representative_image=False,
    )
    return configured, options


# ==================== INDEXING RELIABILITY AND MODEL CHECKS ====================
# Reliability calculations extend indexing, calibration, model selection, and
# evidence reporting while the registered overlay renderer is protected from
# unintended changes that could alter visual q-space registration.

@dataclass
class V9Config(V81Config):
    """Advanced controls for indexing reliability, calibration checks, and model validation."""

    # Renderer-integrity guard. Function signatures identify the expected q-space
    # rendering and registration implementation; indexing stops if executable
    # renderer logic changes unexpectedly.
    freeze_overlay_renderer: bool = True

    # Empirical per-reflection stability under feature/calibration perturbations.
    assignment_stability_trials: int = 96
    assignment_stability_tilt_jitter_deg: float = 0.30
    assignment_stability_scale_jitter: float = 0.0030
    assignment_stability_offset_jitter_q: float = 0.0045
    assignment_stability_axis_jitter_deg: float = 0.18
    assignment_stability_shear_jitter: float = 0.0025
    assignment_stability_min_margin_q: float = 0.004
    bootstrap_exact_primary_search: bool = True
    stability_robust_threshold: float = 0.90
    stability_supported_threshold: float = 0.70
    stability_provisional_threshold: float = 0.50

    # Whole-solution decoy-orientation specificity test.
    solution_decoy_trials: int = 48
    solution_decoy_refine_top: int = 6
    solution_decoy_min_normal_separation_deg: float = 10.0
    solution_decoy_reject_p: float = 0.20
    solution_decoy_pass_p: float = 0.05

    # Cross-series calibration regularization. Shared detector scale and affine
    # terms are accepted only when peak assignments remain identical and fit quality
    # is not materially degraded.
    enable_shared_calibration_refit: bool = True

    # Held-angle multidomain model selection. An additional domain identified from
    # training angles is retained only when it improves prediction of excluded angles.
    enable_multidomain_holdout_selection: bool = True

    # Compare Euclidean and covariance-aware consensus methods using held-angle
    # prediction so the preferred rule is supported by predictive evidence.
    benchmark_consensus_methods: bool = True

    # Full candidate-CIF indexing. Alternative CIFs are searched with the same
    # orientation solver and held-angle checks, rather than radial screening only.
    enable_full_candidate_cif_indexing: bool = True

    # Four-state evidence classification summarizing indexing support and uncertainty.
    final_decision_min_assignment_stability: float = 0.70
    final_decision_reject_assignment_stability: float = 0.45
    final_decision_min_loo_fraction: float = 0.45
    final_decision_reject_loo_fraction: float = 0.25
    final_decision_min_bootstrap_stability: float = 0.80
    final_decision_reject_bootstrap_stability: float = 0.50


# AST-normalized signatures of the protected overlay functions. Formatting and
# comments do not affect these hashes; executable rendering logic does.
_V9_FROZEN_OVERLAY_SIGNATURES = {
    "_warp_registered_image": "068cf52ac2d11b5c8c4eb47f08d6c03efd0c4cc3199c3a20b1c45eea2204747d",
    "registered_composite_image": "0d6b8b461d14ba109260babf17711e49150824c3274518bc4ec60086824b65d5",
    "_plot_binary_overlay_tables": "c62ba4bc91f44625392fd06bb55cbe6198fd948d25d9fdb2c3dc1b4424392d33",
    "plot_indexed_ignored_overlay": "da311336dda3283497f04e2106ffe4cfa1808a705a5a65db1a2eff1030ca7e01",
    "postprocess_registered_series": "a01d91a71cf08f51fc9a374b90588341007c6e664995adfb9097aad437e91f8b",
}


def _v9_function_ast_hash(function) -> str:
    import ast as _ast
    import hashlib as _hashlib
    import inspect as _inspect
    source = _inspect.getsource(function)
    node = _ast.parse(source).body[0]
    normalized = _ast.dump(node, annotate_fields=True, include_attributes=False)
    return _hashlib.sha256(normalized.encode("utf-8")).hexdigest()


def verify_overlay_renderer_frozen(config, output_dir=None):
    """Verify the protected overlay implementation before indexing.

    The renderer is part of the q-space registration workflow. Detecting an
    unexpected executable change before analysis prevents accidental visual
    coordinate changes from being mistaken for scientific indexing differences.
    """
    rows = []
    namespace = globals()
    for name, expected in _V9_FROZEN_OVERLAY_SIGNATURES.items():
        function = namespace.get(name)
        actual = _v9_function_ast_hash(function) if callable(function) else "missing"
        rows.append({
            "function": name,
            "expected_ast_sha256": expected,
            "actual_ast_sha256": actual,
            "unchanged": actual == expected,
        })
    frame = pd.DataFrame(rows)
    if output_dir is not None:
        target = Path(output_dir)
        target.mkdir(parents=True, exist_ok=True)
        frame.to_csv(target / "overlay_renderer_regression_guard.csv", index=False)
    if bool(getattr(config, "freeze_overlay_renderer", True)) and not frame.unchanged.all():
        changed = frame.loc[~frame.unchanged, "function"].tolist()
        raise RuntimeError(
            "Frozen overlay renderer changed unexpectedly: " + ", ".join(changed)
        )
    return frame


    # Run the active multiscale detector and record an audit of feature evidence
    # without modifying the detector's scientific output.


def _high_q_corner_recovery_features(image: dict, config: IndexingConfig) -> pd.DataFrame:
    """Conservatively recover weak peaks in the high-q upper/right field.

    The normal detector uses one global threshold.  On rendered/converted
    GIWAXS images, weak high-q spots sit on a bright, slowly varying background
    and can fall below that global threshold.  This pass uses a local
    background/noise estimate, but recovered points remain low-weight and still
    must survive multi-angle consensus and the normal validation machinery.
    """
    columns = [
        "qr", "qz", "cov_rr", "cov_rz", "cov_zz", "sigma_qr", "sigma_qz",
        "major_width_q", "minor_width_q", "aspect_ratio", "orientation_deg",
        "feature_type", "strength", "snr", "integrated_signal", "source_detector",
        "pixel_x", "pixel_y",
    ]
    if not bool(getattr(config, "enable_high_q_corner_recovery", True)):
        return pd.DataFrame(columns=columns)
    intensity = np.asarray(image.get("intensity"), dtype=float)
    valid = np.asarray(image.get("valid", np.isfinite(intensity)), dtype=bool)
    qr = np.asarray(image.get("qr"), dtype=float)
    qz = np.asarray(image.get("qz"), dtype=float)
    if intensity.ndim != 2 or not qr.size or not qz.size or not valid.any():
        return pd.DataFrame(columns=columns)

    # A tighter DoG filter resolves the small weak spots visible at high q.
    small = gaussian_filter(intensity, 0.90)
    broad = gaussian_filter(intensity, 8.0)
    signal = small - broad
    local_mean = gaussian_filter(signal, 18.0)
    local_variance = gaussian_filter((signal - local_mean) ** 2, 18.0)
    local_sigma = np.sqrt(np.maximum(local_variance, 1e-12))
    local_z = (signal - local_mean) / np.maximum(local_sigma, 1e-12)

    qr_grid, qz_grid = np.meshgrid(qr, qz)
    qr_min, qr_max = float(np.nanmin(qr)), float(np.nanmax(qr))
    qz_min, qz_max = float(np.nanmin(qz)), float(np.nanmax(qz))
    max_radial = math.hypot(max(abs(qr_min), abs(qr_max)), max(abs(qz_min), abs(qz_max)))
    qr_gate = qr_min + float(getattr(config, "high_q_corner_qr_fraction", 0.55)) * (qr_max - qr_min)
    qz_gate = qz_min + float(getattr(config, "high_q_corner_qz_fraction", 0.45)) * (qz_max - qz_min)
    radial_gate = float(getattr(config, "high_q_corner_radial_fraction", 0.62)) * max_radial
    corner = (qr_grid >= qr_gate) & (qz_grid >= qz_gate) & (np.hypot(qr_grid, qz_grid) >= radial_gate)

    snr_gate = float(getattr(config, "high_q_corner_snr_threshold", 2.0))
    spacing = max(5, 2 * max(2, int(config.min_feature_spacing_px) // 2) + 1)
    peaks = (
        valid & corner
        & (signal == maximum_filter(signal, size=spacing))
        & (local_z >= snr_gate)
    )
    candidates = []
    radius = max(2, int(config.subpixel_radius_px))
    for y0, x0 in zip(*np.where(peaks)):
        y1, y2 = max(0, y0 - radius), min(signal.shape[0], y0 + radius + 1)
        x1, x2 = max(0, x0 - radius), min(signal.shape[1], x0 + radius + 1)
        region = np.zeros_like(valid, dtype=bool)
        region[y1:y2, x1:x2] = valid[y1:y2, x1:x2]
        threshold = float(local_mean[y0, x0] + snr_gate * local_sigma[y0, x0])
        feature = _covariance_feature(
            signal, region, image, threshold, "spot",
            "high_q_corner_adaptive", float(x0), float(y0),
        )
        if feature is None:
            continue
        feature["snr"] = float(local_z[y0, x0])
        # Keep this completion evidence below the weight of strong ordinary
        # detections so it cannot take over the primary orientation search.
        feature["integrated_signal"] = 0.55 * float(feature["integrated_signal"])
        candidates.append(feature)

    if not candidates:
        return pd.DataFrame(columns=columns)
    frame = pd.DataFrame(candidates)
    maximum = max(float(frame["integrated_signal"].max()), 1e-12)
    frame["strength"] = np.clip(frame["integrated_signal"] / maximum, 0, 1) * 0.55
    frame["_rank"] = frame["strength"] * np.log1p(np.maximum(frame["snr"], 0.0))
    kept = []
    tolerance = max(0.012, 0.75 * float(config.feature_merge_tolerance_q))
    for row in frame.sort_values("_rank", ascending=False, kind="mergesort").itertuples(index=False):
        if all(math.hypot(float(row.qr) - float(old.qr), float(row.qz) - float(old.qz)) > tolerance for old in kept):
            kept.append(row)
        if len(kept) >= int(getattr(config, "high_q_corner_max_features", 40)):
            break
    return pd.DataFrame([row._asdict() for row in kept]).drop(columns="_rank", errors="ignore").reindex(columns=columns)


def detect_features(image: dict, config: IndexingConfig) -> pd.DataFrame:
    """Active multiscale detector plus conservative high-q edge recovery."""
    result = _v88_detect_features(image, config)
    if "detector_mode" not in result.columns:
        result = result.copy()
        result["detector_mode"] = (
            "multiscale"
            if bool(getattr(config, "enable_multiscale_feature_detection", False))
            else "single_scale"
        )
    recovered = _high_q_corner_recovery_features(image, config)
    if recovered.empty:
        return result

    # Remove recovered points already represented by the ordinary detector.
    if not result.empty:
        tree = cKDTree(result[["qr", "qz"]].to_numpy(float))
        distance, _ = tree.query(recovered[["qr", "qz"]].to_numpy(float), k=1)
        recovered = recovered.loc[
            distance > max(0.012, 0.75 * float(config.feature_merge_tolerance_q))
        ].copy()
    if recovered.empty:
        return result

    for column in result.columns:
        if column not in recovered.columns:
            recovered[column] = np.nan
    recovered["detector_mode"] = "high_q_corner_adaptive"
    combined = pd.concat([result, recovered[result.columns]], ignore_index=True, sort=False)
    return combined.reset_index(drop=True)


def full_search_bootstrap(crystal, members, base_search, config):
    summary = {
        "orientation_stability": np.nan,
        "completed_iterations": 0,
    }
    if bool(getattr(config, "bootstrap_exact_primary_search", False)):
        summary["bootstrap_search_mode"] = "exact_primary_v7_research"
    return pd.DataFrame(), pd.DataFrame(), summary
















def _analyze_single_series(series_id, group, representative_image, image_count, crystal, config, output):
    result = _v88_analyze_single_series(
        series_id, group, representative_image, image_count, crystal, config, output
    )
    series_dir = Path(output) / str(series_id).replace(":", "_series_")

    # These optional model-selection analyses are disabled in the selectable GUI
    # workflows. Empty tables are still written so downstream reporting sees a
    # consistent output schema and can distinguish “not run” from missing output.
    consensus_benchmark = pd.DataFrame()
    consensus_benchmark.to_csv(series_dir / "consensus_method_benchmark.csv", index=False)

    domain_holdout = pd.DataFrame()
    domain_holdout.to_csv(series_dir / "multidomain_exact_holdout_selection.csv", index=False)
    result["multidomain_holdout_selection"] = domain_holdout
    result["summary"]["multidomain_selected_by_holdout"] = True
    return result






def _v9_solution_decoy_test(series_id, result, crystal, config):
    """Estimate whole-solution specificity against physically distinct orientations."""
    search = result["search"]
    primary = search["best"]
    candidates = orientation_candidates(
        crystal, replace(config, max_orientation_candidates=max(
            config.max_orientation_candidates,
            int(config.solution_decoy_trials) * 3,
        ))
    )
    decoys = [
        candidate for candidate in candidates
        if normal_angle(candidate, primary["hkl"], crystal)
        >= float(config.solution_decoy_min_normal_separation_deg)
    ]
    rng = np.random.default_rng(int(config.random_seed) + 2901)
    if len(decoys) > int(config.solution_decoy_trials):
        selected = rng.choice(len(decoys), size=int(config.solution_decoy_trials), replace=False)
        decoys = [decoys[int(index)] for index in selected]
    rows = []
    base_params = np.asarray(primary["params"], float)
    for hkl in decoys:
        best_eval = None
        for tx in config.coarse_tilt_values_deg:
            for ty in config.coarse_tilt_values_deg:
                params = base_params.copy()
                params[0], params[1] = float(tx), float(ty)
                evaluated = _v7_evaluate(
                    hkl, params, crystal, search["fit_features"],
                    search["anchors"], search["validation"], config,
                    materialize=False, raw_features=None,
                )
                if best_eval is None or evaluated["cv_score"] > best_eval["cv_score"]:
                    best_eval = evaluated
        rows.append({
            "normal_h": int(hkl[0]), "normal_k": int(hkl[1]), "normal_l": int(hkl[2]),
            "normal_hkl": str(tuple(hkl)),
            "normal_angle_from_primary_deg": normal_angle(hkl, primary["hkl"], crystal),
            "coarse_decoy_score": float(best_eval["cv_score"]),
            "refined_decoy_score": np.nan,
            "params": best_eval["params"],
            "solution": best_eval,
        })
    rows.sort(key=lambda item: item["coarse_decoy_score"], reverse=True)
    for row in rows[:max(0, int(config.solution_decoy_refine_top))]:
        refined = _v7_refine(
            row["solution"], crystal, search["fit_features"],
            search["anchors"], search["validation"], config,
            raw_features=result.get("members"), fast=False,
        )
        row["refined_decoy_score"] = float(refined["cv_score"])
        row["solution"] = refined
    primary_score = float(primary["cv_score"])
    final_scores = np.array([
        row["refined_decoy_score"] if np.isfinite(row["refined_decoy_score"])
        else row["coarse_decoy_score"] for row in rows
    ], float)
    exceed = int((final_scores >= primary_score).sum()) if len(final_scores) else 0
    p_value = float((1 + exceed) / (1 + len(final_scores))) if len(final_scores) else np.nan
    best_decoy = float(np.max(final_scores)) if len(final_scores) else np.nan
    output_rows = []
    for row, score in zip(rows, final_scores):
        output_rows.append({
            "series_id": series_id,
            "normal_h": row["normal_h"], "normal_k": row["normal_k"], "normal_l": row["normal_l"],
            "normal_hkl": row["normal_hkl"],
            "normal_angle_from_primary_deg": row["normal_angle_from_primary_deg"],
            "coarse_decoy_score": row["coarse_decoy_score"],
            "refined_decoy_score": row["refined_decoy_score"],
            "final_decoy_score": float(score),
            "beats_primary": bool(score >= primary_score),
        })
    summary = {
        "series_id": series_id,
        "primary_score": primary_score,
        "decoy_trials_completed": int(len(final_scores)),
        "decoys_beating_primary": exceed,
        "empirical_decoy_p_value": p_value,
        "best_decoy_score": best_decoy,
        "primary_minus_best_decoy_score": (
            float(primary_score - best_decoy) if np.isfinite(best_decoy) else np.nan
        ),
        "p_value_is_calibrated_statistical_p_value": False,
        "interpretation": "empirical specificity under the configured orientation-decoy generator",
    }
    return pd.DataFrame(output_rows), summary, [row["solution"] for row in rows]


def _v9_prediction_solutions(result):
    solutions = [{
        "domain": "primary",
        "hkl": result["search"]["best"]["hkl"],
        "params": np.asarray(result["search"]["best"]["params"], float),
    }]
    secondary = result.get("secondary") or {}
    if secondary.get("accepted"):
        for item in secondary.get("domain_solutions", []):
            if str(item.get("domain", "primary")) != "primary":
                solutions.append({
                    "domain": str(item.get("domain")),
                    "hkl": tuple(item["hkl"]),
                    "params": np.asarray(item["params"], float),
                })
    return solutions


def _v9_reflection_stability(series_id, result, crystal, config, decoy_solutions):
    assignments = result.get("unified_assignments")
    if assignments is None or assignments.empty:
        assignments = _unified_evidence_table(result)
    if assignments is None or assignments.empty:
        return pd.DataFrame()
    if "display_as_indexed" in assignments:
        assignments = assignments[assignments.display_as_indexed.astype(bool)].copy()
    if assignments.empty:
        return pd.DataFrame()
    for target, fallback in (("qr_exp", "qr"), ("qz_exp", "qz")):
        if target not in assignments and fallback in assignments:
            assignments[target] = assignments[fallback]
    assignments = assignments.dropna(subset=["qr_exp", "qz_exp", "h", "k", "l"]).copy()
    if assignments.empty:
        return pd.DataFrame()
    consensus = result.get("consensus", pd.DataFrame())
    consensus_lookup = consensus.set_index("feature_id") if not consensus.empty and "feature_id" in consensus else None
    bootstrap = result.get("bootstrap_assignments", pd.DataFrame())
    bootstrap_lookup = {}
    if bootstrap is not None and not bootstrap.empty:
        for row in bootstrap.itertuples(index=False):
            bootstrap_lookup[(str(row.feature_id), int(row.h), int(row.k), int(row.l))] = float(row.bootstrap_frequency)
    solutions = _v9_prediction_solutions(result)
    rng = np.random.default_rng(int(config.random_seed) + 4103)
    counts = np.zeros(len(assignments), int)
    trials = max(1, int(config.assignment_stability_trials))
    assigned_keys = []
    for row in assignments.itertuples(index=False):
        domain = str(getattr(row, "orientation_domain", "primary"))
        qr_reference = float(getattr(row, "qr_calc", row.qr_exp))
        assigned_keys.append((domain, int(row.h), int(row.k), int(row.l), 1 if qr_reference >= 0 else -1))
    for _ in range(trials):
        prediction_frames = []
        for solution in solutions:
            params = np.asarray(solution["params"], float).copy()
            params[0:2] += rng.normal(0.0, config.assignment_stability_tilt_jitter_deg, 2)
            params[2:4] += rng.normal(0.0, config.assignment_stability_scale_jitter, 2)
            params[4:6] += rng.normal(0.0, config.assignment_stability_offset_jitter_q, 2)
            params[6] += rng.normal(0.0, config.assignment_stability_axis_jitter_deg)
            params[7] += rng.normal(0.0, config.assignment_stability_shear_jitter)
            array = _v7_prediction_array(
                crystal, solution["hkl"], config, params,
                config.v7_fit_f2_percentile, config.v7_max_fit_predictions,
                solution["domain"],
            )
            prediction_frames.append(_predictions_frame(array))
        predictions = pd.concat(prediction_frames, ignore_index=True, sort=False)
        if predictions.empty:
            continue
        coords = predictions[["qr", "qz"]].to_numpy(float)
        tree = cKDTree(coords)
        exp = assignments[["qr_exp", "qz_exp"]].to_numpy(float).copy()
        for index, row in enumerate(assignments.itertuples(index=False)):
            covariance = np.eye(2) * config.uncertainty_floor_q ** 2
            if consensus_lookup is not None and str(row.feature_id) in consensus_lookup.index:
                feature = consensus_lookup.loc[str(row.feature_id)]
                if isinstance(feature, pd.DataFrame):
                    feature = feature.iloc[0]
                covariance = np.array([
                    [float(feature.get("cov_rr", covariance[0, 0])), float(feature.get("cov_rz", 0.0))],
                    [float(feature.get("cov_rz", 0.0)), float(feature.get("cov_zz", covariance[1, 1]))],
                ]) + np.eye(2) * config.uncertainty_floor_q ** 2
            exp[index] += rng.multivariate_normal([0.0, 0.0], covariance)
        k_query = 2 if len(predictions) > 1 else 1
        distance, nearest = tree.query(exp, k=k_query)
        if k_query == 1:
            distance = distance[:, None]
            nearest = nearest[:, None]
        for index in range(len(assignments)):
            first = int(nearest[index, 0])
            prediction = predictions.iloc[first]
            key = (
                str(prediction.get("orientation_domain", "primary")),
                int(prediction.h), int(prediction.k), int(prediction.l),
                1 if float(prediction.qr) >= 0 else -1,
            )
            margin = (
                float(distance[index, 1] - distance[index, 0])
                if distance.shape[1] > 1 else np.inf
            )
            if (
                key == assigned_keys[index]
                and float(distance[index, 0]) <= float(config.v7_all_feature_tolerance_q)
                and margin >= float(config.assignment_stability_min_margin_q)
            ):
                counts[index] += 1
    perturbation_frequency = counts / float(trials)

    # Feature-level decoy specificity: fraction of decoy solutions unable to
    # explain the feature at the assigned location.
    decoy_match_counts = np.zeros(len(assignments), int)
    usable_decoys = decoy_solutions[:max(1, int(config.solution_decoy_trials))]
    for solution in usable_decoys:
        predictions = solution.get("predictions", pd.DataFrame())
        if predictions is None or predictions.empty:
            continue
        tree = cKDTree(predictions[["qr", "qz"]].to_numpy(float))
        distance, _ = tree.query(assignments[["qr_exp", "qz_exp"]].to_numpy(float))
        decoy_match_counts += distance <= float(config.v7_all_feature_tolerance_q)
    decoy_specificity = (
        1.0 - decoy_match_counts / max(1, len(usable_decoys))
        if usable_decoys else np.full(len(assignments), np.nan)
    )

    rows = []
    for index, row in enumerate(assignments.itertuples(index=False)):
        bootstrap_frequency = bootstrap_lookup.get(
            (str(row.feature_id), int(row.h), int(row.k), int(row.l)), np.nan
        )
        angle_support = np.nan
        if consensus_lookup is not None and str(row.feature_id) in consensus_lookup.index:
            feature = consensus_lookup.loc[str(row.feature_id)]
            if isinstance(feature, pd.DataFrame):
                feature = feature.iloc[0]
            angle_support = float(feature.get("support_fraction", np.nan))
        components = [
            (float(perturbation_frequency[index]), 0.45),
            (float(bootstrap_frequency), 0.25),
            (float(angle_support), 0.15),
            (float(decoy_specificity[index]), 0.15),
        ]
        usable = [(value, weight) for value, weight in components if np.isfinite(value)]
        stability = (
            sum(value * weight for value, weight in usable) / sum(weight for _, weight in usable)
            if usable else np.nan
        )
        if not np.isfinite(stability):
            tier = "unassessed"
        elif stability >= config.stability_robust_threshold:
            tier = "robust"
        elif stability >= config.stability_supported_threshold:
            tier = "supported"
        elif stability >= config.stability_provisional_threshold:
            tier = "provisional"
        else:
            tier = "unstable"
        rows.append({
            "series_id": series_id,
            "feature_id": str(row.feature_id),
            "h": int(row.h), "k": int(row.k), "l": int(row.l),
            "hkl": str(getattr(row, "hkl", f"({int(row.h)} {int(row.k)} {int(row.l)})")),
            "orientation_domain": str(getattr(row, "orientation_domain", "primary")),
            "evidence_stage": str(getattr(row, "evidence_stage", "")),
            "perturbation_assignment_frequency": float(perturbation_frequency[index]),
            "bootstrap_assignment_frequency": bootstrap_frequency,
            "angle_support_fraction": angle_support,
            "decoy_specificity_fraction": float(decoy_specificity[index]),
            "empirical_assignment_stability": stability,
            "stability_tier": tier,
            "stability_is_calibrated_probability": False,
            "stability_conditioning": "configured bootstrap/calibration/feature-jitter and orientation-decoy model",
        })
    return pd.DataFrame(rows)




def _v9_final_decision(summary_row, config):
    flags = []
    loo = float(summary_row.get("leave_one_angle_out_weighted_fraction", np.nan))
    bootstrap = float(summary_row.get("bootstrap_orientation_stability", np.nan))
    stability = float(summary_row.get("median_empirical_assignment_stability", np.nan))
    decoy_p = float(summary_row.get("empirical_decoy_p_value", np.nan))
    boundary = bool(summary_row.get("fit_parameter_at_boundary", False))
    external_status = str(summary_row.get("external_truth_status", "not_supplied"))
    external_passed = summary_row.get("external_truth_passed", np.nan)
    if external_status == "tested" and pd.notna(external_passed) and not bool(external_passed):
        return "REJECTED", ["external_truth_mismatch"]
    reject = (
        (np.isfinite(loo) and loo < config.final_decision_reject_loo_fraction)
        or (np.isfinite(bootstrap) and bootstrap < config.final_decision_reject_bootstrap_stability)
        or (np.isfinite(stability) and stability < config.final_decision_reject_assignment_stability)
        or (np.isfinite(decoy_p) and decoy_p > config.solution_decoy_reject_p)
    )
    if reject:
        if np.isfinite(loo) and loo < config.final_decision_reject_loo_fraction:
            flags.append("held_angle_prediction_failed")
        if np.isfinite(bootstrap) and bootstrap < config.final_decision_reject_bootstrap_stability:
            flags.append("orientation_bootstrap_failed")
        if np.isfinite(stability) and stability < config.final_decision_reject_assignment_stability:
            flags.append("reflection_assignments_unstable")
        if np.isfinite(decoy_p) and decoy_p > config.solution_decoy_reject_p:
            flags.append("decoy_orientations_not_specific")
        return "REJECTED", flags
    required = {
        "held_angle": np.isfinite(loo) and loo >= config.final_decision_min_loo_fraction,
        "bootstrap": np.isfinite(bootstrap) and bootstrap >= config.final_decision_min_bootstrap_stability,
        "assignment_stability": np.isfinite(stability) and stability >= config.final_decision_min_assignment_stability,
        "decoy_specificity": np.isfinite(decoy_p) and decoy_p <= config.solution_decoy_pass_p,
    }
    missing = [name for name, passed in required.items() if not passed]
    if missing:
        return "INCONCLUSIVE", [f"criterion_not_met:{name}" for name in missing]
    if boundary or str(summary_row.get("reliability_flags", "")) not in {"", "[]", "nan"}:
        if boundary:
            flags.append("fit_parameter_at_boundary")
        return "PASS_WITH_WARNINGS", flags
    return "PASS", []


def _v9_write_detection_audit(results):
    features = results.get("features", pd.DataFrame())
    if features is None or features.empty or "detector_mode" not in features:
        return pd.DataFrame()
    columns = [column for column in (
        "series_id", "scan", "angle_deg", "detector_mode",
        "single_scale_feature_count", "multiscale_feature_count",
        "single_multiscale_overlap_fraction",
    ) if column in features]
    if not columns:
        return pd.DataFrame()
    return features[columns].drop_duplicates().sort_values(
        [column for column in ("series_id", "angle_deg", "scan") if column in columns]
    ).reset_index(drop=True)


def run_gixs_indexing(config: IndexingConfig, run_synthetic_test=True):
    """Run the primary indexing analysis, then add reliability diagnostics and guarded refinement."""
    verify_overlay_renderer_frozen(config, config.output_dir)
    results = _v88_run_gixs_indexing(config, run_synthetic_test=run_synthetic_test)
    if not results or results.get("status"):
        return results
    output = Path(results["output_dir"])
    overlay_guard = verify_overlay_renderer_frozen(config, output)
    results["overlay_renderer_regression_guard"] = overlay_guard

    shared = pd.DataFrame()
    shared.to_csv(output / "shared_calibration_refit.csv", index=False)
    results["shared_calibration_refit"] = shared

    summary_updates = {}
    for series_id, result in results["series_results"].items():
        series_dir = output / series_id.replace(":", "_series_")
        decoy_table, decoy_summary, decoy_solutions = _v9_solution_decoy_test(
            series_id, result, results["crystal"], config
        )
        decoy_table.to_csv(series_dir / "solution_orientation_decoy_scores.csv", index=False)
        (series_dir / "solution_orientation_decoy_summary.json").write_text(
            json.dumps(decoy_summary, indent=2, default=str)
        )
        stability = _v9_reflection_stability(
            series_id, result, results["crystal"], config, decoy_solutions
        )
        stability.to_csv(series_dir / "reflection_assignment_stability.csv", index=False)
        if not stability.empty:
            robust_count = int(stability.stability_tier.eq("robust").sum())
            supported_count = int(stability.stability_tier.isin(["robust", "supported"]).sum())
            median_stability = float(stability.empirical_assignment_stability.median())
        else:
            robust_count = supported_count = 0
            median_stability = np.nan
        update = {
            **decoy_summary,
            "median_empirical_assignment_stability": median_stability,
            "robust_stable_reflections": robust_count,
            "supported_or_robust_stable_reflections": supported_count,
            "assignment_stability_trials": int(config.assignment_stability_trials),
        }
        result["reflection_assignment_stability"] = stability
        result["solution_decoy_scores"] = decoy_table
        result["summary"].update(update)
        summary_updates[series_id] = update

    cif_full = pd.DataFrame([{"status": "disabled"}])
    cif_full.to_csv(output / "full_candidate_cif_indexing_validation.csv", index=False)
    results["full_candidate_cif_indexing_validation"] = cif_full

    detector_audit = _v9_write_detection_audit(results)
    detector_audit.to_csv(output / "multiscale_feature_detector_audit.csv", index=False)
    results["multiscale_feature_detector_audit"] = detector_audit

    summary = results.get("summary", pd.DataFrame()).copy()
    if summary.empty:
        summary = pd.DataFrame([result["summary"] for result in results["series_results"].values()])
    else:
        for index, row in summary.iterrows():
            series_id = str(row.series_id)
            for key, value in summary_updates.get(series_id, {}).items():
                summary.loc[index, key] = value
            result_summary = results["series_results"][series_id]["summary"]
            for key in ("shared_calibration_applied", "selected_consensus_method",
                        "multidomain_selected_by_holdout", "multidomain_rejected_by_holdout"):
                if key in result_summary:
                    summary.loc[index, key] = result_summary[key]
    decisions, decision_flags = [], []
    for _, row in summary.iterrows():
        decision, flags = _v9_final_decision(row.to_dict(), config)
        decisions.append(decision)
        decision_flags.append(flags)
    summary["final_decision"] = decisions
    summary["final_decision_flags"] = decision_flags
    results["summary"] = summary
    summary.to_csv(output / "indexing_summary_all_series.csv", index=False)
    for _, row in summary.iterrows():
        series_id = str(row.series_id)
        result = results["series_results"][series_id]
        result["summary"].update(row.to_dict())
        series_dir = output / series_id.replace(":", "_series_")
        (series_dir / "indexing_summary.json").write_text(
            json.dumps(result["summary"], indent=2, default=str)
        )
    external_columns = [column for column in (
        "series_id", "external_truth_status", "external_truth_passed",
        "external_truth_angle_error_deg", "external_truth_source",
        "external_truth_was_blind", "validation_scope", "final_decision",
    ) if column in summary]
    if external_columns:
        summary[external_columns].to_csv(
            output / "external_orientation_truth_validation.csv", index=False
        )

    validation_manifest = {
        "overlay_frozen": bool(overlay_guard.unchanged.all()),
        "workflow_preset": str(getattr(config, "workflow_preset", "unspecified")),
        "bootstrap_iterations": int(config.full_bootstrap_iterations),
        "bootstrap_exact_primary_search": bool(getattr(config, "bootstrap_exact_primary_search", False)),
        "true_image_holdout": bool(getattr(config, "enable_true_image_holdout", False)),
        "exact_primary_holdout": bool(getattr(config, "true_holdout_exact_primary_search", False)),
        "multiscale_detection": bool(getattr(config, "enable_multiscale_feature_detection", False)),
        "covariance_aware_consensus": bool(getattr(config, "enable_covariance_aware_consensus", False)),
        "shared_calibration_refit": bool(getattr(config, "enable_shared_calibration_refit", False)),
        "multidomain_holdout_selection": bool(getattr(config, "enable_multidomain_holdout_selection", False)),
        "candidate_cif_full_indexing": bool(getattr(config, "enable_full_candidate_cif_indexing", False)),
        "external_truth_supplied": bool(getattr(config, "external_truth_normals", ())),
    }
    (output / "v9_validation_manifest.json").write_text(
        json.dumps(validation_manifest, indent=2, default=str)
    )
    return results


def apply_workflow_preset(config: IndexingConfig, name: str):
    requested = str(name).strip().lower().replace("-", "_").replace(" ", "_")

    if requested == "fast_test":
        configured, options = _v88_apply_workflow_preset(config, "recommended")
        configured = replace(
            configured,
            full_leave_one_angle_out=False,
            full_bootstrap_iterations=0,
            bootstrap_exact_primary_search=False,
            test_second_orientation=False,
            true_holdout_exact_primary_search=False,
            enable_true_image_holdout=False,
            enable_multiscale_feature_detection=False,
            enable_covariance_aware_consensus=False,
            benchmark_consensus_methods=False,
            index_ignored_features=False,
            guided_rescue_unclustered_features=False,
            assignment_stability_trials=8,
            solution_decoy_trials=6,
            enable_shared_calibration_refit=False,
            enable_multidomain_holdout_selection=False,
            enable_full_candidate_cif_indexing=False,
        )
        options = dict(options)
        options["workflow_preset"] = "fast_test"
        return configured, options

    if requested == "recommended":
        configured, options = _v88_apply_workflow_preset(config, "recommended")
        configured = replace(
            configured,
            # Standard mode now uses the stronger detector/consensus path.
            enable_multiscale_feature_detection=True,
            enable_covariance_aware_consensus=True,
            benchmark_consensus_methods=False,
            index_ignored_features=True,
            guided_rescue_unclustered_features=True,
            feature_threshold_mad=3.0,
            feature_quantile=0.978,
            ridge_threshold_mad=2.35,
            max_features_per_image=100,
            max_consensus_features=90,
            min_angle_support=2,
            ignored_index_tolerance_q=0.054,
            ignored_index_sigma_limit=2.8,
            v71_completion_max_tolerance_q=0.068,
            guided_rescue_min_angle_support=2,
            assignment_stability_trials=24,
            solution_decoy_trials=12,
            enable_shared_calibration_refit=False,
            enable_multidomain_holdout_selection=False,
            enable_full_candidate_cif_indexing=False,
        )
        options = dict(options)
        options["workflow_preset"] = "recommended"
        return configured, options

    if requested == "improved_coverage":
        configured, options = _v88_apply_workflow_preset(config, "maximum_coverage")
        configured = replace(
            configured,
            # Sensitive but still evidence-constrained feature coverage.
            enable_multiscale_feature_detection=True,
            enable_covariance_aware_consensus=True,
            benchmark_consensus_methods=False,
            feature_threshold_mad=2.8,
            feature_quantile=0.975,
            ridge_threshold_mad=2.2,
            max_features_per_image=120,
            max_consensus_features=110,
            min_angle_support=2,

            # Reconsider ignored consensus peaks after the orientation is fixed.
            index_ignored_features=True,
            ignored_index_tolerance_q=0.058,
            ignored_index_sigma_limit=3.0,
            ignored_index_min_support=2,

            # Search all visible calculated reflections at the fixed solution.
            v71_enable_full_reflection_completion=True,
            v71_completion_base_tolerance_q=0.055,
            v71_completion_max_tolerance_q=0.072,
            v71_completion_sigma_limit=3.2,
            v71_completion_min_support=2,
            v71_completion_min_member_angles=2,

            # Second weak-peak pass around unused calculated reflections.
            guided_rescue_unclustered_features=True,
            guided_rescue_tolerance_q=0.042,
            guided_rescue_sigma_limit=2.4,
            guided_rescue_min_angle_support=2,
            guided_rescue_max_ambiguity=2,
            guided_rescue_min_margin_sigma=0.35,
            guided_rescue_supported_min_angle_support=3,
            guided_rescue_promote_strong_two_angle=True,
            guided_rescue_two_angle_max_normalized_delta=2.10,
            guided_rescue_two_angle_max_delta_q=0.020,
            guided_rescue_two_angle_min_margin_sigma=0.65,

            # Keep improved coverage focused on coverage without expensive full validation.
            full_leave_one_angle_out=False,
            full_bootstrap_iterations=0,
            bootstrap_exact_primary_search=False,
            true_holdout_exact_primary_search=False,
            enable_true_image_holdout=False,
            test_second_orientation=True,
            assignment_stability_trials=16,
            solution_decoy_trials=10,
            enable_shared_calibration_refit=False,
            enable_multidomain_holdout_selection=False,
            enable_full_candidate_cif_indexing=False,
            series_worker_timeout_s=max(int(config.series_worker_timeout_s), 1800),
        )
        options = dict(options)
        options["workflow_preset"] = "improved_coverage"
        return configured, options

    if requested == "preview":
        configured, options = _v88_apply_workflow_preset(config, "preview")
        options = dict(options)
        options["workflow_preset"] = "preview"
        return configured, options

    raise ValueError(
        "workflow preset must be fast_test, recommended, improved_coverage, or preview"
    )

# ================= END INDEXING RELIABILITY AND MODEL CHECKS =================


# GUI/backend module: no top-level configuration and no automatic execution.
# Use gixs_gui_app.backend.service.run_gixs_job().

# Alias used by the service layer after flattening the package into one file.
engine_v9 = sys.modules[__name__]

# ================================ DATA MODELS =================================
@dataclass
class MeasurementInput:
    file: str = ""
    png_file: str = ""
    numerical_file: str = ""
    sample: str = "sample1"
    series: str = "A"
    series_id: str = ""
    angle_deg: float = 0.10
    scan: int = 1
    exposure_s: float = 1.0
    qr_min: float = -1.0
    qr_max: float = 2.2
    qz_min: float = -0.10
    qz_max: float = 2.72
    colormap: str = "jet"
    crop_x0: int | None = None
    crop_y0: int | None = None
    crop_x1: int | None = None
    crop_y1: int | None = None
    enabled: bool = True
    notes: str = ""

    def normalized_series_id(self) -> str:
        return self.series_id.strip() or f"{self.sample.strip()}:{self.series.strip() or 'A'}"

    def to_manifest_row(self) -> dict[str, Any]:
        row = asdict(self)
        row["series_id"] = self.normalized_series_id()
        return row


@dataclass
class GIXSRunRequest:
    cif_path: str
    measurements: list[MeasurementInput]
    output_dir: str
    workflow_preset: str = "improved_coverage"
    alternative_cif_paths: list[str] = field(default_factory=list)
    q_max: float = 2.8
    colormap: str = "jet"
    preview_only: bool = False
    prefer_numerical: bool = True
    run_in_parallel: bool = True
    auto_expand_q_max: bool = True
    high_q_corner_recovery: bool = True

    def to_json_dict(self) -> dict[str, Any]:
        return {
            **asdict(self),
            "measurements": [asdict(item) for item in self.measurements],
        }

    @classmethod
    def from_json_dict(cls, payload: dict[str, Any]) -> "GIXSRunRequest":
        data = dict(payload)
        data["measurements"] = [MeasurementInput(**row) for row in data.get("measurements", [])]
        return cls(**data)


@dataclass
class SeriesResult:
    series_id: str
    final_decision: str = ""
    orientation_hkl: str = ""
    indexed_count: int = 0
    stable_reflection_count: int = 0
    overlay_path: str = ""
    indexed_table_path: str = ""
    stability_table_path: str = ""
    decoy_summary_path: str = ""
    validation_path: str = ""


@dataclass
class GIXSRunResult:
    status: str
    output_dir: str
    summary_csv: str = ""
    html_report: str = ""
    inventory_csv: str = ""
    candidate_cif_csv: str = ""
    validation_manifest_json: str = ""
    warnings: list[str] = field(default_factory=list)
    series_results: list[SeriesResult] = field(default_factory=list)
    error: str = ""

    def to_json_dict(self) -> dict[str, Any]:
        return {
            **asdict(self),
            "series_results": [asdict(item) for item in self.series_results],
        }

    @classmethod
    def from_json_dict(cls, payload: dict[str, Any]) -> "GIXSRunResult":
        data = dict(payload)
        data["series_results"] = [SeriesResult(**row) for row in data.get("series_results", [])]
        return cls(**data)

# ================================ INPUT HELPERS ===============================
_GUI_FILENAME_RE = re.compile(
    _ENGINE_FILENAME_RE.pattern, _ENGINE_FILENAME_RE.flags
)


def parse_measurement_filename(path: str | Path, default_scan: int = 1) -> MeasurementInput:
    p = Path(path)
    compact = p.name.replace(" ", "")
    match = _GUI_FILENAME_RE.search(compact)
    sample, angle, exposure, scan = "sample1", 0.10, 1.0, default_scan
    if match:
        values = match.groupdict()
        sample = f"s{values['sample']}"
        angle = float(values['angle'])
        exposure = float(values['exposure'])
        scan = int(values['scan'])
    suffix = p.suffix.lower()
    numerical = suffix == ".npz"
    return MeasurementInput(
        file=str(p.resolve()),
        png_file=str(p.resolve()) if suffix == ".png" else "",
        numerical_file=str(p.resolve()) if numerical else "",
        sample=sample, series="A", angle_deg=angle, scan=scan, exposure_s=exposure,
    )


def write_manifest(measurements: list[MeasurementInput], path: str | Path) -> Path:
    target = Path(path).expanduser().resolve()
    target.parent.mkdir(parents=True, exist_ok=True)
    rows = [item.to_manifest_row() for item in measurements if item.enabled]
    if not rows:
        raise ValueError("At least one enabled measurement is required.")
    pd.DataFrame(rows).to_csv(target, index=False)
    return target


def validate_measurements(measurements: list[MeasurementInput]) -> list[str]:
    errors: list[str] = []
    enabled = [item for item in measurements if item.enabled]
    if not enabled:
        return ["At least one measurement must be enabled."]
    seen: set[tuple[str, float]] = set()
    for index, item in enumerate(enabled, 1):
        raw = item.numerical_file or item.png_file or item.file
        if not raw or not Path(raw).expanduser().is_file():
            errors.append(f"Measurement row {index} has no existing input file.")
        if not item.sample.strip() and not item.series_id.strip():
            errors.append(f"Measurement row {index} needs a sample or series_id.")
        if not (item.qr_min < item.qr_max and item.qz_min < item.qz_max):
            errors.append(f"Measurement row {index} has invalid q-axis limits.")
        key = (item.normalized_series_id(), round(float(item.angle_deg), 7))
        if key in seen:
            errors.append(f"Duplicate incident angle {item.angle_deg:g}° in {key[0]}.")
        seen.add(key)
    return errors

# =============================== BACKEND SERVICE ==============================
ProgressCallback = Callable[[str, float], None]
LogCallback = Callable[[str], None]


class CancellationToken(Protocol):
    def is_cancelled(self) -> bool: ...


class _CallbackWriter(io.TextIOBase):
    def __init__(self, callback: LogCallback | None):
        self.callback = callback
        self.buffer = ""
    def write(self, text: str) -> int:
        if not text:
            return 0
        self.buffer += text
        while "\n" in self.buffer:
            line, self.buffer = self.buffer.split("\n", 1)
            if self.callback and line.strip():
                current_out, current_err = sys.stdout, sys.stderr
                try:
                    sys.stdout, sys.stderr = sys.__stdout__, sys.__stderr__
                    self.callback(line.rstrip())
                finally:
                    sys.stdout, sys.stderr = current_out, current_err
        return len(text)
    def flush(self) -> None:
        if self.callback and self.buffer.strip():
            current_out, current_err = sys.stdout, sys.stderr
            try:
                sys.stdout, sys.stderr = sys.__stdout__, sys.__stderr__
                self.callback(self.buffer.rstrip())
            finally:
                sys.stdout, sys.stderr = current_out, current_err
        self.buffer = ""


def _progress(callback: ProgressCallback | None, stage: str, fraction: float) -> None:
    if callback:
        callback(stage, max(0.0, min(1.0, fraction)))


def _cancelled(token: CancellationToken | None) -> bool:
    return bool(token and token.is_cancelled())


def validate_request(request: GIXSRunRequest) -> list[str]:
    errors: list[str] = []
    cif = Path(request.cif_path).expanduser()
    if not cif.is_file():
        errors.append("Select an existing primary CIF file.")
    for path in request.alternative_cif_paths:
        if not Path(path).expanduser().is_file():
            errors.append(f"Alternative CIF does not exist: {path}")
    if not request.output_dir.strip():
        errors.append("Select an output directory.")
    errors.extend(validate_measurements(request.measurements))
    return errors



def _finite_scalar(value):
    try:
        array = np.asarray(value)
        if array.size != 1:
            return None
        number = float(array.reshape(-1)[0])
    except Exception:
        return None
    return number if np.isfinite(number) else None


def _npz_auto_physics_metadata(path: str | Path) -> dict[str, float]:
    """Read only small scalar/axis metadata needed by AUTO GIWAXS physics."""
    target = Path(path).expanduser()
    if target.suffix.lower() != ".npz" or not target.is_file():
        return {}
    aliases = {
        "wavelength_A": (
            "xray_wavelength_A", "wavelength_A", "wavelength_angstrom",
            "wavelength_angstroms", "lambda_A", "lambda_angstrom",
        ),
        "energy_keV": (
            "xray_energy_keV", "energy_keV", "beam_energy_keV", "photon_energy_keV",
        ),
        "energy_eV": ("xray_energy_eV", "energy_eV", "beam_energy_eV", "photon_energy_eV"),
        "incidence_angle_deg": (
            "incidence_angle_deg", "incident_angle_deg", "alpha_i_deg", "alpha_i",
        ),
        "film_delta": ("film_delta", "delta_film"),
        "film_beta": ("film_beta", "beta_film"),
        "substrate_delta": ("substrate_delta", "delta_substrate"),
        "substrate_beta": ("substrate_beta", "beta_substrate"),
        "film_thickness_A": ("film_thickness_A", "thickness_A", "film_thickness_angstrom"),
        "film_thickness_nm": ("film_thickness_nm", "thickness_nm"),
    }
    result: dict[str, float] = {}
    try:
        with np.load(target, allow_pickle=False) as archive:
            keys = {str(key).lower(): key for key in archive.files}
            for canonical, names in aliases.items():
                for name in names:
                    key = keys.get(name.lower())
                    if key is None:
                        continue
                    value = _finite_scalar(archive[key])
                    if value is not None:
                        result[canonical] = value
                        break
            # q-axis spacing is used only as diagnostic metadata; accessing these
            # one-dimensional arrays does not load the large intensity image.
            for canonical, names in {
                "qr_spacing": ("qr", "q_r", "qr_axis"),
                "qz_spacing": ("qz", "q_z", "qz_axis"),
            }.items():
                for name in names:
                    key = keys.get(name.lower())
                    if key is None:
                        continue
                    try:
                        axis = np.asarray(archive[key], float).reshape(-1)
                        axis = axis[np.isfinite(axis)]
                        if axis.size > 1:
                            spacing = float(np.median(np.abs(np.diff(axis))))
                            if np.isfinite(spacing) and spacing > 0:
                                result[canonical] = spacing
                                break
                    except Exception:
                        continue
    except Exception:
        return {}
    return result


def _consistent_auto_value(records: list[dict[str, float]], key: str, rtol=1e-5, atol=1e-12):
    values = [float(record[key]) for record in records if key in record and np.isfinite(record[key])]
    if not values:
        return None, "missing"
    reference = float(np.median(values))
    if all(np.isclose(value, reference, rtol=rtol, atol=atol) for value in values):
        return reference, "npz_metadata"
    return None, "inconsistent_across_npz"


def _apply_auto_giwaxs_physics(config: IndexingConfig, request: GIXSRunRequest) -> IndexingConfig:
    """Enable refraction/DWBA only when independently supplied physics supports it.

    AUTO never chooses a model because it improves a diffraction match and never
    invents optical constants.  Missing/invalid information causes a clean
    kinematic fallback, so the normal indexing workflows always remain usable.
    """
    enabled_measurements = [item for item in request.measurements if item.enabled]
    metadata_records = []
    for item in enabled_measurements:
        path = item.numerical_file or item.file
        metadata_records.append(_npz_auto_physics_metadata(path))

    source_notes = []
    fallback = []

    # Prefer explicit NPZ angle metadata when every file agrees. Otherwise use
    # the manifest/filename angles already shown in the measurement table; each
    # series later replaces the global representative value with its own median.
    incidence_angle, angle_source = _consistent_auto_value(
        metadata_records, "incidence_angle_deg", rtol=1e-6, atol=1e-7
    )
    metadata_angles = np.asarray([
        record.get("incidence_angle_deg", np.nan) for record in metadata_records
    ], float)
    metadata_angles = metadata_angles[np.isfinite(metadata_angles) & (metadata_angles > 0)]
    row_angles = np.asarray([float(item.angle_deg) for item in enabled_measurements], float)
    row_angles = row_angles[np.isfinite(row_angles) & (row_angles > 0)]
    if incidence_angle is not None and incidence_angle > 0:
        source_notes.append("incident angle from NPZ metadata")
        physics_angles = metadata_angles if metadata_angles.size else np.asarray([incidence_angle], float)
    else:
        incidence_angle = float(np.median(row_angles)) if row_angles.size else None
        physics_angles = row_angles
        if incidence_angle is not None:
            source_notes.append("incident angle from measurement rows")
        else:
            fallback.append("valid incident angle unavailable")

    wavelength, wavelength_source = _consistent_auto_value(metadata_records, "wavelength_A")
    if wavelength is None:
        energy_kev, energy_source = _consistent_auto_value(metadata_records, "energy_keV")
        if energy_kev is None:
            energy_ev, energy_ev_source = _consistent_auto_value(metadata_records, "energy_eV")
            if energy_ev is not None and energy_ev > 0:
                energy_kev = energy_ev / 1000.0
                energy_source = energy_ev_source
        if energy_kev is not None and energy_kev > 0:
            wavelength = 12.398419843320026 / float(energy_kev)
            wavelength_source = energy_source + " (converted from energy)"
    if wavelength is not None and wavelength > 0:
        source_notes.append("wavelength/energy from NPZ metadata")
    else:
        wavelength = None
        fallback.append("measurement wavelength/energy unavailable")

    parameter_values = {}
    for key in ("film_delta", "film_beta", "substrate_delta", "substrate_beta", "film_thickness_A"):
        value, source = _consistent_auto_value(metadata_records, key)
        if key == "film_thickness_A" and value is None:
            value_nm, source_nm = _consistent_auto_value(metadata_records, "film_thickness_nm")
            if value_nm is not None:
                value, source = 10.0 * value_nm, source_nm + " (nm→Å)"
        parameter_values[key] = value
        if value is not None:
            source_notes.append(f"{key} from NPZ metadata")
        elif source == "inconsistent_across_npz":
            fallback.append(f"{key} inconsistent across NPZ files")

    film_valid = all(
        parameter_values[key] is not None and parameter_values[key] >= 0
        for key in ("film_delta", "film_beta")
    )
    substrate_valid = all(
        parameter_values[key] is not None and parameter_values[key] >= 0
        for key in ("substrate_delta", "substrate_beta")
    )
    thickness_valid = (
        parameter_values["film_thickness_A"] is not None
        and parameter_values["film_thickness_A"] > 0
    )
    wavelength_valid = wavelength is not None and np.isfinite(wavelength) and wavelength > 0
    angle_valid = incidence_angle is not None and np.isfinite(incidence_angle) and incidence_angle > 0

    enable_refraction = bool(wavelength_valid and angle_valid and film_valid)
    enable_dwba = bool(enable_refraction and substrate_valid and thickness_valid)

    if not film_valid:
        fallback.append("film optical constants unavailable")
    if enable_refraction and not enable_dwba:
        if not substrate_valid:
            fallback.append("substrate optical constants unavailable for DWBA")
        if not thickness_valid:
            fallback.append("film thickness unavailable for DWBA")

    if enable_dwba:
        status = "AUTO: refraction position correction + DWBA intensity weighting"
    elif enable_refraction:
        status = "AUTO: refraction position correction; kinematic F² intensity"
    else:
        status = "AUTO: kinematic position + kinematic F² intensity"

    updates = {
        "giwaxs_physics_mode": "auto",
        "giwaxs_physics_status": status,
        "giwaxs_physics_parameter_source": "; ".join(dict.fromkeys(source_notes)),
        "giwaxs_physics_fallback_reason": "; ".join(dict.fromkeys(fallback)),
        "enable_refraction_position_correction": enable_refraction,
        "enable_dwba": enable_dwba,
        "giwaxs_incidence_angles_deg": tuple(
            float(value) for value in np.unique(np.round(physics_angles, 8))
        ) if len(physics_angles) else (),
    }
    if incidence_angle is not None:
        updates["incidence_angle_deg"] = float(incidence_angle)
    if wavelength_valid:
        updates["xray_wavelength_A"] = float(wavelength)
    for key, value in parameter_values.items():
        if value is not None:
            updates[key] = float(value)
    return replace(config, **updates)


def _make_config(request: GIXSRunRequest, manifest_path: Path):
    first = next(item for item in request.measurements if item.enabled)
    base = engine_v9.V9Config(
        cif_path=str(Path(request.cif_path).expanduser().resolve()),
        manifest_path=str(manifest_path),
        search_dirs=(str(manifest_path.parent),),
        output_dir=str(Path(request.output_dir).expanduser().resolve()),
        qr_range=(float(first.qr_min), float(first.qr_max)),
        qz_range=(float(first.qz_min), float(first.qz_max)),
        q_max=float(request.q_max),
        colormap=str(request.colormap),
        preview_only=bool(request.preview_only),
        prefer_numerical=bool(request.prefer_numerical),
        alternative_cif_paths=tuple(str(Path(p).expanduser().resolve()) for p in request.alternative_cif_paths),
        max_parallel_series_workers=0 if request.run_in_parallel else 1,
        freeze_overlay_renderer=True,
    )
    configured, options = engine_v9.apply_workflow_preset(base, request.workflow_preset)
    # These two GUI accuracy options are deliberately applied after the workflow
    # preset so a preset cannot silently turn them back off.
    configured.enable_high_q_corner_recovery = bool(request.high_q_corner_recovery)
    configured.high_q_corner_snr_threshold = 2.0
    configured.high_q_corner_max_features = 40
    configured.high_q_corner_radial_fraction = 0.62
    configured.high_q_corner_qr_fraction = 0.55
    configured.high_q_corner_qz_fraction = 0.45
    if bool(request.high_q_corner_recovery):
        configured.max_features_per_image = max(int(configured.max_features_per_image), 120)
    configured = _apply_auto_giwaxs_physics(configured, request)
    options = dict(options)
    options["giwaxs_physics_mode"] = configured.giwaxs_physics_mode
    options["giwaxs_physics_status"] = configured.giwaxs_physics_status
    return configured, options


def _series_result(output: Path, row: dict) -> SeriesResult:
    series_id = str(row.get("series_id", ""))
    series_dir = output / series_id.replace(":", "_series_")
    orientation = row.get("best_hkl", row.get("normal_hkl", row.get("hkl", "")))
    return SeriesResult(
        series_id=series_id,
        final_decision=str(row.get("final_decision", "")),
        orientation_hkl=str(orientation),
        indexed_count=int(row.get("indexed_count", row.get("n_indexed", 0)) or 0),
        stable_reflection_count=int(row.get("supported_or_robust_stable_reflections", 0) or 0),
        overlay_path=str(series_dir / "indexed_or_ignored_overlay.png"),
        indexed_table_path=str(series_dir / "overlay_indexed_features.csv"),
        stability_table_path=str(series_dir / "reflection_assignment_stability.csv"),
        decoy_summary_path=str(series_dir / "solution_orientation_decoy_summary.json"),
        validation_path=str(series_dir / "indexing_summary.json"),
    )


def _package_result(raw: dict, output: Path) -> GIXSRunResult:
    if raw.get("status") == "preview_only":
        return GIXSRunResult(
            status="preview",
            output_dir=str(output),
            html_report=str(raw.get("html_report") or output / "input_preview_report.html"),
            summary_csv=str(output / "input_preview_summary.csv"),
            validation_manifest_json=str(output / "input_preflight_report.csv"),
        )
    if raw.get("status"):
        return GIXSRunResult(status=str(raw.get("status")), output_dir=str(output), error=str(raw.get("message", "")))
    summary = raw.get("summary", pd.DataFrame())
    rows = summary.to_dict("records") if isinstance(summary, pd.DataFrame) else list(summary or [])
    result = GIXSRunResult(
        status="completed",
        output_dir=str(output),
        summary_csv=str(output / "indexing_summary_all_series.csv"),
        html_report=str(raw.get("html_report") or output / "indexing_report.html"),
        inventory_csv=str(raw.get("output_inventory") or output / "output_inventory.csv"),
        candidate_cif_csv=str(output / "full_candidate_cif_indexing_validation.csv"),
        validation_manifest_json=str(output / "v9_validation_manifest.json"),
        series_results=[_series_result(output, row) for row in rows],
    )
    return result


def run_gixs_job(
    request: GIXSRunRequest,
    progress_callback: ProgressCallback | None = None,
    log_callback: LogCallback | None = None,
    cancellation_token: CancellationToken | None = None,
) -> GIXSRunResult:
    errors = validate_request(request)
    if errors:
        return GIXSRunResult(status="invalid", output_dir=request.output_dir, warnings=errors, error="\n".join(errors))
    output = Path(request.output_dir).expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    request_path = output / "gui_run_request.json"
    request_path.write_text(json.dumps(request.to_json_dict(), indent=2), encoding="utf-8")
    try:
        _progress(progress_callback, "Preparing dataset manifest", 0.05)
        manifest_path = write_manifest(request.measurements, output / "gui_dataset_manifest.csv")
        if _cancelled(cancellation_token):
            return GIXSRunResult(status="cancelled", output_dir=str(output))
        _progress(progress_callback, "Verifying frozen overlay", 0.10)
        config, options = _make_config(request, manifest_path)
        engine_v9.verify_overlay_renderer_frozen(config, output)
        if _cancelled(cancellation_token):
            return GIXSRunResult(status="cancelled", output_dir=str(output))
        _progress(progress_callback, "Running GIXS indexing", 0.15)
        writer = _CallbackWriter(log_callback)
        with contextlib.redirect_stdout(writer), contextlib.redirect_stderr(writer):
            raw = engine_v9.run_user_friendly(config, options)
        writer.flush()
        if _cancelled(cancellation_token):
            return GIXSRunResult(status="cancelled", output_dir=str(output))
        _progress(progress_callback, "Packaging results", 0.95)
        result = _package_result(raw, output)
        (output / "gui_result.json").write_text(json.dumps(result.to_json_dict(), indent=2), encoding="utf-8")
        _progress(progress_callback, "Completed", 1.0)
        return result
    except Exception as exc:
        details = traceback.format_exc()
        (output / "gui_error.txt").write_text(details, encoding="utf-8")
        if log_callback:
            log_callback(details)
        result = GIXSRunResult(status="error", output_dir=str(output), error=str(exc))
        (output / "gui_result.json").write_text(json.dumps(result.to_json_dict(), indent=2), encoding="utf-8")
        return result

# ========================= SINGLE-FILE JOB CONTROLLER =========================
import json as _controller_json
import os as _controller_os
import subprocess as _controller_subprocess
import sys as _controller_sys
import threading as _controller_threading
from pathlib import Path as _ControllerPath
from queue import Queue as _ControllerQueue, Empty as _ControllerEmpty

PREFIX = "GIXS_EVENT "

class _ThreadCancellationToken:
    """Cooperative cancellation token for the in-process fallback worker."""
    def __init__(self):
        self._event = _controller_threading.Event()

    def cancel(self) -> None:
        self._event.set()

    def is_cancelled(self) -> bool:
        return self._event.is_set()


class JobController:
    """Run indexing without freezing the GUI.

    Normal saved ``.py`` execution uses a child Python process. When the code
    is executed as a PyCharm cell/console selection, ``__file__`` is undefined;
    in that case the controller automatically uses a background thread instead
    of failing when the Run indexing button is pressed.
    """
    def __init__(self):
        self.process: _controller_subprocess.Popen | None = None
        self.queue: _ControllerQueue[dict] = _ControllerQueue()
        self._reader: _controller_threading.Thread | None = None
        self._thread: _controller_threading.Thread | None = None
        self._cancel_token: _ThreadCancellationToken | None = None

    @property
    def running(self) -> bool:
        process_running = self.process is not None and self.process.poll() is None
        thread_running = self._thread is not None and self._thread.is_alive()
        return bool(process_running or thread_running)

    @staticmethod
    def _current_script_path() -> _ControllerPath | None:
        """Return this saved application file, or None for cell/console runs."""
        candidates = []
        module_file = globals().get("__file__")
        if module_file:
            candidates.append(module_file)
        if _controller_sys.argv:
            candidates.append(_controller_sys.argv[0])

        seen = set()
        for raw in candidates:
            try:
                path = _ControllerPath(raw).expanduser().resolve()
            except (TypeError, OSError):
                continue
            if path in seen or not path.is_file() or path.suffix.lower() != ".py":
                continue
            seen.add(path)
            try:
                header = path.read_text(encoding="utf-8", errors="ignore")[:12000]
            except OSError:
                continue
            if "class JobController" in header or "GIXS INDEXING WORKBENCH" in header:
                return path
        return None

    def start(self, request: GIXSRunRequest) -> None:
        if self.running:
            raise RuntimeError("A job is already running.")

        errors = validate_request(request)
        if errors:
            raise ValueError("\n".join(errors))

        output = _ControllerPath(request.output_dir).expanduser().resolve()
        output.mkdir(parents=True, exist_ok=True)
        script_path = self._current_script_path()

        if script_path is not None:
            request_path = output / "gui_worker_request.json"
            request_path.write_text(
                _controller_json.dumps(request.to_json_dict(), indent=2),
                encoding="utf-8",
            )
            self.process = _controller_subprocess.Popen(
                [_controller_sys.executable, str(script_path), "--gixs-worker", str(request_path)],
                stdout=_controller_subprocess.PIPE,
                stderr=_controller_subprocess.STDOUT,
                text=True,
                bufsize=1,
                cwd=str(script_path.parent),
                env=_controller_os.environ.copy(),
            )
            self._reader = _controller_threading.Thread(
                target=self._read_output, daemon=True
            )
            self._reader.start()
            self.queue.put({
                "kind": "log",
                "message": f"Started indexing worker process: {script_path.name}",
            })
            return

        # PyCharm's Run Cell / Python Console does not define __file__. Running
        # in a background thread keeps the GUI responsive and avoids requiring
        # a saved neighboring worker script.
        self.process = None
        self._cancel_token = _ThreadCancellationToken()
        self._thread = _controller_threading.Thread(
            target=self._run_in_process,
            args=(request,),
            daemon=True,
        )
        self._thread.start()
        self.queue.put({
            "kind": "log",
            "message": (
                "No saved script path was available, so indexing started in "
                "the PyCharm-safe background-thread mode."
            ),
        })

    def _run_in_process(self, request: GIXSRunRequest) -> None:
        try:
            result = run_gixs_job(
                request,
                progress_callback=lambda stage, fraction: self.queue.put({
                    "kind": "progress", "stage": stage, "fraction": fraction
                }),
                log_callback=lambda message: self.queue.put({
                    "kind": "log", "message": str(message)
                }),
                cancellation_token=self._cancel_token,
            )
            self.queue.put({"kind": "result", "result": result.to_json_dict()})
            code = 0 if result.status in {"completed", "preview", "cancelled"} else 1
        except Exception:
            self.queue.put({"kind": "log", "message": traceback.format_exc()})
            code = 1
        finally:
            self.queue.put({"kind": "exit", "code": code})

    def _read_output(self) -> None:
        assert self.process and self.process.stdout
        for line in self.process.stdout:
            line = line.rstrip("\n")
            if line.startswith(PREFIX):
                try:
                    self.queue.put(_controller_json.loads(line[len(PREFIX):]))
                except _controller_json.JSONDecodeError:
                    self.queue.put({"kind": "log", "message": line})
            elif line.strip():
                self.queue.put({"kind": "log", "message": line})
        code = self.process.wait()
        self.queue.put({"kind": "exit", "code": code})

    def cancel(self) -> None:
        if self.process is not None and self.process.poll() is None:
            self.process.terminate()
            try:
                self.process.wait(timeout=4)
            except _controller_subprocess.TimeoutExpired:
                self.process.kill()
            self.queue.put({"kind": "cancelled"})
            return

        if self._thread is not None and self._thread.is_alive():
            if self._cancel_token is not None:
                self._cancel_token.cancel()
            self.queue.put({
                "kind": "log",
                "message": "Cancellation requested; the active indexing stage will stop at its next cancellation checkpoint.",
            })
            self.queue.put({"kind": "cancelled"})

    def poll(self) -> list[dict]:
        events: list[dict] = []
        while True:
            try:
                events.append(self.queue.get_nowait())
            except _ControllerEmpty:
                break
        return events


def _emit_worker_event(kind: str, **payload) -> None:
    print(PREFIX + _controller_json.dumps({"kind": kind, **payload}, default=str), flush=True)


def _run_worker_request(request_json: str) -> int:
    payload = _controller_json.loads(_ControllerPath(request_json).read_text(encoding="utf-8"))
    request = GIXSRunRequest.from_json_dict(payload)
    result = run_gixs_job(
        request,
        progress_callback=lambda stage, fraction: _emit_worker_event(
            "progress", stage=stage, fraction=fraction
        ),
        log_callback=lambda message: _emit_worker_event("log", message=message),
    )
    _emit_worker_event("result", result=result.to_json_dict())
    return 0 if result.status in {"completed", "preview"} else 1

# ========================== PYQT6 INTEGRATED WORKBENCH ==========================
# Desktop interface for the supported indexing workflows and reciprocal-lattice
# comparison tools. Indexed/ignored classification images are loaded directly from
# backend output so the GUI does not independently redraw or reinterpret those
# scientific classifications.

# Worker processes execute numerical analysis only; GUI construction remains in
# the parent process so multiprocessing cannot open duplicate application windows.
if __name__ == "__main__" and len(sys.argv) >= 3 and sys.argv[1] == "--gixs-worker":
    raise SystemExit(_run_worker_request(sys.argv[2]))

from pathlib import Path as _GuiPath

import matplotlib.image as _mpimg
from matplotlib.figure import Figure as _MplFigure
from matplotlib.backends.backend_qtagg import (
    FigureCanvasQTAgg as _FigureCanvas,
    NavigationToolbar2QT as _NavigationToolbar,
)

from PyQt6.QtCore import Qt, QTimer
from PyQt6.QtGui import QColor
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QTabWidget, QSplitter, QScrollArea,
    QVBoxLayout, QHBoxLayout, QFormLayout, QGroupBox, QLabel,
    QLineEdit, QPushButton, QFileDialog, QMessageBox, QComboBox, QCheckBox,
    QCompleter, QDoubleSpinBox, QSpinBox, QProgressBar, QPlainTextEdit, QTableWidget,
    QTableWidgetItem, QHeaderView, QAbstractItemView,
)


# ---------------------------- manual simulation core ----------------------------

def _manual_reciprocal_basis(a, b, c, alpha_deg, beta_deg, gamma_deg):
    if min(a, b, c) <= 0:
        raise ValueError("a, b, and c must be positive.")
    alpha, beta, gamma = np.deg2rad([alpha_deg, beta_deg, gamma_deg])
    if not all(0 < value < np.pi for value in (alpha, beta, gamma)):
        raise ValueError("Cell angles must be between 0 and 180 degrees.")
    sg = math.sin(gamma)
    if abs(sg) < 1e-10:
        raise ValueError("gamma produces a singular unit cell.")
    avec = np.array([a, 0.0, 0.0], float)
    bvec = np.array([b * math.cos(gamma), b * sg, 0.0], float)
    cx = c * math.cos(beta)
    cy = c * (math.cos(alpha) - math.cos(beta) * math.cos(gamma)) / sg
    cz2 = c * c - cx * cx - cy * cy
    if cz2 <= 1e-12:
        raise ValueError("The lattice parameters produce a zero/invalid cell volume.")
    cvec = np.array([cx, cy, math.sqrt(cz2)], float)
    direct = np.column_stack([avec, bvec, cvec])
    return 2.0 * np.pi * np.linalg.inv(direct).T


def _rotation_align_vector_to_z(vector):
    source = np.asarray(vector, float)
    norm = float(np.linalg.norm(source))
    if norm <= 1e-12:
        raise ValueError("Preferred orientation (H,K,L) cannot be (0,0,0).")
    source = source / norm
    target = np.array([0.0, 0.0, 1.0])
    dot = float(np.clip(np.dot(source, target), -1.0, 1.0))
    if dot > 1.0 - 1e-12:
        return np.eye(3)
    if dot < -1.0 + 1e-12:
        return np.array([[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]])
    axis = np.cross(source, target)
    axis /= np.linalg.norm(axis)
    x, y, z = axis
    skew = np.array([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]])
    angle = math.acos(dot)
    return np.eye(3) + math.sin(angle) * skew + (1.0 - math.cos(angle)) * (skew @ skew)


def _manual_spacegroup_operations(identifier):
    text = str(identifier).strip()
    if not text:
        raise ValueError("Space group cannot be empty.")

    # The editable GUI selector displays entries as, for example,
    # "14 — P 1 21/c 1".  Accept that display form while continuing to accept
    # a bare International Tables number or a Hermann-Mauguin symbol.
    display_match = re.match(r"^\s*(\d{1,3})\s*(?:—|–|-|:|\|)\s*", text)
    if display_match:
        text = display_match.group(1)

    if text.lstrip("+-").isdigit():
        number = int(text)
        if not 1 <= number <= 230:
            raise ValueError("Space-group number must be 1 through 230.")
        group = gemmi.get_spacegroup_reference_setting(number)
    else:
        aliases = {"FCC": 225, "BCC": 229, "DIA": 227, "HCP": 194,
                   "GYR": 230, "SC": 221, "FCT": 123, "BCT": 139, "DBCC": 229}
        group = (gemmi.get_spacegroup_reference_setting(aliases[text.upper()])
                 if text.upper() in aliases else gemmi.find_spacegroup_by_name(text))
    if group is None:
        raise ValueError(f"Gemmi could not resolve space group {identifier!r}.")
    return group.operations()


def _manual_structure_intensity(structure, h, k, l, q_magnitude, calculator=None):
    """Return kinematic X-ray F² for the selected space-group reflection pattern.

    Gemmi's X-ray structure-factor calculator is used so the manual reflection
    pattern and automatic backend are based on the same physical scattering
    quantity. A conservative qualitative fallback is used only when the CIF
    cannot be evaluated by Gemmi.
    """
    if structure is None:
        return 1.0
    try:
        calculator = calculator or gemmi.StructureFactorCalculatorX(structure.cell)
        amplitude = calculator.calculate_sf_from_small_structure(
            structure, (int(h), int(k), int(l))
        )
        value = float(abs(amplitude) ** 2)
        if np.isfinite(value):
            return value
    except Exception:
        pass

    amplitude = 0j
    for site in structure.sites:
        z_number = float(site.element.atomic_number)
        occupancy = float(getattr(site, "occ", 1.0) or 1.0)
        if not np.isfinite(occupancy):
            occupancy = 1.0
        u_iso = float(getattr(site, "u_iso", 0.0) or 0.0)
        if not np.isfinite(u_iso):
            u_iso = 0.0
        damping = math.exp(-max(0.0, u_iso) * q_magnitude * q_magnitude / 2.0)
        phase = 2.0 * np.pi * (h * site.fract.x + k * site.fract.y + l * site.fract.z)
        amplitude += occupancy * z_number * damping * np.exp(1j * phase)
    return float(abs(amplitude) ** 2)


def calculate_manual_reflections(structure, lattice, orientation_hkl, spacegroup,
                                 hkl_max, q_max, dedup_tol=1e-3, physics_config=None):
    basis = _manual_reciprocal_basis(*lattice)
    h0, k0, l0 = map(int, orientation_hkl)
    rotation = _rotation_align_vector_to_z(basis @ np.array([h0, k0, l0], float))
    rotated = rotation @ basis
    operations = _manual_spacegroup_operations(spacegroup)
    sf_calculator = (
        gemmi.StructureFactorCalculatorX(structure.cell)
        if structure is not None else None
    )
    rows = []
    for h in range(-hkl_max, hkl_max + 1):
        for k in range(-hkl_max, hkl_max + 1):
            for l in range(-hkl_max, hkl_max + 1):
                if h == k == l == 0 or operations.is_systematically_absent([h, k, l]):
                    continue
                q = rotated @ np.array([h, k, l], float)
                magnitude = float(np.linalg.norm(q))
                if magnitude > q_max + 1e-9:
                    continue
                qr = float(math.hypot(q[0], q[1]))
                qz = float(q[2])
                rows.append({
                    "h": h, "k": k, "l": l, "hkl": f"({h},{k},{l})",
                    "Qxy": qr, "Qz": qz,
                    "Intensity": _manual_structure_intensity(
                        structure, h, k, l, magnitude, calculator=sf_calculator
                    ),
                })
    if not rows:
        return pd.DataFrame(columns=["h", "k", "l", "hkl", "Qxy", "Qz", "Intensity"])
    frame = pd.DataFrame(rows)
    frame["_qr_key"] = np.rint(frame.Qxy / dedup_tol).astype(np.int64)
    frame["_qz_key"] = np.rint(frame.Qz / dedup_tol).astype(np.int64)
    # Keep all contributing HKLs and sum their qualitative intensities.
    frame = frame.groupby(["_qr_key", "_qz_key"], sort=False, as_index=False).agg(
        Qxy=("Qxy", "mean"), Qz=("Qz", "mean"),
        hkl=("hkl", lambda values: ", ".join(dict.fromkeys(values))),
        Intensity=("Intensity", "sum"),
    )
    frame["QzKinematic"] = frame["Qz"].to_numpy(float)
    frame["KinematicIntensity"] = frame["Intensity"].to_numpy(float)
    frame["DWBAWeight"] = 1.0
    frame["PositionModel"] = "kinematic"
    frame["IntensityModel"] = "kinematic_F2"
    if physics_config is not None:
        if bool(getattr(physics_config, "enable_refraction_position_correction", False)):
            frame["Qz"] = refraction_corrected_qz(
                frame["QzKinematic"].to_numpy(float), physics_config
            )
            frame["PositionModel"] = "GIWAXS_refraction_corrected"
        if bool(getattr(physics_config, "enable_dwba", False)):
            weights = dwba_intensity_envelope(frame["Qz"].to_numpy(float), physics_config)
            frame["DWBAWeight"] = weights
            frame["Intensity"] = frame["KinematicIntensity"].to_numpy(float) * weights
            frame["IntensityModel"] = "DWBA_weighted_F2"
    return frame[[
        "hkl", "Qxy", "Qz", "Intensity", "QzKinematic", "KinematicIntensity",
        "DWBAWeight", "PositionModel", "IntensityModel",
    ]].sort_values(
        ["Qz", "Qxy"], ascending=[False, True], kind="mergesort"
    ).reset_index(drop=True)






class IntegratedGIXSWorkbench(QMainWindow):
    MEASUREMENT_HEADERS = [
        "Use", "File", "Sample", "Series", "Series ID", "Angle (deg)", "Scan",
        "Exposure (s)", "qr min", "qr max", "qz min", "qz max",
    ]

    def __init__(self):
        super().__init__()
        self.setWindowTitle("GIWAXS/GIXS Integrated Indexing Workbench")
        self.resize(1480, 900)
        self.setMinimumSize(1120, 720)
        self.controller = JobController()
        self.current_result = None
        self.alternative_cifs = []
        self.manual_structure = None
        self.manual_cif_spacegroup_number = None
        self.manual_cif_spacegroup_symbol = ""
        self.manual_image = None
        self.manual_image_original = None
        self.manual_image_crop_xyxy = None
        self.manual_image_path = ""
        self.manual_result = pd.DataFrame()
        self._build_ui()
        self.poll_timer = QTimer(self)
        self.poll_timer.timeout.connect(self._poll_controller)
        self.poll_timer.start(150)

    # ------------------------------ generic widgets ------------------------------
    @staticmethod
    def _double(value, minimum=-1e6, maximum=1e6, step=0.1, decimals=5):
        box = QDoubleSpinBox()
        box.setRange(minimum, maximum)
        box.setDecimals(decimals)
        box.setSingleStep(step)
        box.setValue(value)
        box.setKeyboardTracking(False)
        return box

    @staticmethod
    def _integer(value, minimum=-1000, maximum=1000):
        box = QSpinBox()
        box.setRange(minimum, maximum)
        box.setValue(value)
        return box

    @staticmethod
    def _path_row(form, label, line_edit, callback, button_text="Browse…"):
        holder = QWidget()
        layout = QHBoxLayout(holder)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(line_edit, 1)
        button = QPushButton(button_text)
        button.clicked.connect(callback)
        layout.addWidget(button)
        form.addRow(label, holder)
        return button

    def _build_ui(self):
        """Build the combined manual-overlay workspace and its NPZ measurement table.

        The visible controls support CIF reflection generation, numerical NPZ
        indexing, progress reporting, and experimental/calculated peak comparison
        without creating unused interface panels.
        """
        self.main_tabs = QTabWidget()
        self.setCentralWidget(self.main_tabs)
        self.manual_tab = QWidget()
        self.main_tabs.addTab(self.manual_tab, "Manual CIF indexing overlay")

        # Internal CIF-path storage used by the shared structure loader. The visible
        # path control is created separately in the manual indexing panel.
        self.cif_edit = QLineEdit()

        # Measurement table for numerical NPZ inputs used by the automatic calculator.
        # Creating it directly keeps all calculator inputs visible to the analyst.
        self.measurement_table = QTableWidget(0, len(self.MEASUREMENT_HEADERS))
        self.measurement_table.setHorizontalHeaderLabels(self.MEASUREMENT_HEADERS)
        self.measurement_table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.measurement_table.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.measurement_table.setAlternatingRowColors(True)
        self.measurement_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)
        for column in range(len(self.MEASUREMENT_HEADERS)):
            if column != 1:
                self.measurement_table.horizontalHeader().setSectionResizeMode(
                    column, QHeaderView.ResizeMode.ResizeToContents
                )

        self._build_manual_tab()

    # ---------------------------- automatic inputs tab ---------------------------





    def _append_measurement(self, measurement):
        row = self.measurement_table.rowCount()
        self.measurement_table.insertRow(row)
        use_item = QTableWidgetItem()
        use_item.setFlags(use_item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
        use_item.setCheckState(Qt.CheckState.Checked if measurement.enabled else Qt.CheckState.Unchecked)
        self.measurement_table.setItem(row, 0, use_item)
        values = [
            measurement.file or measurement.png_file or measurement.numerical_file,
            measurement.sample, measurement.series, measurement.series_id,
            measurement.angle_deg, measurement.scan, measurement.exposure_s,
            measurement.qr_min, measurement.qr_max, measurement.qz_min, measurement.qz_max,
        ]
        for column, value in enumerate(values, 1):
            item = QTableWidgetItem(str(value))
            if column == 1:
                item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            self.measurement_table.setItem(row, column, item)

    def _remove_measurements(self):
        rows = sorted({index.row() for index in self.measurement_table.selectionModel().selectedRows()}, reverse=True)
        for row in rows:
            self.measurement_table.removeRow(row)

    def _measurement_value(self, row, column, cast, default=None):
        item = self.measurement_table.item(row, column)
        text = item.text().strip() if item else ""
        if not text and default is not None:
            return default
        return cast(text)

    def _measurements_from_table(self):
        measurements = []
        for row in range(self.measurement_table.rowCount()):
            use_item = self.measurement_table.item(row, 0)
            path = self._measurement_value(row, 1, str)
            suffix = _GuiPath(path).suffix.lower()
            measurements.append(MeasurementInput(
                file=path,
                png_file=path if suffix in {".png", ".jpg", ".jpeg"} else "",
                numerical_file=path if suffix == ".npz" else "",
                sample=self._measurement_value(row, 2, str, "sample1"),
                series=self._measurement_value(row, 3, str, "A"),
                series_id=self._measurement_value(row, 4, str, ""),
                angle_deg=self._measurement_value(row, 5, float),
                scan=self._measurement_value(row, 6, int),
                exposure_s=self._measurement_value(row, 7, float),
                qr_min=self._measurement_value(row, 8, float),
                qr_max=self._measurement_value(row, 9, float),
                qz_min=self._measurement_value(row, 10, float),
                qz_max=self._measurement_value(row, 11, float),
                enabled=bool(use_item and use_item.checkState() == Qt.CheckState.Checked),
            ))
        return measurements





    # ------------------------------ manual tab ---------------------------------
    def _build_manual_tab(self):
        splitter = QSplitter(Qt.Orientation.Horizontal)
        layout = QHBoxLayout(self.manual_tab)
        layout.addWidget(splitter)
        controls_container = QWidget()
        controls_layout = QVBoxLayout(controls_container)
        controls_layout.setContentsMargins(6, 6, 6, 6)
        scroll = QScrollArea(); scroll.setWidgetResizable(True); scroll.setWidget(controls_container)
        splitter.addWidget(scroll)

        lattice_group = QGroupBox("Lattice and preferred orientation")
        form = QFormLayout(lattice_group)
        self.manual_cif_edit = QLineEdit()
        self._path_row(form, "CIF:", self.manual_cif_edit, self._choose_manual_cif)
        self.manual_a = self._double(4.0, 0.001, 1000.0, 0.1)
        self.manual_b = self._double(4.0, 0.001, 1000.0, 0.1)
        self.manual_c = self._double(4.0, 0.001, 1000.0, 0.1)
        self.manual_alpha = self._double(90.0, 0.01, 179.99, 1.0)
        self.manual_beta = self._double(90.0, 0.01, 179.99, 1.0)
        self.manual_gamma = self._double(90.0, 0.01, 179.99, 1.0)
        self.manual_h = self._integer(0, -30, 30)
        self.manual_k = self._integer(0, -30, 30)
        self.manual_l = self._integer(1, -30, 30)

        # Searchable selector containing the 230 standard space groups.  Each
        # item shows its International Tables number followed by its
        # Hermann-Mauguin symbol.  The line edit remains editable so typing a
        # number or any part of a symbol filters the completion list.
        self.manual_spacegroup = QComboBox()
        self.manual_spacegroup.setEditable(True)
        self.manual_spacegroup.setInsertPolicy(QComboBox.InsertPolicy.NoInsert)
        self.manual_spacegroup.setMaxVisibleItems(24)
        for spacegroup_number in range(1, 231):
            spacegroup = gemmi.get_spacegroup_reference_setting(spacegroup_number)
            display = f"{spacegroup_number} — {spacegroup.hm}"
            self.manual_spacegroup.addItem(display, spacegroup_number)
        self.manual_spacegroup.setCurrentIndex(0)
        self.manual_spacegroup.lineEdit().setPlaceholderText(
            "Type a number or Hermann-Mauguin symbol"
        )
        spacegroup_completer = self.manual_spacegroup.completer()
        spacegroup_completer.setCaseSensitivity(Qt.CaseSensitivity.CaseInsensitive)
        spacegroup_completer.setFilterMode(Qt.MatchFlag.MatchContains)
        spacegroup_completer.setCompletionMode(QCompleter.CompletionMode.PopupCompletion)
        for label, widget in [
            ("a (Å):", self.manual_a), ("b (Å):", self.manual_b), ("c (Å):", self.manual_c),
            ("alpha (deg):", self.manual_alpha), ("beta (deg):", self.manual_beta),
            ("gamma (deg):", self.manual_gamma), ("Preferred H:", self.manual_h),
            ("Preferred K:", self.manual_k), ("Preferred L:", self.manual_l),
            ("Space group:", self.manual_spacegroup),
        ]:
            form.addRow(label, widget)
        self.manual_spacegroup_guidance = QLabel("")
        self.manual_spacegroup_guidance.setWordWrap(True)
        form.addRow("Space-group consistency:", self.manual_spacegroup_guidance)
        self.manual_spacegroup.currentIndexChanged.connect(
            lambda _=None: self._update_spacegroup_consistency()
        )
        self.manual_spacegroup.lineEdit().textChanged.connect(
            lambda _=None: self._update_spacegroup_consistency()
        )
        controls_layout.addWidget(lattice_group)

        image_group = QGroupBox("Experimental image and display")
        image_form = QFormLayout(image_group)
        self.manual_image_edit = QLineEdit()
        self._path_row(image_form, "Image:", self.manual_image_edit, self._choose_manual_image)
        self.manual_qr_min = self._double(-1.0)
        self.manual_qr_max = self._double(2.2)
        self.manual_qz_min = self._double(-0.10)
        self.manual_qz_max = self._double(2.72)
        self.manual_opacity = self._double(0.85, 0.0, 1.0, 0.05, 2)
        self.manual_hkl_max = self._integer(10, 1, 40)
        self.manual_q_max = self._double(2.8, 0.1, 100.0, 0.1, 4)
        self.manual_mirror = QCheckBox("Mirror simulated points to negative qr")
        self.manual_mirror.setChecked(True)
        self.manual_labels = QCheckBox("Show HKL labels")
        self.manual_labels.setChecked(True)
        self.show_selected_spacegroup_structure = QCheckBox(
            "Show selected space-group reflection pattern"
        )
        self.show_selected_spacegroup_structure.setChecked(True)
        self.show_selected_spacegroup_structure.setToolTip(
            "Shows or hides the cyan reflection pattern generated from the manually selected "
            "space group, lattice parameters, and preferred orientation."
        )
        self.manual_auto_crop = QCheckBox("Automatically crop the q-space panel")
        self.manual_auto_crop.setChecked(True)
        self.manual_fit_frame = QCheckBox("Fit image exactly inside the q-space frame")
        self.manual_fit_frame.setChecked(True)
        self.manual_frame_status = QLabel("No image loaded.")
        self.manual_frame_status.setWordWrap(True)
        self.manual_reframe_button = QPushButton("Re-detect and frame image")
        self.manual_reframe_button.clicked.connect(self._reframe_manual_image)
        for label, widget in [
            ("qr minimum:", self.manual_qr_min), ("qr maximum:", self.manual_qr_max),
            ("qz minimum:", self.manual_qz_min), ("qz maximum:", self.manual_qz_max),
            ("Image opacity:", self.manual_opacity), ("Symmetric HKL limit:", self.manual_hkl_max),
            ("Calculation q max:", self.manual_q_max),
        ]:
            image_form.addRow(label, widget)
        image_form.addRow(self.manual_mirror)
        image_form.addRow(self.manual_labels)
        image_form.addRow(self.show_selected_spacegroup_structure)
        image_form.addRow(self.manual_auto_crop)
        image_form.addRow(self.manual_fit_frame)
        image_form.addRow(self.manual_reframe_button)
        image_form.addRow("Image framing:", self.manual_frame_status)
        controls_layout.addWidget(image_group)

        note = QLabel(
            "Generate the reciprocal-space reflection pattern for the manually selected space "
            "group using the current lattice parameters and preferred orientation. This is a "
            "diffraction-reflection pattern, not a drawing of atoms or real-space symmetry elements. "
            "The full automatic calculator remains a separate NPZ-based indexing calculation."
        )
        note.setWordWrap(True)
        controls_layout.addWidget(note)
        action_row = QHBoxLayout()
        calculate_button = QPushButton("Generate selected space-group reflection pattern")
        calculate_button.setToolTip(
            "Generate the cyan allowed-reflection pattern from the space group selected above."
        )
        calculate_button.clicked.connect(self._calculate_manual)
        action_row.addWidget(calculate_button)
        controls_layout.addLayout(action_row)
        self.manual_status = QLabel("Load a CIF and optional experimental image.")
        self.manual_status.setWordWrap(True)
        controls_layout.addWidget(self.manual_status)
        controls_layout.addStretch(1)

        plot_panel = QWidget(); plot_layout = QVBoxLayout(plot_panel)
        self.manual_figure = _MplFigure()
        self.manual_canvas = _FigureCanvas(self.manual_figure)
        self.manual_axes = self.manual_figure.add_subplot(111)
        plot_layout.addWidget(_NavigationToolbar(self.manual_canvas, self))
        plot_layout.addWidget(self.manual_canvas, 1)
        splitter.addWidget(plot_panel)
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        splitter.setSizes([380, 1050])
        for widget in (self.manual_qr_min, self.manual_qr_max, self.manual_qz_min,
                       self.manual_qz_max, self.manual_opacity):
            widget.valueChanged.connect(lambda _=None: self._redraw_manual())
        self.manual_mirror.stateChanged.connect(lambda _=None: self._redraw_manual())
        self.manual_labels.stateChanged.connect(lambda _=None: self._redraw_manual())
        self.show_selected_spacegroup_structure.stateChanged.connect(
            lambda _=None: self._redraw_manual()
        )
        self.manual_auto_crop.stateChanged.connect(lambda _=None: self._reframe_manual_image())
        self.manual_fit_frame.stateChanged.connect(lambda _=None: self._redraw_manual())
        self._redraw_manual()

    def _selected_manual_spacegroup(self):
        """Return (number, symbol) for the editable manual selector."""
        current_text = self.manual_spacegroup.currentText().strip()
        match = re.match(r"\s*(\d{1,3})", current_text)
        if match:
            number = int(match.group(1))
            if 1 <= number <= 230:
                group = gemmi.get_spacegroup_reference_setting(number)
                return number, str(group.hm)
        group = gemmi.find_spacegroup_by_name(current_text.split("—", 1)[-1].strip())
        if group is not None:
            return int(group.number), str(group.hm)
        try:
            data = self.manual_spacegroup.currentData()
            item_text = self.manual_spacegroup.currentText().strip()
            current_index = self.manual_spacegroup.currentIndex()
            listed_text = self.manual_spacegroup.itemText(current_index).strip() if current_index >= 0 else ""
            if data is not None and item_text == listed_text:
                number = int(data)
                group = gemmi.get_spacegroup_reference_setting(number)
                return number, str(group.hm)
        except Exception:
            pass
        return None, current_text

    def _update_spacegroup_consistency(self):
        label = getattr(self, "manual_spacegroup_guidance", None)
        if label is None:
            return
        selected_number, selected_symbol = self._selected_manual_spacegroup()
        cif_number = getattr(self, "manual_cif_spacegroup_number", None)
        cif_symbol = getattr(self, "manual_cif_spacegroup_symbol", "")
        notes = []
        color = "#555"
        if cif_number is None:
            notes.append("No CIF space group is loaded yet. The manual pattern will use the selected group.")
        elif selected_number == cif_number:
            color = "#207a3c"
            notes.append(f"Matches the loaded CIF: #{cif_number} — {cif_symbol}.")
        else:
            color = "#a35a00"
            notes.append(
                f"Manual override #{selected_number or '?'} — {selected_symbol} differs from the CIF "
                f"#{cif_number} — {cif_symbol}. Allowed reflections follow the manual selection, "
                "while atomic coordinates still come from the CIF."
            )
        if selected_number == 1:
            notes.append(
                "P1 allows nearly every reflection, so a dense visual overlap is not by itself proof of correct indexing."
            )
        label.setText(f"<span style='color:{color}'>{' '.join(notes)}</span>")

    def _choose_manual_cif(self):
        path, _ = QFileDialog.getOpenFileName(self, "Select CIF", "", "CIF files (*.cif);;All files (*)")
        if path:
            self._load_manual_cif(path, show_errors=True)
            if not self.cif_edit.text().strip():
                self.cif_edit.setText(path)

    def _load_manual_cif(self, path, show_errors=True):
        try:
            structure = gemmi.read_small_structure(str(path))
            self.manual_structure = structure
            self.manual_cif_edit.setText(str(path))
            self.manual_a.setValue(float(structure.cell.a))
            self.manual_b.setValue(float(structure.cell.b))
            self.manual_c.setValue(float(structure.cell.c))
            self.manual_alpha.setValue(float(structure.cell.alpha))
            self.manual_beta.setValue(float(structure.cell.beta))
            self.manual_gamma.setValue(float(structure.cell.gamma))
            group_text = str(getattr(structure, "spacegroup_hm", "") or "P 1")
            group = gemmi.find_spacegroup_by_name(group_text)
            if group is not None:
                self.manual_cif_spacegroup_number = int(group.number)
                self.manual_cif_spacegroup_symbol = str(group.hm)
                combo_index = self.manual_spacegroup.findData(int(group.number))
                if combo_index >= 0:
                    self.manual_spacegroup.setCurrentIndex(combo_index)
                else:
                    self.manual_spacegroup.setCurrentText(group_text)
            else:
                self.manual_spacegroup.setCurrentText(group_text)
            self._update_spacegroup_consistency()
            self.manual_status.setText(f"Loaded CIF: {_GuiPath(path).name}")
        except Exception as exc:
            if show_errors:
                QMessageBox.critical(self, "CIF error", str(exc))

    @staticmethod
    def _normalize_manual_image_array(image):
        """Normalize an image for reliable display and crop detection."""
        array = np.asarray(image)
        if array.ndim not in (2, 3):
            raise ValueError(f"Unsupported image shape: {array.shape}")
        if array.ndim == 3 and array.shape[2] == 4:
            rgb = array[..., :3].astype(np.float32, copy=False)
            alpha = array[..., 3].astype(np.float32, copy=False)
            if np.issubdtype(array.dtype, np.integer):
                scale = float(np.iinfo(array.dtype).max)
                rgb = rgb / scale
                alpha = alpha / scale
            elif np.nanmax(rgb) > 1.5 or np.nanmax(alpha) > 1.5:
                rgb = rgb / 255.0
                alpha = alpha / 255.0
            # Composite transparency onto white, avoiding dark/empty borders.
            array = rgb * alpha[..., None] + (1.0 - alpha[..., None])
        elif np.issubdtype(array.dtype, np.integer):
            array = array.astype(np.float32) / float(np.iinfo(array.dtype).max)
        else:
            array = array.astype(np.float32, copy=False)
            finite = np.isfinite(array)
            if finite.any() and float(np.nanmax(array[finite])) > 1.5:
                array = array / 255.0
        return np.clip(array, 0.0, 1.0)

    @staticmethod
    def _manual_detect_panel_crop(image):
        """Find the colored q-space panel while safely retaining raw images."""
        array = np.asarray(image)
        height, width = array.shape[:2]
        if height < 40 or width < 40 or array.ndim != 3 or array.shape[2] < 3:
            return (0, 0, width, height), "full image"

        rgb = np.asarray(array[..., :3], dtype=np.float32)
        try:
            x0, y0, x1, y1 = detect_plot_crop(rgb)
        except Exception:
            return (0, 0, width, height), "full image"

        crop_w, crop_h = x1 - x0, y1 - y0
        area_fraction = (crop_w * crop_h) / max(width * height, 1)
        # Reject tiny colored objects, legends, isolated color bars, or near-empty crops.
        if crop_w < max(80, int(0.25 * width)) or crop_h < max(80, int(0.25 * height)):
            return (0, 0, width, height), "full image"
        if not (0.12 <= area_fraction <= 0.995):
            return (0, 0, width, height), "full image"

        panel = rgb[y0:y1, x0:x1]
        saturation = panel.max(axis=2) - panel.min(axis=2)
        if float(np.mean(saturation > 0.08)) < 0.08:
            return (0, 0, width, height), "full image"
        return (int(x0), int(y0), int(x1), int(y1)), "auto-cropped q-space panel"

    def _apply_manual_image_framing(self):
        if self.manual_image_original is None:
            self.manual_image = None
            self.manual_image_crop_xyxy = None
            if hasattr(self, "manual_frame_status"):
                self.manual_frame_status.setText("No image loaded.")
            return

        original = self.manual_image_original
        height, width = original.shape[:2]
        if self.manual_auto_crop.isChecked():
            crop, description = self._manual_detect_panel_crop(original)
        else:
            crop, description = (0, 0, width, height), "full image (automatic crop disabled)"
        x0, y0, x1, y1 = crop
        self.manual_image_crop_xyxy = crop
        self.manual_image = np.ascontiguousarray(original[y0:y1, x0:x1])
        framed_h, framed_w = self.manual_image.shape[:2]
        percent = 100.0 * framed_w * framed_h / max(width * height, 1)
        self.manual_frame_status.setText(
            f"{description}: {framed_w} × {framed_h} px "
            f"({percent:.1f}% of the loaded image)."
        )

    def _reframe_manual_image(self):
        self._apply_manual_image_framing()
        self._redraw_manual()

    def _choose_manual_image(self):
        path, _ = QFileDialog.getOpenFileName(self, "Select intensity image", "",
                                               "Images (*.png *.jpg *.jpeg *.tif *.tiff);;All files (*)")
        if not path:
            return
        try:
            loaded = _mpimg.imread(path)
            self.manual_image_original = self._normalize_manual_image_array(loaded)
            self.manual_image_path = path
            self.manual_image_edit.setText(path)
            self._apply_manual_image_framing()
            self._redraw_manual()
        except Exception as exc:
            QMessageBox.critical(self, "Image error", str(exc))

    def _current_auto_giwaxs_physics_config(self):
        """Best-effort AUTO physics config for the cyan manual reflection pattern."""
        try:
            request = self._build_request(preview=True)
            base = engine_v9.V9Config(
                cif_path=request.cif_path or "manual.cif",
                output_dir=request.output_dir or ".",
                q_max=float(request.q_max),
            )
            return _apply_auto_giwaxs_physics(base, request)
        except Exception:
            return None

    def _calculate_manual(self):
        try:
            lattice = (
                self.manual_a.value(), self.manual_b.value(), self.manual_c.value(),
                self.manual_alpha.value(), self.manual_beta.value(), self.manual_gamma.value(),
            )
            orientation = (self.manual_h.value(), self.manual_k.value(), self.manual_l.value())
            manual_physics = self._current_auto_giwaxs_physics_config()
            self.manual_result = calculate_manual_reflections(
                self.manual_structure, lattice, orientation, self.manual_spacegroup.currentText(),
                self.manual_hkl_max.value(), self.manual_q_max.value(),
                physics_config=manual_physics,
            )
            if self.manual_result.empty:
                raise ValueError("No allowed reflections were found inside q max.")
            self._redraw_manual()
            selected_group = self.manual_spacegroup.currentText().strip()
            selected_number, _selected_symbol = self._selected_manual_spacegroup()
            caveat = ""
            if selected_number == 1:
                caveat = " P1 is permissive, so use residuals and validation rather than visual density alone."
            elif self.manual_cif_spacegroup_number is not None and selected_number != self.manual_cif_spacegroup_number:
                caveat = " Manual space-group override differs from the CIF; treat this pattern as exploratory."
            physics_note = ""
            if manual_physics is not None:
                physics_note = f" GIWAXS physics: {manual_physics.giwaxs_physics_status}."
            self.manual_status.setText(
                f"Generated {len(self.manual_result)} unique reciprocal-space reflection positions "
                f"for {selected_group}.{caveat}{physics_note}"
            )
            self._update_spacegroup_consistency()
        except Exception as exc:
            QMessageBox.critical(self, "Manual calculation error", str(exc))
            self.manual_status.setText(f"Error: {exc}")

    def _redraw_manual(self):
        if not hasattr(self, "manual_axes"):
            return
        self.manual_axes.clear()
        qr_min, qr_max = self.manual_qr_min.value(), self.manual_qr_max.value()
        qz_min, qz_max = self.manual_qz_min.value(), self.manual_qz_max.value()
        if qr_min >= qr_max or qz_min >= qz_max:
            self.manual_canvas.draw_idle()
            return
        q_span = qr_max - qr_min
        z_span = qz_max - qz_min
        if self.manual_fit_frame.isChecked():
            # Match the physical q-space aspect ratio to the axes box. This makes
            # the image touch all four sides without stretching q coordinates.
            self.manual_axes.set_box_aspect(z_span / q_span)
            image_aspect = "auto"
        else:
            self.manual_axes.set_box_aspect(None)
            image_aspect = "equal"
        self.manual_axes.set_anchor("C")
        self.manual_axes.margins(x=0.0, y=0.0)

        if self.manual_image is not None:
            kwargs = dict(extent=[qr_min, qr_max, qz_min, qz_max], origin="upper",
                          aspect=image_aspect, alpha=self.manual_opacity.value(), zorder=0,
                          interpolation="nearest", resample=False)
            if self.manual_image.ndim == 2:
                kwargs["cmap"] = "turbo"
            self.manual_axes.imshow(self.manual_image, **kwargs)
        if self.manual_result is not None and not self.manual_result.empty:
            qr = self.manual_result.Qxy.to_numpy(float)
            qz = self.manual_result.Qz.to_numpy(float)
            labels = self.manual_result.hkl.astype(str).to_numpy()
            intensity = self.manual_result.Intensity.to_numpy(float)
            if self.manual_mirror.isChecked():
                nonzero = qr > 1e-10
                qr = np.concatenate([qr, -qr[nonzero]])
                qz = np.concatenate([qz, qz[nonzero]])
                labels = np.concatenate([labels, labels[nonzero]])
                intensity = np.concatenate([intensity, intensity[nonzero]])
            visible = (qr >= qr_min) & (qr <= qr_max) & (qz >= qz_min) & (qz <= qz_max)
            qr, qz, labels, intensity = qr[visible], qz[visible], labels[visible], intensity[visible]
            if len(intensity):
                maximum = max(float(np.nanmax(intensity)), 1e-12)
                sizes = 18.0 + 135.0 * np.sqrt(np.clip(intensity / maximum, 0.0, 1.0))
                edge = "white" if self.manual_image is not None else "black"
                self.manual_axes.scatter(qr, qz, s=sizes, marker="o", facecolors="none",
                                         edgecolors=edge, linewidths=1.4, zorder=3,
                                         label="Simulated reflections")
                if self.manual_labels.isChecked():
                    for x, y, text in zip(qr, qz, labels):
                        self.manual_axes.annotate(text, (x, y), xytext=(4, 4),
                                                  textcoords="offset points", fontsize=6.5,
                                                  color=edge, annotation_clip=True)
        self.manual_axes.set_xlim(qr_min, qr_max)
        self.manual_axes.set_ylim(qz_min, qz_max)
        self.manual_axes.set_xlabel(r"$q_r$ ($\AA^{-1}$)")
        self.manual_axes.set_ylabel(r"$q_z$ ($\AA^{-1}$)")
        self.manual_axes.set_title("Manual experimental–simulation overlay")
        self.manual_axes.grid(alpha=0.18)
        # Tight framing leaves room for labels while keeping the image flush with
        # the four data boundaries. The axes resize cleanly with the GUI window.
        self.manual_figure.tight_layout(pad=1.1)
        self.manual_canvas.draw_idle()

    # ------------------------------- preview tab --------------------------------

    # ------------------------------- progress tab -------------------------------



    # -------------------------------- results tab -------------------------------






# ================= MANUAL-STYLE GUI WITH FULL INDEXING BACKEND =================

class ManualFullIndexingWorkbench(IntegratedGIXSWorkbench):
    """Manual CIF overlay interface backed by the full automatic indexing engine.

    Indexed/ignored classifications, HKL labels, calculated coordinates,
    calibration, multidomain analysis, and validation are derived from the
    automatic backend. The GUI therefore displays backend scientific results
    rather than assigning reflections from visual proximity alone.
    """

    def __init__(self):
        # The base constructor builds the manual tab and calls self._redraw_manual().
        # Because Python dispatches that call to this subclass override, the
        # subclass-only checkboxes do not exist yet. Keep a construction guard so
        # the base redraw is used until the full-indexing controls are ready.
        self._manual_full_ui_ready = False
        self.backend_indexed = pd.DataFrame()
        self.backend_ignored = pd.DataFrame()
        self.backend_table = pd.DataFrame()
        self.backend_series_item = None
        super().__init__()
        self.setWindowTitle("GIWAXS/GIXS Manual CIF Overlay — Full Indexing Backend")
        self.resize(1510, 980)
        self._convert_to_manual_only_workspace()
        self._manual_full_ui_ready = True
        self._redraw_manual()

    def _convert_to_manual_only_workspace(self):
        # Keep the interface focused on the combined manual-overlay and automatic
        # NPZ-calculator workspace.
        for index in range(self.main_tabs.count()):
            self.main_tabs.setTabVisible(index, self.main_tabs.widget(index) is self.manual_tab)
        self.main_tabs.setCurrentWidget(self.manual_tab)
        self.main_tabs.setTabText(self.main_tabs.indexOf(self.manual_tab), "Manual CIF indexing overlay")

        root_layout = self.manual_tab.layout()
        horizontal = root_layout.itemAt(0).widget()
        root_layout.removeWidget(horizontal)
        vertical = QSplitter(Qt.Orientation.Vertical)
        vertical.addWidget(horizontal)

        table_panel = QGroupBox("Indexed HKL, experimental/calculated q values, and correlated intensities")
        table_layout = QVBoxLayout(table_panel)
        self.full_correlation_label = QLabel(
            "Run full indexing to populate the experimental/calculated coordinate and intensity table."
        )
        self.full_correlation_label.setWordWrap(True)
        table_layout.addWidget(self.full_correlation_label)

        # Simple view presents the primary experimental/calculated comparison
        # quantities for quick interpretation. Advanced view exposes the complete
        # diagnostic table. Switching views changes presentation only and never
        # modifies indexing or assignment results.
        table_view_row = QHBoxLayout()
        table_view_label = QLabel("Table view:")
        self.table_view_combo = QComboBox()
        self.table_view_combo.addItems(["Simple view", "Advanced view"])
        self.table_view_combo.setCurrentText("Simple view")
        self.table_view_combo.setToolTip(
            "Simple view shows the main peak-identification, q-position, d-spacing, "
            "intensity, and match-quality results. Advanced view shows every available field."
        )
        table_view_row.addWidget(table_view_label)
        table_view_row.addWidget(self.table_view_combo)
        table_view_row.addStretch(1)
        table_layout.addLayout(table_view_row)
        self.table_view_combo.currentTextChanged.connect(
            lambda _=None: self._refresh_visible_data_table()
        )

        self.full_index_table = QTableWidget()
        self.full_index_table.setAlternatingRowColors(True)
        self.full_index_table.setSortingEnabled(True)
        self.full_index_table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.full_index_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        table_layout.addWidget(self.full_index_table, 1)
        bottom_buttons = QHBoxLayout()
        export_button = QPushButton("Export indexed correlation table…")
        export_button.clicked.connect(self._export_full_index_table)
        bottom_buttons.addWidget(export_button)
        bottom_buttons.addStretch(1)
        table_layout.addLayout(bottom_buttons)
        vertical.addWidget(table_panel)
        vertical.setSizes([665, 285])
        root_layout.addWidget(vertical)

        # Add the full-engine controls into the existing left scroll panel.
        scroll = horizontal.widget(0)
        controls_layout = scroll.widget().layout()
        engine_group = QGroupBox("Full automatic indexing calculations")
        engine_layout = QVBoxLayout(engine_group)

        note = QLabel(
            "The interface stays manual-style, but this section runs the complete indexing engine: "
            "multiscale feature detection, multi-angle registration/consensus, orientation search, "
            "calibration refinement, ignored/weak-peak recovery, multidomain testing, and validation."
        )
        note.setWordWrap(True)
        engine_layout.addWidget(note)

        form = QFormLayout()
        self.manual_workflow_combo = QComboBox()
        self.manual_workflow_combo.addItems([
            "Fast test indexing",
            "Standard indexing",
            "Improved coverage indexing",
        ])
        self.manual_workflow_combo.setCurrentText("Improved coverage indexing")
        self.manual_output_edit = QLineEdit(str(_GuiPath.cwd() / "gixs_manual_full_indexing_results"))
        output_holder = QWidget()
        output_holder_layout = QHBoxLayout(output_holder)
        output_holder_layout.setContentsMargins(0, 0, 0, 0)
        output_holder_layout.addWidget(self.manual_output_edit, 1)
        output_button = QPushButton("Browse…")
        output_button.clicked.connect(self._choose_manual_output)
        output_holder_layout.addWidget(output_button)
        form.addRow("Workflow:", self.manual_workflow_combo)
        form.addRow("Output folder:", output_holder)
        engine_layout.addLayout(form)

        measurement_buttons = QHBoxLayout()
        add_images = QPushButton("Add NPZ calculation files…")
        add_images.clicked.connect(self._add_npz_calculation_files)
        remove_images = QPushButton("Remove selected")
        remove_images.clicked.connect(self._remove_npz_calculation_files)
        measurement_buttons.addWidget(add_images)
        measurement_buttons.addWidget(remove_images)
        engine_layout.addLayout(measurement_buttons)

        # Display the editable NPZ measurement manifest in a compact form so scan
        # metadata can be reviewed without reducing the reciprocal-space plot area.
        self.measurement_table.setParent(engine_group)
        self.measurement_table.setMaximumHeight(165)
        for column in (4, 7, 8, 9, 10, 11):
            self.measurement_table.setColumnHidden(column, True)
        self.measurement_table.setColumnWidth(1, 155)
        engine_layout.addWidget(self.measurement_table)
        hint = QLabel("Edit Sample, Series, Angle, and Scan directly in the table. q ranges come from the image controls above.")
        hint.setWordWrap(True)
        engine_layout.addWidget(hint)

        display_row = QHBoxLayout()
        self.show_full_simulation = QCheckBox("Show all simulated CIF points")
        self.show_full_simulation.setChecked(False)
        self.show_indexed_points = QCheckBox("Show indexed points")
        self.show_indexed_points.setChecked(True)
        self.show_not_indexed_points = QCheckBox("Show not-indexed points")
        self.show_not_indexed_points.setChecked(True)
        display_row.addWidget(self.show_full_simulation)
        display_row.addWidget(self.show_indexed_points)
        display_row.addWidget(self.show_not_indexed_points)
        engine_layout.addLayout(display_row)
        for box in (self.show_full_simulation, self.show_indexed_points, self.show_not_indexed_points):
            box.stateChanged.connect(lambda _=None: self._redraw_manual())

        run_row = QHBoxLayout()
        run_button = QPushButton("Run full indexing onto this overlay")
        run_button.setMinimumHeight(34)
        run_button.clicked.connect(self._start_full_indexing)
        cancel_button = QPushButton("Cancel")
        cancel_button.clicked.connect(self._cancel_full_indexing)
        cancel_button.setEnabled(False)
        self.manual_cancel_button = cancel_button
        run_row.addWidget(run_button, 1)
        run_row.addWidget(cancel_button)
        engine_layout.addLayout(run_row)

        self.manual_progress_label = QLabel("Ready")
        self.manual_progress_bar = QProgressBar()
        self.manual_progress_bar.setRange(0, 100)
        engine_layout.addWidget(self.manual_progress_label)
        engine_layout.addWidget(self.manual_progress_bar)
        self.manual_progress_log = QPlainTextEdit()
        self.manual_progress_log.setReadOnly(True)
        self.manual_progress_log.setMaximumHeight(115)
        engine_layout.addWidget(self.manual_progress_log)

        series_row = QHBoxLayout()
        series_row.addWidget(QLabel("Displayed series:"))
        self.manual_series_combo = QComboBox()
        self.manual_series_combo.currentTextChanged.connect(self._load_full_series)
        series_row.addWidget(self.manual_series_combo, 1)
        self.manual_decision_label = QLabel("")
        self.manual_decision_label.setWordWrap(True)
        engine_layout.addLayout(series_row)
        engine_layout.addWidget(self.manual_decision_label)

        # Place automatic-calculator controls before manual reflection-generation
        # actions so the two analysis paths remain visually separated.
        controls_layout.insertWidget(max(0, controls_layout.count() - 4), engine_group)

        # The manual reflection generator is a visual crystallographic simulation
        # based on selected structure parameters and preferred orientation. Its
        # explanatory text distinguishes it from automatic NPZ indexing.
        for child in scroll.widget().findChildren(QLabel):
            if "Manual mode is exploratory" in child.text():
                child.setText(
                    "Generate selected space-group reflection pattern uses the manually chosen "
                    "space group, lattice parameters, and preferred orientation. Run full automatic "
                    "indexing separately to generate the NPZ-calculator reflection model."
                )

    def _choose_manual_output(self):
        path = QFileDialog.getExistingDirectory(self, "Select output folder", self.manual_output_edit.text())
        if path:
            self.manual_output_edit.setText(path)

    def _load_manual_cif(self, path, show_errors=True):
        super()._load_manual_cif(path, show_errors=show_errors)
        if path:
            self.cif_edit.setText(str(path))



    def _workflow_key_manual(self):
        return {
            "Fast test indexing": "fast_test",
            "Standard indexing": "recommended",
            "Improved coverage indexing": "improved_coverage",
        }[self.manual_workflow_combo.currentText()]


    def _start_full_indexing(self):
        if self.controller.running:
            QMessageBox.information(self, "Job running", "An indexing job is already running.")
            return
        try:
            request = self._build_request(preview=False)
            errors = validate_request(request)
            if errors:
                raise ValueError("\n".join(errors))
            self.controller.start(request)
        except Exception as exc:
            QMessageBox.critical(self, "Cannot start", str(exc))
            return
        self.backend_indexed = pd.DataFrame()
        self.backend_ignored = pd.DataFrame()
        self.backend_table = pd.DataFrame()
        self.manual_progress_bar.setValue(0)
        self.manual_progress_label.setText("Starting full indexing")
        self.manual_progress_log.clear()
        self.manual_cancel_button.setEnabled(True)
        self.manual_status.setText("Full indexing is running. The GUI remains responsive.")
        if getattr(self, "giwaxs_physics_status_label", None) is not None:
            self.giwaxs_physics_status_label.setText(
                "AUTO — evaluating measurement metadata; safe kinematic fallback remains available."
            )

    def _cancel_full_indexing(self):
        self.controller.cancel()
        self.manual_cancel_button.setEnabled(False)
        self.manual_progress_label.setText("Cancellation requested")

    def _poll_controller(self):
        for event in self.controller.poll():
            kind = event.get("kind")
            if kind == "progress":
                stage = str(event.get("stage", ""))
                fraction = float(event.get("fraction", 0.0))
                self.manual_progress_label.setText(stage)
                self.manual_progress_bar.setValue(int(round(100 * fraction)))
            elif kind == "log":
                self.manual_progress_log.appendPlainText(str(event.get("message", "")))
            elif kind == "result":
                self.current_result = GIXSRunResult.from_json_dict(event["result"])
                self.manual_cancel_button.setEnabled(False)
                self._populate_full_manual_result()
            elif kind == "cancelled":
                self.manual_progress_log.appendPlainText("Cancellation requested by user.")
            elif kind == "exit":
                self.manual_cancel_button.setEnabled(False)

    def _populate_full_manual_result(self):
        result = self.current_result
        if result is None:
            return
        if result.status != "completed":
            QMessageBox.critical(self, "Run failed", result.error or "The run did not complete.")
            self.manual_status.setText(f"Full indexing failed: {result.error}")
            return
        self.manual_series_combo.blockSignals(True)
        self.manual_series_combo.clear()
        self.manual_series_combo.addItems([item.series_id for item in result.series_results])
        self.manual_series_combo.blockSignals(False)
        if self.manual_series_combo.count():
            self.manual_series_combo.setCurrentIndex(0)
            self._load_full_series(self.manual_series_combo.currentText())
        self.manual_progress_label.setText("Completed")
        self.manual_progress_bar.setValue(100)
        self._update_giwaxs_physics_status_from_output()

    def _update_giwaxs_physics_status_from_output(self):
        label = getattr(self, "giwaxs_physics_status_label", None)
        result = getattr(self, "current_result", None)
        if label is None or result is None:
            return
        metadata_path = _GuiPath(result.output_dir) / "dwba_model_metadata.json"
        try:
            payload = json.loads(metadata_path.read_text())
        except Exception:
            label.setText("AUTO — physics metadata unavailable; indexing result remains valid.")
            return
        status = str(payload.get("status") or "AUTO")
        reason = str(payload.get("fallback_reason") or "").strip()
        angle = payload.get("incidence_angle_deg")
        wavelength = payload.get("xray_wavelength_A")
        details = []
        try:
            if np.isfinite(float(angle)):
                details.append(f"αi≈{float(angle):.4g}°")
        except Exception:
            pass
        try:
            physics_active = bool(payload.get("enabled") or payload.get("position_correction_enabled"))
            if physics_active and np.isfinite(float(wavelength)) and float(wavelength) > 0:
                details.append(f"λ={float(wavelength):.5g} Å")
        except Exception:
            pass
        text = status
        if details:
            text += " | " + ", ".join(details)
        if reason:
            text += " | fallback notes: " + reason
        label.setText(text)

    @staticmethod
    def _pick_column(frame, names):
        for name in names:
            if name in frame.columns:
                return name
        return None

    def _load_full_series(self, series_id):
        if not self.current_result or not series_id:
            return
        item = next((entry for entry in self.current_result.series_results if entry.series_id == series_id), None)
        if item is None:
            return
        self.backend_series_item = item
        indexed_path = _GuiPath(item.indexed_table_path)
        series_dir = _GuiPath(item.overlay_path).parent
        ignored_path = series_dir / "overlay_ignored_features.csv"
        coverage_path = series_dir / "consensus_feature_classification_exhaustive.csv"
        try:
            self.backend_indexed = pd.read_csv(indexed_path) if indexed_path.is_file() else pd.DataFrame()
        except Exception:
            self.backend_indexed = pd.DataFrame()
        try:
            self.backend_ignored = pd.read_csv(ignored_path) if ignored_path.is_file() else pd.DataFrame()
        except Exception:
            self.backend_ignored = pd.DataFrame()
        if self.backend_ignored.empty and coverage_path.is_file():
            try:
                coverage = pd.read_csv(coverage_path)
                class_col = self._pick_column(coverage, ["overlay_class", "classification", "status"])
                if class_col:
                    self.backend_ignored = coverage[
                        coverage[class_col].astype(str).str.lower().isin(["ignored", "not_indexed", "not indexed"])
                    ].copy()
            except Exception:
                pass
        self.backend_table = self._build_backend_correlation_table(self.backend_indexed)
        self._fill_backend_table_widget(self.backend_table)
        self._update_backend_correlation_label()
        self.manual_decision_label.setText(
            f"Decision: {item.final_decision or 'not reported'} | orientation: {item.orientation_hkl or 'not reported'} | "
            f"indexed: {len(self.backend_indexed)} | not indexed: {len(self.backend_ignored)}"
        )
        self.manual_status.setText(
            f"Full indexing loaded for {series_id}: {len(self.backend_indexed)} indexed and "
            f"{len(self.backend_ignored)} not-indexed experimental points."
        )
        self._redraw_manual()

    def _build_backend_correlation_table(self, indexed):
        if indexed is None or indexed.empty:
            return pd.DataFrame()
        frame = indexed.copy().reset_index(drop=True)
        qr_exp = self._pick_column(frame, ["qr_exp", "qr", "experimental_qr"])
        qz_exp = self._pick_column(frame, ["qz_exp", "qz", "experimental_qz"])
        qr_calc = self._pick_column(frame, ["qr_calc", "predicted_qr", "qr_pred"])
        qz_calc = self._pick_column(frame, ["qz_calc", "predicted_qz", "qz_pred"])
        exp_i = self._pick_column(frame, ["experimental_integrated_intensity", "integrated_intensity", "strength", "experimental_intensity", "peak_intensity", "local_strength"])
        calc_i = self._pick_column(
            frame, ["effective_intensity", "calculated_intensity", "f2", "prediction_weight"]
        )

        def values(column, default=np.nan):
            if column is None:
                return np.full(len(frame), default, dtype=float)
            return pd.to_numeric(frame[column], errors="coerce").to_numpy(float)

        exp_qr = values(qr_exp)
        exp_qz = values(qz_exp)
        calc_qr = values(qr_calc)
        calc_qz = values(qz_calc)
        experimental_intensity = values(exp_i)
        calculated_intensity = values(calc_i)

        def normalize(array):
            result = np.full(len(array), np.nan, dtype=float)
            finite = np.isfinite(array)
            if finite.any():
                maximum = float(np.nanmax(np.maximum(array[finite], 0.0)))
                if maximum > 0:
                    result[finite] = np.maximum(array[finite], 0.0) / maximum
                else:
                    result[finite] = 0.0
            return result

        exp_norm = normalize(experimental_intensity)
        calc_norm = normalize(calculated_intensity)
        agreement = 1.0 - np.abs(exp_norm - calc_norm)
        hkl_col = self._pick_column(frame, ["hkl", "overlay_hkl_text"])
        if hkl_col:
            hkl = frame[hkl_col].astype(str).to_numpy()
        elif all(name in frame.columns for name in ("h", "k", "l")):
            hkl = np.array([f"({int(h)} {int(k)} {int(l)})" for h, k, l in frame[["h", "k", "l"]].itertuples(index=False, name=None)])
        else:
            hkl = np.array(["" for _ in range(len(frame))])
        delta_qr = exp_qr - calc_qr
        delta_qz = exp_qz - calc_qz
        delta_q = np.hypot(delta_qr, delta_qz)
        q_total_exp = np.hypot(exp_qr, exp_qz)
        q_total_calc = np.hypot(calc_qr, calc_qz)
        d_exp = np.divide(2.0 * np.pi, q_total_exp, out=np.full(len(frame), np.nan), where=q_total_exp > 0)
        d_calc = np.divide(2.0 * np.pi, q_total_calc, out=np.full(len(frame), np.nan), where=q_total_calc > 0)
        table = pd.DataFrame({
            "ID": np.arange(1, len(frame) + 1),
            "HKL": hkl,
            "q_r experimental": exp_qr,
            "q_z experimental": exp_qz,
            "q_r calculated": calc_qr,
            "q_z calculated": calc_qz,
            "d experimental": d_exp,
            "d calculated": d_calc,
            "delta q_r": delta_qr,
            "delta q_z": delta_qz,
            "delta q": delta_q,
            "experimental intensity": experimental_intensity,
            "calculated intensity": calculated_intensity,
            "normalized experimental intensity": exp_norm,
            "normalized calculated intensity": calc_norm,
            "intensity agreement": agreement,
        })
        for source, target in [
            ("experimental_intensity_sigma", "experimental intensity sigma"),
            ("experimental_integrated_snr", "experimental integrated SNR"),
            ("experimental_intensity_quality_score", "experimental intensity quality"),
        ]:
            if source in frame.columns:
                table[target] = pd.to_numeric(frame[source], errors="coerce").to_numpy(float)
        for source, target in [
            ("orientation_domain", "orientation domain"),
            ("index_source", "index source"),
            ("salvage_evidence_tier", "evidence tier"),
            ("assignment_support_score", "assignment support"),
            ("support", "angle support"),
            ("feature_id", "feature ID"),
        ]:
            if source in frame.columns:
                table[target] = frame[source].to_numpy()
        return table

    def _table_view_is_simple(self):
        combo = getattr(self, "table_view_combo", None)
        return combo is None or combo.currentText() == "Simple view"

    @staticmethod
    def _match_quality_label(delta_q):
        try:
            value = float(delta_q)
        except (TypeError, ValueError):
            return ""
        if not np.isfinite(value):
            return ""
        if value <= 0.020:
            return "Very Good"
        if value <= 0.040:
            return "Good"
        if value <= 0.060:
            return "Fair"
        return "Poor"

    def _refresh_visible_data_table(self):
        # When manual click assignments are available, display those experimental
        # versus calculated pairs. Otherwise show the automatic correlation table.
        manual_filler = getattr(self, "_fill_manual_assignment_table", None)
        if hasattr(self, "manual_click_assignments") and callable(manual_filler):
            manual_filler()
        else:
            self._fill_backend_table_widget(getattr(self, "backend_table", pd.DataFrame()))

    def _simple_backend_table(self, frame):
        if frame is None or frame.empty:
            return pd.DataFrame()
        intensity_match = 100.0 * pd.to_numeric(frame.get("intensity agreement"), errors="coerce")
        delta_q = pd.to_numeric(frame.get("delta q"), errors="coerce")
        return pd.DataFrame({
            "Peak": frame.get("ID"),
            "Crystal Plane (HKL)": frame.get("HKL"),
            "Measured qᵣ (Å⁻¹)": frame.get("q_r experimental"),
            "Measured qz (Å⁻¹)": frame.get("q_z experimental"),
            "Predicted qᵣ (Å⁻¹)": frame.get("q_r calculated"),
            "Predicted qz (Å⁻¹)": frame.get("q_z calculated"),
            "Measured d-spacing (Å)": frame.get("d experimental"),
            "Predicted d-spacing (Å)": frame.get("d calculated"),
            "Match Error Δq (Å⁻¹)": delta_q,
            "Measured Intensity": frame.get("experimental intensity"),
            "Predicted Intensity": frame.get("calculated intensity"),
            "Intensity Match (%)": intensity_match,
            "Match Quality": [self._match_quality_label(value) for value in delta_q],
        })

    def _apply_simple_table_header_tooltips(self):
        tooltips = {
            "Peak": "Row number for the indexed peak.",
            "Crystal Plane (HKL)": "Miller indices assigned to this diffraction peak.",
            "Measured qᵣ (Å⁻¹)": "Experimental in-plane reciprocal-space position.",
            "Measured qz (Å⁻¹)": "Experimental out-of-plane reciprocal-space position.",
            "Predicted qᵣ (Å⁻¹)": "Calculated in-plane reciprocal-space position.",
            "Predicted qz (Å⁻¹)": "Calculated out-of-plane reciprocal-space position.",
            "Measured d-spacing (Å)": "Experimental lattice spacing calculated as 2π/|q|.",
            "Predicted d-spacing (Å)": "Calculated lattice spacing calculated as 2π/|q|.",
            "Match Error Δq (Å⁻¹)": "Overall distance between measured and predicted q positions; smaller is better.",
            "Measured Intensity": "Experimental peak intensity.",
            "Predicted Intensity": "Calculated reflection intensity used for comparison.",
            "Intensity Match (%)": "Relative measured/calculated intensity agreement; closer to 100% is better.",
            "Match Quality": "Plain-language position match: Very Good ≤0.02, Good ≤0.04, Fair ≤0.06, Poor >0.06 Å⁻¹.",
        }
        header = self.full_index_table.horizontalHeader()
        for column in range(self.full_index_table.columnCount()):
            item = self.full_index_table.horizontalHeaderItem(column)
            if item is not None:
                item.setToolTip(tooltips.get(item.text(), ""))

    def _fill_backend_table_widget(self, frame):
        self.full_index_table.setSortingEnabled(False)
        self.full_index_table.clear()
        display_frame = self._simple_backend_table(frame) if self._table_view_is_simple() else frame
        if display_frame is None or display_frame.empty:
            self.full_index_table.setRowCount(0)
            self.full_index_table.setColumnCount(1)
            self.full_index_table.setHorizontalHeaderLabels(["Status"])
            self.full_index_table.setRowCount(1)
            self.full_index_table.setItem(0, 0, QTableWidgetItem("No indexed points available."))
            self.full_index_table.setSortingEnabled(True)
            return
        self.full_index_table.setColumnCount(len(display_frame.columns))
        self.full_index_table.setHorizontalHeaderLabels([str(column) for column in display_frame.columns])
        self.full_index_table.setRowCount(len(display_frame))
        for row_index, row in enumerate(display_frame.itertuples(index=False, name=None)):
            for column_index, value in enumerate(row):
                if pd.isna(value):
                    text = ""
                elif isinstance(value, (float, np.floating)):
                    text = f"{float(value):.6g}"
                else:
                    text = str(value)
                self.full_index_table.setItem(row_index, column_index, QTableWidgetItem(text))
        self.full_index_table.resizeColumnsToContents()
        if self._table_view_is_simple():
            self._apply_simple_table_header_tooltips()
        self.full_index_table.setSortingEnabled(True)

    def _update_backend_correlation_label(self):
        table = self.backend_table
        if table is None or table.empty:
            self.full_correlation_label.setText("No indexed points are available for correlation.")
            return
        x = pd.to_numeric(table["experimental intensity"], errors="coerce")
        y = pd.to_numeric(table["calculated intensity"], errors="coerce")
        valid = np.isfinite(x) & np.isfinite(y)
        pearson = np.nan
        spearman = np.nan
        if int(valid.sum()) >= 3 and float(np.nanstd(x[valid])) > 0 and float(np.nanstd(y[valid])) > 0:
            pearson = float(pd.Series(x[valid]).corr(pd.Series(y[valid]), method="pearson"))
            spearman = float(pd.Series(x[valid]).corr(pd.Series(y[valid]), method="spearman"))
        p_text = "not available" if not np.isfinite(pearson) else f"{pearson:.4f}"
        s_text = "not available" if not np.isfinite(spearman) else f"{spearman:.4f}"
        median_delta = float(np.nanmedian(pd.to_numeric(table["delta q"], errors="coerce")))
        self.full_correlation_label.setText(
            f"Indexed points: {len(table)} | median coordinate residual: {median_delta:.5f} Å⁻¹ | "
            f"Pearson intensity correlation: {p_text} | Spearman intensity correlation: {s_text}. "
            "Calculated intensity uses the backend effective reflection intensity (DWBA-weighted when AUTO physics is available); "
            "experimental intensity uses robust integrated numerical intensity when propagated, otherwise the legacy detected strength."
        )

    @staticmethod
    def _xy_columns(frame, experimental=True):
        if frame is None or frame.empty:
            return None, None
        if experimental:
            qr_names = ["qr_exp", "qr", "experimental_qr"]
            qz_names = ["qz_exp", "qz", "experimental_qz"]
        else:
            qr_names = ["qr_calc", "predicted_qr", "qr_pred"]
            qz_names = ["qz_calc", "predicted_qz", "qz_pred"]
        qr = next((name for name in qr_names if name in frame.columns), None)
        qz = next((name for name in qz_names if name in frame.columns), None)
        return qr, qz

    def _redraw_manual(self):
        # During initial GUI construction, subclass display controls may not yet
        # exist. Use the base redraw path until the combined indexing workspace is
        # fully initialized, then enable the richer overlay controls.
        if not getattr(self, "_manual_full_ui_ready", False):
            return IntegratedGIXSWorkbench._redraw_manual(self)

        # When automatic results are available, hide the dense full-CIF simulation
        # by default so measured/calculated assignments remain readable. The complete
        # calculated reflection pattern can still be shown on request.
        saved_result = self.manual_result
        show_full_simulation = bool(self.show_full_simulation.isChecked())
        if not show_full_simulation and not self.backend_indexed.empty:
            self.manual_result = pd.DataFrame()
        IntegratedGIXSWorkbench._redraw_manual(self)
        self.manual_result = saved_result

        if not hasattr(self, "manual_axes"):
            return
        if self.show_not_indexed_points.isChecked() and not self.backend_ignored.empty:
            qr_col, qz_col = self._xy_columns(self.backend_ignored, experimental=True)
            if qr_col and qz_col:
                qr = pd.to_numeric(self.backend_ignored[qr_col], errors="coerce").to_numpy(float)
                qz = pd.to_numeric(self.backend_ignored[qz_col], errors="coerce").to_numpy(float)
                valid = np.isfinite(qr) & np.isfinite(qz)
                self.manual_axes.scatter(
                    qr[valid], qz[valid], s=28, marker="x", c="0.72", alpha=0.78,
                    linewidths=0.9, label=f"not indexed ({int(valid.sum())})", zorder=5,
                )

        if self.show_indexed_points.isChecked() and not self.backend_indexed.empty:
            exp_qr_col, exp_qz_col = self._xy_columns(self.backend_indexed, experimental=True)
            calc_qr_col, calc_qz_col = self._xy_columns(self.backend_indexed, experimental=False)
            if exp_qr_col and exp_qz_col:
                exp_qr = pd.to_numeric(self.backend_indexed[exp_qr_col], errors="coerce").to_numpy(float)
                exp_qz = pd.to_numeric(self.backend_indexed[exp_qz_col], errors="coerce").to_numpy(float)
                valid = np.isfinite(exp_qr) & np.isfinite(exp_qz)
                if calc_qr_col and calc_qz_col:
                    calc_qr = pd.to_numeric(self.backend_indexed[calc_qr_col], errors="coerce").to_numpy(float)
                    calc_qz = pd.to_numeric(self.backend_indexed[calc_qz_col], errors="coerce").to_numpy(float)
                    connected = valid & np.isfinite(calc_qr) & np.isfinite(calc_qz)
                    for x1, y1, x2, y2 in zip(calc_qr[connected], calc_qz[connected], exp_qr[connected], exp_qz[connected]):
                        self.manual_axes.plot([x1, x2], [y1, y2], color="white", alpha=0.42, linewidth=0.65, zorder=5)
                    self.manual_axes.scatter(
                        calc_qr[connected], calc_qz[connected], s=20, marker="o",
                        facecolors="none", edgecolors="white", linewidths=0.8,
                        alpha=0.75, label="calculated matched positions", zorder=6,
                    )
                self.manual_axes.scatter(
                    exp_qr[valid], exp_qz[valid], s=58, marker="o", facecolors="none",
                    edgecolors="lime", linewidths=1.45,
                    label=f"indexed experimental points ({int(valid.sum())})", zorder=7,
                )
                if self.manual_labels.isChecked():
                    hkl_col = self._pick_column(self.backend_indexed, ["hkl", "overlay_hkl_text"])
                    for index in np.flatnonzero(valid):
                        if hkl_col:
                            label = str(self.backend_indexed.iloc[index][hkl_col])
                        elif all(name in self.backend_indexed.columns for name in ("h", "k", "l")):
                            row = self.backend_indexed.iloc[index]
                            label = f"({int(row.h)} {int(row.k)} {int(row.l)})"
                        else:
                            label = str(index + 1)
                        annotation = self.manual_axes.annotate(
                            label, (exp_qr[index], exp_qz[index]), xytext=(4, 4),
                            textcoords="offset points", fontsize=6.2, color="white",
                            annotation_clip=True, zorder=8,
                        )
                        annotation.set_path_effects([
                            matplotlib_patheffects.Stroke(linewidth=1.5, foreground="black"),
                            matplotlib_patheffects.Normal(),
                        ])
        if not self.backend_indexed.empty or not self.backend_ignored.empty:
            self.manual_axes.set_title("Experimental image with full-backend CIF indexing")
            handles, labels = self.manual_axes.get_legend_handles_labels()
            if handles:
                self.manual_axes.legend(fontsize=7.5, loc="upper right")
            self.manual_figure.tight_layout(pad=1.1)
            self.manual_canvas.draw_idle()




class _NonVisualToggle:
    __slots__ = ("_checked",)
    def __init__(self, checked):
        self._checked = bool(checked)
    def isChecked(self):
        return bool(self._checked)
    def setChecked(self, checked):
        self._checked = bool(checked)


class _NonVisualNumber:
    __slots__ = ("_value",)
    def __init__(self, value):
        self._value = value
    def value(self):
        return self._value
    def setValue(self, value):
        self._value = value


class _NullStatusLabel:
    __slots__ = ("_text",)
    def __init__(self):
        self._text = ""
    def setText(self, text):
        self._text = str(text)
    def text(self):
        return self._text
    def setWordWrap(self, _enabled):
        pass

class ManualCalculatorClickWorkbench(ManualFullIndexingWorkbench):
    """Full automatic NPZ indexing calculator plus manual PNG peak comparison.

    The automatic engine consumes one or more numerical NPZ measurements and
    determines the orientation, calibration, domains, and calculated reflection
    positions. The PNG is display-only. Experimental points are selected by the
    user on that PNG and paired with the automatic calculator's reflection list.
    """

    _ASSIGNMENT_COLUMNS = [
        "ID", "NormalizationReference", "HKL", "OrientationDomain", "CalculatorSource",
        "IndexingEvidenceSource", "DetectionSource", "InputProvenance",
        "QrExp", "QzExp", "QrCalc", "QzCalc", "QrCalcOriginal", "QzCalcOriginal",
        "QTotalExp", "QTotalCalc", "DSpacingExp", "DSpacingCalc", "DeltaD",
        "DeltaQr", "DeltaQz", "DeltaQ", "ResidualDirectionDeg",
        "ExpIntensity", "SigmaExpIntensity", "RelativeIntensityUncertainty",
        "CalcIntensity", "CalcSingleIntensity", "CalcComparisonIntensity", "CalcF2", "DWBAWeight",
        "CalcOverlapCount", "CalcOverlapRadius", "CalcOverlapHKLs", "CalcIntensityModel",
        "ExpRelative", "CalcRelative",
        "CalcScaledToExperiment", "ScaledResidual", "ScaledRelativeResidual",
        "IntensityAgreement", "NormalizationScaleFactor", "ReferenceExperimentalIntensity",
        "ReferenceCalculatedIntensity", "RobustScaleFactor", "CalcRobustScaledToExperiment",
        "LogIntensityResidual", "RobustIntensityAgreement",
        "EmpiricalOrientationScale", "CalcOrientationScaledToExperiment",
        "OrientationLogIntensityResidual", "OrientationIntensityAgreement",
        "OrientationCorrectionApplied", "OrientationCorrectionCVBaseline", "OrientationCorrectionCVModel",
        "IntensityScalePairCount", "ReferenceNPZ", "PixelX", "PixelY",
        "SubpixelX", "SubpixelY", "PeakHeight", "LocalBackground", "BackgroundNoise",
        "BackgroundGradient", "PeakSNR", "IntegratedSNR", "PeakAreaPixels",
        "PeakSigmaQrWidth", "PeakSigmaQzWidth", "OverlapPeakCount", "Deblended",
        "SaturationFraction", "ValidPixelFraction", "EdgeTruncated",
        "IntensityQualityScore", "IntensityQuality", "IntensityNormalizationFactor",
        "IntensityNormalization", "IntensitySource", "IntensityUncertaintySource", "PoissonTermUsed",
        "SigmaQrExp", "SigmaQzExp", "SigmaQrSystematic", "SigmaQzSystematic",
        "SigmaQrTotal", "SigmaQzTotal", "SigmaQExp", "PositionQuality", "UncertaintyRatio",
        "AssignmentScore", "BackendSupportScore", "SuggestedHKL", "SuggestedDeltaQ",
        "SuggestionImprovement", "ReassignmentRecommended", "SuggestionStatus", "AssignmentStability", "StabilityTier",
    ]

    # Human-facing table order. Keep internal/export field names stable while
    # placing physically meaningful numerical quantities before metadata.
    _ASSIGNMENT_DISPLAY_COLUMNS = [
        "ID", "HKL",
        # Reciprocal-space / crystallographic quantities
        "QrExp", "QzExp", "QrCalc", "QzCalc", "QrCalcOriginal", "QzCalcOriginal",
        "QTotalExp", "QTotalCalc", "DSpacingExp", "DSpacingCalc", "DeltaD",
        "DeltaQr", "DeltaQz", "DeltaQ", "ResidualDirectionDeg",
        "SigmaQrExp", "SigmaQzExp", "SigmaQrSystematic", "SigmaQzSystematic",
        "SigmaQrTotal", "SigmaQzTotal", "SigmaQExp",
        # Intensity-like quantities
        "ExpIntensity", "SigmaExpIntensity", "CalcComparisonIntensity", "CalcSingleIntensity",
        "CalcF2", "DWBAWeight", "ReferenceExperimentalIntensity",
        "ReferenceCalculatedIntensity", "CalcScaledToExperiment", "ScaledResidual",
        "CalcRobustScaledToExperiment", "CalcOrientationScaledToExperiment",
        "PeakHeight", "LocalBackground", "BackgroundNoise", "BackgroundGradient",
        # Image-coordinate quantities
        "PixelX", "PixelY", "SubpixelX", "SubpixelY",
        # Dimensionless diagnostics
        "ExpRelative", "CalcRelative", "ScaledRelativeResidual", "IntensityAgreement",
        "NormalizationScaleFactor", "RobustScaleFactor", "LogIntensityResidual",
        "RobustIntensityAgreement", "EmpiricalOrientationScale", "OrientationLogIntensityResidual",
        "OrientationIntensityAgreement", "OrientationCorrectionCVBaseline", "OrientationCorrectionCVModel",
        "PeakSNR", "IntegratedSNR", "RelativeIntensityUncertainty", "PeakAreaPixels",
        "PeakSigmaQrWidth", "PeakSigmaQzWidth", "SaturationFraction", "ValidPixelFraction",
        "IntensityQualityScore", "IntensityNormalizationFactor", "OverlapPeakCount",
        "CalcOverlapCount", "CalcOverlapRadius", "IntensityScalePairCount",
        "UncertaintyRatio", "AssignmentScore", "BackendSupportScore", "AssignmentStability",
        "SuggestedDeltaQ", "SuggestionImprovement",
        # Categorical / provenance information
        "NormalizationReference", "PositionQuality", "IntensityQuality", "Deblended", "EdgeTruncated",
        "PoissonTermUsed", "CalcOverlapHKLs", "CalcIntensityModel", "IntensityNormalization",
        "IntensitySource", "IntensityUncertaintySource", "OrientationCorrectionApplied", "StabilityTier", "SuggestedHKL",
        "ReassignmentRecommended", "SuggestionStatus", "OrientationDomain", "CalculatorSource",
        "IndexingEvidenceSource", "DetectionSource", "ReferenceNPZ", "InputProvenance",
    ]

    _ASSIGNMENT_DISPLAY_HEADERS = {
        "ID": "ID",
        "HKL": "HKL",
        "QrExp": "qᵣ exp (Å⁻¹)",
        "QzExp": "q_z exp (Å⁻¹)",
        "QrCalc": "qᵣ calc/refined (Å⁻¹)",
        "QzCalc": "q_z calc/refined (Å⁻¹)",
        "QrCalcOriginal": "qᵣ calc original (Å⁻¹)",
        "QzCalcOriginal": "q_z calc original (Å⁻¹)",
        "QTotalExp": "q total exp (Å⁻¹)",
        "QTotalCalc": "q total calc (Å⁻¹)",
        "DSpacingExp": "d exp (Å)",
        "DSpacingCalc": "d calc (Å)",
        "DeltaD": "Δd (Å)",
        "DeltaQr": "Δqᵣ (Å⁻¹)",
        "DeltaQz": "Δq_z (Å⁻¹)",
        "DeltaQ": "Δq total (Å⁻¹)",
        "ResidualDirectionDeg": "Residual direction (deg)",
        "SigmaQrExp": "σ qᵣ centroid (Å⁻¹)",
        "SigmaQzExp": "σ q_z centroid (Å⁻¹)",
        "SigmaQrSystematic": "σ qᵣ mapping (Å⁻¹)",
        "SigmaQzSystematic": "σ q_z mapping (Å⁻¹)",
        "SigmaQrTotal": "σ qᵣ total (Å⁻¹)",
        "SigmaQzTotal": "σ q_z total (Å⁻¹)",
        "SigmaQExp": "σ q total exp (Å⁻¹)",
        "ExpIntensity": "Exp integrated intensity (exp a.u.)",
        "SigmaExpIntensity": "σ exp intensity (exp a.u.)",
        "RelativeIntensityUncertainty": "Relative intensity uncertainty (unitless)",
        "CalcIntensity": "Calc selected-reflection intensity (calc a.u.)",
        "CalcSingleIntensity": "Calc single-reflection intensity (calc a.u.)",
        "CalcComparisonIntensity": "Calc unresolved-spot intensity (calc a.u.)",
        "CalcF2": "Kinematic F² (calc a.u.)",
        "DWBAWeight": "DWBA weight (unitless)",
        "CalcOverlapCount": "Calc unresolved reflection count",
        "CalcOverlapRadius": "Calc overlap radius (Å⁻¹)",
        "CalcOverlapHKLs": "Calc unresolved HKLs",
        "CalcIntensityModel": "Calculated intensity model",
        "ReferenceExperimentalIntensity": "Reference exp intensity (exp a.u.)",
        "ReferenceCalculatedIntensity": "Reference calc intensity (calc a.u.)",
        "CalcScaledToExperiment": "Calc reference-scaled (exp a.u.)",
        "ScaledResidual": "Reference-scaled residual (exp a.u.)",
        "CalcRobustScaledToExperiment": "Calc robust-scaled (exp a.u.)",
        "PeakHeight": "Peak height (exp a.u.)",
        "LocalBackground": "Local background (exp a.u.)",
        "PixelX": "Pixel X (px)",
        "PixelY": "Pixel Y (px)",
        "SubpixelX": "Subpixel X (px)",
        "SubpixelY": "Subpixel Y (px)",
        "ExpRelative": "Exp relative (unitless)",
        "CalcRelative": "Calc relative (unitless)",
        "ScaledRelativeResidual": "Reference-scaled relative residual (unitless)",
        "IntensityAgreement": "Reference intensity agreement (unitless)",
        "NormalizationScaleFactor": "Reference scale factor (exp/calc)",
        "RobustScaleFactor": "Robust multi-peak scale (exp/calc)",
        "LogIntensityResidual": "log10 intensity residual (unitless)",
        "RobustIntensityAgreement": "Robust intensity agreement (unitless)",
        "EmpiricalOrientationScale": "Empirical orientation scale (exp/calc)",
        "CalcOrientationScaledToExperiment": "Calc orientation-corrected (exp a.u.)",
        "OrientationLogIntensityResidual": "Orientation log10 residual (unitless)",
        "OrientationIntensityAgreement": "Orientation intensity agreement (unitless)",
        "OrientationCorrectionApplied": "Empirical orientation correction applied",
        "OrientationCorrectionCVBaseline": "CV baseline log-RMSE (dex)",
        "OrientationCorrectionCVModel": "CV orientation log-RMSE (dex)",
        "IntensityScalePairCount": "Intensity scale trusted-pair count",
        "PeakSNR": "Peak-height SNR (unitless)",
        "IntegratedSNR": "Integrated intensity SNR (unitless)",
        "PeakAreaPixels": "Integrated peak area (px)",
        "PeakSigmaQrWidth": "Peak σ width qᵣ (Å⁻¹)",
        "PeakSigmaQzWidth": "Peak σ width q_z (Å⁻¹)",
        "BackgroundNoise": "Local background σ (exp a.u.)",
        "BackgroundGradient": "Background gradient (exp a.u./px)",
        "OverlapPeakCount": "Nearby experimental peak count",
        "Deblended": "Local deblending applied",
        "SaturationFraction": "Clipping/saturation fraction (unitless)",
        "ValidPixelFraction": "Valid aperture fraction (unitless)",
        "EdgeTruncated": "Peak aperture edge-truncated",
        "IntensityQualityScore": "Intensity quality score (unitless)",
        "IntensityQuality": "Intensity quality",
        "IntensityNormalizationFactor": "Experimental normalization factor",
        "IntensityNormalization": "Experimental normalization",
        "IntensitySource": "Experimental intensity source",
        "IntensityUncertaintySource": "Experimental intensity uncertainty source",
        "PoissonTermUsed": "Poisson uncertainty term used",
        "UncertaintyRatio": "Δq / σ total (unitless)",
        "AssignmentScore": "Candidate assignment score (unitless)",
        "BackendSupportScore": "Backend support score (unitless)",
        "AssignmentStability": "Assignment stability (unitless)",
        "SuggestedDeltaQ": "Suggested Δq (Å⁻¹)",
        "SuggestionImprovement": "Suggested improvement (Å⁻¹)",
        "NormalizationReference": "Normalization reference",
        "PositionQuality": "Position quality",
        "StabilityTier": "Stability tier",
        "SuggestedHKL": "Suggested one-to-one HKL",
        "ReassignmentRecommended": "Reassignment recommended",
        "SuggestionStatus": "One-to-one review status",
        "OrientationDomain": "Orientation domain",
        "CalculatorSource": "Calculator source",
        "IndexingEvidenceSource": "Indexing evidence source",
        "DetectionSource": "Detection source",
        "ReferenceNPZ": "Reference NPZ",
        "InputProvenance": "Input provenance",
    }

    # Positional-audit thresholds are centralized so the table classifications,
    # summary percentages, and exported audit use the same definitions.
    _VERY_GOOD_DELTA_Q = 0.020
    _ACCEPTABLE_DELTA_Q = 0.040
    _BORDERLINE_DELTA_Q = 0.060
    _WITHIN_UNCERTAINTY_RATIO = 2.0

    @classmethod
    def _position_quality_metrics(cls, delta_q, sigma_qr, sigma_qz):
        """Return uncertainty ratios and categorical positional-quality labels."""
        delta_q = np.asarray(delta_q, dtype=float)
        sigma_total = np.hypot(
            np.asarray(sigma_qr, dtype=float),
            np.asarray(sigma_qz, dtype=float),
        )
        uncertainty_ratio = np.full(delta_q.shape, np.nan, dtype=float)
        usable_sigma = np.isfinite(sigma_total) & (sigma_total > 1e-8)
        uncertainty_ratio[usable_sigma] = delta_q[usable_sigma] / sigma_total[usable_sigma]

        quality = np.full(delta_q.shape, "unassessed", dtype=object)
        finite = np.isfinite(delta_q)
        quality[finite & (delta_q <= cls._VERY_GOOD_DELTA_Q)] = "very_good"
        quality[
            finite
            & (delta_q > cls._VERY_GOOD_DELTA_Q)
            & (delta_q <= cls._ACCEPTABLE_DELTA_Q)
        ] = "acceptable"
        quality[
            finite
            & (delta_q > cls._ACCEPTABLE_DELTA_Q)
            & (delta_q <= cls._BORDERLINE_DELTA_Q)
        ] = "borderline"
        quality[finite & (delta_q > cls._BORDERLINE_DELTA_Q)] = "poor"
        quality[
            usable_sigma
            & (uncertainty_ratio <= cls._WITHIN_UNCERTAINTY_RATIO)
            & (delta_q <= cls._ACCEPTABLE_DELTA_Q)
        ] = "within_uncertainty"
        return uncertainty_ratio, quality

    @staticmethod
    def _relative_intensity_statistics(experimental, calculated):
        """Return Pearson, Spearman, and RMSE for paired relative intensities."""
        experimental = np.asarray(experimental, dtype=float)
        calculated = np.asarray(calculated, dtype=float)
        valid = np.isfinite(experimental) & np.isfinite(calculated)
        pearson = spearman = rmse = np.nan
        if int(valid.sum()) >= 2:
            rmse = float(
                np.sqrt(np.mean((experimental[valid] - calculated[valid]) ** 2))
            )
        if (
            int(valid.sum()) >= 3
            and np.nanstd(experimental[valid]) > 0
            and np.nanstd(calculated[valid]) > 0
        ):
            exp_series = pd.Series(experimental[valid])
            calc_series = pd.Series(calculated[valid])
            pearson = float(exp_series.corr(calc_series, method="pearson"))
            spearman = float(exp_series.corr(calc_series, method="spearman"))
        return pearson, spearman, rmse

    @staticmethod
    def _safe_d_spacing(q_values):
        q_values = np.asarray(q_values, dtype=float)
        result = np.full(q_values.shape, np.nan, dtype=float)
        valid = np.isfinite(q_values) & (np.abs(q_values) > 1e-12)
        result[valid] = (2.0 * np.pi) / np.abs(q_values[valid])
        return result

    @staticmethod
    def _prediction_support_score_values(calculator_source, evidence_source, stability):
        """Map backend evidence/stability to a conservative 0..1 support score."""
        source = str(calculator_source or "").lower()
        evidence = str(evidence_source or "").lower()
        if "indexed solution" in source:
            base = 0.82
        elif "unused" in source:
            base = 0.28
        else:
            base = 0.50
        if "joint_multidomain" in evidence or "held" in evidence or "validated" in evidence:
            base = max(base, 0.92)
        elif "completion" in evidence:
            base = max(base, 0.70)
        elif "unused_calculated_prediction" in evidence:
            base = min(base, 0.30)
        try:
            stability_value = float(stability)
        except Exception:
            stability_value = np.nan
        if np.isfinite(stability_value):
            stability_value = float(np.clip(stability_value, 0.0, 1.0))
            base = 0.65 * base + 0.35 * stability_value
        return float(np.clip(base, 0.0, 1.0))

    @classmethod
    def _robust_intensity_fit(cls, experimental, calculated, trusted_mask=None,
                              sigma_experimental=None, quality_score=None):
        """Fit a robust multiplicative scale in log space across trusted peaks.

        Experimental intensity uncertainty and the automatically derived
        intensity-quality score are used as weights when available.  The fit is
        still only a relative scale between different arbitrary units.
        """
        experimental = np.asarray(experimental, dtype=float)
        calculated = np.asarray(calculated, dtype=float)
        valid = (
            np.isfinite(experimental) & np.isfinite(calculated)
            & (experimental > 0) & (calculated > 0)
        )
        if trusted_mask is not None:
            trusted_mask = np.asarray(trusted_mask, dtype=bool)
            trusted_valid = valid & trusted_mask
            if int(trusted_valid.sum()) >= 3:
                valid = trusted_valid
        if int(valid.sum()) < 2:
            return np.nan, np.full(experimental.shape, np.nan), np.full(experimental.shape, np.nan), np.nan
        log_ratio_all = np.full(experimental.shape, np.nan, dtype=float)
        log_ratio_all[valid] = np.log(experimental[valid]) - np.log(calculated[valid])
        log_ratio = log_ratio_all[valid]
        weights = np.ones(len(log_ratio), dtype=float)
        if sigma_experimental is not None:
            sigma = np.asarray(sigma_experimental, dtype=float)[valid]
            relative = sigma / np.maximum(experimental[valid], 1e-12)
            relative = np.where(np.isfinite(relative) & (relative > 0), np.clip(relative, 0.03, 2.0), 0.35)
            weights *= 1.0 / (relative * relative)
        if quality_score is not None:
            quality = np.asarray(quality_score, dtype=float)[valid]
            quality = np.where(np.isfinite(quality), np.clip(quality, 0.10, 1.0), 0.50)
            weights *= quality * quality
        finite_positive = np.isfinite(weights) & (weights > 0)
        if finite_positive.any():
            weights /= max(float(np.nanmedian(weights[finite_positive])), 1e-12)
        weights = np.clip(weights, 0.05, 20.0)
        initial = float(np.nanmedian(log_ratio))
        try:
            fit = least_squares(
                lambda x: (log_ratio - float(x[0])) * np.sqrt(weights),
                x0=np.array([initial], dtype=float),
                loss="soft_l1",
                f_scale=0.35,
            )
            log_scale = float(fit.x[0])
        except Exception:
            log_scale = initial
        scale = float(np.exp(log_scale))
        scaled = calculated * scale
        log_residual = np.full(experimental.shape, np.nan, dtype=float)
        compare = (
            np.isfinite(experimental) & np.isfinite(scaled)
            & (experimental > 0) & (scaled > 0)
        )
        log_residual[compare] = np.log10(experimental[compare] / scaled[compare])
        rmse = (
            float(np.sqrt(np.average(log_residual[valid] ** 2, weights=weights)))
            if int(valid.sum()) else np.nan
        )
        return scale, scaled, log_residual, rmse

    @staticmethod
    def _orientation_intensity_basis(qr, qz):
        """Smooth low-order basis for empirical thin-film orientation weighting."""
        qr = np.asarray(qr, dtype=float)
        qz = np.asarray(qz, dtype=float)
        chi = np.arctan2(qr, qz)
        return np.column_stack([
            np.ones(len(chi)),
            np.sin(2.0 * chi), np.cos(2.0 * chi),
            np.sin(4.0 * chi), np.cos(4.0 * chi),
        ])

    @classmethod
    def _cross_validated_orientation_intensity_correction(
            cls, experimental, calculated, qr, qz, trusted_mask,
            sigma_experimental=None, quality_score=None):
        """Fit an empirical angular intensity envelope only when LOO-CV supports it.

        This is not presented as a crystallographic structure factor correction.
        It is an empirical preferred-orientation/texture envelope applied only to
        the manual intensity comparison.  It never participates in HKL indexing.
        """
        experimental = np.asarray(experimental, dtype=float)
        calculated = np.asarray(calculated, dtype=float)
        qr = np.asarray(qr, dtype=float)
        qz = np.asarray(qz, dtype=float)
        trusted = np.asarray(trusted_mask, dtype=bool)
        valid = trusted & np.isfinite(experimental) & np.isfinite(calculated) & (experimental > 0) & (calculated > 0)
        if int(valid.sum()) < 8:
            n = len(experimental)
            return False, np.full(n, np.nan), np.full(n, np.nan), np.full(n, np.nan), np.nan, np.nan
        chi = np.degrees(np.arctan2(qr[valid], qz[valid]))
        if float(np.nanmax(chi) - np.nanmin(chi)) < 12.0:
            n = len(experimental)
            return False, np.full(n, np.nan), np.full(n, np.nan), np.full(n, np.nan), np.nan, np.nan
        X_all = cls._orientation_intensity_basis(qr, qz)
        indices = np.where(valid)[0]
        y_all = np.log(experimental / calculated)
        weights_all = np.ones(len(experimental), dtype=float)
        if sigma_experimental is not None:
            sigma = np.asarray(sigma_experimental, dtype=float)
            rel = sigma / np.maximum(experimental, 1e-12)
            rel = np.where(np.isfinite(rel) & (rel > 0), np.clip(rel, 0.03, 2.0), 0.35)
            weights_all *= 1.0 / (rel * rel)
        if quality_score is not None:
            quality = np.asarray(quality_score, dtype=float)
            quality = np.where(np.isfinite(quality), np.clip(quality, 0.10, 1.0), 0.50)
            weights_all *= quality * quality
        weights_all = np.clip(weights_all / max(float(np.nanmedian(weights_all[valid])), 1e-12), 0.05, 20.0)

        def fit_beta(train_indices):
            X = X_all[train_indices]
            y = y_all[train_indices]
            w = weights_all[train_indices]
            root_w = np.sqrt(w)
            ridge = np.diag([0.0, 0.35, 0.35, 0.60, 0.60])
            lhs = (X * root_w[:, None]).T @ (X * root_w[:, None]) + ridge
            rhs = (X * root_w[:, None]).T @ (y * root_w)
            try:
                return np.linalg.solve(lhs, rhs)
            except Exception:
                return np.linalg.lstsq(lhs, rhs, rcond=None)[0]

        base_errors = []
        model_errors = []
        for held in indices:
            train = indices[indices != held]
            if len(train) < 7:
                continue
            baseline = float(np.nanmedian(y_all[train]))
            beta = fit_beta(train)
            base_errors.append((y_all[held] - baseline) / math.log(10.0))
            model_errors.append((y_all[held] - float(X_all[held] @ beta)) / math.log(10.0))
        if len(model_errors) < 6:
            n = len(experimental)
            return False, np.full(n, np.nan), np.full(n, np.nan), np.full(n, np.nan), np.nan, np.nan
        cv_baseline = float(np.sqrt(np.mean(np.square(base_errors))))
        cv_model = float(np.sqrt(np.mean(np.square(model_errors))))
        # Require both a relative and absolute improvement to avoid activating a
        # fitted texture envelope for negligible numerical gains.
        applied = bool(
            np.isfinite(cv_baseline) and np.isfinite(cv_model)
            and cv_model <= 0.90 * cv_baseline
            and (cv_baseline - cv_model) >= 0.025
        )
        n = len(experimental)
        if not applied:
            return False, np.full(n, np.nan), np.full(n, np.nan), np.full(n, np.nan), cv_baseline, cv_model
        beta = fit_beta(indices)
        log_factor = X_all @ beta
        # The intercept is the arbitrary experimental/calculated unit conversion
        # and may legitimately span many orders of magnitude. Bound only the
        # angular texture modulation around that intercept, not the global scale.
        # This prevents extrapolated orientation weights from becoming absurd
        # while preserving a physically necessary arbitrary-unit normalization.
        global_log_scale = float(beta[0])
        relative_log_factor = np.clip(
            log_factor - global_log_scale, math.log(0.08), math.log(12.5)
        )
        factor = np.exp(global_log_scale + relative_log_factor)
        scaled = calculated * factor
        residual = np.full(n, np.nan, dtype=float)
        compare = np.isfinite(experimental) & np.isfinite(scaled) & (experimental > 0) & (scaled > 0)
        residual[compare] = np.log10(experimental[compare] / scaled[compare])
        return True, factor, scaled, residual, cv_baseline, cv_model

    @staticmethod
    def _transform_q_points(qr, qz, parameters):
        """Apply the bounded 2-D q-map refinement used only for manual comparison."""
        theta, scale_qr, scale_qz, offset_qr, offset_qz = [float(v) for v in parameters]
        c = math.cos(theta)
        s = math.sin(theta)
        qr = np.asarray(qr, dtype=float)
        qz = np.asarray(qz, dtype=float)
        rotated_qr = c * qr - s * qz
        rotated_qz = s * qr + c * qz
        return (
            scale_qr * rotated_qr + offset_qr,
            scale_qz * rotated_qz + offset_qz,
        )

    @classmethod
    def _fit_q_mapping_parameters(cls, qr_calc, qz_calc, qr_exp, qz_exp, sigma_qr=None, sigma_qz=None):
        """Robustly fit a small rotation, anisotropic scale, and q offsets.

        This is a 2-D reciprocal-space mapping refinement for the displayed
        calculator predictions. It is intentionally bounded and regularized; it
        is not a replacement for the backend's full crystallographic orientation
        search.
        """
        qr_calc = np.asarray(qr_calc, dtype=float)
        qz_calc = np.asarray(qz_calc, dtype=float)
        qr_exp = np.asarray(qr_exp, dtype=float)
        qz_exp = np.asarray(qz_exp, dtype=float)
        n = len(qr_calc)
        if sigma_qr is None:
            sigma_qr = np.full(n, 0.01, dtype=float)
        if sigma_qz is None:
            sigma_qz = np.full(n, 0.01, dtype=float)
        sigma_qr = np.asarray(sigma_qr, dtype=float)
        sigma_qz = np.asarray(sigma_qz, dtype=float)
        sigma_qr = np.where(np.isfinite(sigma_qr) & (sigma_qr > 0), np.maximum(sigma_qr, 0.006), 0.01)
        sigma_qz = np.where(np.isfinite(sigma_qz) & (sigma_qz > 0), np.maximum(sigma_qz, 0.006), 0.01)

        def residual(parameters):
            pred_qr, pred_qz = cls._transform_q_points(qr_calc, qz_calc, parameters)
            data = np.concatenate([
                (pred_qr - qr_exp) / sigma_qr,
                (pred_qz - qz_exp) / sigma_qz,
            ])
            theta, scale_r, scale_z, off_r, off_z = parameters
            # Mild regularization discourages a visually good but physically
            # implausibly large correction when only a few trusted peaks exist.
            reg = 0.20 * np.array([
                theta / np.deg2rad(1.0),
                (scale_r - 1.0) / 0.02,
                (scale_z - 1.0) / 0.02,
                off_r / 0.02,
                off_z / 0.02,
            ])
            return np.concatenate([data, reg])

        lower = np.array([np.deg2rad(-3.0), 0.95, 0.95, -0.08, -0.08], dtype=float)
        upper = np.array([np.deg2rad(3.0), 1.05, 1.05, 0.08, 0.08], dtype=float)
        fit = least_squares(
            residual,
            x0=np.array([0.0, 1.0, 1.0, 0.0, 0.0], dtype=float),
            bounds=(lower, upper),
            loss="soft_l1",
            f_scale=1.0,
            max_nfev=3000,
        )
        return np.asarray(fit.x, dtype=float)

    def _mapping_systematic_q_for_path(self, path):
        """Estimate PNG-to-NPZ coordinate uncertainty from the configured pixel-position uncertainty."""
        try:
            resolved = str(_GuiPath(path).expanduser().resolve())
            image = self._reference_npz_cache.get(resolved)
            if image is None:
                image = self._probe_npz(resolved)
                self._reference_npz_cache[resolved] = image
            qr = np.asarray(image["qr"], dtype=float)
            qz = np.asarray(image["qz"], dtype=float)
            qr_step = float(np.nanmedian(np.abs(np.diff(qr)))) if len(qr) > 1 else np.nan
            qz_step = float(np.nanmedian(np.abs(np.diff(qz)))) if len(qz) > 1 else np.nan
            pixels = float(getattr(self, "click_mapping_uncertainty_pixels", None).value()) if getattr(self, "click_mapping_uncertainty_pixels", None) is not None else 1.0
            return abs(qr_step) * pixels, abs(qz_step) * pixels
        except Exception:
            return np.nan, np.nan

    def _augment_experimental_uncertainty(self, result, qr_axis, qz_axis):
        qr_axis = np.asarray(qr_axis, dtype=float)
        qz_axis = np.asarray(qz_axis, dtype=float)
        qr_step = float(np.nanmedian(np.abs(np.diff(qr_axis)))) if len(qr_axis) > 1 else np.nan
        qz_step = float(np.nanmedian(np.abs(np.diff(qz_axis)))) if len(qz_axis) > 1 else np.nan
        pixels = float(self.click_mapping_uncertainty_pixels.value()) if getattr(self, "click_mapping_uncertainty_pixels", None) is not None else 1.0
        sys_qr = abs(qr_step) * pixels if np.isfinite(qr_step) else np.nan
        sys_qz = abs(qz_step) * pixels if np.isfinite(qz_step) else np.nan
        stat_qr = float(result.get("SigmaQrExp", np.nan))
        stat_qz = float(result.get("SigmaQzExp", np.nan))
        total_qr = math.hypot(stat_qr if np.isfinite(stat_qr) else 0.0, sys_qr if np.isfinite(sys_qr) else 0.0)
        total_qz = math.hypot(stat_qz if np.isfinite(stat_qz) else 0.0, sys_qz if np.isfinite(sys_qz) else 0.0)
        if not np.isfinite(stat_qr) and not np.isfinite(sys_qr):
            total_qr = np.nan
        if not np.isfinite(stat_qz) and not np.isfinite(sys_qz):
            total_qz = np.nan
        result["SigmaQrSystematic"] = sys_qr
        result["SigmaQzSystematic"] = sys_qz
        result["SigmaQrTotal"] = total_qr
        result["SigmaQzTotal"] = total_qz
        q_total = math.hypot(float(result["QrExp"]), float(result["QzExp"]))
        if q_total > 1e-12 and np.isfinite(total_qr) and np.isfinite(total_qz):
            result["SigmaQExp"] = math.sqrt(
                (float(result["QrExp"]) / q_total * total_qr) ** 2
                + (float(result["QzExp"]) / q_total * total_qz) ** 2
            )
        else:
            result["SigmaQExp"] = np.nan
        return result

    def _candidate_score_arrays(self, dr, dz, predictions, exp_point):
        """Rank candidates primarily by position, with backend support as a tie-breaker."""
        dr = np.asarray(dr, dtype=float)
        dz = np.asarray(dz, dtype=float)
        distance = np.hypot(dr, dz)
        sigma_r = float(exp_point.get("SigmaQrTotal", np.nan))
        sigma_z = float(exp_point.get("SigmaQzTotal", np.nan))
        if not np.isfinite(sigma_r) or sigma_r <= 0:
            sigma_r = 0.010
        if not np.isfinite(sigma_z) or sigma_z <= 0:
            sigma_z = 0.010
        sigma_r = max(sigma_r, 0.006)
        sigma_z = max(sigma_z, 0.006)
        normalized = np.hypot(dr / sigma_r, dz / sigma_z)
        support_series = predictions["BackendSupportScore"] if "BackendSupportScore" in predictions.columns else pd.Series(0.5, index=predictions.index)
        support = pd.to_numeric(support_series, errors="coerce").to_numpy(float)
        support = np.where(np.isfinite(support), np.clip(support, 0.0, 1.0), 0.5)
        stability_series = predictions["AssignmentStability"] if "AssignmentStability" in predictions.columns else pd.Series(np.nan, index=predictions.index)
        stability = pd.to_numeric(stability_series, errors="coerce").to_numpy(float)
        stability_penalty = np.where(np.isfinite(stability), 0.35 * (1.0 - np.clip(stability, 0.0, 1.0)), 0.18)
        support_penalty = 0.75 * (1.0 - support)
        score = normalized + support_penalty + stability_penalty
        return score, distance, support

    def __init__(self):
        # Internal defaults for click uncertainty, candidate ranking, one-to-one
        # assignment, rejection, and refinement. These safeguards run automatically
        # without adding extra controls to the primary analysis workflow.
        self.click_mapping_uncertainty_pixels = _NonVisualNumber(1.0)
        self.use_support_candidate_ranking = _NonVisualToggle(True)
        self.enforce_one_to_one_assignments = _NonVisualToggle(True)
        self.block_very_poor_assignments = _NonVisualToggle(True)
        self.assignment_reject_threshold = _NonVisualNumber(0.10)
        self.refinement_trusted_threshold = _NonVisualNumber(0.04)
        self.refinement_min_pairs = _NonVisualNumber(5)
        self.refinement_status_label = _NullStatusLabel()
        # GIWAXS physics selection is evaluated automatically from measurement
        # metadata even though no dedicated status-row widget is displayed.
        self.giwaxs_physics_status_label = None
        self._hybrid_ui_ready = False
        self.calculator_predictions = pd.DataFrame()
        self.manual_click_assignments = pd.DataFrame(columns=self._ASSIGNMENT_COLUMNS)
        self.manual_click_assignments_by_series = {}
        self.pending_clicked_experimental = None
        self._candidate_prediction_indices = []
        self._reference_npz_cache = {}
        self._manual_intensity_physics_cache = {}
        self._displayed_calculator_prediction_ids = set()
        self._q_refinement = None
        self._last_suggestion_summary = ""
        super().__init__()
        self.setWindowTitle(
            "GIWAXS/GIXS — Automatic NPZ Calculator + Manual PNG Peak Indexing"
        )
        self._configure_hybrid_workspace()
        self._hybrid_ui_ready = True
        self._redraw_manual()

    # -------------------------- display-only PNG loading -----------------------
    def _choose_manual_image(self):
        """Load the overlay image without adding it to the NPZ calculator list.

        This workspace intentionally separates the two roles: PNG is display-only
        and NPZ is calculator-only. The base image chooser is called directly so
        loading a display image never changes the NPZ measurement table.
        """
        previous = str(self.manual_image_path or "")
        IntegratedGIXSWorkbench._choose_manual_image(self)
        current = str(self.manual_image_path or "")
        if not current or current == previous:
            return

        # Enforce the numerical-input contract of the automatic calculator:
        # measurement rows must reference NPZ data, while PNG remains display-only.
        removed = 0
        if hasattr(self, "measurement_table"):
            for row in range(self.measurement_table.rowCount() - 1, -1, -1):
                item = self.measurement_table.item(row, 1)
                candidate = item.text().strip() if item else ""
                if candidate and _GuiPath(candidate).suffix.lower() != ".npz":
                    self.measurement_table.removeRow(row)
                    removed += 1

        if hasattr(self, "_refresh_reference_npz_combo"):
            self._refresh_reference_npz_combo()
        if hasattr(self, "_auto_select_reference_npz_for_png"):
            self._auto_select_reference_npz_for_png()

        message = (
            f"Loaded overlay image: {_GuiPath(current).name}. "
            "It is display-only and was not added to the automatic-indexing inputs. "
            "Use 'Add NPZ calculation files…' to add calculator data."
        )
        if removed:
            message += f" Removed {removed} non-NPZ row(s) from the calculator list."
        self.manual_status.setText(message)

    # ----------------------------- GUI conversion -----------------------------
    def _configure_hybrid_workspace(self):
        # Configure the combined workspace to use the supported indexing workflows
        # while keeping the interface focused on primary analysis controls.
        if hasattr(self, "manual_workflow_combo"):
            self.manual_workflow_combo.setCurrentText("Improved coverage indexing")

        # Hide automatic experimental markers in the manual-comparison plot.
        # Calculated reflections come from the backend model; experimental points
        # are selected directly from the displayed PNG by the analyst.
        for name in ("show_full_simulation", "show_indexed_points", "show_not_indexed_points"):
            widget = getattr(self, name, None)
            if widget is not None:
                widget.setVisible(False)

        # Restrict automatic-calculator inputs to numerical NPZ measurements.
        # The PNG chooser is reserved for the display-only experimental overlay.
        for button in self.findChildren(QPushButton):
            text = button.text()
            if text == "Add angle images…":
                button.setText("Add NPZ calculation files…")
                try:
                    button.clicked.disconnect()
                except TypeError:
                    pass
                button.clicked.connect(self._add_npz_calculation_files)
            elif text == "Remove selected":
                try:
                    button.clicked.disconnect()
                except TypeError:
                    pass
                button.clicked.connect(self._remove_npz_calculation_files)
            elif text == "Use current image":
                button.setVisible(False)
            elif text == "Run full indexing onto this overlay":
                button.setText("Run full automatic indexing calculator")
            elif text in {
                "Calculate manual overlay",
                "Optional quick CIF simulation",
            }:
                button.setText("Generate selected space-group reflection pattern")
            elif text == "Open approved backend overlay":
                button.setVisible(False)

        # The auto calculator should not receive the PNG as a measurement.
        for label in self.findChildren(QLabel):
            if "Edit Sample, Series, Angle, and Scan directly" in label.text():
                label.setText(
                    "Add one or more NPZ measurements. Edit Sample, Series, Angle, and Scan here. "
                    "The PNG above is display-only and is never sent to the automatic calculator. "
                    "NPZ files without qr/qz axes use the q limits shown above."
                )
            elif "The interface stays manual-style" in label.text():
                label.setText(
                    "The full automatic indexing calculator runs on the NPZ files: multiscale feature "
                    "detection, registration, consensus, orientation search, calibration, recovery, "
                    "multidomain testing, and validation. Its calculated reflection positions are then "
                    "shown on the PNG so you can select the experimental peaks manually."
                )

        # Rename the existing table section and repurpose it for manual comparisons.
        for group in self.findChildren(QGroupBox):
            if group.title().startswith("Indexed HKL, experimental/calculated"):
                group.setTitle(
                    "Manually clicked experimental peaks vs automatic-calculator reflections"
                )

        self.full_correlation_label.setText(
            "Run the automatic NPZ calculator, then click an experimental peak on the PNG and its "
            "calculated reflection. The strongest matched experimental peak remains the reference-pair "
            "normalization, while robust trusted-peak scaling is also calculated automatically."
        )

        # Add controls for pairing a selected experimental peak with a calculated
        # reflection adjacent to the reciprocal-space display.
        root_layout = self.manual_tab.layout()
        vertical = root_layout.itemAt(0).widget()
        horizontal = vertical.widget(0)
        scroll = horizontal.widget(0)
        controls_layout = scroll.widget().layout()

        group = QGroupBox("Manual experimental peak selection and comparison")
        layout = QVBoxLayout(group)
        explanation = QLabel(
            "1) Load a PNG above. 2) Add NPZ files and run the full automatic calculator. "
            "3) Left-click an experimental peak on the PNG. 4) Left-click its calculated white circle "
            "or choose an HKL below. The click position is selected on the PNG; numerical q and "
            "intensity are sampled from the chosen reference NPZ at that location."
        )
        explanation.setWordWrap(True)
        layout.addWidget(explanation)

        # Basic workflow controls
        basic_group = QGroupBox("Basic workflow")
        basic_form = QFormLayout(basic_group)
        self.click_reference_npz_combo = QComboBox()
        self.click_reference_npz_combo.currentIndexChanged.connect(self._reference_npz_changed)
        self.click_use_npz_axes = QCheckBox("Frame the PNG using the selected NPZ q axes")
        self.click_use_npz_axes.setChecked(True)
        self.click_use_npz_axes.stateChanged.connect(lambda _=None: self._reference_npz_changed())
        self.click_snap_to_max = QCheckBox("Snap click to the nearby experimental peak center")
        self.click_snap_to_max.setChecked(True)
        self.click_snap_radius = self._integer(10, 0, 50)
        self.click_snap_radius.setToolTip(
            "Search radius used to move an edge click toward the center of the nearby experimental peak."
        )
        self.click_calc_tolerance = self._double(0.10, 0.001, 1.0, 0.005, 4)
        self.click_candidate_count = self._integer(12, 1, 100)
        basic_form.addRow("Reference NPZ for clicked intensity:", self.click_reference_npz_combo)
        basic_form.addRow(self.click_use_npz_axes)
        basic_form.addRow(self.click_snap_to_max)
        basic_form.addRow("Local snap radius (pixels):", self.click_snap_radius)
        basic_form.addRow("Calculated-point click tolerance (Å⁻¹):", self.click_calc_tolerance)
        basic_form.addRow("Nearest HKL candidates:", self.click_candidate_count)
        layout.addWidget(basic_group)

        # Display controls
        display_group = QGroupBox("Display controls")
        display_form = QFormLayout(display_group)
        self.show_calculator_points = QCheckBox("Show automatic-calculator simulation points")
        self.show_calculator_points.setChecked(True)
        self.show_calculator_points.stateChanged.connect(lambda _=None: self._redraw_manual())
        self.calculator_strongest_only = QCheckBox(
            "Show only the strongest automatic-calculator reflection positions"
        )
        self.calculator_strongest_only.setChecked(True)
        self.calculator_strongest_only.stateChanged.connect(lambda _=None: self._redraw_manual())
        self.calculator_max_hkls = self._integer(120, 1, 10000)
        self.calculator_max_hkls.setToolTip(
            "Number of strongest distinct automatic-calculator reflection positions to display. "
            "The limit is applied after q-frame and NPZ-mask filtering."
        )
        self.calculator_max_hkls.valueChanged.connect(lambda _=None: self._redraw_manual())
        self.calculator_hide_masked = QCheckBox(
            "Hide automatic-calculator reflections in masked/invalid NPZ regions"
        )
        self.calculator_hide_masked.setChecked(True)
        self.calculator_hide_masked.setToolTip(
            "Uses the selected reference NPZ validity mask before selecting strongest points. "
            "HKL text remains controlled only by 'Show HKL labels'."
        )
        self.calculator_hide_masked.stateChanged.connect(lambda _=None: self._redraw_manual())
        display_form.addRow(self.show_calculator_points)
        display_form.addRow(self.calculator_strongest_only)
        display_form.addRow("Number of strongest automatic-calculator points:", self.calculator_max_hkls)
        display_form.addRow(self.calculator_hide_masked)
        layout.addWidget(display_group)

        # Advanced accuracy controls
        advanced_group = QGroupBox("Advanced accuracy settings")
        advanced_form = QFormLayout(advanced_group)
        self.click_subpixel_fit = QCheckBox(
            "Use subpixel centroid and background-integrated experimental intensity"
        )
        self.click_subpixel_fit.setChecked(True)
        self.click_peak_radius = self._integer(4, 1, 20)
        self.click_background_radius = self._integer(9, 3, 40)
        self.click_peak_radius.setToolTip("Seed radius for the automatic adaptive peak-integration model.")
        self.click_background_radius.setToolTip("Minimum local radius used by the robust background model.")
        self.click_require_explicit_axes = QCheckBox("Require explicit qr and qz axes inside each NPZ")
        self.click_require_explicit_axes.setChecked(True)
        self.auto_expand_q_max = QCheckBox("Automatically expand q max to cover the entire NPZ frame")
        self.auto_expand_q_max.setChecked(True)
        self.high_q_corner_recovery = QCheckBox("Recover weak high-q edge/corner peaks")
        self.high_q_corner_recovery.setChecked(True)
        self.input_provenance_label = QLabel("No reference NPZ selected.")
        self.input_provenance_label.setWordWrap(True)
        advanced_form.addRow(self.click_subpixel_fit)
        advanced_form.addRow("Peak integration radius (pixels):", self.click_peak_radius)
        advanced_form.addRow("Background radius (pixels):", self.click_background_radius)
        advanced_form.addRow(self.click_require_explicit_axes)
        advanced_form.addRow(self.auto_expand_q_max)
        advanced_form.addRow(self.high_q_corner_recovery)
        # Retain provenance metadata internally for auditing while keeping the
        # manual-click interface focused on quantities needed for peak comparison.
        self.input_provenance_label.setVisible(False)
        layout.addWidget(advanced_group)

        audit_group = QGroupBox("Accuracy and validation summary")
        audit_layout = QVBoxLayout(audit_group)
        self.accuracy_summary_label = QLabel(
            "Run the automatic calculator and add manual peak pairs to populate the audit."
        )
        self.accuracy_summary_label.setWordWrap(True)
        audit_layout.addWidget(self.accuracy_summary_label)
        self.export_accuracy_audit_button = QPushButton("Export accuracy audit…")
        self.export_accuracy_audit_button.clicked.connect(self._export_accuracy_audit)
        audit_layout.addWidget(self.export_accuracy_audit_button)
        self.accuracy_audit_group = audit_group
        layout.addWidget(audit_group)
        audit_group.setVisible(False)

        self.click_pending_label = QLabel(
            "No pending point. First click an experimental peak on the PNG."
        )
        self.click_pending_label.setWordWrap(True)
        layout.addWidget(self.click_pending_label)

        candidate_row = QHBoxLayout()
        self.click_candidate_combo = QComboBox()
        self.click_candidate_combo.setEnabled(False)
        assign_button = QPushButton("Assign selected HKL")
        assign_button.clicked.connect(self._assign_selected_calculator_candidate)
        candidate_row.addWidget(self.click_candidate_combo, 1)
        candidate_row.addWidget(assign_button)
        layout.addLayout(candidate_row)

        actions = QHBoxLayout()
        cancel_button = QPushButton("Cancel pending")
        cancel_button.clicked.connect(self._cancel_pending_click)
        undo_button = QPushButton("Undo last")
        undo_button.clicked.connect(self._undo_manual_assignment)
        clear_button = QPushButton("Clear assignments")
        clear_button.clicked.connect(self._clear_manual_assignments)
        actions.addWidget(cancel_button)
        actions.addWidget(undo_button)
        actions.addWidget(clear_button)
        layout.addLayout(actions)

        transfer = QHBoxLayout()
        import_button = QPushButton("Import assignment CSV…")
        import_button.clicked.connect(self._import_manual_assignments)
        export_button = QPushButton("Export assignment CSV…")
        export_button.clicked.connect(self._export_full_index_table)
        transfer.addWidget(import_button)
        transfer.addWidget(export_button)
        layout.addLayout(transfer)

        # Place immediately after the automatic calculator group.
        engine_index = -1
        for index in range(controls_layout.count()):
            widget = controls_layout.itemAt(index).widget()
            if isinstance(widget, QGroupBox) and widget.title() == "Full automatic indexing calculations":
                engine_index = index
                break
        controls_layout.insertWidget(engine_index + 1 if engine_index >= 0 else controls_layout.count() - 1, group)

        self.manual_canvas.mpl_connect("button_press_event", self._on_manual_plot_click)
        self.measurement_table.cellChanged.connect(lambda _row, _col: self._refresh_reference_npz_combo())
        self._refresh_reference_npz_combo()
        self._fill_manual_assignment_table()
        self._update_spacegroup_consistency()
        self._update_accuracy_summary()

    # ---------------------------- robust NPZ input ----------------------------
    def _npz_probe_config(self, qr_min=None, qr_max=None, qz_min=None, qz_max=None):
        config_class = globals().get("V9Config", globals().get("V8Config", IndexingConfig))
        return config_class(
            cif_path=self.manual_cif_edit.text().strip() or self.cif_edit.text().strip() or "placeholder.cif",
            qr_range=(self.manual_qr_min.value() if qr_min is None else float(qr_min),
                      self.manual_qr_max.value() if qr_max is None else float(qr_max)),
            qz_range=(self.manual_qz_min.value() if qz_min is None else float(qz_min),
                      self.manual_qz_max.value() if qz_max is None else float(qz_max)),
            prefer_numerical=True,
        )

    def _probe_npz(self, path, measurement=None):
        if measurement is None:
            config = self._npz_probe_config()
        else:
            config = self._npz_probe_config(
                measurement.qr_min, measurement.qr_max, measurement.qz_min, measurement.qz_max
            )
        return load_numerical_qspace(str(path), config)

    def _enforce_explicit_npz_axes(self, path, image):
        require = getattr(self, "click_require_explicit_axes", None)
        if require is not None and require.isChecked() and not bool(image.get("explicit_q_axes", False)):
            raise ValueError(
                f"{_GuiPath(path).name} does not contain usable explicit qr and qz axes. "
                "For high-accuracy indexing, save intensity together with qr and qz in the NPZ. "
                "Uncheck 'Require explicit qr and qz axes inside each NPZ' only for an approximate run."
            )

    def _add_npz_calculation_files(self):
        paths, _ = QFileDialog.getOpenFileNames(
            self, "Select numerical NPZ measurements", "",
            "NumPy q-space archives (*.npz);;All files (*)",
        )
        if not paths:
            return
        failures = []
        added = 0
        for path in paths:
            try:
                measurement = parse_measurement_filename(path, self.measurement_table.rowCount() + 1)
                measurement.file = str(_GuiPath(path).resolve())
                measurement.numerical_file = measurement.file
                measurement.png_file = ""
                measurement.qr_min = self.manual_qr_min.value()
                measurement.qr_max = self.manual_qr_max.value()
                measurement.qz_min = self.manual_qz_min.value()
                measurement.qz_max = self.manual_qz_max.value()
                image = self._probe_npz(path, measurement)
                self._enforce_explicit_npz_axes(path, image)
                measurement.qr_min = float(np.nanmin(image["qr"]))
                measurement.qr_max = float(np.nanmax(image["qr"]))
                measurement.qz_min = float(np.nanmin(image["qz"]))
                measurement.qz_max = float(np.nanmax(image["qz"]))
                self._reference_npz_cache[str(_GuiPath(path).resolve())] = image
                self._append_measurement(measurement)
                added += 1
            except Exception as exc:
                failures.append(f"{_GuiPath(path).name}: {exc}")
        self._refresh_reference_npz_combo()
        if added:
            self._auto_select_reference_npz_for_png()
            self._update_qmax_from_loaded_npz()
            self.manual_status.setText(
                f"Added {added} NPZ calculation file(s). Edit incident angles if the filenames did not contain them. "
                f"Current radial q max: {self.manual_q_max.value():.3f} Å⁻¹."
            )
        if failures:
            QMessageBox.warning(
                self, "Some NPZ files could not be loaded", "\n\n".join(failures)
            )

    def _remove_npz_calculation_files(self):
        self._remove_measurements()
        self._refresh_reference_npz_combo()

    def _refresh_reference_npz_combo(self):
        if not hasattr(self, "click_reference_npz_combo"):
            return
        previous = self.click_reference_npz_combo.currentData()
        paths = []
        for row in range(self.measurement_table.rowCount()):
            item = self.measurement_table.item(row, 1)
            path = item.text().strip() if item else ""
            if path and _GuiPath(path).suffix.lower() == ".npz":
                resolved = str(_GuiPath(path).expanduser().resolve())
                if resolved not in paths:
                    paths.append(resolved)
        self.click_reference_npz_combo.blockSignals(True)
        self.click_reference_npz_combo.clear()
        for path in paths:
            self.click_reference_npz_combo.addItem(_GuiPath(path).name, path)
        if previous in paths:
            self.click_reference_npz_combo.setCurrentIndex(paths.index(previous))
        self.click_reference_npz_combo.blockSignals(False)
        if paths:
            self._reference_npz_changed()

    @staticmethod
    def _normalized_overlay_measurement_stem(path):
        """Normalize only known converter suffixes for exact PNG↔NPZ pairing."""
        stem = _GuiPath(str(path)).stem.lower().strip()
        for suffix in ("_with_q_axes", "_axes", "_converted"):
            if stem.endswith(suffix):
                stem = stem[:-len(suffix)]
                break
        return stem

    def _auto_select_reference_npz_for_png(self):
        """Select the unique NPZ whose measurement stem exactly matches the display PNG.

        Exact matching prevents intensity from being measured from a different
        scan when several NPZ files are loaded. If there is not one unique exact
        stem match, the current reference selection is left unchanged.
        """
        if not getattr(self, "manual_image_path", None) or not hasattr(self, "click_reference_npz_combo"):
            return False
        png_stem = self._normalized_overlay_measurement_stem(self.manual_image_path)
        matches = []
        for index in range(self.click_reference_npz_combo.count()):
            path = self.click_reference_npz_combo.itemData(index)
            if path and self._normalized_overlay_measurement_stem(path) == png_stem:
                matches.append(index)
        if len(matches) != 1:
            return False
        target = matches[0]
        if self.click_reference_npz_combo.currentIndex() != target:
            self.click_reference_npz_combo.setCurrentIndex(target)
        else:
            self._reference_npz_changed()
        return True

    def _reference_npz_changed(self):
        if not hasattr(self, "click_reference_npz_combo"):
            return
        path = self.click_reference_npz_combo.currentData()
        if not path:
            return
        try:
            image = self._reference_npz_cache.get(path)
            if image is None:
                image = self._probe_npz(path)
                self._reference_npz_cache[path] = image
            if self.click_use_npz_axes.isChecked():
                self.manual_qr_min.setValue(float(np.nanmin(image["qr"])))
                self.manual_qr_max.setValue(float(np.nanmax(image["qr"])))
                self.manual_qz_min.setValue(float(np.nanmin(image["qz"])))
                self.manual_qz_max.setValue(float(np.nanmax(image["qz"])))
            self._update_qmax_from_loaded_npz()
            note = image.get("numerical_loader_note", "recognized numerical intensity")
            provenance = image.get("input_provenance", self._npz_provenance_for_path(path))
            normalization_note = str(image.get("intensity_normalization_note", "no explicit normalization metadata"))
            self.manual_status.setText(
                f"Reference NPZ: {_GuiPath(path).name} | shape {np.asarray(image['raw_intensity']).shape} | {note}. "
                f"Radial q max: {self.manual_q_max.value():.3f} Å⁻¹. {provenance} "
                f"Intensity: {normalization_note}."
            )
            self._update_accuracy_summary()
        except Exception as exc:
            QMessageBox.critical(self, "Reference NPZ error", str(exc))

    # ------------------------ automatic calculator request ---------------------
    def _required_qmax_for_npz_frame(self, measurements=None):
        """Return a radial q limit that covers every enabled NPZ reciprocal-space frame.

        The radial magnitude is ``sqrt(qr**2 + qz**2)``, so the required limit is
        determined by the furthest corner of the q_r/q_z rectangle rather than by
        the largest value on either individual axis.
        """
        if measurements is None:
            measurements = self._measurements_from_table()
        required = float(self.manual_q_max.value())
        for measurement in measurements:
            if not bool(getattr(measurement, "enabled", True)):
                continue
            path = measurement.numerical_file or measurement.file
            image = self._reference_npz_cache.get(str(_GuiPath(path).expanduser().resolve()))
            if image is None:
                try:
                    image = self._probe_npz(path, measurement)
                    self._reference_npz_cache[str(_GuiPath(path).expanduser().resolve())] = image
                except Exception:
                    image = None
            if image is not None:
                qr = np.asarray(image.get("qr", []), dtype=float)
                qz = np.asarray(image.get("qz", []), dtype=float)
                if qr.size and qz.size:
                    max_qr = float(np.nanmax(np.abs(qr)))
                    max_qz = float(np.nanmax(np.abs(qz)))
                    required = max(required, math.hypot(max_qr, max_qz) + 0.08)
                    continue
            max_qr = max(abs(float(measurement.qr_min)), abs(float(measurement.qr_max)))
            max_qz = max(abs(float(measurement.qz_min)), abs(float(measurement.qz_max)))
            required = max(required, math.hypot(max_qr, max_qz) + 0.08)
        return float(required)

    def _update_qmax_from_loaded_npz(self):
        if not getattr(self, "auto_expand_q_max", None) or not self.auto_expand_q_max.isChecked():
            return
        try:
            required = self._required_qmax_for_npz_frame()
            if required > self.manual_q_max.value() + 1e-6:
                self.manual_q_max.setValue(required)
                self.manual_status.setText(
                    f"Calculation q max expanded to {required:.3f} Å⁻¹ so the full NPZ frame, including its corners, is covered."
                )
        except Exception:
            pass

    def _build_request(self, preview=False):
        measurements = self._measurements_from_table()
        for measurement in measurements:
            path = measurement.numerical_file or measurement.file
            if _GuiPath(path).suffix.lower() != ".npz":
                raise ValueError(
                    f"Automatic calculator input must be NPZ, not {_GuiPath(path).name}. "
                    "Load the PNG only in the Experimental image field above."
                )
            image = self._reference_npz_cache.get(str(_GuiPath(path).expanduser().resolve()))
            if image is None:
                image = self._probe_npz(path, measurement)
                self._reference_npz_cache[str(_GuiPath(path).expanduser().resolve())] = image
            self._enforce_explicit_npz_axes(path, image)
            measurement.file = path
            measurement.numerical_file = path
            measurement.png_file = ""
            measurement.colormap = "jet"
        if not measurements:
            raise ValueError(
                "Add at least one NPZ calculation file. The overlay PNG is display-only and is not an automatic-calculator input."
            )
        q_max = float(self.manual_q_max.value())
        if getattr(self, "auto_expand_q_max", None) is not None and self.auto_expand_q_max.isChecked():
            q_max = self._required_qmax_for_npz_frame(measurements)
            self.manual_q_max.setValue(q_max)
        return GIXSRunRequest(
            cif_path=self.manual_cif_edit.text().strip() or self.cif_edit.text().strip(),
            measurements=measurements,
            output_dir=self.manual_output_edit.text().strip(),
            workflow_preset="preview" if preview else self._workflow_key_manual(),
            alternative_cif_paths=list(self.alternative_cifs),
            q_max=q_max,
            colormap="jet",
            preview_only=preview,
            prefer_numerical=True,
            run_in_parallel=True,
            auto_expand_q_max=(
                getattr(self, "auto_expand_q_max", None) is not None
                and self.auto_expand_q_max.isChecked()
            ),
            high_q_corner_recovery=(
                getattr(self, "high_q_corner_recovery", None) is not None
                and self.high_q_corner_recovery.isChecked()
            ),
        )

    def _start_full_indexing(self):
        self.calculator_predictions = pd.DataFrame()
        self._manual_intensity_physics_cache.clear()
        self._q_refinement = None
        if getattr(self, "refinement_status_label", None) is not None:
            self.refinement_status_label.setText("No manual q-space refinement is active.")
        self.pending_clicked_experimental = None
        self._candidate_prediction_indices = []
        self.click_candidate_combo.clear()
        self.click_candidate_combo.setEnabled(False)
        super()._start_full_indexing()

    # --------------------- calculator output and predictions -------------------
    @staticmethod
    def _flex_column(frame, names):
        if frame is None or frame.empty:
            return None
        lower = {str(column).lower(): column for column in frame.columns}
        for name in names:
            if name in frame.columns:
                return name
            if str(name).lower() in lower:
                return lower[str(name).lower()]
        return None

    def _prediction_rows_from_frame(self, frame, source):
        if frame is None or frame.empty:
            return []
        qr_col = self._flex_column(frame, ["qr_calc", "qr", "predicted_qr", "qr_pred"])
        qz_col = self._flex_column(frame, ["qz_calc", "qz", "predicted_qz", "qz_pred"])
        if qr_col is None or qz_col is None:
            return []
        hkl_col = self._flex_column(frame, ["hkl", "overlay_hkl_text", "HKL"])
        intensity_col = self._flex_column(
            frame, ["effective_intensity", "calculated_intensity", "f2", "prediction_weight", "CalcIntensity"]
        )
        f2_col = self._flex_column(frame, ["f2", "F2", "structure_factor_intensity"])
        dwba_col = self._flex_column(frame, ["dwba_weight", "DWBAWeight"])
        intensity_model_col = self._flex_column(frame, ["intensity_model", "IntensityModel"])
        domain_col = self._flex_column(frame, ["orientation_domain", "domain"])
        evidence_col = self._flex_column(frame, ["index_source", "assignment_source", "role"])
        detector_col = self._flex_column(
            frame, ["detection_source", "source_detector", "detection_source_mix"]
        )
        rows = []
        for index, row in frame.reset_index(drop=True).iterrows():
            try:
                qr = float(row[qr_col])
                qz = float(row[qz_col])
            except Exception:
                continue
            if not np.isfinite(qr) or not np.isfinite(qz):
                continue
            if hkl_col is not None and str(row[hkl_col]).strip() not in ("", "nan"):
                hkl = str(row[hkl_col])
            elif all(name in frame.columns for name in ("h", "k", "l")):
                hkl = f"({int(row['h'])} {int(row['k'])} {int(row['l'])})"
            else:
                hkl = f"reflection {index + 1}"
            try:
                calc_i = float(row[intensity_col]) if intensity_col is not None else np.nan
            except Exception:
                calc_i = np.nan
            try:
                calc_f2 = float(row[f2_col]) if f2_col is not None else calc_i
            except Exception:
                calc_f2 = calc_i
            try:
                dwba_weight = float(row[dwba_col]) if dwba_col is not None else 1.0
            except Exception:
                dwba_weight = 1.0
            if intensity_model_col is not None and str(row[intensity_model_col]) not in ("", "nan"):
                intensity_model = str(row[intensity_model_col])
            else:
                intensity_model = "DWBA-weighted F²" if np.isfinite(dwba_weight) and abs(dwba_weight - 1.0) > 1e-8 else "kinematic F²"
            domain = str(row[domain_col]) if domain_col is not None else "primary"
            stability_col = self._flex_column(
                frame, ["empirical_assignment_stability", "assignment_stability"]
            )
            tier_col = self._flex_column(frame, ["stability_tier", "assignment_stability_tier"])
            try:
                stability = float(row[stability_col]) if stability_col is not None else np.nan
            except Exception:
                stability = np.nan
            tier = str(row[tier_col]) if tier_col is not None and str(row[tier_col]) != "nan" else "unassessed"
            evidence_source = (
                str(row[evidence_col])
                if evidence_col is not None and str(row[evidence_col]) not in ("", "nan")
                else ("unused_calculated_prediction" if "unused" in source else "unreported")
            )
            detection_source = (
                str(row[detector_col])
                if detector_col is not None and str(row[detector_col]) not in ("", "nan")
                else ("not_applicable" if "unused" in source else "unreported")
            )
            rows.append({
                "HKL": hkl,
                "QrCalc": qr,
                "QzCalc": qz,
                "CalcIntensity": calc_i,
                "CalcF2": calc_f2,
                "DWBAWeight": dwba_weight,
                "CalcIntensityModel": intensity_model,
                "OrientationDomain": domain,
                "CalculatorSource": source,
                "IndexingEvidenceSource": evidence_source,
                "DetectionSource": detection_source,
                "AssignmentStability": stability,
                "StabilityTier": tier,
            })
        return rows

    def _build_calculator_prediction_inventory(self, item):
        rows = self._prediction_rows_from_frame(self.backend_indexed, "automatic indexed solution")
        series_dir = _GuiPath(item.overlay_path).parent
        unused_path = series_dir / "unused_calculated_reflections.csv"
        if unused_path.is_file():
            try:
                unused = pd.read_csv(unused_path)
                rows.extend(self._prediction_rows_from_frame(unused, "automatic unused prediction"))
            except Exception as exc:
                self.manual_progress_log.appendPlainText(
                    f"Could not read unused calculated reflections: {exc}"
                )
        if not rows:
            return pd.DataFrame(columns=[
                "CalcID", "HKL", "QrCalc", "QzCalc", "CalcIntensity",
                "OrientationDomain", "CalculatorSource", "IndexingEvidenceSource", "DetectionSource",
            ])
        frame = pd.DataFrame(rows)
        frame["_key"] = (
            frame.HKL.astype(str) + "|" + frame.OrientationDomain.astype(str) + "|" +
            frame.QrCalc.round(6).astype(str) + "|" + frame.QzCalc.round(6).astype(str)
        )
        frame = frame.drop_duplicates("_key", keep="first").drop(columns="_key")
        frame = frame.sort_values(
            ["CalcIntensity", "QzCalc", "QrCalc"], ascending=[False, True, True],
            na_position="last", kind="mergesort",
        ).reset_index(drop=True)
        frame.insert(0, "CalcID", np.arange(1, len(frame) + 1, dtype=int))

        stability_path = _GuiPath(getattr(item, "stability_table_path", ""))
        if stability_path.is_file():
            try:
                stability = pd.read_csv(stability_path)
                if not stability.empty and "hkl" in stability.columns:
                    stability = stability.copy()
                    stability["_hkl_key"] = stability["hkl"].astype(str).str.findall(
                        r"[-+]?\d+"
                    ).str.join(",")
                    domain_name = self._flex_column(stability, ["orientation_domain", "domain"])
                    stability["_domain_key"] = (
                        stability[domain_name].astype(str) if domain_name is not None else "primary"
                    )
                    frame["_hkl_key"] = frame["HKL"].astype(str).str.findall(
                        r"[-+]?\d+"
                    ).str.join(",")
                    frame["_domain_key"] = frame["OrientationDomain"].astype(str)
                    keep = ["_hkl_key", "_domain_key"]
                    for column in ("empirical_assignment_stability", "stability_tier"):
                        if column in stability.columns:
                            keep.append(column)
                    frame = frame.merge(
                        stability[keep].drop_duplicates(["_hkl_key", "_domain_key"]),
                        on=["_hkl_key", "_domain_key"], how="left", suffixes=("", "_validated")
                    )
                    if "empirical_assignment_stability" in frame:
                        frame["AssignmentStability"] = pd.to_numeric(
                            frame["empirical_assignment_stability"], errors="coerce"
                        ).combine_first(pd.to_numeric(frame.get("AssignmentStability"), errors="coerce"))
                    if "stability_tier" in frame:
                        validated_tier = frame["stability_tier"].astype(str)
                        frame["StabilityTier"] = np.where(
                            validated_tier.ne("nan") & validated_tier.ne(""),
                            validated_tier, frame.get("StabilityTier", "unassessed")
                        )
                    frame = frame.drop(columns=[
                        column for column in (
                            "_hkl_key", "_domain_key", "empirical_assignment_stability", "stability_tier"
                        ) if column in frame.columns
                    ])
            except Exception as exc:
                self.manual_progress_log.appendPlainText(
                    f"Could not merge reflection-stability results: {exc}"
                )
        frame["QrCalcOriginal"] = pd.to_numeric(frame["QrCalc"], errors="coerce")
        frame["QzCalcOriginal"] = pd.to_numeric(frame["QzCalc"], errors="coerce")
        frame["BackendSupportScore"] = [
            self._prediction_support_score_values(source, evidence, stability)
            for source, evidence, stability in zip(
                frame.get("CalculatorSource", pd.Series("", index=frame.index)),
                frame.get("IndexingEvidenceSource", pd.Series("", index=frame.index)),
                frame.get("AssignmentStability", pd.Series(np.nan, index=frame.index)),
            )
        ]
        return frame

    def _reference_specific_giwaxs_physics_config(self, reference_npz):
        """Return cached AUTO physics with the clicked measurement's incidence angle."""
        if not reference_npz:
            return self._current_auto_giwaxs_physics_config()
        try:
            cache_key = str(_GuiPath(reference_npz).expanduser().resolve())
        except Exception:
            cache_key = str(reference_npz)
        if cache_key in self._manual_intensity_physics_cache:
            return self._manual_intensity_physics_cache[cache_key]
        base = self._current_auto_giwaxs_physics_config()
        if base is None:
            self._manual_intensity_physics_cache[cache_key] = None
            return None
        angle = None
        try:
            metadata = _npz_auto_physics_metadata(reference_npz)
            value = metadata.get("incidence_angle_deg")
            if value is not None and np.isfinite(float(value)) and float(value) > 0:
                angle = float(value)
        except Exception:
            pass
        if angle is None:
            try:
                parsed = _parse_measurement_name(_GuiPath(reference_npz))
                if parsed is not None:
                    value = float(parsed.get("angle_deg", np.nan))
                    if np.isfinite(value) and value > 0:
                        angle = value
            except Exception:
                pass
        if angle is None:
            self._manual_intensity_physics_cache[cache_key] = base
            return base
        configured = replace(
            base,
            incidence_angle_deg=float(angle),
            giwaxs_incidence_angles_deg=(float(angle),),
        )
        self._manual_intensity_physics_cache[cache_key] = configured
        return configured

    def _calculated_unresolved_intensity(self, prediction_index, experimental_point):
        """Return measurement-specific single/cluster calculated intensity.

        When AUTO DWBA is physically supported, the clicked reference NPZ's own
        incident angle is used instead of the series-mean DWBA envelope.  Nearby
        calculated reflections that fall inside the experimental resolution are
        summed because the detector observes them as one unresolved spot.
        """
        if self.calculator_predictions is None or self.calculator_predictions.empty:
            return np.nan, np.nan, 1, np.nan, "", np.nan, 1.0, "unreported"
        predictions = self.calculator_predictions.reset_index(drop=True)
        prediction_index = int(prediction_index)
        center = predictions.iloc[prediction_index]
        path = experimental_point.get("ReferenceNPZ", "")
        dqr = dqz = np.nan
        try:
            resolved = str(_GuiPath(path).expanduser().resolve())
            image = self._reference_npz_cache.get(resolved)
            if image is None:
                image = self._probe_npz(path)
                self._reference_npz_cache[resolved] = image
            qr_axis = np.asarray(image["qr"], float)
            qz_axis = np.asarray(image["qz"], float)
            dqr = float(np.nanmedian(np.abs(np.diff(qr_axis)))) if len(qr_axis) > 1 else np.nan
            dqz = float(np.nanmedian(np.abs(np.diff(qz_axis)))) if len(qz_axis) > 1 else np.nan
        except Exception:
            pass
        width_r = float(experimental_point.get("PeakSigmaQrWidth", np.nan))
        width_z = float(experimental_point.get("PeakSigmaQzWidth", np.nan))
        axis_r = max(0.004, 2.0 * dqr if np.isfinite(dqr) else 0.0, 1.5 * width_r if np.isfinite(width_r) else 0.0)
        axis_z = max(0.004, 2.0 * dqz if np.isfinite(dqz) else 0.0, 1.5 * width_z if np.isfinite(width_z) else 0.0)
        axis_r = min(axis_r, 0.025)
        axis_z = min(axis_z, 0.025)
        qr = pd.to_numeric(predictions["QrCalc"], errors="coerce").to_numpy(float)
        qz = pd.to_numeric(predictions["QzCalc"], errors="coerce").to_numpy(float)
        metric = np.sqrt(((qr - float(center.QrCalc)) / axis_r) ** 2 + ((qz - float(center.QzCalc)) / axis_z) ** 2)
        members = np.isfinite(metric) & (metric <= 1.0)

        f2 = pd.to_numeric(
            predictions.get("CalcF2", predictions.get("CalcIntensity")), errors="coerce"
        ).to_numpy(float)
        backend_effective = pd.to_numeric(predictions.get("CalcIntensity"), errors="coerce").to_numpy(float)
        physics = self._reference_specific_giwaxs_physics_config(path)
        weights = np.ones(len(predictions), dtype=float)
        model = "kinematic F²"
        if physics is not None and bool(getattr(physics, "enable_dwba", False)):
            qz_physics = pd.to_numeric(
                predictions.get("QzCalcOriginal", predictions.get("QzCalc")), errors="coerce"
            ).to_numpy(float)
            try:
                weights = dwba_intensity_envelope(qz_physics, physics)
                values = f2 * weights
                model = f"single-angle DWBA-weighted F² (αi={float(physics.incidence_angle_deg):.5g}°)"
            except Exception:
                values = backend_effective
                weights = pd.to_numeric(
                    predictions.get("DWBAWeight", pd.Series(1.0, index=predictions.index)), errors="coerce"
                ).fillna(1.0).to_numpy(float)
                model = str(getattr(center, "CalcIntensityModel", "backend effective intensity"))
        else:
            # If AUTO physics is not available, use F² when present.  This avoids
            # accidentally treating a series-mean DWBA value as a single-angle
            # prediction when the per-measurement physics cannot be reconstructed.
            values = np.where(np.isfinite(f2), f2, backend_effective)
            model = "kinematic F²"
        positive = members & np.isfinite(values) & (values > 0)
        selected_single = float(values[prediction_index]) if np.isfinite(values[prediction_index]) else float(center.CalcIntensity)
        selected_f2 = float(f2[prediction_index]) if np.isfinite(f2[prediction_index]) else np.nan
        selected_weight = float(weights[prediction_index]) if np.isfinite(weights[prediction_index]) else 1.0
        if not positive.any():
            return (
                selected_single, selected_single, 1, float(math.hypot(axis_r, axis_z)),
                str(center.HKL), selected_f2, selected_weight, model,
            )
        total = float(np.sum(values[positive]))
        hkls = ", ".join(predictions.loc[positive, "HKL"].astype(str).tolist())
        return (
            selected_single, total, int(positive.sum()), float(math.hypot(axis_r, axis_z)),
            hkls, selected_f2, selected_weight, model,
        )

    def _load_full_series(self, series_id):
        # Parent loads the complete automatic result and validation metadata.
        super()._load_full_series(series_id)
        if self.backend_series_item is None:
            return
        self._q_refinement = None
        if getattr(self, "refinement_status_label", None) is not None:
            self.refinement_status_label.setText("No manual q-space refinement is active.")
        self.calculator_predictions = self._build_calculator_prediction_inventory(
            self.backend_series_item
        )
        self.manual_click_assignments = self.manual_click_assignments_by_series.get(
            series_id, pd.DataFrame(columns=self._ASSIGNMENT_COLUMNS)
        ).copy()
        # Loading a series restores the backend's unrefined coordinates. Manual
        # q-map refinement is intentionally session/series-local and must be run
        # again if desired after switching series.
        if not self.manual_click_assignments.empty and "_CalcID" in self.manual_click_assignments.columns:
            lookup = self.calculator_predictions.set_index("CalcID")
            for row_index, calc_id in enumerate(
                pd.to_numeric(self.manual_click_assignments["_CalcID"], errors="coerce").to_numpy(float)
            ):
                if np.isfinite(calc_id) and int(calc_id) in lookup.index:
                    pred = lookup.loc[int(calc_id)]
                    self.manual_click_assignments.loc[row_index, "QrCalc"] = float(pred["QrCalc"])
                    self.manual_click_assignments.loc[row_index, "QzCalc"] = float(pred["QzCalc"])
                    self.manual_click_assignments.loc[row_index, "QrCalcOriginal"] = float(pred["QrCalcOriginal"])
                    self.manual_click_assignments.loc[row_index, "QzCalcOriginal"] = float(pred["QzCalcOriginal"])
        self.pending_clicked_experimental = None
        self._candidate_prediction_indices = []
        self.click_candidate_combo.clear()
        self.click_candidate_combo.setEnabled(False)
        self._recalculate_manual_assignment_results()
        self.manual_decision_label.setText(
            f"Automatic calculator decision: {self.backend_series_item.final_decision or 'not reported'} | "
            f"orientation: {self.backend_series_item.orientation_hkl or 'not reported'} | "
            f"calculated visible reflections: {len(self.calculator_predictions)}. "
            "Experimental comparison rows are added only by your clicks."
        )
        self.manual_status.setText(
            f"Automatic calculator loaded for {series_id}. Click a PNG experimental peak, then its calculated circle."
        )
        self._update_accuracy_summary()
        self._redraw_manual()

    # -------------------------- click-data extraction --------------------------
    def _reference_npz_data(self):
        path = self.click_reference_npz_combo.currentData()
        if not path:
            raise ValueError("Select a reference NPZ file for clicked q/intensity values.")
        image = self._reference_npz_cache.get(path)
        if image is None:
            image = self._probe_npz(path)
            self._reference_npz_cache[path] = image
        return path, image

    @staticmethod
    def _nearest_axis(axis, value):
        return int(np.argmin(np.abs(np.asarray(axis, dtype=float) - float(value))))

    def _refine_manual_click_position(self, image, row, column):
        """Return a stable subpixel center for a manually selected experimental peak.

        A click can fall near the edge of a diffraction spot. The local NPZ signal
        is smoothed to identify a credible nearby maximum, followed by a compact
        intensity-weighted centroid around that maximum. This refines only the
        recorded experimental q coordinate; robust intensity integration remains
        a separate calculation.
        """
        intensity = np.asarray(
            image.get("quantitative_intensity", image.get("raw_intensity", image["intensity"])),
            dtype=float,
        )
        valid = np.asarray(image.get("valid", np.isfinite(intensity)), dtype=bool)
        height, width = intensity.shape
        row = int(np.clip(int(row), 0, height - 1))
        column = int(np.clip(int(column), 0, width - 1))

        radius = max(int(self.click_snap_radius.value()), 1)
        r0, r1 = max(0, row - radius), min(height, row + radius + 1)
        c0, c1 = max(0, column - radius), min(width, column + radius + 1)
        patch = intensity[r0:r1, c0:c1].astype(float, copy=True)
        patch_valid = valid[r0:r1, c0:c1] & np.isfinite(patch)
        if not patch_valid.any():
            return float(row), float(column)

        # Fill invalid pixels with the local median only for smoothing; invalid
        # pixels are never allowed to become peak candidates or centroid weight.
        local_values = patch[patch_valid]
        local_median = float(np.nanmedian(local_values))
        work = np.where(patch_valid, patch, local_median)
        smooth = gaussian_filter(work, sigma=1.0, mode="nearest")

        # Estimate local contrast/noise so tiny pixel-scale bumps do not steal a
        # click from the real peak center in crowded or textured regions.
        smooth_values = smooth[patch_valid]
        background = float(np.nanmedian(smooth_values))
        residual_values = smooth_values - background
        mad = float(np.nanmedian(np.abs(residual_values - np.nanmedian(residual_values))))
        noise = max(1.4826 * mad, float(np.nanstd(residual_values)) * 0.25, 1e-12)
        signal = smooth - background
        max_signal = float(np.nanmax(signal[patch_valid]))
        if not np.isfinite(max_signal) or max_signal <= 0:
            return float(row), float(column)

        yy, xx = np.indices(patch.shape, dtype=float)
        click_y = float(row - r0)
        click_x = float(column - c0)
        distance = np.hypot(yy - click_y, xx - click_x)

        local_maxima = (
            patch_valid
            & (smooth == maximum_filter(np.where(patch_valid, smooth, -np.inf), size=3, mode="nearest"))
            & (signal >= max(2.0 * noise, 0.18 * max_signal))
            & (distance <= float(radius))
        )

        if local_maxima.any():
            my, mx = np.where(local_maxima)
            strengths = np.maximum(signal[my, mx], 0.0) / max(max_signal, 1e-12)
            distances = np.hypot(my - click_y, mx - click_x) / max(float(radius), 1.0)
        # Prefer a credible local intensity maximum near the click while limiting
        # displacement toward a brighter but more distant feature. This keeps the
        # refined coordinate representative of the peak selected by the analyst.
            score = strengths - 0.32 * distances
            best = int(np.argmax(score))
            center_y = float(my[best])
            center_x = float(mx[best])
        else:
            candidate = np.where(patch_valid, signal, -np.inf)
            flat = int(np.nanargmax(candidate))
            center_y, center_x = map(float, np.unravel_index(flat, patch.shape))

        # Compact subpixel centroid around the chosen maximum.  Keeping this
        # aperture small makes the marker center stable on broad peaks/streaks.
        core_radius = max(2.0, min(4.5, radius * 0.45))
        core_distance = np.hypot(yy - center_y, xx - center_x)
        core = patch_valid & (core_distance <= core_radius)
        core_values = smooth[core]
        if core_values.size:
            core_background = float(np.nanpercentile(core_values, 25.0))
            weights = np.where(core, np.maximum(smooth - core_background, 0.0), 0.0)
            # Emphasize the peak core so a sloping tail/background does not pull
            # the centroid toward the edge of the spot.
            weights *= np.exp(-0.5 * (core_distance / max(core_radius * 0.75, 1.0)) ** 2)
            total = float(np.sum(weights))
            if total > 0:
                sub_y = float(np.sum((yy + r0) * weights) / total)
                sub_x = float(np.sum((xx + c0) * weights) / total)
            else:
                sub_y, sub_x = center_y + r0, center_x + c0
        else:
            sub_y, sub_x = center_y + r0, center_x + c0

        # Constrain subpixel refinement to the local neighborhood around the
        # selected feature so the refined coordinate cannot migrate to another peak.
        delta_y = sub_y - float(row)
        delta_x = sub_x - float(column)
        shift = math.hypot(delta_x, delta_y)
        if shift > float(radius) and shift > 0:
            scale = float(radius) / shift
            sub_y = float(row) + delta_y * scale
            sub_x = float(column) + delta_x * scale
        return sub_y, sub_x

    def _experimental_from_png_click(self, qr_click, qz_click):
        path, image = self._reference_npz_data()
        self._enforce_explicit_npz_axes(path, image)
        qr = np.asarray(image["qr"], dtype=float)
        qz = np.asarray(image["qz"], dtype=float)
        intensity = np.asarray(image.get("quantitative_intensity", image.get("raw_intensity", image["intensity"])), dtype=float)
        valid_mask = np.asarray(image.get("valid", np.isfinite(intensity)), dtype=bool)
        column = self._nearest_axis(qr, qr_click)
        row = self._nearest_axis(qz, qz_click)

        snapped_sub_y = float(row)
        snapped_sub_x = float(column)
        if self.click_snap_to_max.isChecked() and self.click_snap_radius.value() > 0:
            snapped_sub_y, snapped_sub_x = self._refine_manual_click_position(
                image, row, column
            )
            row = int(np.clip(round(snapped_sub_y), 0, intensity.shape[0] - 1))
            column = int(np.clip(round(snapped_sub_x), 0, intensity.shape[1] - 1))

        value = float(intensity[row, column])
        if not np.isfinite(value):
            raise ValueError("The selected experimental point has a non-finite numerical intensity.")

        result = {
            "QrExp": float(qr[column]), "QzExp": float(qz[row]),
            "ExpIntensity": value, "ReferenceNPZ": str(path),
            "PixelX": int(column), "PixelY": int(row),
            "SubpixelX": float(column), "SubpixelY": float(row),
            "PeakHeight": value, "LocalBackground": np.nan, "BackgroundNoise": np.nan,
            "BackgroundGradient": np.nan, "PeakSNR": np.nan, "IntegratedSNR": np.nan,
            "SigmaExpIntensity": np.nan, "RelativeIntensityUncertainty": np.nan,
            "PeakAreaPixels": 1, "PeakSigmaQrWidth": np.nan, "PeakSigmaQzWidth": np.nan,
            "OverlapPeakCount": 0, "Deblended": False, "SaturationFraction": np.nan,
            "ValidPixelFraction": 1.0, "EdgeTruncated": False,
            "IntensityQualityScore": np.nan, "IntensityQuality": "unassessed",
            "IntensityNormalizationFactor": float(image.get("intensity_normalization_factor", 1.0)),
            "IntensityNormalization": str(image.get("intensity_normalization_note", "none")),
            "IntensitySource": "PNG-reconstructed NPZ" if image.get("png_reconstructed", False) else "native numerical NPZ",
            "PoissonTermUsed": False,
            "SigmaQrExp": np.nan, "SigmaQzExp": np.nan,
            "SigmaQrSystematic": np.nan, "SigmaQzSystematic": np.nan,
            "SigmaQrTotal": np.nan, "SigmaQzTotal": np.nan, "SigmaQExp": np.nan,
        }
        if self.click_subpixel_fit.isChecked():
            measured = _robust_local_peak_intensity(
                image, row, column,
                base_peak_radius=int(self.click_peak_radius.value()),
                background_radius=max(int(self.click_background_radius.value()), int(self.click_peak_radius.value()) + 3),
            )
            result.update(measured)
            result["ReferenceNPZ"] = str(path)

        # Position snapping is deliberately separate from intensity integration.
        # Keep the compact click-centered centroid as the plotted/manual q point,
        # even when the wider robust integration aperture is enabled.
        if self.click_snap_to_max.isChecked() and self.click_snap_radius.value() > 0:
            result["PixelX"] = int(column)
            result["PixelY"] = int(row)
            result["SubpixelX"] = float(snapped_sub_x)
            result["SubpixelY"] = float(snapped_sub_y)
            result["QrExp"] = float(np.interp(snapped_sub_x, np.arange(len(qr), dtype=float), qr))
            result["QzExp"] = float(np.interp(snapped_sub_y, np.arange(len(qz), dtype=float), qz))
        return self._augment_experimental_uncertainty(result, qr, qz)

    def _populate_calculator_candidates(self):
        self.click_candidate_combo.clear()
        self._candidate_prediction_indices = []
        if self.pending_clicked_experimental is None or self.calculator_predictions.empty:
            self.click_candidate_combo.setEnabled(False)
            return
        qr = float(self.pending_clicked_experimental["QrExp"])
        qz = float(self.pending_clicked_experimental["QzExp"])
        predictions = self.calculator_predictions
        dr = pd.to_numeric(predictions.QrCalc, errors="coerce").to_numpy(float) - qr
        dz = pd.to_numeric(predictions.QzCalc, errors="coerce").to_numpy(float) - qz
        score, distance, support = self._candidate_score_arrays(
            dr, dz, predictions, self.pending_clicked_experimental
        )

        if not getattr(self, "use_support_candidate_ranking", None) or not self.use_support_candidate_ranking.isChecked():
            score = distance.copy()

        # In one-to-one mode, remove calculated reflections that are already assigned
        # so one predicted reflection cannot be silently reused for multiple
        # experimental peaks.
        if (
            getattr(self, "enforce_one_to_one_assignments", None) is not None
            and self.enforce_one_to_one_assignments.isChecked()
            and not self.manual_click_assignments.empty
            and "_CalcID" in self.manual_click_assignments.columns
        ):
            used = set(
                pd.to_numeric(self.manual_click_assignments["_CalcID"], errors="coerce")
                .dropna().astype(int)
            )
            calc_ids = pd.to_numeric(predictions["CalcID"], errors="coerce").to_numpy(float)
            for i, calc_id in enumerate(calc_ids):
                if np.isfinite(calc_id) and int(calc_id) in used:
                    score[i] = np.inf

        finite = np.isfinite(score) & np.isfinite(distance)
        candidate_indices = np.where(finite)[0]
        if not len(candidate_indices):
            self.click_candidate_combo.setEnabled(False)
            return
        candidate_indices = candidate_indices[
            np.lexsort((distance[candidate_indices], score[candidate_indices]))
        ]
        order = candidate_indices[: int(self.click_candidate_count.value())]
        for index in order:
            row = predictions.iloc[int(index)]
            self._candidate_prediction_indices.append(int(index))
            stability_value = getattr(row, "AssignmentStability", np.nan)
            try:
                stability_text = f"{float(stability_value):.3f}" if np.isfinite(float(stability_value)) else "n/a"
            except Exception:
                stability_text = "n/a"
            self.click_candidate_combo.addItem(
                f"{row.HKL} | qcalc=({row.QrCalc:.5f}, {row.QzCalc:.5f}) | "
                f"Δq={distance[index]:.5f} | score={score[index]:.3f} | "
                f"support={support[index]:.2f} | Icalc={row.CalcIntensity:.6g} | "
                f"stability={getattr(row, 'StabilityTier', 'unassessed')} ({stability_text}) | "
                f"{row.OrientationDomain} | evidence={getattr(row, 'IndexingEvidenceSource', 'unreported')}"
            )
        self.click_candidate_combo.setEnabled(bool(self._candidate_prediction_indices))

    def _on_manual_plot_click(self, event):
        if not self._hybrid_ui_ready or event.inaxes is not self.manual_axes:
            return
        if event.xdata is None or event.ydata is None:
            return
        toolbar = getattr(self.manual_canvas, "toolbar", None)
        if toolbar is not None and str(getattr(toolbar, "mode", "")):
            return
        if event.button == 3:
            if self.pending_clicked_experimental is not None:
                self._cancel_pending_click()
            else:
                self._remove_nearest_manual_assignment(float(event.xdata), float(event.ydata))
            return
        if event.button != 1:
            return
        try:
            if self.pending_clicked_experimental is None:
                if self.manual_image is None:
                    raise ValueError("Load the PNG overlay image first.")
                if self.calculator_predictions.empty:
                    raise ValueError("Run the full automatic NPZ indexing calculator first.")
                self.pending_clicked_experimental = self._experimental_from_png_click(
                    float(event.xdata), float(event.ydata)
                )
                self._populate_calculator_candidates()
                point = self.pending_clicked_experimental
                self.click_pending_label.setText(
                    f"Experimental peak selected at q=({point['QrExp']:.5f}, {point['QzExp']:.5f}), "
                    f"Iexp={point['ExpIntensity']:.6g}, SNR={point.get('PeakSNR', np.nan):.2f}. "
                    "Now click the matching calculated white circle."
                )
                self._redraw_manual()
                return

            distance = np.hypot(
                self.calculator_predictions.QrCalc.to_numpy(float) - float(event.xdata),
                self.calculator_predictions.QzCalc.to_numpy(float) - float(event.ydata),
            )
            index = int(np.argmin(distance))
            if float(distance[index]) > float(self.click_calc_tolerance.value()):
                raise ValueError(
                    f"No calculated point is within {self.click_calc_tolerance.value():.4f} Å⁻¹ of that click. "
                    "Click closer to a white calculated circle or choose an HKL from the list."
                )
            self._complete_manual_assignment(index)
        except Exception as exc:
            QMessageBox.warning(self, "Manual comparison", str(exc))

    def _assign_selected_calculator_candidate(self):
        if self.pending_clicked_experimental is None:
            QMessageBox.information(self, "Assignment", "Click an experimental PNG peak first.")
            return
        combo_index = self.click_candidate_combo.currentIndex()
        if combo_index < 0 or combo_index >= len(self._candidate_prediction_indices):
            QMessageBox.information(self, "Assignment", "Choose a calculated candidate.")
            return
        self._complete_manual_assignment(self._candidate_prediction_indices[combo_index])

    def _complete_manual_assignment(self, prediction_index):
        if self.pending_clicked_experimental is None:
            return
        prediction = self.calculator_predictions.iloc[int(prediction_index)]
        exp = dict(self.pending_clicked_experimental)
        calc_key = int(prediction.CalcID)
        delta_q = math.hypot(
            float(exp["QrExp"]) - float(prediction.QrCalc),
            float(exp["QzExp"]) - float(prediction.QzCalc),
        )

        if (
            getattr(self, "block_very_poor_assignments", None) is not None
            and self.block_very_poor_assignments.isChecked()
            and delta_q > float(self.assignment_reject_threshold.value())
        ):
            QMessageBox.warning(
                self,
                "Assignment blocked",
                f"This pairing has Δq={delta_q:.5f} Å⁻¹, above the current "
                f"{self.assignment_reject_threshold.value():.5f} Å⁻¹ limit.\n\n"
                "Choose a closer HKL, run the trusted-match refinement, or explicitly uncheck "
                "the blocking option if you intentionally need to keep this outlier.",
            )
            return

        if (
            getattr(self, "enforce_one_to_one_assignments", None) is not None
            and self.enforce_one_to_one_assignments.isChecked()
            and not self.manual_click_assignments.empty
        ):
            existing_exp_distance = np.hypot(
                pd.to_numeric(self.manual_click_assignments["QrExp"], errors="coerce").to_numpy(float) - float(exp["QrExp"]),
                pd.to_numeric(self.manual_click_assignments["QzExp"], errors="coerce").to_numpy(float) - float(exp["QzExp"]),
            )
            if np.isfinite(existing_exp_distance).any() and float(np.nanmin(existing_exp_distance)) <= 0.004:
                QMessageBox.warning(
                    self,
                    "Experimental peak already paired",
                    "This experimental peak is already represented by an existing manual pair. "
                    "One-to-one assignment enforcement prevents adding the same peak twice.",
                )
                return

        duplicate = (
            not self.manual_click_assignments.empty
            and "_CalcID" in self.manual_click_assignments.columns
            and calc_key in set(
                pd.to_numeric(self.manual_click_assignments["_CalcID"], errors="coerce")
                .dropna().astype(int)
            )
        )
        if duplicate:
            if (
                getattr(self, "enforce_one_to_one_assignments", None) is not None
                and self.enforce_one_to_one_assignments.isChecked()
            ):
                QMessageBox.warning(
                    self,
                    "Calculated reflection already assigned",
                    f"{prediction.HKL} is already paired with another experimental peak. "
                    "One-to-one assignment enforcement is enabled, so choose another reflection.",
                )
                return
            answer = QMessageBox.question(
                self,
                "Calculated reflection already assigned",
                f"{prediction.HKL} is already paired with another experimental peak. Add another pairing anyway?",
            )
            if answer != QMessageBox.StandardButton.Yes:
                return

        dr = np.array([float(prediction.QrCalc) - float(exp["QrExp"])])
        dz = np.array([float(prediction.QzCalc) - float(exp["QzExp"])])
        assignment_score, _distance, support = self._candidate_score_arrays(
            dr, dz, self.calculator_predictions.iloc[[int(prediction_index)]], exp
        )
        (
            selected_calc_intensity, cluster_intensity, cluster_count, cluster_radius,
            cluster_hkls, selected_f2, selected_dwba_weight, selected_intensity_model,
        ) = self._calculated_unresolved_intensity(prediction_index, exp)
        row = {
            "ID": len(self.manual_click_assignments) + 1,
            "NormalizationReference": False,
            "HKL": str(prediction.HKL),
            "OrientationDomain": str(prediction.OrientationDomain),
            "CalculatorSource": str(prediction.CalculatorSource),
            "IndexingEvidenceSource": str(getattr(prediction, "IndexingEvidenceSource", "unreported")),
            "DetectionSource": str(getattr(prediction, "DetectionSource", "unreported")),
            "InputProvenance": self._npz_provenance_for_path(exp["ReferenceNPZ"]),
            "QrExp": exp["QrExp"], "QzExp": exp["QzExp"],
            "QrCalc": float(prediction.QrCalc), "QzCalc": float(prediction.QzCalc),
            "QrCalcOriginal": float(getattr(prediction, "QrCalcOriginal", prediction.QrCalc)),
            "QzCalcOriginal": float(getattr(prediction, "QzCalcOriginal", prediction.QzCalc)),
            "DeltaQr": exp["QrExp"] - float(prediction.QrCalc),
            "DeltaQz": exp["QzExp"] - float(prediction.QzCalc),
            "DeltaQ": delta_q,
            "ExpIntensity": exp["ExpIntensity"],
            "SigmaExpIntensity": exp.get("SigmaExpIntensity", np.nan),
            "RelativeIntensityUncertainty": exp.get("RelativeIntensityUncertainty", np.nan),
            "CalcIntensity": selected_calc_intensity,
            "CalcSingleIntensity": selected_calc_intensity,
            "CalcComparisonIntensity": cluster_intensity,
            "CalcF2": selected_f2,
            "DWBAWeight": selected_dwba_weight,
            "CalcOverlapCount": cluster_count,
            "CalcOverlapRadius": cluster_radius,
            "CalcOverlapHKLs": cluster_hkls,
            "CalcIntensityModel": selected_intensity_model,
            "ReferenceNPZ": exp["ReferenceNPZ"],
            "PixelX": exp["PixelX"], "PixelY": exp["PixelY"],
            "SubpixelX": exp.get("SubpixelX", exp["PixelX"]),
            "SubpixelY": exp.get("SubpixelY", exp["PixelY"]),
            "PeakHeight": exp.get("PeakHeight", np.nan),
            "LocalBackground": exp.get("LocalBackground", np.nan),
            "BackgroundNoise": exp.get("BackgroundNoise", np.nan),
            "BackgroundGradient": exp.get("BackgroundGradient", np.nan),
            "PeakSNR": exp.get("PeakSNR", np.nan),
            "IntegratedSNR": exp.get("IntegratedSNR", np.nan),
            "PeakAreaPixels": exp.get("PeakAreaPixels", np.nan),
            "PeakSigmaQrWidth": exp.get("PeakSigmaQrWidth", np.nan),
            "PeakSigmaQzWidth": exp.get("PeakSigmaQzWidth", np.nan),
            "OverlapPeakCount": exp.get("OverlapPeakCount", 0),
            "Deblended": exp.get("Deblended", False),
            "SaturationFraction": exp.get("SaturationFraction", np.nan),
            "ValidPixelFraction": exp.get("ValidPixelFraction", np.nan),
            "EdgeTruncated": exp.get("EdgeTruncated", False),
            "IntensityQualityScore": exp.get("IntensityQualityScore", np.nan),
            "IntensityQuality": exp.get("IntensityQuality", "unassessed"),
            "IntensityNormalizationFactor": exp.get("IntensityNormalizationFactor", 1.0),
            "IntensityNormalization": exp.get("IntensityNormalization", "none"),
            "IntensitySource": exp.get("IntensitySource", "unreported"),
            "IntensityUncertaintySource": exp.get("IntensityUncertaintySource", "unreported"),
            "PoissonTermUsed": exp.get("PoissonTermUsed", False),
            "SigmaQrExp": exp.get("SigmaQrExp", np.nan),
            "SigmaQzExp": exp.get("SigmaQzExp", np.nan),
            "SigmaQrSystematic": exp.get("SigmaQrSystematic", np.nan),
            "SigmaQzSystematic": exp.get("SigmaQzSystematic", np.nan),
            "SigmaQrTotal": exp.get("SigmaQrTotal", np.nan),
            "SigmaQzTotal": exp.get("SigmaQzTotal", np.nan),
            "SigmaQExp": exp.get("SigmaQExp", np.nan),
            "AssignmentScore": float(assignment_score[0]),
            "BackendSupportScore": float(support[0]),
            "AssignmentStability": float(getattr(prediction, "AssignmentStability", np.nan)),
            "StabilityTier": str(getattr(prediction, "StabilityTier", "unassessed")),
            "SuggestedHKL": "",
            "SuggestedDeltaQ": np.nan,
            "SuggestionImprovement": np.nan,
            "ReassignmentRecommended": False,
            "SuggestionStatus": "not_reviewed",
            "_CalcID": calc_key,
            "_SuggestedCalcID": np.nan,
        }
        self.manual_click_assignments = pd.concat(
            [self.manual_click_assignments, pd.DataFrame([row])], ignore_index=True, sort=False
        )
        self.pending_clicked_experimental = None
        self._candidate_prediction_indices = []
        self.click_candidate_combo.clear()
        self.click_candidate_combo.setEnabled(False)
        self.click_pending_label.setText(
            "Assignment added. Click another experimental PNG peak to continue."
        )
        self._recalculate_manual_assignment_results()
        self._save_manual_assignments_for_current_series()
        self._redraw_manual()

    # --------------------- strongest-peak normalization/table ------------------
    def _save_manual_assignments_for_current_series(self):
        series = self.manual_series_combo.currentText() if hasattr(self, "manual_series_combo") else ""
        if series:
            self.manual_click_assignments_by_series[series] = self.manual_click_assignments.copy()

    def _recalculate_manual_assignment_results(self):
        frame = self.manual_click_assignments.copy()
        if frame.empty:
            self.manual_click_assignments = pd.DataFrame(columns=self._ASSIGNMENT_COLUMNS)
            self._fill_manual_assignment_table()
            self.full_correlation_label.setText(
                "No manual comparison rows yet. Click an experimental PNG peak and pair it with an automatic calculated reflection."
            )
            self._update_accuracy_summary()
            return
        frame = frame.reset_index(drop=True)
        frame["ID"] = np.arange(1, len(frame) + 1, dtype=int)

        # Store backend-predicted q coordinates separately from optional bounded
        # manual mapping refinement so original and adjusted coordinates remain
        # independently auditable.
        for column in ("QrExp", "QzExp", "QrCalc", "QzCalc"):
            if column not in frame.columns:
                frame[column] = np.nan
        if "QrCalcOriginal" not in frame.columns:
            frame["QrCalcOriginal"] = pd.to_numeric(frame["QrCalc"], errors="coerce")
        if "QzCalcOriginal" not in frame.columns:
            frame["QzCalcOriginal"] = pd.to_numeric(frame["QzCalc"], errors="coerce")
        frame["QrCalcOriginal"] = pd.to_numeric(frame["QrCalcOriginal"], errors="coerce").combine_first(
            pd.to_numeric(frame["QrCalc"], errors="coerce")
        )
        frame["QzCalcOriginal"] = pd.to_numeric(frame["QzCalcOriginal"], errors="coerce").combine_first(
            pd.to_numeric(frame["QzCalc"], errors="coerce")
        )

        # Recompute systematic PNG-to-NPZ mapping uncertainty from the configured
        # pixel uncertainty so every assignment uses the same uncertainty model.
        for column in ("SigmaQrExp", "SigmaQzExp"):
            if column not in frame.columns:
                frame[column] = np.nan
        sys_qr = np.full(len(frame), np.nan, dtype=float)
        sys_qz = np.full(len(frame), np.nan, dtype=float)
        if "ReferenceNPZ" in frame.columns:
            for i, path in enumerate(frame["ReferenceNPZ"].astype(str)):
                sys_qr[i], sys_qz[i] = self._mapping_systematic_q_for_path(path)
        frame["SigmaQrSystematic"] = sys_qr
        frame["SigmaQzSystematic"] = sys_qz
        stat_qr = pd.to_numeric(frame["SigmaQrExp"], errors="coerce").to_numpy(float)
        stat_qz = pd.to_numeric(frame["SigmaQzExp"], errors="coerce").to_numpy(float)
        stat_qr_zero = np.where(np.isfinite(stat_qr), stat_qr, 0.0)
        stat_qz_zero = np.where(np.isfinite(stat_qz), stat_qz, 0.0)
        sys_qr_zero = np.where(np.isfinite(sys_qr), sys_qr, 0.0)
        sys_qz_zero = np.where(np.isfinite(sys_qz), sys_qz, 0.0)
        total_qr = np.hypot(stat_qr_zero, sys_qr_zero)
        total_qz = np.hypot(stat_qz_zero, sys_qz_zero)
        total_qr[~(np.isfinite(stat_qr) | np.isfinite(sys_qr))] = np.nan
        total_qz[~(np.isfinite(stat_qz) | np.isfinite(sys_qz))] = np.nan
        frame["SigmaQrTotal"] = total_qr
        frame["SigmaQzTotal"] = total_qz

        qr_exp = pd.to_numeric(frame["QrExp"], errors="coerce").to_numpy(float)
        qz_exp = pd.to_numeric(frame["QzExp"], errors="coerce").to_numpy(float)
        qr_calc = pd.to_numeric(frame["QrCalc"], errors="coerce").to_numpy(float)
        qz_calc = pd.to_numeric(frame["QzCalc"], errors="coerce").to_numpy(float)
        frame["DeltaQr"] = qr_exp - qr_calc
        frame["DeltaQz"] = qz_exp - qz_calc
        frame["DeltaQ"] = np.hypot(frame["DeltaQr"].to_numpy(float), frame["DeltaQz"].to_numpy(float))
        frame["ResidualDirectionDeg"] = np.degrees(
            np.arctan2(frame["DeltaQz"].to_numpy(float), frame["DeltaQr"].to_numpy(float))
        )

        frame["QTotalExp"] = np.hypot(qr_exp, qz_exp)
        frame["QTotalCalc"] = np.hypot(qr_calc, qz_calc)
        frame["DSpacingExp"] = self._safe_d_spacing(frame["QTotalExp"].to_numpy(float))
        frame["DSpacingCalc"] = self._safe_d_spacing(frame["QTotalCalc"].to_numpy(float))
        frame["DeltaD"] = frame["DSpacingExp"].to_numpy(float) - frame["DSpacingCalc"].to_numpy(float)
        q_total_exp = frame["QTotalExp"].to_numpy(float)
        sigma_q = np.full(len(frame), np.nan, dtype=float)
        valid_q = (
            np.isfinite(q_total_exp) & (q_total_exp > 1e-12)
            & np.isfinite(total_qr) & np.isfinite(total_qz)
        )
        sigma_q[valid_q] = np.sqrt(
            (qr_exp[valid_q] / q_total_exp[valid_q] * total_qr[valid_q]) ** 2
            + (qz_exp[valid_q] / q_total_exp[valid_q] * total_qz[valid_q]) ** 2
        )
        frame["SigmaQExp"] = sigma_q

        # Recompute measurement-specific DWBA weighting and unresolved calculated
        # spot intensity after q-map refinement or series restoration. The internal
        # reflection identifier preserves deterministic linkage to the selected HKL.
        if (
            "_CalcID" in frame.columns
            and self.calculator_predictions is not None
            and not self.calculator_predictions.empty
        ):
            calc_id_values = pd.to_numeric(self.calculator_predictions["CalcID"], errors="coerce").to_numpy(float)
            id_to_index = {int(value): index for index, value in enumerate(calc_id_values) if np.isfinite(value)}
            for row_index in frame.index:
                try:
                    calc_id = int(float(frame.at[row_index, "_CalcID"]))
                except Exception:
                    continue
                prediction_index = id_to_index.get(calc_id)
                if prediction_index is None:
                    continue
                experimental_point = frame.loc[row_index].to_dict()
                try:
                    (
                        selected_calc, cluster_calc, cluster_count, cluster_radius, cluster_hkls,
                        selected_f2, selected_weight, intensity_model,
                    ) = self._calculated_unresolved_intensity(prediction_index, experimental_point)
                except Exception:
                    continue
                frame.at[row_index, "CalcIntensity"] = selected_calc
                frame.at[row_index, "CalcSingleIntensity"] = selected_calc
                frame.at[row_index, "CalcComparisonIntensity"] = cluster_calc
                frame.at[row_index, "CalcF2"] = selected_f2
                frame.at[row_index, "DWBAWeight"] = selected_weight
                frame.at[row_index, "CalcOverlapCount"] = cluster_count
                frame.at[row_index, "CalcOverlapRadius"] = cluster_radius
                frame.at[row_index, "CalcOverlapHKLs"] = cluster_hkls
                frame.at[row_index, "CalcIntensityModel"] = intensity_model

        # Retain the theoretical intensity of the selected reflection, while using
        # the summed intensity of unresolved calculated reflections for experimental
        # spot comparison when multiple predictions fall within one resolvable peak.
        for column in ("ExpIntensity", "CalcIntensity"):
            if column not in frame.columns:
                frame[column] = np.nan
        if "CalcSingleIntensity" not in frame.columns:
            frame["CalcSingleIntensity"] = pd.to_numeric(frame["CalcIntensity"], errors="coerce")
        if "CalcComparisonIntensity" not in frame.columns:
            frame["CalcComparisonIntensity"] = pd.to_numeric(frame["CalcIntensity"], errors="coerce")
        frame["CalcComparisonIntensity"] = pd.to_numeric(frame["CalcComparisonIntensity"], errors="coerce").combine_first(
            pd.to_numeric(frame["CalcIntensity"], errors="coerce")
        )
        exp = pd.to_numeric(frame["ExpIntensity"], errors="coerce").to_numpy(float)
        calc = pd.to_numeric(frame["CalcComparisonIntensity"], errors="coerce").to_numpy(float)
        valid_exp = np.isfinite(exp)
        reference_index = int(np.nanargmax(np.where(valid_exp, exp, -np.inf))) if valid_exp.any() else 0
        reference_exp = float(exp[reference_index])
        reference_calc = float(calc[reference_index])
        frame["NormalizationReference"] = False
        frame.loc[reference_index, "NormalizationReference"] = True
        if not np.isfinite(reference_exp) or abs(reference_exp) < 1e-15:
            reference_exp = np.nan
        if not np.isfinite(reference_calc) or abs(reference_calc) < 1e-15:
            reference_calc = np.nan
        frame["ReferenceExperimentalIntensity"] = reference_exp
        frame["ReferenceCalculatedIntensity"] = reference_calc
        frame["ExpRelative"] = exp / reference_exp if np.isfinite(reference_exp) else np.nan
        frame["CalcRelative"] = calc / reference_calc if np.isfinite(reference_calc) else np.nan
        scale = reference_exp / reference_calc if np.isfinite(reference_exp) and np.isfinite(reference_calc) else np.nan
        frame["NormalizationScaleFactor"] = scale
        frame["CalcScaledToExperiment"] = calc * scale
        frame["ScaledResidual"] = exp - frame["CalcScaledToExperiment"].to_numpy(float)
        denominator = np.maximum(
            np.abs(exp),
            max(abs(reference_exp) if np.isfinite(reference_exp) else 1.0, 1e-12) * 1e-12,
        )
        frame["ScaledRelativeResidual"] = frame["ScaledResidual"].to_numpy(float) / denominator
        exp_rel = pd.to_numeric(frame["ExpRelative"], errors="coerce").to_numpy(float)
        calc_rel = pd.to_numeric(frame["CalcRelative"], errors="coerce").to_numpy(float)
        frame["IntensityAgreement"] = np.clip(
            1.0 - np.abs(exp_rel - calc_rel) / (np.abs(exp_rel) + np.abs(calc_rel) + 1e-12),
            0.0, 1.0,
        )

        # Robust multi-peak scaling now also uses experimental intensity uncertainty
        # and local quality.  Very low-quality, saturated, edge-truncated, or
        # positionally questionable peaks do not define the scale when enough good
        # alternatives exist.
        delta_q = pd.to_numeric(frame["DeltaQ"], errors="coerce").to_numpy(float)
        sigma_i = pd.to_numeric(frame.get("SigmaExpIntensity", pd.Series(np.nan, index=frame.index)), errors="coerce").to_numpy(float)
        quality_i = pd.to_numeric(frame.get("IntensityQualityScore", pd.Series(np.nan, index=frame.index)), errors="coerce").to_numpy(float)
        integrated_snr = pd.to_numeric(frame.get("IntegratedSNR", pd.Series(np.nan, index=frame.index)), errors="coerce").to_numpy(float)
        saturation = pd.to_numeric(frame.get("SaturationFraction", pd.Series(np.nan, index=frame.index)), errors="coerce").to_numpy(float)
        edge = frame.get("EdgeTruncated", pd.Series(False, index=frame.index)).fillna(False).astype(bool).to_numpy()
        trusted_intensity = np.isfinite(delta_q) & (delta_q <= self._ACCEPTABLE_DELTA_Q)
        quality_gate = (~np.isfinite(quality_i)) | (quality_i >= 0.40)
        snr_gate = (~np.isfinite(integrated_snr)) | (integrated_snr >= 3.0)
        saturation_gate = (~np.isfinite(saturation)) | (saturation <= 0.20)
        strict_trusted = trusted_intensity & quality_gate & snr_gate & saturation_gate & (~edge)
        if int(strict_trusted.sum()) >= 3:
            trusted_intensity = strict_trusted
        robust_scale, robust_scaled, log_residual, robust_log_rmse = self._robust_intensity_fit(
            exp, calc, trusted_intensity, sigma_i, quality_i
        )
        frame["IntensityScalePairCount"] = int(np.sum(
            trusted_intensity & np.isfinite(exp) & np.isfinite(calc) & (exp > 0) & (calc > 0)
        ))
        frame["RobustScaleFactor"] = robust_scale
        frame["CalcRobustScaledToExperiment"] = robust_scaled
        frame["LogIntensityResidual"] = log_residual
        robust_agreement = np.full(len(frame), np.nan, dtype=float)
        valid_log = np.isfinite(log_residual)
        robust_agreement[valid_log] = np.power(10.0, -np.abs(log_residual[valid_log]))
        frame["RobustIntensityAgreement"] = robust_agreement

        # If enough trusted peaks span a useful angular range, test a smooth empirical
        # preferred-orientation intensity envelope by leave-one-out validation. Apply
        # it only when held-out log-intensity prediction error clearly improves.
        orientation_applied, orientation_scale, orientation_scaled, orientation_residual, cv_base, cv_model = (
            self._cross_validated_orientation_intensity_correction(
                exp, calc, qr_exp, qz_exp, trusted_intensity, sigma_i, quality_i
            )
        )
        frame["EmpiricalOrientationScale"] = orientation_scale
        frame["CalcOrientationScaledToExperiment"] = orientation_scaled
        frame["OrientationLogIntensityResidual"] = orientation_residual
        orientation_agreement = np.full(len(frame), np.nan, dtype=float)
        valid_orientation = np.isfinite(orientation_residual)
        orientation_agreement[valid_orientation] = np.power(10.0, -np.abs(orientation_residual[valid_orientation]))
        frame["OrientationIntensityAgreement"] = orientation_agreement
        frame["OrientationCorrectionApplied"] = bool(orientation_applied)
        frame["OrientationCorrectionCVBaseline"] = cv_base
        frame["OrientationCorrectionCVModel"] = cv_model

        # Positional quality now uses centroid + PNG-to-NPZ mapping uncertainty,
        # not the centroid precision alone.
        uncertainty_ratio, quality = self._position_quality_metrics(
            delta_q, total_qr, total_qz
        )
        frame["UncertaintyRatio"] = uncertainty_ratio
        frame["PositionQuality"] = quality

        # Backend support is a secondary diagnostic used for candidate ranking;
        # q position remains the dominant assignment criterion.
        frame["BackendSupportScore"] = [
            self._prediction_support_score_values(source, evidence, stability)
            for source, evidence, stability in zip(
                frame.get("CalculatorSource", pd.Series("", index=frame.index)),
                frame.get("IndexingEvidenceSource", pd.Series("", index=frame.index)),
                frame.get("AssignmentStability", pd.Series(np.nan, index=frame.index)),
            )
        ]
        sigma_r_score = np.where(np.isfinite(total_qr) & (total_qr > 0), np.maximum(total_qr, 0.006), 0.010)
        sigma_z_score = np.where(np.isfinite(total_qz) & (total_qz > 0), np.maximum(total_qz, 0.006), 0.010)
        normalized_position = np.hypot(
            frame["DeltaQr"].to_numpy(float) / sigma_r_score,
            frame["DeltaQz"].to_numpy(float) / sigma_z_score,
        )
        stability_series = frame["AssignmentStability"] if "AssignmentStability" in frame.columns else pd.Series(np.nan, index=frame.index)
        stability = pd.to_numeric(stability_series, errors="coerce").to_numpy(float)
        stability_penalty = np.where(np.isfinite(stability), 0.35 * (1.0 - np.clip(stability, 0.0, 1.0)), 0.18)
        support_penalty = 0.75 * (1.0 - frame["BackendSupportScore"].to_numpy(float))
        frame["AssignmentScore"] = normalized_position + stability_penalty + support_penalty

        # Keep the most recent suggestion diagnostics until a new suggestion analysis
        # is computed so unrelated table refreshes do not erase interpretive metadata.
        defaults = {
            "SuggestedHKL": "",
            "SuggestedDeltaQ": np.nan,
            "SuggestionImprovement": np.nan,
            "ReassignmentRecommended": False,
            "SuggestionStatus": "not_reviewed",
        }
        for column, default in defaults.items():
            if column not in frame.columns:
                frame[column] = default

        for column in self._ASSIGNMENT_COLUMNS:
            if column not in frame.columns:
                frame[column] = np.nan
        private = [column for column in frame.columns if str(column).startswith("_")]
        self.manual_click_assignments = frame[self._ASSIGNMENT_COLUMNS + private]
        self._fill_manual_assignment_table()

        pearson, spearman, rmse = self._relative_intensity_statistics(exp_rel, calc_rel)
        reference_row = frame.iloc[reference_index]
        self.full_correlation_label.setText(
            f"Manual pairs: {len(frame)} | reference: ID {int(reference_row.ID)}, {reference_row.HKL}, "
            f"Iexp={reference_exp:.6g}, Icalc={reference_calc:.6g}, ref scale={scale:.6g} | "
            f"robust multi-peak scale={robust_scale:.6g}, robust log-RMSE={robust_log_rmse:.4f} | "
            f"orientation envelope={'applied' if orientation_applied else 'not applied'} "
            f"(LOO {cv_base:.4f}→{cv_model:.4f} dex) | "
            f"Pearson={pearson:.4f} | Spearman={spearman:.4f} | reference-normalized RMSE={rmse:.4f}"
        )
        self._update_accuracy_summary()

    def _refinement_source_coordinates_for_assignments(self, frame):
        """Return original backend q coordinates for each assignment when available."""
        qr = pd.to_numeric(frame.get("QrCalcOriginal", frame.get("QrCalc")), errors="coerce").to_numpy(float)
        qz = pd.to_numeric(frame.get("QzCalcOriginal", frame.get("QzCalc")), errors="coerce").to_numpy(float)
        if "_CalcID" in frame.columns and not self.calculator_predictions.empty:
            lookup = self.calculator_predictions.set_index("CalcID")
            for i, calc_id in enumerate(pd.to_numeric(frame["_CalcID"], errors="coerce").to_numpy(float)):
                if not np.isfinite(calc_id) or int(calc_id) not in lookup.index:
                    continue
                row = lookup.loc[int(calc_id)]
                try:
                    qr[i] = float(row.get("QrCalcOriginal", row.get("QrCalc")))
                    qz[i] = float(row.get("QzCalcOriginal", row.get("QzCalc")))
                except Exception:
                    pass
        return qr, qz

    def _refine_q_mapping_from_trusted_matches(self):
        """Fit and apply a conservative 2-D q-map refinement from trusted manual pairs."""
        if self.manual_click_assignments is None or self.manual_click_assignments.empty:
            QMessageBox.information(self, "q refinement", "Add manual experimental/calculated peak pairs first.")
            return
        if self.calculator_predictions is None or self.calculator_predictions.empty:
            QMessageBox.information(self, "q refinement", "Run and load the automatic calculator first.")
            return

        frame = self.manual_click_assignments.copy().reset_index(drop=True)
        dq = pd.to_numeric(frame.get("DeltaQ"), errors="coerce").to_numpy(float)
        threshold = float(self.refinement_trusted_threshold.value())
        trusted = np.isfinite(dq) & (dq <= threshold)
        # Avoid allowing an extremely low-SNR point to drive a refinement when
        # SNR is available, but do not discard rows whose SNR is simply unknown.
        if "PeakSNR" in frame.columns:
            snr = pd.to_numeric(frame["PeakSNR"], errors="coerce").to_numpy(float)
            trusted &= (~np.isfinite(snr)) | (snr >= 2.0)
        min_pairs = int(self.refinement_min_pairs.value())
        if int(trusted.sum()) < min_pairs:
            QMessageBox.warning(
                self,
                "Not enough trusted matches",
                f"Only {int(trusted.sum())} pairs satisfy Δq ≤ {threshold:.4f} Å⁻¹ (and SNR ≥ 2 when available). "
                f"At least {min_pairs} are required to avoid an unstable refinement.",
            )
            return

        qr_source, qz_source = self._refinement_source_coordinates_for_assignments(frame)
        qr_exp = pd.to_numeric(frame["QrExp"], errors="coerce").to_numpy(float)
        qz_exp = pd.to_numeric(frame["QzExp"], errors="coerce").to_numpy(float)
        sigma_qr_series = frame["SigmaQrTotal"] if "SigmaQrTotal" in frame.columns else pd.Series(np.nan, index=frame.index)
        sigma_qz_series = frame["SigmaQzTotal"] if "SigmaQzTotal" in frame.columns else pd.Series(np.nan, index=frame.index)
        sigma_qr = pd.to_numeric(sigma_qr_series, errors="coerce").to_numpy(float)
        sigma_qz = pd.to_numeric(sigma_qz_series, errors="coerce").to_numpy(float)
        valid = trusted & np.isfinite(qr_source) & np.isfinite(qz_source) & np.isfinite(qr_exp) & np.isfinite(qz_exp)
        if int(valid.sum()) < min_pairs:
            QMessageBox.warning(self, "q refinement", "Too few finite trusted points remain after validation.")
            return

        src_r = qr_source[valid]
        src_z = qz_source[valid]
        dst_r = qr_exp[valid]
        dst_z = qz_exp[valid]
        sig_r = sigma_qr[valid]
        sig_z = sigma_qz[valid]
        before = np.hypot(src_r - dst_r, src_z - dst_z)
        try:
            params = self._fit_q_mapping_parameters(src_r, src_z, dst_r, dst_z, sig_r, sig_z)
        except Exception as exc:
            QMessageBox.critical(self, "q refinement failed", str(exc))
            return
        fit_r, fit_z = self._transform_q_points(src_r, src_z, params)
        after = np.hypot(fit_r - dst_r, fit_z - dst_z)

        cv_before = np.nan
        cv_after = np.nan
        if int(valid.sum()) >= 6:
            held_before = []
            held_after = []
            for hold in range(len(src_r)):
                train = np.ones(len(src_r), dtype=bool)
                train[hold] = False
                try:
                    p_cv = self._fit_q_mapping_parameters(
                        src_r[train], src_z[train], dst_r[train], dst_z[train], sig_r[train], sig_z[train]
                    )
                    pred_r, pred_z = self._transform_q_points(
                        np.array([src_r[hold]]), np.array([src_z[hold]]), p_cv
                    )
                    held_before.append(math.hypot(src_r[hold] - dst_r[hold], src_z[hold] - dst_z[hold]))
                    held_after.append(math.hypot(float(pred_r[0]) - dst_r[hold], float(pred_z[0]) - dst_z[hold]))
                except Exception:
                    continue
            if held_before and held_after:
                cv_before = float(np.median(held_before))
                cv_after = float(np.median(held_after))

        median_before = float(np.median(before))
        median_after = float(np.median(after))
        if median_after >= median_before - 1e-5:
            QMessageBox.warning(
                self,
                "q refinement not applied",
                f"The robust fit did not improve the trusted-match median Δq ({median_before:.5f} → {median_after:.5f} Å⁻¹).",
            )
            return
        if np.isfinite(cv_before) and np.isfinite(cv_after) and cv_after > cv_before * 1.05:
            QMessageBox.warning(
                self,
                "q refinement rejected by holdout check",
                f"Training matches improved, but leave-one-out median Δq worsened ({cv_before:.5f} → {cv_after:.5f} Å⁻¹). "
                "The correction was not applied to avoid overfitting.",
            )
            return

        # Always transform from the unmodified backend coordinates so repeated
        # refinement clicks do not compound prior transforms.
        pred_qr_original = pd.to_numeric(
            self.calculator_predictions.get("QrCalcOriginal", self.calculator_predictions["QrCalc"]), errors="coerce"
        ).to_numpy(float)
        pred_qz_original = pd.to_numeric(
            self.calculator_predictions.get("QzCalcOriginal", self.calculator_predictions["QzCalc"]), errors="coerce"
        ).to_numpy(float)
        refined_qr, refined_qz = self._transform_q_points(pred_qr_original, pred_qz_original, params)
        self.calculator_predictions["QrCalc"] = refined_qr
        self.calculator_predictions["QzCalc"] = refined_qz

        lookup = self.calculator_predictions.set_index("CalcID")
        assignments = self.manual_click_assignments.copy()
        if "_CalcID" in assignments.columns:
            for row_index, calc_id in enumerate(pd.to_numeric(assignments["_CalcID"], errors="coerce").to_numpy(float)):
                if np.isfinite(calc_id) and int(calc_id) in lookup.index:
                    pred = lookup.loc[int(calc_id)]
                    assignments.loc[row_index, "QrCalc"] = float(pred["QrCalc"])
                    assignments.loc[row_index, "QzCalc"] = float(pred["QzCalc"])
                    assignments.loc[row_index, "QrCalcOriginal"] = float(pred["QrCalcOriginal"])
                    assignments.loc[row_index, "QzCalcOriginal"] = float(pred["QzCalcOriginal"])
        self.manual_click_assignments = assignments
        theta, scale_r, scale_z, off_r, off_z = params
        self._q_refinement = {
            "parameters": params.tolist(),
            "trusted_pairs": int(valid.sum()),
            "median_before": median_before,
            "median_after": median_after,
            "cv_before": cv_before,
            "cv_after": cv_after,
        }
        status = (
            f"Active: rotation {math.degrees(theta):+.3f}°, scale qr {scale_r:.5f}, scale qz {scale_z:.5f}, "
            f"offsets ({off_r:+.5f}, {off_z:+.5f}) Å⁻¹; trusted median Δq {median_before:.5f} → {median_after:.5f} Å⁻¹"
        )
        if np.isfinite(cv_before) and np.isfinite(cv_after):
            status += f"; leave-one-out {cv_before:.5f} → {cv_after:.5f} Å⁻¹"
        self.refinement_status_label.setText(status)
        self._recalculate_manual_assignment_results()
        self._suggest_one_to_one_reassignments(show_message=False)
        self._save_manual_assignments_for_current_series()
        self._redraw_manual()
        QMessageBox.information(
            self,
            "q refinement applied",
            status + "\n\nThis is a bounded 2-D comparison/mapping refinement. It does not rewrite the backend crystallographic orientation solution or files on disk.",
        )

    def _reset_q_mapping_refinement(self):
        if self.calculator_predictions is None or self.calculator_predictions.empty:
            return
        if "QrCalcOriginal" in self.calculator_predictions.columns:
            self.calculator_predictions["QrCalc"] = pd.to_numeric(
                self.calculator_predictions["QrCalcOriginal"], errors="coerce"
            )
        if "QzCalcOriginal" in self.calculator_predictions.columns:
            self.calculator_predictions["QzCalc"] = pd.to_numeric(
                self.calculator_predictions["QzCalcOriginal"], errors="coerce"
            )
        lookup = self.calculator_predictions.set_index("CalcID")
        assignments = self.manual_click_assignments.copy()
        if not assignments.empty and "_CalcID" in assignments.columns:
            for row_index, calc_id in enumerate(pd.to_numeric(assignments["_CalcID"], errors="coerce").to_numpy(float)):
                if np.isfinite(calc_id) and int(calc_id) in lookup.index:
                    pred = lookup.loc[int(calc_id)]
                    assignments.loc[row_index, "QrCalc"] = float(pred["QrCalc"])
                    assignments.loc[row_index, "QzCalc"] = float(pred["QzCalc"])
        self.manual_click_assignments = assignments
        self._q_refinement = None
        self.refinement_status_label.setText("No manual q-space refinement is active.")
        self._recalculate_manual_assignment_results()
        if not self.manual_click_assignments.empty:
            self._suggest_one_to_one_reassignments(show_message=False)
        self._save_manual_assignments_for_current_series()
        self._redraw_manual()

    def _suggest_one_to_one_reassignments(self, show_message=True):
        """Find a global one-to-one set of better HKL suggestions without silently changing rows."""
        if self.manual_click_assignments is None or self.manual_click_assignments.empty:
            if show_message:
                QMessageBox.information(self, "HKL suggestions", "Add manual peak pairs first.")
            return
        if self.calculator_predictions is None or self.calculator_predictions.empty:
            if show_message:
                QMessageBox.information(self, "HKL suggestions", "Load automatic-calculator predictions first.")
            return

        frame = self.manual_click_assignments.copy().reset_index(drop=True)
        predictions = self.calculator_predictions.reset_index(drop=True)
        exp_r = pd.to_numeric(frame["QrExp"], errors="coerce").to_numpy(float)
        exp_z = pd.to_numeric(frame["QzExp"], errors="coerce").to_numpy(float)
        pred_r = pd.to_numeric(predictions["QrCalc"], errors="coerce").to_numpy(float)
        pred_z = pd.to_numeric(predictions["QzCalc"], errors="coerce").to_numpy(float)
        support_series = predictions["BackendSupportScore"] if "BackendSupportScore" in predictions.columns else pd.Series(0.5, index=predictions.index)
        support = pd.to_numeric(support_series, errors="coerce").to_numpy(float)
        support = np.where(np.isfinite(support), np.clip(support, 0.0, 1.0), 0.5)
        stability_series = predictions["AssignmentStability"] if "AssignmentStability" in predictions.columns else pd.Series(np.nan, index=predictions.index)
        stability = pd.to_numeric(stability_series, errors="coerce").to_numpy(float)
        stability_penalty = np.where(np.isfinite(stability), 0.006 * (1.0 - np.clip(stability, 0.0, 1.0)), 0.003)
        support_penalty = 0.012 * (1.0 - support)
        max_delta = float(self.assignment_reject_threshold.value())
        n_rows = len(frame)
        n_pred = len(predictions)
        real_cost = np.full((n_rows, n_pred), 1e4, dtype=float)
        distance_matrix = np.full((n_rows, n_pred), np.inf, dtype=float)
        sigma_r_series = frame["SigmaQrTotal"] if "SigmaQrTotal" in frame.columns else pd.Series(np.nan, index=frame.index)
        sigma_z_series = frame["SigmaQzTotal"] if "SigmaQzTotal" in frame.columns else pd.Series(np.nan, index=frame.index)
        sigma_r_rows = pd.to_numeric(sigma_r_series, errors="coerce").to_numpy(float)
        sigma_z_rows = pd.to_numeric(sigma_z_series, errors="coerce").to_numpy(float)
        for i in range(n_rows):
            if not (np.isfinite(exp_r[i]) and np.isfinite(exp_z[i])):
                continue
            dr = pred_r - exp_r[i]
            dz = pred_z - exp_z[i]
            d = np.hypot(dr, dz)
            distance_matrix[i] = d
            allowed = np.isfinite(d) & (d <= max_delta)
            sigma_r = max(float(sigma_r_rows[i]), 0.006) if np.isfinite(sigma_r_rows[i]) else 0.010
            sigma_z = max(float(sigma_z_rows[i]), 0.006) if np.isfinite(sigma_z_rows[i]) else 0.010
            uncertainty_position_cost = 0.010 * np.hypot(dr / sigma_r, dz / sigma_z)
            real_cost[i, allowed] = (
                uncertainty_position_cost[allowed]
                + support_penalty[allowed]
                + stability_penalty[allowed]
            )
        # Dummy columns allow a peak to remain unassigned rather than forcing a
        # scientifically implausible match just to complete the assignment matrix.
        dummy_cost = np.full((n_rows, n_rows), max_delta + 0.025, dtype=float)
        for i in range(n_rows):
            dummy_cost[i, i] = max_delta + 0.015
        rows, cols = linear_sum_assignment(np.hstack([real_cost, dummy_cost]))

        frame["SuggestedHKL"] = ""
        frame["SuggestedDeltaQ"] = np.nan
        frame["SuggestionImprovement"] = np.nan
        frame["ReassignmentRecommended"] = False
        frame["SuggestionStatus"] = "not_reviewed"
        frame["_SuggestedCalcID"] = np.nan
        recommended = 0
        unassigned = 0
        for row_i, col in zip(rows, cols):
            current_dq = float(pd.to_numeric(pd.Series([frame.loc[row_i, "DeltaQ"]]), errors="coerce").iloc[0])
            if col >= n_pred or not np.isfinite(distance_matrix[row_i, col]) or distance_matrix[row_i, col] > max_delta:
                frame.loc[row_i, "SuggestionStatus"] = "no_supported_candidate"
                unassigned += 1
                continue
            pred = predictions.iloc[int(col)]
            suggested_dq = float(distance_matrix[row_i, col])
            improvement = current_dq - suggested_dq if np.isfinite(current_dq) else np.nan
            current_calc_id = np.nan
            if "_CalcID" in frame.columns:
                try:
                    current_calc_id = float(frame.loc[row_i, "_CalcID"])
                except Exception:
                    current_calc_id = np.nan
            new_calc_id = int(pred.CalcID)
            changed = not (np.isfinite(current_calc_id) and int(current_calc_id) == new_calc_id)
            recommend = bool(
                changed
                and np.isfinite(improvement)
                and improvement >= 0.003
                and suggested_dq < current_dq
                and (current_dq > self._ACCEPTABLE_DELTA_Q or improvement >= 0.008)
            )
            frame.loc[row_i, "SuggestedHKL"] = str(pred.HKL)
            frame.loc[row_i, "SuggestedDeltaQ"] = suggested_dq
            frame.loc[row_i, "SuggestionImprovement"] = improvement
            frame.loc[row_i, "ReassignmentRecommended"] = recommend
            frame.loc[row_i, "SuggestionStatus"] = "better_candidate" if recommend else "keep_current"
            frame.loc[row_i, "_SuggestedCalcID"] = new_calc_id
            if recommend:
                recommended += 1

        self.manual_click_assignments = frame
        self._last_suggestion_summary = (
            f"One-to-one review: {recommended} reassignment(s) recommended; "
            f"{unassigned} point(s) had no supported candidate within {max_delta:.4f} Å⁻¹."
        )
        self._recalculate_manual_assignment_results()
        self._save_manual_assignments_for_current_series()
        if show_message:
            QMessageBox.information(self, "HKL suggestions", self._last_suggestion_summary)

    def _apply_recommended_reassignments(self):
        if self.manual_click_assignments is None or self.manual_click_assignments.empty:
            QMessageBox.information(self, "Apply suggestions", "No manual assignments are available.")
            return
        frame = self.manual_click_assignments.copy().reset_index(drop=True)
        if "ReassignmentRecommended" not in frame.columns or not frame["ReassignmentRecommended"].fillna(False).astype(bool).any():
            self._suggest_one_to_one_reassignments(show_message=False)
            frame = self.manual_click_assignments.copy().reset_index(drop=True)
        recommended_mask = frame.get("ReassignmentRecommended", False)
        recommended_mask = pd.Series(recommended_mask).fillna(False).astype(bool).to_numpy()
        count = int(recommended_mask.sum())
        if not count:
            QMessageBox.information(self, "Apply suggestions", "No better one-to-one reassignment is currently recommended.")
            return
        answer = QMessageBox.question(
            self,
            "Apply HKL suggestions",
            f"Apply {count} recommended one-to-one HKL reassignment(s)?\n\n"
            "Only rows with a clearly smaller q residual are changed. The experimental click coordinates remain untouched.",
        )
        if answer != QMessageBox.StandardButton.Yes:
            return
        lookup = self.calculator_predictions.set_index("CalcID")
        changed = 0
        for i in np.where(recommended_mask)[0]:
            try:
                suggested_id = int(float(frame.loc[i, "_SuggestedCalcID"]))
            except Exception:
                continue
            if suggested_id not in lookup.index:
                continue
            pred = lookup.loc[suggested_id]
            for column in (
                "HKL", "OrientationDomain", "CalculatorSource", "IndexingEvidenceSource",
                "DetectionSource", "CalcIntensity", "AssignmentStability", "StabilityTier",
                "BackendSupportScore", "QrCalc", "QzCalc", "QrCalcOriginal", "QzCalcOriginal",
            ):
                if column in pred.index:
                    frame.loc[i, column] = pred[column]
            frame.loc[i, "_CalcID"] = suggested_id
            changed += 1
        frame["SuggestedHKL"] = ""
        frame["SuggestedDeltaQ"] = np.nan
        frame["SuggestionImprovement"] = np.nan
        frame["ReassignmentRecommended"] = False
        frame["SuggestionStatus"] = "applied_or_keep"
        frame["_SuggestedCalcID"] = np.nan
        self.manual_click_assignments = frame
        self._recalculate_manual_assignment_results()
        self._save_manual_assignments_for_current_series()
        self._redraw_manual()
        QMessageBox.information(self, "HKL suggestions applied", f"Applied {changed} reassignment(s).")

    def _remove_unsupported_assignments(self):
        if self.manual_click_assignments is None or self.manual_click_assignments.empty:
            return
        frame = self.manual_click_assignments.copy().reset_index(drop=True)
        dq = pd.to_numeric(frame.get("DeltaQ"), errors="coerce").to_numpy(float)
        threshold = float(self.assignment_reject_threshold.value())
        remove = np.isfinite(dq) & (dq > threshold)
        if "SuggestionStatus" in frame.columns:
            remove |= frame["SuggestionStatus"].astype(str).eq("no_supported_candidate").to_numpy()
        count = int(remove.sum())
        if not count:
            QMessageBox.information(
                self, "Remove unsupported pairs", "No assignments exceed the current threshold or lack a supported one-to-one candidate."
            )
            return
        answer = QMessageBox.question(
            self,
            "Remove unsupported pairs",
            f"Remove {count} assignment(s) that exceed Δq={threshold:.4f} Å⁻¹ or have no supported one-to-one candidate?\n\n"
            "This removes only the manual pair rows; it does not alter the automatic calculator output.",
        )
        if answer != QMessageBox.StandardButton.Yes:
            return
        self.manual_click_assignments = frame.loc[~remove].reset_index(drop=True)
        self._recalculate_manual_assignment_results()
        self._save_manual_assignments_for_current_series()
        self._redraw_manual()

    def _residual_diagnostics(self, frame):
        result = {
            "mean_delta_qr": np.nan,
            "mean_delta_qz": np.nan,
            "residual_bias_magnitude": np.nan,
            "corr_delta_qr_with_qr": np.nan,
            "corr_delta_qz_with_qz": np.nan,
            "poor_count": 0,
            "very_poor_count": 0,
            "domain_count": 0,
            "domain_summary": "",
        }
        if frame is None or frame.empty:
            return result
        dr = pd.to_numeric(frame.get("DeltaQr"), errors="coerce").to_numpy(float)
        dz = pd.to_numeric(frame.get("DeltaQz"), errors="coerce").to_numpy(float)
        dq = pd.to_numeric(frame.get("DeltaQ"), errors="coerce").to_numpy(float)
        qr = pd.to_numeric(frame.get("QrCalc"), errors="coerce").to_numpy(float)
        qz = pd.to_numeric(frame.get("QzCalc"), errors="coerce").to_numpy(float)
        finite_r = np.isfinite(dr)
        finite_z = np.isfinite(dz)
        if finite_r.any():
            result["mean_delta_qr"] = float(np.mean(dr[finite_r]))
        if finite_z.any():
            result["mean_delta_qz"] = float(np.mean(dz[finite_z]))
        if np.isfinite(result["mean_delta_qr"]) and np.isfinite(result["mean_delta_qz"]):
            result["residual_bias_magnitude"] = math.hypot(result["mean_delta_qr"], result["mean_delta_qz"])
        both = np.isfinite(dr) & np.isfinite(qr)
        if int(both.sum()) >= 4 and np.nanstd(dr[both]) > 0 and np.nanstd(qr[both]) > 0:
            result["corr_delta_qr_with_qr"] = float(np.corrcoef(dr[both], qr[both])[0, 1])
        both = np.isfinite(dz) & np.isfinite(qz)
        if int(both.sum()) >= 4 and np.nanstd(dz[both]) > 0 and np.nanstd(qz[both]) > 0:
            result["corr_delta_qz_with_qz"] = float(np.corrcoef(dz[both], qz[both])[0, 1])
        finite_dq = dq[np.isfinite(dq)]
        result["poor_count"] = int(np.sum(finite_dq > self._BORDERLINE_DELTA_Q))
        result["very_poor_count"] = int(np.sum(finite_dq > float(self.assignment_reject_threshold.value())))
        if "OrientationDomain" in frame.columns:
            counts = frame["OrientationDomain"].astype(str).replace("nan", "").value_counts()
            counts = counts[counts.index != ""]
            result["domain_count"] = int(len(counts))
            result["domain_summary"] = ", ".join(f"{name}:{int(count)}" for name, count in counts.items())
        return result

    def _npz_provenance_for_path(self, path):
        if not path:
            return "No reference NPZ"
        resolved = str(_GuiPath(path).expanduser().resolve())
        image = self._reference_npz_cache.get(resolved)
        if image is not None and image.get("input_provenance"):
            provenance = str(image["input_provenance"])
            normalization = str(image.get("intensity_normalization_note", ""))
            if normalization:
                provenance += f" Experimental intensity normalization: {normalization}."
            return provenance
        try:
            with np.load(resolved, allow_pickle=False) as archive:
                keys = set(archive.files)
                conversion_note = ""
                if "conversion_note" in keys:
                    conversion_note = str(np.asarray(archive["conversion_note"]).item())
                if (
                    "source_png" in keys
                    or "rendered png" in conversion_note.lower()
                    or "reconstructed" in conversion_note.lower()
                ):
                    return "PNG-reconstructed NPZ: explicit axes, approximate color-derived intensities"
                lower_keys = {key.lower() for key in keys}
                if {"qr", "qz"}.issubset(lower_keys):
                    return "Numerical NPZ with explicit reciprocal-space axes"
        except Exception:
            pass
        return "NPZ provenance unknown"

    def _accuracy_audit_record(self):
        record = {
            "audit_scope": "internal consistency only; no external crystallographic ground truth",
            "reference_npz": "",
            "input_provenance": "No reference NPZ",
            "explicit_q_axes": False,
            "automatic_decision": "not_run",
            "orientation_hkl": "",
            "indexed_count": 0,
            "stable_reflection_count": 0,
            "automatic_prediction_count": int(len(self.calculator_predictions)),
            "high_q_corner_indexed_count": 0,
            "indexing_evidence_sources": "",
            "manual_pair_count": int(len(self.manual_click_assignments)),
            "median_delta_q": np.nan,
            "mean_delta_q": np.nan,
            "max_delta_q": np.nan,
            "fraction_delta_q_le_0p02": np.nan,
            "fraction_delta_q_le_0p04": np.nan,
            "poor_pair_count": 0,
            "very_poor_pair_count": 0,
            "mean_delta_qr": np.nan,
            "mean_delta_qz": np.nan,
            "residual_bias_magnitude": np.nan,
            "corr_delta_qr_with_qr": np.nan,
            "corr_delta_qz_with_qz": np.nan,
            "median_peak_snr": np.nan,
            "median_sigma_q_total": np.nan,
            "median_backend_support": np.nan,
            "pearson_relative_intensity": np.nan,
            "spearman_relative_intensity": np.nan,
            "normalized_intensity_rmse": np.nan,
            "robust_multi_peak_scale": np.nan,
            "robust_log_intensity_rmse": np.nan,
            "median_robust_intensity_agreement": np.nan,
            "reassignment_recommended_count": 0,
            "orientation_domain_count": 0,
            "orientation_domain_summary": "",
            "q_refinement_active": bool(self._q_refinement is not None),
            "q_refinement_summary": "",
            "manual_spacegroup_number": np.nan,
            "manual_spacegroup_symbol": "",
            "cif_spacegroup_number": self.manual_cif_spacegroup_number,
            "cif_spacegroup_symbol": self.manual_cif_spacegroup_symbol,
            "spacegroup_consistent": False,
            "warnings": "",
        }
        path = (
            self.click_reference_npz_combo.currentData()
            if hasattr(self, "click_reference_npz_combo") else None
        )
        if path:
            resolved = str(_GuiPath(path).expanduser().resolve())
            record["reference_npz"] = resolved
            record["input_provenance"] = self._npz_provenance_for_path(resolved)
            image = self._reference_npz_cache.get(resolved)
            if image is not None:
                record["explicit_q_axes"] = bool(image.get("explicit_q_axes", False))
        item = getattr(self, "backend_series_item", None)
        if item is not None:
            record["automatic_decision"] = str(
                getattr(item, "final_decision", "") or "unreported"
            )
            record["orientation_hkl"] = str(getattr(item, "orientation_hkl", ""))
            record["indexed_count"] = int(getattr(item, "indexed_count", 0) or 0)
            record["stable_reflection_count"] = int(
                getattr(item, "stable_reflection_count", 0) or 0
            )
        predictions = self.calculator_predictions
        if predictions is not None and not predictions.empty:
            if "DetectionSource" in predictions.columns:
                detector = predictions["DetectionSource"].astype(str)
                record["high_q_corner_indexed_count"] = int(
                    detector.str.contains("high_q_corner", case=False, na=False).sum()
                )
            if "IndexingEvidenceSource" in predictions.columns:
                sources = sorted(
                    value for value in set(predictions["IndexingEvidenceSource"].astype(str))
                    if value not in ("", "nan")
                )
                record["indexing_evidence_sources"] = ", ".join(sources)

        frame = self.manual_click_assignments
        if frame is not None and not frame.empty:
            dq = pd.to_numeric(frame.get("DeltaQ"), errors="coerce").to_numpy(float)
            finite_dq = dq[np.isfinite(dq)]
            if len(finite_dq):
                record["median_delta_q"] = float(np.median(finite_dq))
                record["mean_delta_q"] = float(np.mean(finite_dq))
                record["max_delta_q"] = float(np.max(finite_dq))
                record["fraction_delta_q_le_0p02"] = float(
                    np.mean(finite_dq <= self._VERY_GOOD_DELTA_Q)
                )
                record["fraction_delta_q_le_0p04"] = float(
                    np.mean(finite_dq <= self._ACCEPTABLE_DELTA_Q)
                )
            residual = self._residual_diagnostics(frame)
            record.update({
                "poor_pair_count": residual["poor_count"],
                "very_poor_pair_count": residual["very_poor_count"],
                "mean_delta_qr": residual["mean_delta_qr"],
                "mean_delta_qz": residual["mean_delta_qz"],
                "residual_bias_magnitude": residual["residual_bias_magnitude"],
                "corr_delta_qr_with_qr": residual["corr_delta_qr_with_qr"],
                "corr_delta_qz_with_qz": residual["corr_delta_qz_with_qz"],
                "orientation_domain_count": residual["domain_count"],
                "orientation_domain_summary": residual["domain_summary"],
            })
            snr = pd.to_numeric(frame.get("PeakSNR"), errors="coerce").to_numpy(float)
            snr = snr[np.isfinite(snr)]
            if len(snr):
                record["median_peak_snr"] = float(np.median(snr))
            sigma_q = pd.to_numeric(frame.get("SigmaQExp"), errors="coerce").to_numpy(float)
            sigma_q = sigma_q[np.isfinite(sigma_q)]
            if len(sigma_q):
                record["median_sigma_q_total"] = float(np.median(sigma_q))
            support = pd.to_numeric(frame.get("BackendSupportScore"), errors="coerce").to_numpy(float)
            support = support[np.isfinite(support)]
            if len(support):
                record["median_backend_support"] = float(np.median(support))
            exp_rel = pd.to_numeric(frame.get("ExpRelative"), errors="coerce").to_numpy(float)
            calc_rel = pd.to_numeric(frame.get("CalcRelative"), errors="coerce").to_numpy(float)
            pearson, spearman, rmse = self._relative_intensity_statistics(exp_rel, calc_rel)
            record["pearson_relative_intensity"] = pearson
            record["spearman_relative_intensity"] = spearman
            record["normalized_intensity_rmse"] = rmse
            robust_scale_values = pd.to_numeric(frame.get("RobustScaleFactor"), errors="coerce").to_numpy(float)
            robust_scale_values = robust_scale_values[np.isfinite(robust_scale_values)]
            if len(robust_scale_values):
                record["robust_multi_peak_scale"] = float(np.median(robust_scale_values))
            log_res = pd.to_numeric(frame.get("LogIntensityResidual"), errors="coerce").to_numpy(float)
            trusted = np.isfinite(dq) & (dq <= self._ACCEPTABLE_DELTA_Q)
            valid_log = np.isfinite(log_res) & trusted
            if int(valid_log.sum()) < 2:
                valid_log = np.isfinite(log_res)
            if int(valid_log.sum()):
                record["robust_log_intensity_rmse"] = float(
                    np.sqrt(np.mean(log_res[valid_log] ** 2))
                )
            robust_agreement = pd.to_numeric(frame.get("RobustIntensityAgreement"), errors="coerce").to_numpy(float)
            robust_agreement = robust_agreement[np.isfinite(robust_agreement)]
            if len(robust_agreement):
                record["median_robust_intensity_agreement"] = float(np.median(robust_agreement))
            if "ReassignmentRecommended" in frame.columns:
                record["reassignment_recommended_count"] = int(
                    frame["ReassignmentRecommended"].fillna(False).astype(bool).sum()
                )
            sigma_i = pd.to_numeric(frame.get("SigmaExpIntensity"), errors="coerce").to_numpy(float)
            rel_i = pd.to_numeric(frame.get("RelativeIntensityUncertainty"), errors="coerce").to_numpy(float)
            quality_i = pd.to_numeric(frame.get("IntensityQualityScore"), errors="coerce").to_numpy(float)
            integrated_snr = pd.to_numeric(frame.get("IntegratedSNR"), errors="coerce").to_numpy(float)
            for key, values in (
                ("median_sigma_exp_intensity", sigma_i),
                ("median_relative_intensity_uncertainty", rel_i),
                ("median_intensity_quality_score", quality_i),
                ("median_integrated_snr", integrated_snr),
            ):
                finite_values = values[np.isfinite(values)]
                if len(finite_values):
                    record[key] = float(np.median(finite_values))
            quality_labels = frame.get("IntensityQuality", pd.Series("unassessed", index=frame.index)).astype(str).str.lower()
            record["high_intensity_quality_count"] = int((quality_labels == "high").sum())
            record["low_intensity_quality_count"] = int((quality_labels == "low").sum())
            record["deblended_peak_count"] = int(frame.get("Deblended", pd.Series(False, index=frame.index)).fillna(False).astype(bool).sum())
            record["edge_truncated_intensity_count"] = int(frame.get("EdgeTruncated", pd.Series(False, index=frame.index)).fillna(False).astype(bool).sum())
            saturation_values = pd.to_numeric(frame.get("SaturationFraction"), errors="coerce").to_numpy(float)
            record["saturation_flagged_count"] = int(np.sum(np.isfinite(saturation_values) & (saturation_values > 0.10)))
            overlap_values = pd.to_numeric(frame.get("CalcOverlapCount"), errors="coerce").to_numpy(float)
            record["calculated_unresolved_overlap_count"] = int(np.sum(np.isfinite(overlap_values) & (overlap_values > 1)))
            orientation_applied_values = frame.get("OrientationCorrectionApplied", pd.Series(False, index=frame.index)).fillna(False).astype(bool)
            record["empirical_orientation_intensity_correction_applied"] = bool(orientation_applied_values.any())
            cv_base_values = pd.to_numeric(frame.get("OrientationCorrectionCVBaseline"), errors="coerce").to_numpy(float)
            cv_model_values = pd.to_numeric(frame.get("OrientationCorrectionCVModel"), errors="coerce").to_numpy(float)
            if np.isfinite(cv_base_values).any():
                record["orientation_intensity_cv_baseline"] = float(np.nanmedian(cv_base_values))
            if np.isfinite(cv_model_values).any():
                record["orientation_intensity_cv_model"] = float(np.nanmedian(cv_model_values))
            scale_pair_values = pd.to_numeric(frame.get("IntensityScalePairCount"), errors="coerce").to_numpy(float)
            if np.isfinite(scale_pair_values).any():
                record["intensity_scale_pair_count"] = int(np.nanmax(scale_pair_values))

        if self._q_refinement is not None:
            p = self._q_refinement
            params = p.get("parameters", [np.nan] * 5)
            record["q_refinement_summary"] = (
                f"rotation_deg={math.degrees(float(params[0])):.5f}; "
                f"scale_qr={float(params[1]):.7g}; scale_qz={float(params[2]):.7g}; "
                f"offset_qr={float(params[3]):.7g}; offset_qz={float(params[4]):.7g}; "
                f"trusted_pairs={p.get('trusted_pairs', 0)}; "
                f"median_before={p.get('median_before', np.nan):.7g}; "
                f"median_after={p.get('median_after', np.nan):.7g}; "
                f"cv_before={p.get('cv_before', np.nan):.7g}; cv_after={p.get('cv_after', np.nan):.7g}"
            )

        manual_number, manual_symbol = self._selected_manual_spacegroup()
        record["manual_spacegroup_number"] = manual_number
        record["manual_spacegroup_symbol"] = manual_symbol
        record["spacegroup_consistent"] = bool(
            manual_number is not None
            and self.manual_cif_spacegroup_number is not None
            and manual_number == self.manual_cif_spacegroup_number
        )
        warnings = []
        if "PNG-reconstructed" in record["input_provenance"]:
            warnings.append("experimental intensities are reconstructed from PNG colors")
        if not record["explicit_q_axes"]:
            warnings.append("reference NPZ lacks explicit reciprocal-space axes")
        if manual_number == 1:
            warnings.append("P1 is highly permissive; chance overlaps are common")
        if (
            self.manual_cif_spacegroup_number is not None
            and manual_number != self.manual_cif_spacegroup_number
        ):
            warnings.append("manual space group differs from CIF space group")
        if record["automatic_decision"] in ("not_run", "", "INCONCLUSIVE", "REJECTED"):
            warnings.append(f"automatic validation decision is {record['automatic_decision']}")
        if record["manual_pair_count"] < 6:
            warnings.append("fewer than 6 manual peak pairs; positional statistics are preliminary")
        if record.get("low_intensity_quality_count", 0):
            warnings.append(
                f"{record['low_intensity_quality_count']} manual peak(s) have low experimental intensity quality"
            )
        if record.get("saturation_flagged_count", 0):
            warnings.append(
                f"{record['saturation_flagged_count']} manual peak(s) show possible clipping/saturation"
            )
        if record.get("edge_truncated_intensity_count", 0):
            warnings.append(
                f"{record['edge_truncated_intensity_count']} manual peak aperture(s) are truncated by an edge/mask"
            )
        if record.get("empirical_orientation_intensity_correction_applied", False):
            warnings.append(
                "an empirical preferred-orientation intensity envelope passed leave-one-out validation; "
                "it improves comparison but is not an independent structure-factor prediction"
            )
        if record["very_poor_pair_count"]:
            warnings.append(
                f"{record['very_poor_pair_count']} manual pair(s) exceed the current very-poor Δq threshold"
            )
        if np.isfinite(record["residual_bias_magnitude"]) and record["residual_bias_magnitude"] > 0.015:
            warnings.append("manual residuals show a systematic q-space offset; refinement/calibration should be checked")
        if np.isfinite(record["corr_delta_qr_with_qr"]) and abs(record["corr_delta_qr_with_qr"]) > 0.55:
            warnings.append("Δqr changes systematically with qr, suggesting scale/rotation mismatch")
        if np.isfinite(record["corr_delta_qz_with_qz"]) and abs(record["corr_delta_qz_with_qz"]) > 0.55:
            warnings.append("Δqz changes systematically with qz, suggesting scale/rotation mismatch")
        if record["orientation_domain_count"] > 1:
            warnings.append(
                f"manual pairs span {record['orientation_domain_count']} orientation domains; verify the multidomain interpretation"
            )
        if record["q_refinement_active"]:
            warnings.append(
                "manual 2-D q-map refinement is active; the automatic validation decision still refers to the original backend solution"
            )
        record["warnings"] = "; ".join(warnings)
        return record

    def _update_accuracy_summary(self):
        label = getattr(self, "accuracy_summary_label", None)
        if label is None:
            return
        record = self._accuracy_audit_record()
        provenance_label = getattr(self, "input_provenance_label", None)
        if provenance_label is not None:
            provenance_label.setText(record["input_provenance"])

        def fmt(value, digits=4):
            try:
                return f"{float(value):.{digits}f}" if np.isfinite(float(value)) else "n/a"
            except Exception:
                return "n/a"

        decision = record["automatic_decision"]
        color = (
            "#207a3c" if decision == "PASS"
            else "#a35a00" if decision in ("PASS_WITH_WARNINGS", "INCONCLUSIVE")
            else "#9b2226" if decision == "REJECTED"
            else "#555"
        )
        refinement_text = "inactive"
        if record["q_refinement_active"] and self._q_refinement is not None:
            p = self._q_refinement
            refinement_text = (
                f"active; trusted median Δq {fmt(p.get('median_before'), 5)} → {fmt(p.get('median_after'), 5)} Å⁻¹"
            )
            if np.isfinite(float(p.get("cv_after", np.nan))):
                refinement_text += (
                    f", holdout {fmt(p.get('cv_before'), 5)} → {fmt(p.get('cv_after'), 5)} Å⁻¹"
                )

        label.setText(
            f"<b>Automatic validation:</b> <span style='color:{color}'>{decision}</span> | "
            f"orientation {record['orientation_hkl'] or 'n/a'} | indexed {record['indexed_count']} | "
            f"stable {record['stable_reflection_count']} | display inventory {record['automatic_prediction_count']} | "
            f"high-q adaptive indexed {record['high_q_corner_indexed_count']}<br>"
            f"<b>Manual positional check:</b> {record['manual_pair_count']} pairs | "
            f"median Δq {fmt(record['median_delta_q'], 5)} Å⁻¹ | mean Δq {fmt(record['mean_delta_q'], 5)} Å⁻¹ | "
            f"≤{self._VERY_GOOD_DELTA_Q:.2f}: {fmt(100.0 * record['fraction_delta_q_le_0p02'], 1)}% | "
            f"≤{self._ACCEPTABLE_DELTA_Q:.2f}: {fmt(100.0 * record['fraction_delta_q_le_0p04'], 1)}% | "
            f"poor {record['poor_pair_count']} | very poor {record['very_poor_pair_count']}<br>"
            f"<b>Residual diagnostics:</b> mean (Δqᵣ, Δq_z)=({fmt(record['mean_delta_qr'], 5)}, "
            f"{fmt(record['mean_delta_qz'], 5)}) Å⁻¹ | bias {fmt(record['residual_bias_magnitude'], 5)} Å⁻¹ | "
            f"corr(Δqᵣ,qᵣ) {fmt(record['corr_delta_qr_with_qr'], 3)} | corr(Δq_z,q_z) {fmt(record['corr_delta_qz_with_qz'], 3)}<br>"
            f"<b>Uncertainty/support:</b> median σq {fmt(record['median_sigma_q_total'], 5)} Å⁻¹ | "
            f"median backend support {fmt(record['median_backend_support'], 3)} | "
            f"recommended one-to-one changes {record['reassignment_recommended_count']} | "
            f"domains {record['orientation_domain_summary'] or 'n/a'}<br>"
            f"<b>q-map refinement:</b> {refinement_text}<br>"
            f"<b>Intensity check:</b> Pearson {fmt(record['pearson_relative_intensity'])} | "
            f"Spearman {fmt(record['spearman_relative_intensity'])} | reference-normalized RMSE {fmt(record['normalized_intensity_rmse'])} | "
            f"robust log-RMSE {fmt(record['robust_log_intensity_rmse'])} | "
            f"median robust agreement {fmt(record['median_robust_intensity_agreement'])}<br>"
            f"<b>Experimental intensity quality:</b> high {record.get('high_intensity_quality_count', 0)} | "
            f"low {record.get('low_intensity_quality_count', 0)} | deblended {record.get('deblended_peak_count', 0)} | "
            f"possible clipping {record.get('saturation_flagged_count', 0)} | edge-truncated {record.get('edge_truncated_intensity_count', 0)} | "
            f"median integrated SNR {fmt(record.get('median_integrated_snr'))} | "
            f"median relative σI {fmt(record.get('median_relative_intensity_uncertainty'))}<br>"
            f"<b>Calculated intensity model:</b> unresolved-overlap spots {record.get('calculated_unresolved_overlap_count', 0)} | "
            f"scale pairs {record.get('intensity_scale_pair_count', 0)} | empirical orientation envelope "
            f"{'applied' if record.get('empirical_orientation_intensity_correction_applied', False) else 'not applied'} "
            f"(LOO {fmt(record.get('orientation_intensity_cv_baseline'))} → {fmt(record.get('orientation_intensity_cv_model'))} dex)<br>"
            f"<b>Input:</b> {record['input_provenance']}<br>"
            f"<b>Warnings:</b> {record['warnings'] or 'None from the internal audit.'}"
        )

    def _export_accuracy_audit(self):
        record = self._accuracy_audit_record()
        path, _ = QFileDialog.getSaveFileName(
            self,
            "Export accuracy audit",
            "giwaxs_internal_accuracy_audit.csv",
            "CSV files (*.csv)",
        )
        if not path:
            return
        if not path.lower().endswith(".csv"):
            path += ".csv"
        pd.DataFrame([record]).to_csv(path, index=False)
        details_path = str(
            _GuiPath(path).with_name(_GuiPath(path).stem + "_manual_pairs.csv")
        )
        saved_details = False
        if (
            self.manual_click_assignments is not None
            and not self.manual_click_assignments.empty
        ):
            self.manual_click_assignments.to_csv(details_path, index=False)
            saved_details = True
        message = f"Saved internal audit to:\n{path}"
        if saved_details:
            message += f"\n\nSaved manual-pair details to:\n{details_path}"
        QMessageBox.information(self, "Accuracy audit exported", message)

    def _simple_manual_assignment_table(self, frame):
        if frame is None or frame.empty:
            return pd.DataFrame()

        def numeric(name):
            if name not in frame.columns:
                return pd.Series(np.nan, index=frame.index, dtype=float)
            return pd.to_numeric(frame[name], errors="coerce")

        # Use the unresolved-spot comparison intensity when available because it
        # is the same calculated quantity used by the manual intensity comparison.
        predicted_intensity = numeric("CalcComparisonIntensity")
        if "CalcIntensity" in frame.columns:
            predicted_intensity = predicted_intensity.combine_first(numeric("CalcIntensity"))
        intensity_match = 100.0 * numeric("IntensityAgreement")
        delta_q = numeric("DeltaQ")

        return pd.DataFrame({
            "Peak": frame.get("ID", pd.Series(np.arange(1, len(frame) + 1), index=frame.index)),
            "Crystal Plane (HKL)": frame.get("HKL", pd.Series("", index=frame.index)),
            "Measured qᵣ (Å⁻¹)": numeric("QrExp"),
            "Measured qz (Å⁻¹)": numeric("QzExp"),
            "Predicted qᵣ (Å⁻¹)": numeric("QrCalc"),
            "Predicted qz (Å⁻¹)": numeric("QzCalc"),
            "Measured d-spacing (Å)": numeric("DSpacingExp"),
            "Predicted d-spacing (Å)": numeric("DSpacingCalc"),
            "Match Error Δq (Å⁻¹)": delta_q,
            "Measured Intensity": numeric("ExpIntensity"),
            "Predicted Intensity": predicted_intensity,
            "Intensity Match (%)": intensity_match,
            "Match Quality": [self._match_quality_label(value) for value in delta_q],
        })

    def _fill_manual_assignment_table(self):
        frame = self.manual_click_assignments.copy()
        simple_view = self._table_view_is_simple()
        if simple_view:
            display_frame = self._simple_manual_assignment_table(frame)
            public_columns = list(display_frame.columns)
            display_headers = public_columns
        else:
            public_columns = [
                column for column in self._ASSIGNMENT_DISPLAY_COLUMNS if column in frame.columns
            ]
            display_headers = [
                self._ASSIGNMENT_DISPLAY_HEADERS.get(column, column) for column in public_columns
            ]
            display_frame = frame[public_columns].copy() if public_columns else pd.DataFrame(index=frame.index)

        self.full_index_table.setSortingEnabled(False)
        self.full_index_table.clear()
        if frame.empty:
            self.full_index_table.setColumnCount(1)
            self.full_index_table.setHorizontalHeaderLabels(["Status"])
            self.full_index_table.setRowCount(1)
            self.full_index_table.setItem(
                0, 0, QTableWidgetItem("No manually paired experimental/calculated points yet.")
            )
            self.full_index_table.setSortingEnabled(True)
            return

        self.full_index_table.setColumnCount(len(public_columns))
        self.full_index_table.setHorizontalHeaderLabels(display_headers)
        self.full_index_table.setRowCount(len(display_frame))
        for row_index, row in display_frame.iterrows():
            try:
                dq_value = float(frame.loc[row_index, "DeltaQ"])
            except Exception:
                dq_value = np.nan
            recommendation_value = (
                frame.loc[row_index, "ReassignmentRecommended"]
                if "ReassignmentRecommended" in frame.columns else False
            )
            recommended = bool(pd.notna(recommendation_value) and bool(recommendation_value))
            tooltip = ""
            if np.isfinite(dq_value):
                tooltip = f"Experimental-to-calculated Δq = {dq_value:.5f} Å⁻¹."
            if recommended:
                tooltip += " A better one-to-one HKL reassignment is currently recommended."
            for column_index, value in enumerate(row):
                if pd.isna(value):
                    text = ""
                elif isinstance(value, (bool, np.bool_)):
                    text = "YES" if bool(value) else ""
                elif isinstance(value, (float, np.floating)):
                    if simple_view and display_headers[column_index] == "Intensity Match (%)":
                        text = f"{float(value):.1f}"
                    elif simple_view and "d-spacing" in display_headers[column_index]:
                        text = f"{float(value):.4f}"
                    elif simple_view and ("qᵣ" in display_headers[column_index] or "qz" in display_headers[column_index] or "Δq" in display_headers[column_index]):
                        text = f"{float(value):.5f}"
                    else:
                        text = f"{float(value):.7g}"
                else:
                    text = str(value)
                item = QTableWidgetItem(text)
                if tooltip:
                    item.setToolTip(tooltip)
                self.full_index_table.setItem(row_index, column_index, item)
        self.full_index_table.resizeColumnsToContents()
        if simple_view:
            self._apply_simple_table_header_tooltips()
        self.full_index_table.setSortingEnabled(True)

    def _cancel_pending_click(self):
        self.pending_clicked_experimental = None
        self._candidate_prediction_indices = []
        self.click_candidate_combo.clear()
        self.click_candidate_combo.setEnabled(False)
        self.click_pending_label.setText(
            "No pending point. First click an experimental peak on the PNG."
        )
        self._redraw_manual()

    def _undo_manual_assignment(self):
        if self.manual_click_assignments.empty:
            return
        self.manual_click_assignments = self.manual_click_assignments.iloc[:-1].copy()
        self._recalculate_manual_assignment_results()
        self._save_manual_assignments_for_current_series()
        self._redraw_manual()

    def _clear_manual_assignments(self):
        if self.manual_click_assignments.empty:
            return
        answer = QMessageBox.question(
            self, "Clear assignments", "Remove every manual experimental/calculated pairing?"
        )
        if answer != QMessageBox.StandardButton.Yes:
            return
        self.manual_click_assignments = pd.DataFrame(columns=self._ASSIGNMENT_COLUMNS)
        self._recalculate_manual_assignment_results()
        self._save_manual_assignments_for_current_series()
        self._redraw_manual()

    def _remove_nearest_manual_assignment(self, qr, qz):
        if self.manual_click_assignments.empty:
            return
        distance = np.hypot(
            pd.to_numeric(self.manual_click_assignments.QrExp, errors="coerce").to_numpy(float) - qr,
            pd.to_numeric(self.manual_click_assignments.QzExp, errors="coerce").to_numpy(float) - qz,
        )
        index = int(np.nanargmin(distance))
        if float(distance[index]) <= max(float(self.click_calc_tolerance.value()), 0.05):
            self.manual_click_assignments = self.manual_click_assignments.drop(
                self.manual_click_assignments.index[index]
            ).reset_index(drop=True)
            self._recalculate_manual_assignment_results()
            self._save_manual_assignments_for_current_series()
            self._redraw_manual()

    def _export_full_index_table(self):
        if self.manual_click_assignments.empty:
            QMessageBox.information(self, "Export", "No manual comparison rows are available.")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Save manual experimental/calculated comparison",
            "manual_png_peaks_vs_automatic_calculator.csv", "CSV files (*.csv)"
        )
        if path:
            self.manual_click_assignments[
                [column for column in self._ASSIGNMENT_COLUMNS if column in self.manual_click_assignments]
            ].to_csv(path, index=False)

    def _import_manual_assignments(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Import manual assignment CSV", "", "CSV files (*.csv);;All files (*)"
        )
        if not path:
            return
        try:
            frame = pd.read_csv(path)
            required = {"HKL", "QrExp", "QzExp", "QrCalc", "QzCalc", "ExpIntensity", "CalcIntensity"}
            missing = sorted(required - set(frame.columns))
            if missing:
                raise ValueError("Missing required columns: " + ", ".join(missing))
            self.manual_click_assignments = frame.copy()
            self._recalculate_manual_assignment_results()
            self._save_manual_assignments_for_current_series()
            self._redraw_manual()
        except Exception as exc:
            QMessageBox.critical(self, "Import error", str(exc))

    # ------------------------------- rendering --------------------------------
    def _redraw_manual(self):
        if not getattr(self, "_hybrid_ui_ready", False):
            return super()._redraw_manual()
        # Redraw only the display PNG and manual overlay elements. Automatic indexing
        # result files remain separate scientific outputs and are not modified by
        # visualization refreshes.
        saved_result = self.manual_result
        self.manual_result = pd.DataFrame()
        IntegratedGIXSWorkbench._redraw_manual(self)
        self.manual_result = saved_result
        axes = self.manual_axes

        # Draw the reciprocal-space reflection pattern generated from the
        # manually selected space group. This is independent of the automatic
        # NPZ-calculator reflection model and remains controlled by its own
        # visibility checkbox. HKL text is controlled only by Show HKL labels.
        if (
            getattr(self, "show_selected_spacegroup_structure", None) is not None
            and self.show_selected_spacegroup_structure.isChecked()
            and saved_result is not None
            and not saved_result.empty
        ):
            structure_qr = pd.to_numeric(
                saved_result["Qxy"], errors="coerce"
            ).to_numpy(float)
            structure_qz = pd.to_numeric(
                saved_result["Qz"], errors="coerce"
            ).to_numpy(float)
            structure_labels = saved_result["hkl"].astype(str).to_numpy()
            structure_intensity = pd.to_numeric(
                saved_result["Intensity"], errors="coerce"
            ).to_numpy(float)

            if self.manual_mirror.isChecked():
                nonzero = structure_qr > 1e-10
                structure_qr = np.concatenate([
                    structure_qr, -structure_qr[nonzero]
                ])
                structure_qz = np.concatenate([
                    structure_qz, structure_qz[nonzero]
                ])
                structure_labels = np.concatenate([
                    structure_labels, structure_labels[nonzero]
                ])
                structure_intensity = np.concatenate([
                    structure_intensity, structure_intensity[nonzero]
                ])

            visible_structure = (
                np.isfinite(structure_qr)
                & np.isfinite(structure_qz)
                & (structure_qr >= self.manual_qr_min.value())
                & (structure_qr <= self.manual_qr_max.value())
                & (structure_qz >= self.manual_qz_min.value())
                & (structure_qz <= self.manual_qz_max.value())
            )
            structure_qr = structure_qr[visible_structure]
            structure_qz = structure_qz[visible_structure]
            structure_labels = structure_labels[visible_structure]
            structure_intensity = structure_intensity[visible_structure]

            if len(structure_qr):
                finite_i = np.isfinite(structure_intensity) & (structure_intensity >= 0)
                structure_sizes = np.full(len(structure_qr), 28.0)
                if finite_i.any() and float(np.nanmax(structure_intensity[finite_i])) > 0:
                    structure_sizes[finite_i] = 24.0 + 72.0 * np.sqrt(
                        structure_intensity[finite_i]
                        / float(np.nanmax(structure_intensity[finite_i]))
                    )
                axes.scatter(
                    structure_qr,
                    structure_qz,
                    s=structure_sizes,
                    marker="D",
                    facecolors="none",
                    edgecolors="cyan",
                    linewidths=1.0,
                    alpha=0.78,
                    label=(
                        "selected space-group reflection pattern "
                        f"({len(structure_qr)})"
                    ),
                    zorder=3,
                )
                if self.manual_labels.isChecked():
                    for x, y, label in zip(
                        structure_qr, structure_qz, structure_labels
                    ):
                        annotation = axes.annotate(
                            str(label),
                            (float(x), float(y)),
                            xytext=(4, 4),
                            textcoords="offset points",
                            fontsize=6.2,
                            color="cyan",
                            annotation_clip=True,
                            zorder=4,
                        )
                        annotation.set_path_effects([
                            matplotlib_patheffects.Stroke(
                                linewidth=1.4, foreground="black"
                            ),
                            matplotlib_patheffects.Normal(),
                        ])

        # Display the automatic calculator's reflection model. All eligibility
        # filters are applied before intensity ranking. The strongest-point count
        # therefore equals the requested value whenever enough distinct valid
        # reciprocal-space positions exist.
        if self.show_calculator_points.isChecked() and not self.calculator_predictions.empty:
            frame = self.calculator_predictions.copy()
            frame["QrCalc"] = pd.to_numeric(frame["QrCalc"], errors="coerce")
            frame["QzCalc"] = pd.to_numeric(frame["QzCalc"], errors="coerce")
            frame["CalcIntensity"] = pd.to_numeric(
                frame["CalcIntensity"], errors="coerce"
            )
            visible = (
                np.isfinite(frame["QrCalc"])
                & np.isfinite(frame["QzCalc"])
                & frame.QrCalc.between(
                    self.manual_qr_min.value(), self.manual_qr_max.value()
                )
                & frame.QzCalc.between(
                    self.manual_qz_min.value(), self.manual_qz_max.value()
                )
            )
            shown = frame.loc[visible].copy()

            # Apply the selected numerical NPZ mask before ranking so masked
            # points do not consume slots in the requested strongest-point count.
            if (
                not shown.empty
                and getattr(self, "calculator_hide_masked", None) is not None
                and self.calculator_hide_masked.isChecked()
            ):
                try:
                    _mask_path, mask_image = self._reference_npz_data()
                    mask_qr = np.asarray(mask_image["qr"], dtype=float)
                    mask_qz = np.asarray(mask_image["qz"], dtype=float)
                    mask_array = np.asarray(
                        mask_image.get(
                            "valid",
                            np.isfinite(
                                mask_image.get(
                                    "raw_intensity", mask_image["intensity"]
                                )
                            ),
                        ),
                        dtype=bool,
                    )
                    mask_keep = np.zeros(len(shown), dtype=bool)
                    qr_low = float(np.nanmin(mask_qr))
                    qr_high = float(np.nanmax(mask_qr))
                    qz_low = float(np.nanmin(mask_qz))
                    qz_high = float(np.nanmax(mask_qz))
                    for point_index, row in enumerate(shown.itertuples(index=False)):
                        point_qr = float(row.QrCalc)
                        point_qz = float(row.QzCalc)
                        if not (
                            qr_low <= point_qr <= qr_high
                            and qz_low <= point_qz <= qz_high
                        ):
                            continue
                        column = self._nearest_axis(mask_qr, point_qr)
                        mask_row = self._nearest_axis(mask_qz, point_qz)
                        if (
                            0 <= mask_row < mask_array.shape[0]
                            and 0 <= column < mask_array.shape[1]
                        ):
                            mask_keep[point_index] = bool(
                                mask_array[mask_row, column]
                            )
                    shown = shown.loc[mask_keep].copy()
                except Exception:
                    # A missing reference mask should not prevent the automatic
                    # reflection model from being displayed.
                    pass

            # Rank deterministically and count distinct plotted positions rather
            # than overlapping HKL rows. This prevents a request for 120 markers
            # from looking like only a few dozen because several HKLs share the
            # same reciprocal-space coordinate.
            if not shown.empty:
                shown["_rank_intensity"] = shown["CalcIntensity"].where(
                    np.isfinite(shown["CalcIntensity"]), -np.inf
                )
                shown["_qr_key"] = shown["QrCalc"].round(6)
                shown["_qz_key"] = shown["QzCalc"].round(6)
                shown["_hkl_sort"] = shown.get(
                    "HKL", pd.Series("", index=shown.index)
                ).astype(str)
                shown = shown.sort_values(
                    ["_rank_intensity", "QzCalc", "QrCalc", "_hkl_sort"],
                    ascending=[False, True, True, True],
                    na_position="last",
                    kind="mergesort",
                )
                shown = shown.drop_duplicates(
                    ["_qr_key", "_qz_key"], keep="first"
                )

            available_positions = int(len(shown))
            strongest_only = bool(
                getattr(self, "calculator_strongest_only", None) is not None
                and self.calculator_strongest_only.isChecked()
            )
            requested_positions = int(
                self.calculator_max_hkls.value()
                if getattr(self, "calculator_max_hkls", None) is not None
                else available_positions
            )
            if strongest_only:
                shown = shown.head(requested_positions).copy()

            values = pd.to_numeric(
                shown["CalcIntensity"], errors="coerce"
            ).to_numpy(float)
            finite = np.isfinite(values) & (values >= 0)
            sizes = np.full(len(shown), 34.0)
            if finite.any() and float(np.nanmax(values[finite])) > 0:
                sizes[finite] = 22.0 + 80.0 * np.sqrt(
                    values[finite] / np.nanmax(values[finite])
                )

            displayed_positions = int(len(shown))
            if strongest_only and available_positions >= requested_positions:
                automatic_label = (
                    f"strongest automatic calculated reflections "
                    f"({displayed_positions})"
                )
            elif strongest_only:
                automatic_label = (
                    f"automatic calculated reflections "
                    f"({displayed_positions} available; "
                    f"{requested_positions} requested)"
                )
            else:
                automatic_label = (
                    f"automatic calculated reflections ({displayed_positions})"
                )

            if displayed_positions:
                axes.scatter(
                    shown.QrCalc,
                    shown.QzCalc,
                    s=sizes,
                    marker="o",
                    facecolors="none",
                    edgecolors="white",
                    linewidths=0.9,
                    alpha=0.72,
                    label=automatic_label,
                    zorder=4,
                )
                if self.manual_labels.isChecked() and "HKL" in shown.columns:
                    for row in shown.itertuples(index=False):
                        annotation = axes.annotate(
                            str(row.HKL),
                            (float(row.QrCalc), float(row.QzCalc)),
                            xytext=(4, 4),
                            textcoords="offset points",
                            fontsize=6.2,
                            color="white",
                            annotation_clip=True,
                            zorder=5,
                        )
                        annotation.set_path_effects([
                            matplotlib_patheffects.Stroke(
                                linewidth=1.4, foreground="black"
                            ),
                            matplotlib_patheffects.Normal(),
                        ])

        if not self.manual_click_assignments.empty:
            frame = self.manual_click_assignments
            for row in frame.itertuples(index=False):
                axes.plot(
                    [float(row.QrCalc), float(row.QrExp)],
                    [float(row.QzCalc), float(row.QzExp)],
                    color="white", alpha=0.55, linewidth=0.75, zorder=5,
                )
            axes.scatter(
                frame.QrCalc, frame.QzCalc, s=48, marker="o", facecolors="none",
                edgecolors="white", linewidths=1.25, label="selected calculated points", zorder=6,
            )
            axes.scatter(
                frame.QrExp, frame.QzExp, s=68, marker="o", facecolors="none",
                edgecolors="lime", linewidths=1.65,
                label=f"manually clicked experimental peaks ({len(frame)})", zorder=7,
            )
            reference = frame[frame.NormalizationReference.astype(bool)]
            if not reference.empty:
                axes.scatter(
                    reference.QrExp, reference.QzExp, s=125, marker="*", c="white",
                    edgecolors="black", linewidths=0.65, label="normalization reference", zorder=9,
                )
            if self.manual_labels.isChecked():
                for row in frame.itertuples(index=False):
                    annotation = axes.annotate(
                        str(row.HKL), (float(row.QrExp), float(row.QzExp)), xytext=(4, 4),
                        textcoords="offset points", fontsize=6.5, color="white",
                        annotation_clip=True, zorder=9,
                    )
                    annotation.set_path_effects([
                        matplotlib_patheffects.Stroke(linewidth=1.5, foreground="black"),
                        matplotlib_patheffects.Normal(),
                    ])

        if self.pending_clicked_experimental is not None:
            point = self.pending_clicked_experimental
            axes.scatter(
                [point["QrExp"]], [point["QzExp"]], s=90, marker="x", c="yellow",
                linewidths=1.8, label="pending experimental peak", zorder=10,
            )

        axes.set_title("PNG overlay: manual experimental peaks vs automatic NPZ calculator")
        handles, labels = axes.get_legend_handles_labels()
        if handles:
            # Deduplicate labels while preserving drawing order.
            unique = dict(zip(labels, handles))
            axes.legend(unique.values(), unique.keys(), fontsize=7.3, loc="upper right")
        self.manual_figure.tight_layout(pad=1.1)
        self.manual_canvas.draw_idle()

_MANUAL_FULL_INDEXING_WINDOW = None


def _running_inside_ipython() -> bool:
    """Return True when the program is running inside IPython/Jupyter."""
    try:
        return get_ipython() is not None
    except NameError:
        return False


def _launch_manual_full_indexing_application():
    """Launch the workbench correctly in scripts *and* Jupyter.

    In a normal Python process we use Qt's standard blocking event loop.
    In Jupyter/IPython, IPython integrates Qt's event loop with the notebook
    (the same idea as running ``%gui qt``).  We then create/show the window
    without calling ``app.exec()``, so the notebook remains usable.

    This function changes only GUI startup behavior; the indexing engine,
    manual-click calculations, table data, and plotting logic are untouched.
    """
    global _MANUAL_FULL_INDEXING_WINDOW

    app = QApplication.instance() or QApplication(sys.argv)
    app.setQuitOnLastWindowClosed(True)

    if _running_inside_ipython():
        shell = get_ipython()
        try:
            # Equivalent to `%gui qt`; required for a responsive Qt GUI in Jupyter.
            shell.run_line_magic("gui", "qt")
        except Exception as exc:
            print(
                "Jupyter could not enable Qt GUI integration automatically.\n"
                "If the window does not appear, run `%gui qt` in a cell and then "
                "run this file again.\n"
                f"Details: {exc}"
            )

        window = ManualCalculatorClickWorkbench()
        _MANUAL_FULL_INDEXING_WINDOW = window
        window.show()
        window.raise_()
        window.activateWindow()
        app.processEvents()

        print(
            "GIWAXS/GIXS workbench launched. "
            "The GUI should appear as a separate desktop window."
        )
        print(
            "Jupyter remains available while the GUI is open. "
            "Close the GUI normally when finished."
        )
        return window

    window = ManualCalculatorClickWorkbench()
    _MANUAL_FULL_INDEXING_WINDOW = window
    window.show()
    window.raise_()
    window.activateWindow()
    return app.exec()


if __name__ == "__main__":
    _launch_result = _launch_manual_full_indexing_application()
    if isinstance(_launch_result, int):
        raise SystemExit(_launch_result)
