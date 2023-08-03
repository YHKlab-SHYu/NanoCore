#
# VASP
#
### supply POTCAR directory
vaspini = '/TGM/Apps/VASP/POTCAR/'

### supply VASP executable directory
vasp_odir="/TGM/Apps/VASP/OLD_BIN/5.4.4/O2/NORMAL"
vasp_ndir="/TGM/Apps/VASP/5.4.4.pl2"
vasp_dirv6="/TGM/Apps/VASP/bin/6.3.1"

vasp_dir = vasp_odir

vasp_calculator  = vasp_dir + '/vasp.5.4.4.pl2.O2.NORMAL.std.x' 
vasp_POTCAR_LDA  = vaspini + '1.POTPAW.LDA.54.RECOMMEND'
vasp_POTCAR_PBE  = vaspini + '2.POTPAW.PBE.54.RECOMMEND'
vasp_POTCAR_PW91 = vaspini + '2.POTPAW.PBE.54.RECOMMEND'        # connect to correct directory

#
# SIESTA
#

siesta_dir = '/home2/jwlee/siesta-v4.1-b4'
# Siesta calculator location
siesta_calculator = siesta_dir + '/Obj/siesta'

# Pseudopentiential files location
siesta_psf_location = '/home2/jwlee/pseudopotential'

# Pyprojection files location
import os
temporary_location  = os.getcwd()
siesta_pyprojection = temporary_location + '/pyprojection/pyprojection.py'

# Siesta utilities location
siesta_util_tbtrans = siesta_dir + '/Util/TS/TBtrans/tbtrans'
siesta_util_band = 'Bands/gnubands'
siesta_util_dos  = 'Eig2DOS/Eig2DOS'
siesta_util_pdos = 'Contrib/APostnikov/fmpdos'
siesta_util_rho  = 'Contrib/APostnikov/rho2xsf'
siesta_util_vh   = 'Macroave/Src/macroave'

