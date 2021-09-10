'multi_generate_networks.py' is designed to be used from terminal (bash, etc.).

Usage:
python multi_generate_networks.py [argumentName] [argument]


Included are 3 pre-generated networks used for simulations.
Network used for DC dynamics and chaos sections is in folder
'2016-09-08-155153_asn_nw_00100_nj_00261_seed_042_avl_100.00_disp_10.00.mat'

Network used for Supplementary Figure 5:
'2020-04-02-121953_asn_nw_00540_nj_01079_seed_001_avl_07.00_disp_01.60_gns_05.00_cdisp_700.00.mat'

Network used for Supplementary Figure 20:
'asn_nw_00992_nj_03008_seed_2821_avl_10.00_disp_01.00_lx_100.00_ly_100.00.mat'



To generate networks used in avalanche analysis run the following from commandline (bash):

'python multi_generate_networks.py --Lx 50 --LxMax 200 --numSims 4 --seedMax 3100 --density 0.06 --folder Density0.06ChangeSize'

'python multi_generate_networks.py --Lx 50 --LxMax 200 --numSims 4 --seedMax 1050 --density 0.10 --folder Density0.10ChangeSize'

'python multi_generate_networks.py --Lx 50 --LxMax 200 --numSims 4 --seedMax 1010 --density 0.14 --folder Density0.14ChangeSize'


Dependencies:
- networkX
- numpy
- scipy


Code written by past and present members (and collaborators) of Kuncic group at University of Sydney:
- Joel Hochstetter
- Miro Astore
- Paula Sanz-Leon
