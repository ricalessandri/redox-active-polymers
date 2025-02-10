
# Outer reorganization energy

Scripts to reproduce Figure S27 of [JACS Au, **2024**, DOI: 10.1021/jacsau.4c00276](https://doi.org/10.1021/jacsau.4c00276).

Execute:
```
PROC_OLAM_step01_get_cavity_radius.bash
PROC_OLAM_step02_compute_outer_lambda.py       # *note*: cavity values are hard-coded here!!
PROC_OLAM_step03_plot_lambdaout_vs_epsilon.py  # produces Figure S27 
```


## Appendix
1. Example Gaussian input for cavity calculation: 
   ```
   g16<<EOF
   %nprocshared=28
   %mem=16GB
   %Chk=NMPHTH-cavity-in-DMF-wB97XD.chk
   
   #P wB97XD/Def2SVPP opt=(MaxCycles=500, RecalcFC=50) freq scrf=(pcm,solvent=n,n-DiMethylFormamide) nosymm int(grid=ultrafine) scf=(xqc,tight)
   
   Cavity NMPHTH solvent
   
   0 1
   C   -3.166316  0.265612  0.220018
   N   -1.732181  0.156507  0.123598
   C   -1.087004  -0.994520  -0.237962
   O   -1.554294  -2.078564  -0.529826
   C   0.342924  -0.691548  -0.215206
   C   1.432805  -1.492129  -0.495705
   C   2.698957  -0.908399  -0.384467
   C   2.842488  0.441794  -0.001098
   C   1.722344  1.231674  0.277682
   C   0.483454  0.632442  0.160733
   C   -0.856757  1.175643  0.378246
   O   -1.081533  2.321012  0.719144
   H   -3.603908  0.023477  -0.753041
   H   -3.526930  -0.467577  0.947585
   H   -3.464967  1.271762  0.524061
   H   1.307301  -2.529221  -0.789483
   H   3.586730  -1.502142  -0.596094
   H   3.839186  0.872791  0.078236
   H   1.817701  2.271385  0.573580
   
   EOF
   ```

2. List of solvents (Gaussian):
   ```
   https://web.archive.org/web/20160922154735/http://www.gaussian.com/g_tech/g_ur/k_scrf.htm
   ```

