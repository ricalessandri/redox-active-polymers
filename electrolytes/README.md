# Electrolytes

## TBA & TEA
Merged topology from LigParGen with nonbonded parameters from [2014JNCanongiaLopes-JPCB](https://doi.org/10.1021/jp0476545).
Charged are subsequently scaled by 0.8 following [2017BDoherty-JCTC](https://doi.org/10.1021/acs.jctc.7b00520).

## PF6, BF4, & ClO4
From [https://github.com/orlandoacevedo/IL](https://github.com/orlandoacevedo/IL):
```
mkdir -p PF6_08scaled
cp IL/0.8\*2009IL/ITP/PF6_scale0.8.itp           PF6_08scaled/PF6_08scaled.itp
cp IL/0.8\*2009IL/ITP/PF6_atomtypes_scale0.8.itp PF6_08scaled/PF6_atomtypes_08scaled.itp
cp IL/PDB/PF6.pdb                                PF6_08scaled/.

mkdir -p BF4_08scaled
cp IL/0.8\*2009IL/ITP/BF4_scale0.8.itp           BF4_08scaled/BF4_08scaled.itp
cp IL/0.8\*2009IL/ITP/BF4_atomtypes_scale0.8.itp BF4_08scaled/BF4_atomtypes_08scaled.itp
cp IL/PDB/BF4.pdb                                BF4_08scaled/.

mkdir -p ClO4_08scaled
cp IL/0.8\*2009IL/ITP/ClO4_scale0.8.itp           ClO4_08scaled/ClO4_08scaled.itp
cp IL/0.8\*2009IL/ITP/ClO4_atomtypes_scale0.8.itp ClO4_08scaled/ClO4_atomtypes_08scaled.itp
cp IL/PDB/ClO4.pdb                                ClO4_08scaled/.
```

## Li
```
echo "[ atomtypes ]"                                              >  Li_08scaled/Li_atomtypes_08scaled.itp
echo "; scaled by 0.8"                                            >> Li_08scaled/Li_atomtypes_08scaled.itp
echo "; name  mass  charge  ptype   sigma   epsilon"              >> Li_08scaled/Li_atomtypes_08scaled.itp
grep opls_406 /usr/local/gromacs/share/gromacs/top/oplsaa.ff/*itp >> Li_08scaled/Li_atomtypes_08scaled.itp
```

