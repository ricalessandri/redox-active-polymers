# Redox Active Polymers

If you use or wish to cite the models and data contained in this repo,
please refer to the following manuscript:

- R. Alessandri, C.-H. Li, S. Keating, K. T. Mohanty, A. Peng, J. L. Lutkenhaus, S. J. Rowan, D. P. Tabor, J. J. de Pablo,
  "Structural, Ionic, and Electronic Properties of Solid-State Phthalimide-Containing Polymers for All-Organic Batteries",
  *Chemrxiv*, **2024**, [10.26434/chemrxiv-2024-qkjr8](https://doi.org/10.26434/chemrxiv-2024-qkjr8).


## Atomistic models

### Polymers

- The PMAP AA model can be found in the [library of Polyply](https://github.com/marrink-lab/polyply_1.0/blob/master/LIBRARY.md).
- The PEPP AA model can be found in the [library of Polyply](https://github.com/marrink-lab/polyply_1.0/blob/master/LIBRARY.md).
- The PVBP AA model can be found in the [library of Polyply](https://github.com/marrink-lab/polyply_1.0/blob/master/LIBRARY.md).
- The PTMA AA model can be found in the [library of Polyply](https://github.com/marrink-lab/polyply_1.0/blob/master/LIBRARY.md).

### Electrolytes
- [electrolytes](./electrolytes)


## System configurations

- PMAP:DME:TBAPF6 system configurations relaxed at 300K (in Gromacs format):

| swelling % / polymer SoC |                0% |               20% |               60% |
|--------------------------|-------------------|-------------------|-------------------|
|  5%                      | [.gro](PMAP05000) | [.gro](PMAP05020) | [.gro](PMAP05060) |
| 10%                      | [.gro](PMAP10000) | [.gro](PMAP10020) | [.gro](PMAP10060) |
| 20%                      | [.gro](PMAP20000) | [.gro](PMAP20020) | [.gro](PMAP20060) |


[PMAP05000]: https://github.com/ricalessandri/redox-active-polymers/tree/main/configurations/PMAP000charge_DME_TBAPF6_05percent/relax-30mer-300K-D/1-relax-100ns-whole.gro
[PMAP05020]: https://github.com/ricalessandri/redox-active-polymers/tree/main/configurations/PMAP020charge_DME_TBAPF6_05percent/relax-30mer-300K-D/1-relax-100ns-whole.gro
[PMAP05060]: https://github.com/ricalessandri/redox-active-polymers/tree/main/configurations/PMAP060charge_DME_TBAPF6_05percent/relax-30mer-300K-D/1-relax-100ns-whole.gro
[PMAP10000]: https://github.com/ricalessandri/redox-active-polymers/tree/main/configurations/PMAP000charge_DME_TBAPF6_10percent/relax-30mer-300K-D/1-relax-100ns-whole.gro
[PMAP10020]: https://github.com/ricalessandri/redox-active-polymers/tree/main/configurations/PMAP020charge_DME_TBAPF6_10percent/relax-30mer-300K-D/1-relax-100ns-whole.gro
[PMAP10060]: https://github.com/ricalessandri/redox-active-polymers/tree/main/configurations/PMAP060charge_DME_TBAPF6_10percent/relax-30mer-300K-D/1-relax-100ns-whole.gro
[PMAP20000]: https://github.com/ricalessandri/redox-active-polymers/tree/main/configurations/PMAP000charge_DME_TBAPF6_20percent/relax-30mer-300K-D/1-relax-100ns-whole.gro
[PMAP20020]: https://github.com/ricalessandri/redox-active-polymers/tree/main/configurations/PMAP020charge_DME_TBAPF6_20percent/relax-30mer-300K-D/1-relax-100ns-whole.gro
[PMAP20060]: https://github.com/ricalessandri/redox-active-polymers/tree/main/configurations/PMAP060charge_DME_TBAPF6_20percent/relax-30mer-300K-D/1-relax-100ns-whole.gro

