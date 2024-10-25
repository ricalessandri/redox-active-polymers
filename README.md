# Redox-Active Polymers

If you use or wish to cite the models and data contained in this repo,
please refer to the following manuscript:

- R. Alessandri, C.-H. Li, S. Keating, K. T. Mohanty, A. Peng, J. L. Lutkenhaus, S. J. Rowan, D. P. Tabor, J. J. de Pablo,
  "Structural, Ionic, and Electronic Properties of Solid-State Phthalimide-Containing Polymers for All-Organic Batteries",
  *Chemrxiv*, **2024**, [10.26434/chemrxiv-2024-qkjr8](https://doi.org/10.26434/chemrxiv-2024-qkjr8).

Feel free to ask questions about the files in this repo by opening an [issue](https://github.com/ricalessandri/redox-active-polymers/issues).


## Atomistic models

### Polymers

- The PMAP AA model can be found in the [library of Polyply](https://github.com/marrink-lab/polyply_1.0/blob/master/LIBRARY.md);
  the topology for a 30-mer (with state of charge = 0%) can be found in [`polymers/`](./polymers).
- The PEPP AA model can be found in the [library of Polyply](https://github.com/marrink-lab/polyply_1.0/blob/master/LIBRARY.md);
  the topology for a 30-mer (with state of charge = 0%) can be found in [`polymers/`](./polymers).
- The PVBP AA model can be found in the [library of Polyply](https://github.com/marrink-lab/polyply_1.0/blob/master/LIBRARY.md);
  the topology for a 30-mer (with state of charge = 0%) can be found in [`polymers/`](./polymers).

### Electrolytes and solvents

- [electrolytes](./electrolytes)
- [solvents](./solvents)


## System configurations

- PMAP:DME:TBAPF6 system configurations relaxed at 300K (in Gromacs format):

| swelling % / polymer SoC |                0%                    |               20% |               60% |
|--------------------------|--------------------------------------|-------------------|-------------------|
|  5%                      | [.gro][PMAP05000]; [.top][tPMAP05000]| [.gro][PMAP05020] | [.gro][PMAP05060] |
| 10%                      | [.gro][PMAP10000]                    | [.gro][PMAP10020] | [.gro][PMAP10060] |
| 20%                      | [.gro][PMAP20000]                    | [.gro][PMAP20020] | [.gro][PMAP20060] |

- PEPP:DME:TBAPF6 system configurations relaxed at 300K (in Gromacs format):

| swelling % / polymer SoC |                0% |               20% |               60% |
|--------------------------|-------------------|-------------------|-------------------|
|  5%                      | [.gro][PEPP05000] | [.gro][PEPP05020] | [.gro][PEPP05060] |
| 10%                      | [.gro][PEPP10000] | [.gro][PEPP10020] | [.gro][PEPP10060] |
| 20%                      | [.gro][PEPP20000] | [.gro][PEPP20020] | [.gro][PEPP20060] |

- PVBP:DME:TBAPF6 system configurations relaxed at 300K (in Gromacs format):

| swelling % / polymer SoC |                0% |               20% |               60% |
|--------------------------|-------------------|-------------------|-------------------|
|  5%                      | [.gro][PVBP05000] | [.gro][PVBP05020] | [.gro][PVBP05060] |
| 10%                      | [.gro][PVBP10000] | [.gro][PVBP10020] | [.gro][PVBP10060] |
| 20%                      | [.gro][PVBP20000] | [.gro][PVBP20020] | [.gro][PVBP20060] |


## MD simulation protocol

- **equilibration**: [step1 - minimization](./MD_settings/eq_step1_min.mdp), [step2 - NVT](./MD_settings/eq_step2_NVT.mdp), [step3 - 100 bar, 900K](./MD_settings/eq_step3_NPT_highP.mdp), [step4 - 900K](./MD_settings/eq_step4_NPT_highT.mdp) 
- **cooling**: [NPT (900K to 100K)](./MD_settings/cool_NPT.mdp) 
- **relaxation**: [NPT (300K)](./MD_settings/relax_NPT.mdp) 

<h1 align="center">
  <div>
    <a href="https://github.com/ricalessandri/redox-active-polymers" target="_blank">
      <img src="docs/assets/MD_protocol.png" alt="MD" height="250">
    </a>
  </div>
</h1>


## Use the ML surrogate model to infer the electronic couplings
``` 
./PROC_MLVij_step01_predictions.bash                   # infer orbital overlaps
./PROC_MLVij_step02_plot_the_inferred_overlaps.bash    # plot the inferred orbital overlaps
./PROC_MLVij_step03_convert_to_EC_and_plot.bash        # convert to electronic couplings and plot
``` 
Example result can be found [here](./predictions/MODELMONO007-monomers-2x600K595K585K575K475K465K455K445K435K345K335K325K315K305K/PMAP000charge_DME_TBAPF6_05percent_30mer_300KD100ns/pair-predictions-300K-D-100ns/0_data_couplings_cutoff09A.pdf).



[PMAP05000]: https://github.com/ricalessandri/redox-active-polymers/tree/main/configurations/PMAP000charge_DME_TBAPF6_05percent/relax-30mer-300K-D/1-relax-100ns-whole.gro
[tPMAP05000]: https://github.com/ricalessandri/redox-active-polymers/tree/main/configurations/PMAP000charge_DME_TBAPF6_05percent/relax-30mer-300K-D/system_melt.top
[PMAP05020]: https://github.com/ricalessandri/redox-active-polymers/tree/main/configurations/PMAP020charge_DME_TBAPF6_05percent/relax-30mer-300K-D/1-relax-100ns-whole.gro
[PMAP05060]: https://github.com/ricalessandri/redox-active-polymers/tree/main/configurations/PMAP060charge_DME_TBAPF6_05percent/relax-30mer-300K-D/1-relax-100ns-whole.gro
[PMAP10000]: https://github.com/ricalessandri/redox-active-polymers/tree/main/configurations/PMAP000charge_DME_TBAPF6_10percent/relax-30mer-300K-D/1-relax-100ns-whole.gro
[PMAP10020]: https://github.com/ricalessandri/redox-active-polymers/tree/main/configurations/PMAP020charge_DME_TBAPF6_10percent/relax-30mer-300K-D/1-relax-100ns-whole.gro
[PMAP10060]: https://github.com/ricalessandri/redox-active-polymers/tree/main/configurations/PMAP060charge_DME_TBAPF6_10percent/relax-30mer-300K-D/1-relax-100ns-whole.gro
[PMAP20000]: https://github.com/ricalessandri/redox-active-polymers/tree/main/configurations/PMAP000charge_DME_TBAPF6_20percent/relax-30mer-300K-D/1-relax-100ns-whole.gro
[PMAP20020]: https://github.com/ricalessandri/redox-active-polymers/tree/main/configurations/PMAP020charge_DME_TBAPF6_20percent/relax-30mer-300K-D/1-relax-100ns-whole.gro
[PMAP20060]: https://github.com/ricalessandri/redox-active-polymers/tree/main/configurations/PMAP060charge_DME_TBAPF6_20percent/relax-30mer-300K-D/1-relax-100ns-whole.gro
[PEPP05000]: https://github.com/ricalessandri/redox-active-polymers/tree/main/configurations/PEPP000charge_DME_TBAPF6_05percent/relax-30mer-300K-D/1-relax-100ns-whole.gro
[PEPP05020]: https://github.com/ricalessandri/redox-active-polymers/tree/main/configurations/PEPP020charge_DME_TBAPF6_05percent/relax-30mer-300K-D/1-relax-100ns-whole.gro
[PEPP05060]: https://github.com/ricalessandri/redox-active-polymers/tree/main/configurations/PEPP060charge_DME_TBAPF6_05percent/relax-30mer-300K-D/1-relax-100ns-whole.gro
[PEPP10000]: https://github.com/ricalessandri/redox-active-polymers/tree/main/configurations/PEPP000charge_DME_TBAPF6_10percent/relax-30mer-300K-D/1-relax-100ns-whole.gro
[PEPP10020]: https://github.com/ricalessandri/redox-active-polymers/tree/main/configurations/PEPP020charge_DME_TBAPF6_10percent/relax-30mer-300K-D/1-relax-100ns-whole.gro
[PEPP10060]: https://github.com/ricalessandri/redox-active-polymers/tree/main/configurations/PEPP060charge_DME_TBAPF6_10percent/relax-30mer-300K-D/1-relax-100ns-whole.gro
[PEPP20000]: https://github.com/ricalessandri/redox-active-polymers/tree/main/configurations/PEPP000charge_DME_TBAPF6_20percent/relax-30mer-300K-D/1-relax-100ns-whole.gro
[PEPP20020]: https://github.com/ricalessandri/redox-active-polymers/tree/main/configurations/PEPP020charge_DME_TBAPF6_20percent/relax-30mer-300K-D/1-relax-100ns-whole.gro
[PEPP20060]: https://github.com/ricalessandri/redox-active-polymers/tree/main/configurations/PEPP060charge_DME_TBAPF6_20percent/relax-30mer-300K-D/1-relax-100ns-whole.gro
[PVBP05000]: https://github.com/ricalessandri/redox-active-polymers/tree/main/configurations/PVBP000charge_DME_TBAPF6_05percent/relax-30mer-300K-D/1-relax-100ns-whole.gro
[PVBP05020]: https://github.com/ricalessandri/redox-active-polymers/tree/main/configurations/PVBP020charge_DME_TBAPF6_05percent/relax-30mer-300K-D/1-relax-100ns-whole.gro
[PVBP05060]: https://github.com/ricalessandri/redox-active-polymers/tree/main/configurations/PVBP060charge_DME_TBAPF6_05percent/relax-30mer-300K-D/1-relax-100ns-whole.gro
[PVBP10000]: https://github.com/ricalessandri/redox-active-polymers/tree/main/configurations/PVBP000charge_DME_TBAPF6_10percent/relax-30mer-300K-D/1-relax-100ns-whole.gro
[PVBP10020]: https://github.com/ricalessandri/redox-active-polymers/tree/main/configurations/PVBP020charge_DME_TBAPF6_10percent/relax-30mer-300K-D/1-relax-100ns-whole.gro
[PVBP10060]: https://github.com/ricalessandri/redox-active-polymers/tree/main/configurations/PVBP060charge_DME_TBAPF6_10percent/relax-30mer-300K-D/1-relax-100ns-whole.gro
[PVBP20000]: https://github.com/ricalessandri/redox-active-polymers/tree/main/configurations/PVBP000charge_DME_TBAPF6_20percent/relax-30mer-300K-D/1-relax-100ns-whole.gro
[PVBP20020]: https://github.com/ricalessandri/redox-active-polymers/tree/main/configurations/PVBP020charge_DME_TBAPF6_20percent/relax-30mer-300K-D/1-relax-100ns-whole.gro
[PVBP20060]: https://github.com/ricalessandri/redox-active-polymers/tree/main/configurations/PVBP060charge_DME_TBAPF6_20percent/relax-30mer-300K-D/1-relax-100ns-whole.gro

