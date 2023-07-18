# Synthetic-NEM-2000bus-Data
This repository contains the Matpower files of a synthetic network model for National Electricity Market (NEM) of Australia.

- `snem2000.m'  contains synthetic network data for all states.
- `snem1803.m'  contains synthetic network data for the mainland sub-network.
- `snem197.m'   contains synthetic network data for Tamania sub-network.

The data set development algorithm and technical report can be found at https://arxiv.org/abs/2306.08176,
with CSIRO Data Access Portal release at https://doi.org/10.25919/esxz-q276.

## MatPower
MATPOWER solves the cases using either IPOPT or Artelys Knitro with the default parameters.
 - r = runopf('snem2000', mpoption('verbose', 2, 'opf.ac.solver', 'IPOPT'));
 - r = runopf('snem2000', mpoption('verbose', 2, 'opf.ac.solver', 'KNITRO'));

The default MIPS solver can also solve the problem with on step-size control and increase the number of iterations.
 - e.g. r = runopf('snem2000', mpoption('verbose', 2, 'mips.step_control', 1, 'mips.max_it', 400));


## License
This data is licensed under CC-BY (https://creativecommons.org/licenses/by/4.0/).

## Contributors
- Rahmat Heidari
- Matthew Amos
- Frederik Geth

## Acknowledgement
We express our gratitude to Dr. Felipe Arraño-Vargas and Dr. Georgios Konstantinou who released the original data used in this study, and for their ongoing support and advice on gaining insights from the data. The original S-NEM2300bus benchmark data and sotftware are licensed under BSD https://opensource.org/license/BSD-3-clause/ and as part of the following work:
- F. Arraño-Vargas and G. Konstantinou, "Modular Design and Real-Time Simulators Toward Power System Digital Twins Implementation," 
 	  in IEEE Transactions on Industrial Informatics, doi: 10.1109/TII.2022.3178713 (https://doi.org/10.1109/TII.2022.3178713
- F. Arraño-Vargas and G. Konstantinou, "Synthetic Grid Modeling for Real-Time Simulations," 2021 IEEE PES Innovative Smart Grid Technologies - Asia (ISGT Asia),
    2021, pp. 1-5, doi: 10.1109/ISGTAsia49270.2021.9715654 (https://doi.org/10.1109/ISGTAsia49270.2021.9715654_
