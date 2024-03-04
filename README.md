# net90

Nuclear network for hydrodynamical simulations. Specially suited for Type Ia SN simulations.
It includes 89 nuclear species + electrons
Electron/positron captures are included only on protons/neutrons.
Temperature is coupled and solved jointly with the nuclear species.

<p align="center">
  <img src="https://raw.githubusercontent.com/realnewton/net90/main/Net90architecture.png" alt="Network architecture" width="1000"/>
</p>

# Compilation and usage
We provide a makefile that compiles `net90` and `testnuclear.f90`, along with all needed EOS. `testnuclear.f90` calls `net90` iteratively and mimicking different user-defined scenarios. This is equivalent to the calls to `net90` that a hydrodynamic code will do on a single fluid element. The objective of `testnuclear.f90` is just to show how `net90` should be called and perform toy scenarios. If the objective is to use `net90` coupled with a hydrodynamics code, only the folder `src/` is needed.
To compile simply do `make` and this will generate an executable `runnet90`.
Only a fortran compiler is needed. Our makefile assumes that `gfortran` is available.

## Calling `net90`
From the main code:
(D) Double precision, (I) Integer, (L) Logical
```
call net90(tempin,xin,rho,delta,screen,cv,val1,val2,dUedYein,theta,tempout,xout,sumtot,nucenergy,k_iter)
```

INPUT:
* tempin    (D): temperature
* xin       (D): mass fractions
* rho       (D): density (g/cm^3) of the current time-step
* delta     (D): time-step in seconds
* screen    (L): Turns the screening corrections to the rates on/off
* cv        (D): heat capacity (from EOS)
* val1      (D): (dp/dTemp)*(deltarho/rho**2)
* val2      (D): (dp/dTemp)*(deltarho/rho**2)*Temp
* dUedYein  (D): dUe/dYe (from EOS)
* theta     (D): Chooses between explicit (0), implicit (1), or mixed. Recommended value is 0.7d0

OUTPUT:
* tempout   (D): temperature
* xout      (D): mass fractions
* sumtot    (D): total mass fraction (=1 if conserved)
* nucenergy (D): total nuclear energy generation rate (erg/s)
* k_iter    (I): Total number of Newton-Raphson iterations performed

## Description of `parameters.in` for toy scenarios

- Mass fractions of each species (D)
- Initial_temp  (D): initial temperature in K
- Initial_dens  (D): initial density in g/cm^3
- Screening     (L): Activates the screening corrections
- Adiabatic_exp (L): (T) an adiabatic expansion follows after NSE; (F) the process is isochoric
- Isotherm_test (L): (T) temperature is kept constant via an artificially high cv; (F) temperature is updated by net90
- electron_capt (L): Activates e-/e+ reactions
- tabul_rates   (L): (T) e-/e+ rates interpolating on a table; (F) e-/e+ rates analytically
- Theta         (D): Chooses between explicit (0), implicit (1), or mixed. Default value is 0.7d0
- C(c_p)Na_prob (D): Probability of 12C+12C->23Na+p branch. Default value is 0.5d0
- C(c_a)Ne_prob (D): Probability of 12C+12C->20Ne+4He branch. Default value is 0.5d0
- O(O_p)P_prob  (D): Probability of 16O+16O->31P+p branch. Default value is 0.6d0
- O(O_a)S_prob  (D): Probability of 16O+16O->28S+a branch. Default value is 0.4d0
- NR_max_iter   (I): Maximum number of Newton-Raphson iterations in net90 per timestep. Default value is 10
- NR_limit      (D): Convergence limit for NR iterations. Default value is 1.d-7  
- type_eos      (I): Selects the EOS. (1) relativistic electrons (Nadyozhin 1974), ions (Bravo & Garcia-Senz 1992), and radiation; (2) helmholtz EOS (Cox & Giuli chapter 24; Timmes & Swesty ApJ 1999)

# Nuclear species included 
... and their corresponding index in the mass fraction array

| index | Element | index | Element | index | Element | index | Element | index | Element |
|--|---|--|---|--|---|--|---|--|---|
|1|p|21|Al26|41|Ca39|61|V46|81|Zn59|
|2|n|22|Mg26|42|Ar37|62|Ti46|82|Ni57|
|3|He4|23|P29|43|K39|63|Mn49|83|Cu59|
|4|C12|24|S32|44|Ca38|64|Fe52|84|Zn58|
|5|O16|25|S31|45|K38|65|Fe51|85|Cu58|
|6|Ne20|26|Si29|46|Ar38|66|Cr49|86|Ni58|
|7|Na21|27|P31|47|Sc41|67|Mn51|87|Co57|
|8|Mg24|28|S30|48|Ti44|68|Fe50|88|Co58|
|9|Mg23|29|P30|49|Ti43|69|Mn50|89|N13
|10|Ne21|30|Si30|50|Ca41|70|Cr50|
|11|Na23|31|Cl33|51|Sc43|71|Co53|
|12|Mg22|32|Ar36|52|Ti42|72|Ni56|
|13|Na22|33|Ar35|53|Sc42|73|Ni55|
|14|Ne22|34|S33|54|Ca42|74|Fe53|
|15|Al25|35|Cl35|55|V45|75|Co55|
|16|Si28|36|Ar34|56|Cr48|76|Ni54|
|17|Si27|37|Cl34|57|Cr47|77|Co54|
|18|Mg25|38|S34|58|Ti45|78|Fe54|
|19|Al27|39|K37|59|V47|79|Cu57|
|20|Si26|40|Ca40|60|Cr46|80|Zn60|

# References
Cabezon et al. ApJS 151 (2004)
Garcia-Senz el al. (2024) Upcoming

# Authors
Ruben M. Cabezon & D. García-Senz (2023)
