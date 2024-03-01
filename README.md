# net90

Nuclear network for hydrodynamical simulations. Specially suited for Type Ia SN simulations.
It includes 89 nuclear species + electrons
Electron/positron captures are included only on protons/neutrons.
Temperature is coupled and solved jointly with the nuclear species.

<p align="center">
  <img src="[https://raw.githubusercontent.com/unibas-dmi-hpc/SPH-EXA/develop/docs/artwork/SPH-EXA_logo.png](https://raw.githubusercontent.com/realnewton/net90/main/Net90architecture.png)" alt="Network architecture" width="300"/>
</p>

# Nuclear species included 
and its corresponding index in the mass fraction array
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
