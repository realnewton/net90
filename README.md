# net90

Nuclear network for hydrodynamical simulations. Specially suited for Type Ia SN simulations. It includes 89 nuclear species + electrons.</BR>
Electron/positron captures are included only on free protons/neutrons. Temperature is coupled and solved jointly with the nuclear species with an implicit scheme.

<p align="center">
  <img src="https://raw.githubusercontent.com/realnewton/net90/main/Net90architecture.png" alt="Network architecture" width="1000"/>
</p>

# Compilation and usage
We provide a makefile that compiles `net90` and `testnuclear.f90`, along with all needed EOS. `testnuclear.f90` calls `net90` iteratively and mimicking different user-defined scenarios. This is equivalent to the calls to `net90` that a hydrodynamic code will do on a single fluid element. The objective of `testnuclear.f90` is just to show how `net90` should be called and perform toy scenarios. If the objective is to use `net90` coupled with a hydrodynamics code, only the folder `src/` is needed.
To compile simply do `make` and this will generate an executable `runnet90`.
Only a fortran compiler is needed. Our makefile assumes that `gfortran` is available.

## Calling `net90`
Only abundance-related magnitudes need to be imported from the main code:</br>
(D) Double precision, (I) Integer, (L) Logical
```
USE nuclear90_module, only:niso,ye,a,z
```
* niso      (I): nuclear species (=89)
* ye        (D): electron abundance
* a[1:niso] (D): array of atomic mass of each species
* z[1:niso] (D): array of atomic number of each species

Then, to call `net90`:
```
call net90(tempin,xin,rho,delta,screen,ecapture,tabulated,cv,val1,val2,dUedYein,theta,&
         & tempout,xout,sumtot,nucenergy,k_iter)
```

INPUT:
* tempin    (D): temperature
* xin       (D): mass fractions
* rho       (D): density (g/cm^3) of the current time-step
* delta     (D): time-step in seconds
* screen    (L): Activates/deactivates the screening corrections to the rates
* ecapture  (L): Activates/deactivates e-/e+ reactions
* tabulated (L): (T) e-/e+ rates interpolating on a table; (F) e-/e+ rates analytically
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
- Screening     (L): Activates/deactivates the screening corrections
- Adiabatic_exp (L): (T) an adiabatic expansion follows after NSE; (F) the process is isochoric
- Isotherm_test (L): (T) temperature is kept constant via an artificially high cv; (F) temperature is updated by net90
- electron_capt (L): Activates/deactivates e-/e+ reactions
- tabul_rates   (L): (T) e-/e+ rates interpolating on a table; (F) e-/e+ rates analytically
- Theta         (D): Chooses between explicit (0), implicit (1), or mixed. Default value is 0.7d0
- type_eos      (I): Selects the EOS. (1) relativistic electrons (Nadyozhin 1974), ions (Bravo & Garcia-Senz 1992), and radiation; (2) helmholtz EOS (Cox & Giuli chapter 24; Timmes & Swesty ApJ 1999)

# Nuclear species included 
Corresponding index in the mass fraction array

| index | Element | index | Element | index | Element | index | Element | index | Element |
|--|---|--|---|--|---|--|---|--|---|
|1|p|21|26Al|41|39Ca|61|46V|81|59Zn|
|2|n|22|26Mg|42|37Ar|62|46Ti|82|57Ni|
|3|4He|23|29P|43|39K|63|49Mn|83|59Cu|
|4|12C|24|32S|44|38Ca|64|52Fe|84|58Zn|
|5|16O|25|31S|45|38K|65|51Fe|85|58Cu|
|6|20Ne|26|29Si|46|38Ar|66|49Cr|86|58Ni|
|7|21Na|27|31P|47|41Sc|67|51Mn|87|57Co|
|8|24Mg|28|30S|48|44Ti|68|50Fe|88|59Co|
|9|23Mg|29|30P|49|43Ti|69|50Mn|89|13N|
|10|21Ne|30|30Si|50|41Ca|70|50Cr|
|11|23Na|31|33Cl|51|43Sc|71|53Co|
|12|22Mg|32|36Ar|52|42Ti|72|56Ni|
|13|22Na|33|35Ar|53|42Sc|73|55Ni|
|14|22Ne|34|33S|54|42Ca|74|53Fe|
|15|25Al|35|35Cl|55|45V|75|55Co|
|16|28Si|36|34Ar|56|48Cr|76|54Ni|
|17|27Si|37|34Cl|57|47Cr|77|54Co|
|18|25Mg|38|34S|58|45Ti|78|54Fe|
|19|27Al|39|37K|59|47V|79|57Cu|
|20|26Si|40|40Ca|60|46Cr|80|60Zn|

Corresponding index in the reaction rate arrays

| index | Reaction | index | Reaction | index | Reaction | index | Reaction 
|--|---|--|---|--|---|--|---|
|1|12C(12C,a)20Ne [^1]|46|37K(n,g)38K|91|23Na(p,g)24Mg|136|52Fe(a,g)56Ni
|2|12C(16O,a)24Mg|47|37Ar(n,g)38Ar|92|25Mg(p,g)26Al|137|56Ni(a,g)60Zn|
|3|16O(16O,a)28Si [^2]|48|40Ca(p,g)41Sc|93|26Al(p,g)27Si|138|20Ne(a,n)23Mg|
|4|-|49|43Ti(n,g)44Ti|94|26Mg(p,g)27Al|139|24Mg(a,n)27Si|
|5|a(2a,g)12C [^3]|50|42Ti(n,g)43Ti|95|27Al(p,g)28Si|140|28Si(a,n)31S|
|6|12C(a,g)16O|51|40Ca(n,g)41Ca|96|29Si(p,g)30P|141|32S(a,n)35Ar|
|7|16O(a,g)20Ne|52|42Sc(n,g)43Sc|97|30P(p,g)31S|142|36Ar(a,n)39Ca|
|8|20Ne(p,g)21Na|53|41Sc(p,g)42Ti|98|30Si(p,g)31P|143|40Ca(a,n)43Ti|
|9|23Mg(n,g)24Mg|54|41Sc(n,g)42Sc|99|31P(p,g)32S|144|44Ti(a,n)47Cr|
|10|22Mg(n,g)23Mg|55|41Ca(n,g)42Ca|100|33S(p,g)34Cl|145|48Cr(a,n)51Fe|
|11|20Ne(n,g)21Ne|56|44Ti(p,g)45V|101|34Cl(p,g)35Ar|146|52Fe(a,n)55Ni|
|12|22Na(n,g)23Na|57|47Cr(n,g)48Cr|102|34S(p,g)35Cl|147|56Ni(a,n)59Zn|
|13|21Na(p,g)22Mg|58|46Cr(n,g)47Cr|103|35Cl(p,g)36Ar|148|20Ne(a,p)23Na|
|14|21Na(n,g)22Na|59|44Ti(n,g)45Ti|104|37Ar(p,g)38K|149|24Mg(a,p)27Al|
|15|21Ne(n,g)22Ne|60|46V(n,g)47V|105|38K(p,g)39Ca|150|28Si(a,p)31P|
|16|24Mg(p,g)25Al|61|45V(p,g)46Cr|106|38Ar(p,g)39K|151|32S(a,p)35Cl|
|17|27Si(n,g)28Si|62|45V(n,g)46V|107|39K(p,g)40Ca|152|36Ar(a,p)39K|
|18|26Si(n,g)27Si|63|45Ti(n,g)46Ti|108|41Ca(p,g)42Sc|153|40ca(a,p)43Sc|
|19|24Mg(n,g)25Mg|64|48Cr(p,g)49Mn|109|42Sc(p,g)43Ti|154|44Ti(a,p)47V|
|20|26Al(n,g)27Al|65|51Fe(n,g)52Fe|110|42Ca(p,g)43Sc|155|48Cr(a,p)51Mn|
|21|25Al(p,g)26Si|66|50Fe(n,g)51Fe|111|43Sc(p,g)44Ti|156|52Fe(a,p)55Co|
|22|25Al(n,g)26Al|67|48Cr(n,g)49Cr|112|45Ti(p,g)46V|157|56Ni(a,p)59Cu|
|23|25Mg(n,g)26Mg|68|50Mn(n,g)51Mn|113|46V(p,g)47Cr|158|54Fe(a,p)57Co|
|24|28Si(p,g)29P|69|49Mn(p,g)50Fe|114|46Ti(p,g)47V|159|57Co(n,g)58Co|
|25|31S(n,g)32S|70|49Mn(n,g)50Mn|115|47V(p,g)48Cr|160|12C(p,g)13N|
|26|30S(n,g)31S|71|49Cr(n,g)50Cr|116|49Cr(p,g)50Mn|161|13N(a,p)16O|
|27|28Si(n,g)29Si|72|52Fe(p,g)53Co|117|50Mn(p,g)51Fe|-|p(e-,nu)n [^4]|
|28|30P(n,g)31P|73|55Ni(n,g)56Ni|118|50Cr(p,g)51Mn|-|n(e+,anu)p [^4]|
|29|29P(p,g)30S|74|54Ni(n,g)55Ni|119|51Mn(p,g)52Fe|
|30|29P(n,g)30P|75|52Fe(n,g)53Fe|120|53Fe(p,g)54Co|
|31|29Si(n,g)30Si|76|54Co(n,g)55Co|121|54Co(p,g)55Ni|
|32|32S(p,g)33Cl|77|53Co(p,g)54Ni|122|54Fe(p,g)55Co|
|33|35Ar(n,g)36Ar|78|53Co(n,g)54Co|123|55Co(p,g)56Ni|
|34|34Ar(n,g)35Ar|79|53Fe(n,g)54Fe|124|57Ni(p,g)58Cu|
|35|32S(n,g)33S|80|56Ni(p,g)57Cu|125|58Cu(p,g)59Zn|
|36|34Cl(n,g)35Cl|81|59Zn(n,g)60Zn|126|58Ni(p,g)59Cu|
|37|33Cl(p,g)34Ar|82|58Zn(n,g)59Zn|127|59Cu(p,g)60Zn|
|38|33Cl(n,g)34Cl|83|56Ni(n,g)57Ni|128|20Ne(a,g)24Mg|
|39|33S(n,g)34S|84|58Cu(n,g)59Cu|129|24Mg(a,g)28Si|
|40|36Ar(p,g)37K|85|57Cu(p,g)58Zn|130|28Si(a,g)32S|
|41|39Ca(n,g)40Ca|86|57Cu(n,g)58Cu|131|32S(a,g)36Ar|
|42|38Ca(n,g)39Ca|87|57Ni(n,g)58Ni|132|36Ar(a,g)40ca|
|43|36Ar(n,g)37Ar|88|21Ne(p,g)22Na|133|40ca(a,g)44Ti|
|44|38K(n,g)39K|89|22Na(p,g)23Mg|134|44Ti(a,g)48Cr|
|45|37K(p,g)38Ca|90|22Ne(p,g)23Na|135|48Cr(a,g)52Fe|

[^1]: This reaction is split into two branches: 50% 12C(12C,a)20Ne, 50% 12C(12C,p)23Na
[^2]: This reaction is split into two branches: 40% 16O(16O,a)28Si, 60% 16O(16O,p)31P
[^3]: The triple alpha reaction is calculated in two steps via 8Be
[^4]: electron and positron reactions are stored independently, hence they don't have an associated index

# References
- Garcia-Senz, Cabezon, Reichert, Sanz, Escartin, Psaltis, Arcones, Thielemann; (2024) Upcoming
- Sanz, Cabezon, Garcia-Senz; [NIC-XVI, 260 (2022)](https://ui.adsabs.harvard.edu/abs/2022EPJWC.26011036S/abstract)
- Cabezon, Garcia-Senz, Bravo; [ApJS 151 (2004)](https://ui.adsabs.harvard.edu/abs/2004ApJS..151..345C/abstract)
- Garcia-Senz, Cabezon; [NPA, 718 (2003)](https://ui.adsabs.harvard.edu/abs/2003NuPhA.718..566G/abstract)

# Authors (in alphabetical order)
Ruben M. Cabezon: (sciCORE) University of Basel</br>
Domingo Garcia-Senz: Polytechnic University of Catalonia</br>
Axel Sanz: Polytechnic University of Catalonia
