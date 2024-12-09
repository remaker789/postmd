# Common unit transformation

## Electric

$$
\mathrm{elementary \ charge}: \quad 1 \, e = 1.602176634\times 10^{-19} \, \mathrm{C}
$$

$$
1 \mathrm{S} = 1 \mathrm{\Omega^{-1}}= 1 \mathrm{A/V} = 1 \mathrm{s\cdot C^2/(kg \cdot m^2)}
$$

$$
\mathrm{unit \ of \ dipole :} \quad 1 \, \mathrm{Debye} = 3.33564 \times 10^{-30} \, \mathrm{C\cdot m}
$$

## Energy

$$
\mathrm{Avogadro \ constant:} \quad N_A = 6.02214076 \times10^{23} \, \mathrm{mol^{-1}}
$$

$$
\mathrm{Boltzmann \ constant:} \quad k_\mathrm{B} = 1.380649 \times 10^{-23} \, \mathrm{J/K}
$$

$$
\mathrm{molar \ gas \ constant:} \quad R= 8.314462618 \, \mathrm{J\cdot K^{-1}\cdot mol^{-1}}
$$

$$
k_\mathrm{B} = R/N_A
$$

$$
1 \, \mathrm{cal} = 4.184 \, \mathrm{J}
$$

$$
1\, \mathrm{kcal/mol} = 4.184\, \mathrm{kJ/mol} = 0.0434\ \mathrm{eV}=1.688\, \mathrm{k_B T}(\mathrm{at \ 298.15\, K})
$$

$$
1 \ \mathrm{Ry} = 1/2 \ \mathrm{Ha}= 13.6057 \ \mathrm{eV}
$$

$\mathrm{Ry}$的全称为Rydberg(里德伯)。$\mathrm{Ha}$的全称为Hartree，是atomic units中的能量单位。

## Length

$$
1 \, \mathrm{Bohr} = 0.529177210909 \, \mathrm{Angstrom}
$$

## Time

$$
1 \, a.u. = 0.024188843265857 \, \mathrm{fs}
$$

## Spectrum

在光谱学中，波数 $\tilde{\nu}$ 定义为:

$$
\tilde{\nu} = \frac{1}{\lambda}
$$

其中$\lambda$为电磁辐射在真空中的波长。波数的单位为$\mathrm{m}^{-1}$，但是在光谱学中，经常使用$\mathrm{cm^{-1}}$。

$$
E = \hbar \omega = h \nu = h \frac{c}{\lambda} = hc\tilde{\nu}
$$

其中:

- $E$: 能量， in $\mathrm{J}$
- $h$: 普朗克常数，in $\mathrm{J \cdot s}$
- $\omega=2\pi \nu$: 电磁波角频率，in $\mathrm{rad/s}$
- $\nu$: 电磁波频率，in $\mathrm{Hz}$(or $\mathrm{s^{-1}}$)
- $c$: 真空光速，in $\mathrm{m/s}$
- $\lambda$: 电磁波波长，in $\mathrm{m}$
- $\tilde{\nu}$: 电磁波波数，in $\mathrm{m^{-1}}$

$$
k_\mathrm{B}T= 208.51 \, \mathrm{cm^{-1}} \quad \mathrm{for} \  T=300K
$$

```python
import scipy.constants as C
k = C.Boltzmann
h = C.h
T = 300
c = C.speed_of_light
print(f"The wave number corresponding to temperature of {T} K is {k*T/h/c/100} cm^-1")
```
