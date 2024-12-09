# Atomic Units

## 简介

a.u.是一套广泛应用于**原子物理学**中的单位制，在研究**电子**的相关性质时，应用尤为广泛。
常用的一套a.u.是Hartree atomic units。在Hartree a.u.中，以下4个基本物理常量被定义为1:

- 约化普朗克常数(reduced Planck constant) $\hbar$ , 也被称为 atomic unit of action
- 元电荷(elementary charge) $e$, 也被称为 atomic unit of charge
- 电子质量(electron mass) $m_e$, 也被称为 atomic unit of mass
- 玻尔半径(Bohr radius) $a_0$, 也被称为 atomic unit of length
- $4\pi \epsilon_0$, 这个与玻尔半径被定义为1是等效的。
  同时，可以推导出一个Hartree energy: $E_\mathrm{h}=\frac{\hbar^2}{m_e a_0^2}=\frac{e^2}{4 \pi \epsilon_0 a_0}$在 Hartree a.u.下为1!
  1 Hartree = 	27.211386245988 eV，正好等于H原子电子能量的中13.6eV的两倍。

## 方程表示中实际应用

- 无量纲的距离$\tilde{r}=r/a_0$, 其中$a_0$为玻尔半径。
- 两个电子间的库仑力$\frac{e^2}{4 \pi \epsilon_0 r}$简化为$\frac{1}{\tilde{r}}$，需要注意的是，这里的$\tilde{r}$是一个无量纲的距离
- 动能算符$-\frac{\hbar^2}{2m}\nabla^2$ 简化为$-\frac{1}{2\tilde{m}}\tilde{\nabla}^2$，其中$\tilde{m}=m/m_e$是一个无量纲质量，$m_e$为电子质量; $\tilde{\nabla}^2=a_0^2 \nabla^2$。

## Hartree atomic units

|        Atomic units        |  Name  |           Expression           |                        Value in SI units                        |                                    Other equivalents                                    |
| :------------------------: | :-----: | :----------------------------: | :--------------------------------------------------------------: | :-------------------------------------------------------------------------------------: |
|      **action**      |        |           $\hbar$           |         1.054571817...×10$^{−34}$ $\mathrm{J⋅s}$         |                                                                                        |
|      **charge**      |        |             $e$             |                  1.602176634×10$^{−19}$ C                  |                                                                                        |
|      **length**      |  bohr  |            $a_0$            |               5.29177210903(80)×10$^{−11}$ m               |                     $\hbar/m_e c \alpha$, 0.529177210903(80) Ang                     |
|      **energy**      | hartree |        $E_\mathrm{h}$        |              4.3597447222071(85)×10$^{−18}$ J              |              $2R_{\infty}hc$, $\alpha^2m_ec^2$, 27.211386245988(53) eV              |
|       **mass**       |        |            $m_e$            |               9.1093837015(28)×10$^{−31}$ kg               |                                                                                        |
|       **time**       |        |     $\hbar/E_\mathrm{h}$     |              2.4188843265857(47)×10$^{−17}$ s              |                                                                                        |
|     **velocity**     |        |   $a_0 E_\mathrm{h}/\hbar$   |       2.18769126364(33)×10$^6$ $\mathrm{m⋅s^{−1}}$       |                                      $\alpha c$                                      |
|           force           |        |      $E_\mathrm{h}/a_0$      |                8.2387234983(12)×10$^{−8}$ N                |                       82.387 nN, 51.421$\mathrm{eV·Ang^{−1}}$                       |
|   electric dipole moment   |        |            $ea_0$            |        8.4783536255(13)×10$^{−30}$ $\mathrm{C⋅m}$        |                                $\simeq$ 2.541746473 D                                |
|          momentum          |        |         $\hbar /a_0$         |  1.99285191410(30)×10$^{−24}$ $\mathrm{kg·m·s^{−1}}$  |                                                                                        |
|  1st hyperpolarizability  |        |       $e^3a_0^3/E_h^2$       | 3.2063613061(15)×10$^{−53}$ $\mathrm{C^3⋅m^3⋅J^{−2}}$ |                                                                                        |
|  2nd hyperpolarizability  |        |       $e^4a_0^4/E_h^3$       | 6.2353799905(38)×10$^{−65}$ $\mathrm{C^4⋅m^4⋅J^{−3}}$ |                                                                                        |
|       charge density       |        |          $e/a_0^3$          |     1.08120238457(49)×10$^{12}$ $\mathrm{C⋅m^{−3}}$     |                                                                                        |
|          current          |        |    $eE_\mathrm{h} /\hbar$    |               6.623618237510(13)×10$^{−3}$ A               |                                                                                        |
|       electric field       |        |     $E_\mathrm{h}/ea_0$     |     5.14220674763(78)×10$^{11}$ $\mathrm{V⋅m^{−1}}$     | 5.14220674763(78)$\mathrm{GV⋅cm^{−1}}$, 51.4220674763(78) $\mathrm{V⋅ Ang^{-1}}$ |
|  electric field gradient  |        |    $E_\mathrm{h}/ea_0^2$    |      9.7173624292(29)×10$^{21}$ $\mathrm{V⋅m^{−2}}$      |                                                                                        |
|  electric polarizability  |        |   $e^2a_0^2/E_\mathrm{h}$   | 1.64877727436(50)×10$^{−41}$ $\mathrm{C^2⋅m^2⋅J^{−1}}$ |                                                                                        |
|     electric potential     |        |       $E_\mathrm{h}/e$       |                      27.211386245988(53) V                      |                                                                                        |
| electric quadrupole moment |        |           $ea_0^2$           |       4.4865515246(14)×10$^{−40}$ $\mathrm{C⋅m^2}$       |                                                                                        |
|   magnetic dipole moment   |        |        $\hbar e/m_e$        |    1.85480201566(56)×10$^{−23}$ $\mathrm{J⋅T^{−1}}$    |                                   $2\mu_\mathrm{B}$                                   |
|   magnetic flux density   |        |        $\hbar/ea_0^2$        |                  2.35051756758(71)×10$^5$ T                  |                        $\simeq$ 2.35051756758(71)×10$^9$ G                        |
|      magnetizability      |        |       $e^2 a_0^2/m_e$       |     7.8910366008(48)×10$^{−29}$ $\mathrm{J⋅T^{−2}}$     |                                                                                        |
|        permittivity        |        |    $e^2/a_0E_\mathrm{h}$    |    1.11265005545(17)×10$^{−10}$ $\mathrm{F⋅m^{−1}}$    |                                   $4\pi \epsilon_0$                                   |
|          pressure          |        |     $E_\mathrm{h}/a_0^3$     |                2.9421015697(13)×10$^{13}$ Pa                |                                                                                        |
|         irradiance         |        | $E_\mathrm{h}^2/\hbar a_0^2$ |      6.4364099007(19)×10$^{19}$ $\mathrm{W⋅m^{−2}}$      |                                                                                        |

- $c$ is the [speed of light](https://en.wikipedia.org/wiki/Speed_of_light "Speed of light")
- $\epsilon_0$ is the [vacuum permittivity](https://en.wikipedia.org/wiki/Vacuum_permittivity "Vacuum permittivity")
- $R_{\infty}$ is the [Rydberg constant](https://en.wikipedia.org/wiki/Rydberg_constant "Rydberg constant") $=\frac{m_e e^4}{8\epsilon_0^2 h^3 c}=10973731.568160(21) \mathrm{m^{−1}}$
- $h$ is the [Planck constant](https://en.wikipedia.org/wiki/Planck_constant "Planck constant")
- $\alpha$ is the [fine-structure constant](https://en.wikipedia.org/wiki/Fine-structure_constant "Fine-structure constant") $=\frac{e^2}{2 \epsilon_0 hc}=\frac{e^2}{4 \pi \epsilon_0 \hbar c}=  0.0072973525693(11) \simeq \frac{1}{137.036}$
- $\mu_\mathrm{B}$ is the [Bohr magneton](https://en.wikipedia.org/wiki/Bohr_magneton "Bohr magneton")
