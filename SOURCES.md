# PVT Lab — textbook basis, equations and data provenance

This project is an educational **pure-component** real-fluid calculator. Its primary course reference is:

J. M. Smith, H. C. Van Ness, M. M. Abbott and M. T. Swihart, *Introduction to Chemical Engineering Thermodynamics*, 8th edition, McGraw-Hill Education (2018).

The implementation follows the organization and notation of the textbook. The equation numbers and printed page references below make it easier to connect the program to the course material.

## 1. Cubic equations of state

The backend uses the generic cubic equation (textbook Eq. 3.41):

\[
P=\frac{RT}{V-b}-\frac{a(T)}{(V+\epsilon b)(V+\sigma b)}
=\frac{RT}{V-b}-\frac{a(T)}{V^2+ubV+wb^2}
\]

where \(u=\epsilon+\sigma\), \(w=\epsilon\sigma\), and the attraction term is \(a(T)=a_c\alpha(T_r,\omega)\). The critical-property parameterization is (Eqs. 3.44–3.45):

\[
a_c=\Psi\frac{R^2T_c^2}{P_c},\qquad b=\Omega\frac{RT_c}{P_c},\qquad T_r=T/T_c
\]

The exact constants and alpha functions implemented are from **Table 3.1, p. 100**:

| EOS | \(u\) | \(w\) | \(\Omega\) | \(\Psi\) | \(\alpha(T_r,\omega)\) |
|---|---:|---:|---:|---:|---|
| van der Waals | 0 | 0 | 0.12500 | 0.421875 | 1 |
| Redlich–Kwong | 1 | 0 | 0.08664 | 0.42748 | \(T_r^{-1/2}\) |
| Soave–Redlich–Kwong | 1 | 0 | 0.08664 | 0.42748 | \([1+m(1-\sqrt{T_r})]^2\), \(m=0.480+1.574\omega-0.176\omega^2\) |
| Peng–Robinson | 2 | −1 | 0.07780 | 0.45724 | \([1+m(1-\sqrt{T_r})]^2\), \(m=0.37464+1.54226\omega-0.26992\omega^2\) |

The solver forms the cubic polynomial in compressibility factor \(Z=PV/(RT)\). It keeps real roots satisfying \(Z>B\), where \(B=bP/(RT)\), equivalent to a real molar volume greater than \(b\). The smallest and largest roots are labelled liquid-like and vapour-like. If three real roots exist, the intermediate mathematical root is displayed for transparency but is excluded from the stable-phase comparison because the textbook identifies the middle segment as physically unstable in the two-phase region. The lowest-fugacity option compares the EOS fugacity coefficients of the outer roots (Sec. 3.6, pp. 96–97; Example 3.9, pp. 100–101).

## 2. Residual enthalpy and residual entropy

The textbook defines a residual property as the difference between the real-fluid property and the ideal-gas-state property at the same temperature and pressure (Eq. 6.41, p. 221):

\[
M^R=M-M^{ig},\qquad M=M^{ig}+M^R
\]

For cubic equations of state, the backend uses the generic-cubic expressions in **Sec. 13.6, Eqs. 13.74–13.76, pp. 488–489**. Define:

\[
B=\frac{bP}{RT},\qquad q=\frac{a(T)}{bRT},\qquad
I=\frac{1}{\sigma-\epsilon}\ln\left(\frac{Z+\sigma B}{Z+\epsilon B}\right)
\]

For the van der Waals limiting case, \(\epsilon=\sigma=0\), the limit is \(I=B/Z\). Then:

\[
\frac{H^R}{RT}=Z-1+\left(\frac{d\ln\alpha}{d\ln T_r}-1\right)qI
\]

\[
\frac{S^R}{R}=\ln(Z-B)+\frac{d\ln\alpha}{d\ln T_r}qI
\]

The implementation evaluates these expressions for each selected EOS, including the temperature derivative of its alpha function. These are molar residual/departure properties.

## 3. Ideal-gas enthalpy and entropy contributions

The textbook combines ideal-gas and residual contributions in **Eqs. 6.50 and 6.51, p. 222**:

\[
H=H_0^{ig}+\int_{T_0}^{T}C_P^{ig}\,dT+H^R
\]

\[
S=S_0^{ig}+\int_{T_0}^{T}\frac{C_P^{ig}}{T}\,dT-R\ln\left(\frac{P}{P_0}\right)+S^R
\]

Unlike the earlier draft's constant-\(C_p\) approximation, this version uses the book's temperature-dependent ideal-gas heat-capacity correlation from **Appendix C, Table C.1, p. 656**:

\[
\frac{C_P^{ig}}{R}=A+BT+CT^2+DT^{-2}
\]

The coefficients are stored in the gas database in SI-compatible polynomial form. The ideal-gas integrals are evaluated analytically:

\[
\frac{\Delta H^{ig}}{R}=A(T-T_0)+\frac{B}{2}(T^2-T_0^2)+\frac{C}{3}(T^3-T_0^3)-D\left(\frac{1}{T}-\frac{1}{T_0}\right)
\]

\[
\frac{\Delta S^{ig}}{R}=A\ln\frac{T}{T_0}+B(T-T_0)+\frac{C}{2}(T^2-T_0^2)-\frac{D}{2}\left(\frac{1}{T^2}-\frac{1}{T_0^2}\right)-\ln\frac{P}{P_0}
\]

The table states its correlation range for each gas: from 298 K to the listed maximum temperature \(T_{max}\). The app enforces the reference temperature to be within that range. If a calculation temperature falls outside it, the app still evaluates the polynomial but displays a warning that the ideal-gas contribution is extrapolated.

The defaults \(H_0^{ig}=0\) and \(S_0^{ig}=0\) at \(T_0=298.15\,K\), \(P_0=1\,bar\) make the displayed values **reference-relative**. They are not standard enthalpies of formation or absolute entropies. Users may enter other consistent reference values.

## 4. Pure-component fugacity

For a pure species, the textbook defines fugacity through Eq. 10.31 and defines the fugacity coefficient by (Eq. 10.33, pp. 366–367):

\[
\phi=\frac{f}{P},\qquad \ln\phi=\frac{G^R}{RT}
\]

The backend uses the cubic-EOS residual Gibbs-energy relation in **Eq. 13.85, p. 492**:

\[
\ln\phi=Z-1-\ln(Z-B)-qI,\qquad f=\phi P
\]

The displayed fugacity is in bar and the fugacity coefficient is dimensionless.

## 5. Selected pure components and property-data source

The ten gases are all listed in **Appendix B, Table B.1, pp. 650–652**. Their molecular weights, acentric factors, critical temperatures, critical pressures, critical compressibility factors, critical molar volumes and normal boiling points are taken from that table (rounded only as displayed there).

The final set is:

1. Methane (CH₄)
2. Ethane (C₂H₆)
3. Propane (C₃H₈)
4. n-Butane (C₄H₁₀)
5. n-Pentane (C₅H₁₂)
6. Nitrogen (N₂)
7. Oxygen (O₂)
8. Carbon dioxide (CO₂)
9. Ammonia (NH₃)
10. Water vapour (H₂O)

The five consecutive alkanes let the user compare how critical properties and non-ideality change with molecular size. Nitrogen and oxygen provide small atmospheric-gas examples; carbon dioxide adds a common industrial fluid; ammonia and water add polar/associating examples that are useful for discussing model limitations.

**Hydrogen was removed from the original component list.** Although Appendix B lists hydrogen's critical properties, Sec. 3.7 (pp. 104–105) cautions that quantum gases such as hydrogen, helium and neon do not follow ordinary corresponding-states behaviour and recommends temperature-dependent effective critical parameters for hydrogen. The simple fixed-critical-constant cubic implementation in this project did not implement that correction, so n-pentane is a more defensible textbook-aligned selection for this version.

The ideal-gas Cp polynomial coefficients for the selected components are from **Appendix C, Table C.1, p. 656**. The program converts the table's scaled columns before storing them: table entries for \(10^3B\) are multiplied by \(10^{-3}\), entries for \(10^6C\) by \(10^{-6}\), and entries for \(10^{-5}D\) by \(10^5\).

## 6. Model scope and limitations

- Pure-component calculations only. Mixture mixing rules, binary interaction parameters, and mixture fugacity are not implemented.
- Cubic EOS methods are approximate. Water and ammonia are polar/associating fluids; the project does not include an association term or a specialized water/ammonia model. Treat those results as educational estimates rather than high-accuracy property predictions.
- The Cp polynomial is taken from the textbook's Table C.1 and is stated for 298 K to each gas's \(T_{max}\). Calculations outside that range are extrapolations and are flagged.
- The P–V curve is the selected EOS isotherm. It is not an experimental phase envelope and does not apply a Maxwell equal-area construction or solve the saturation pressure.
- Lowest-fugacity root selection is a cubic-EOS phase preference at the specified state, not a complete phase-equilibrium calculation.

## 7. Cross-check references

These can be used as independent checks, not as substitutes for the course textbook:

- NIST Chemistry WebBook, Standard Reference Database 69: https://webbook.nist.gov/chemistry/
- NIST teqp documentation, generalized cubic EOS: https://pages.nist.gov/teqp-docs/en/latest/models/cubics.html
- Bell, I. H. & Jäger, A. (2016), “Helmholtz Energy Transformations of Common Cubic Equations of State for Use with Pure Fluids and Mixtures,” *Journal of Research of the NIST*, 121, 238–263. https://doi.org/10.6028/jres.121.011
