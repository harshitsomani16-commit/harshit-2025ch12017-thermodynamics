# PVT Lab — Real-Gas Thermodynamics

**Project by Harshit Somani**

A Flask + HTML/CSS/JavaScript application for pure-component PVT calculations using four cubic equations of state. It calculates compressibility factor, molar volume, pure-component fugacity coefficient and fugacity, plus molar enthalpy and entropy using textbook Cp(T) ideal-gas contributions and EOS residual functions.

## Features

- Ten pure components from Smith et al. (2018), Appendix B, Table B.1: methane, ethane, propane, n-butane, n-pentane, nitrogen, oxygen, carbon dioxide, ammonia and water vapour.
- Four cubic EOS models: van der Waals, Redlich–Kwong, Soave–Redlich–Kwong and Peng–Robinson.
- Three calculation modes: solve V from P,T; solve P from T,V; solve T from P,V.
- Cubic compressibility-factor root calculation, with vapour-like, liquid-like and lowest-fugacity stable-outer-root selections when solving volume. If three real roots occur, the intermediate mathematical root is shown but excluded from stable-phase comparison.
- Pure-component fugacity: f = φP.
- Molar enthalpy and entropy: ideal-gas contribution integrated from the Cp(T) polynomial in Appendix C, Table C.1, plus cubic-EOS residual/departure functions.
- Dynamic component colour palette, side information panels, selected-EOS parameter panel and interactive P–V isotherm.

## Primary textbook basis

J. M. Smith, H. C. Van Ness, M. M. Abbott and M. T. Swihart, *Introduction to Chemical Engineering Thermodynamics*, 8th ed., McGraw-Hill Education (2018):

- Cubic EOS: Sec. 3.6, Table 3.1, pp. 95–101.
- Residual properties: Sec. 6.2, Eqs. 6.50–6.51, pp. 221–222.
- Pure-component fugacity: Sec. 10.5, Eqs. 10.31–10.33, pp. 366–367.
- Cubic-EOS residual enthalpy and entropy: Sec. 13.6, Eqs. 13.74–13.76, pp. 488–489.
- Critical properties: Appendix B, Table B.1, pp. 650–652.
- Ideal-gas heat capacities: Appendix C, Table C.1, p. 656.

See [SOURCES.md](SOURCES.md) for derivations, coefficient conversions, assumptions and limitations.

## Project structure

```text
pvt-thermo-lab/
├── app.py
├── requirements.txt
├── README.md
├── SOURCES.md
├── test_calculations.py
├── templates/
│   └── index.html
└── static/
    ├── styles.css
    └── app.js
```

## Run locally in PyCharm (Windows)

1. Extract the ZIP file.
2. Open the extracted `pvt-thermo-lab` folder in PyCharm.
3. Configure a Python virtual environment (Python 3.10+ recommended).
4. Open PyCharm's Terminal and run:

   ```powershell
   pip install -r requirements.txt
   ```

5. Right-click `app.py` and choose **Run 'app'** (or run `python app.py` in the terminal).
6. Open http://127.0.0.1:5000 in Chrome or another browser.
7. Optional regression checks: run `python -m unittest test_calculations.py`. The tests include the textbook n-butane/RK worked example (Example 13.5).

Keep the Flask run process active while using the website. Press Ctrl+C in the terminal to stop it.

## Calculation assumptions

- The backend uses SI units internally; the interface uses bar, K and L/mol.
- Default ideal-gas reference state: 298.15 K and 1 bar; reference enthalpy and entropy default to zero. Results are relative to these user-selected reference values, not standard formation enthalpy or absolute entropy.
- Cp/R is evaluated as A + BT + CT² + D/T² from Appendix C, Table C.1. The book states a range of 298 K to each component's Tmax. The app warns when a state-temperature calculation extrapolates outside that range.
- The pure-component fugacity coefficient follows ln φ = GR/(RT) and fugacity is f = φP.
- Water and ammonia are included for teaching comparison, but their polarity/association can make simple cubic-EOS estimates less accurate.
- Hydrogen was not retained because the textbook notes special effective critical-parameter treatment for hydrogen among quantum gases; that correction is not implemented in this app.

## API endpoints

- `GET /api/gases` — component database and model names.
- `POST /api/calculate` — solve state and return PVT, φ, f, H and S.
- `POST /api/curve` — generate a P–V isotherm for the selected component and EOS.
