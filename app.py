"""PVT & Thermodynamic Properties Lab — Harshit Somani.

Pure-component PVT, fugacity, molar enthalpy and molar entropy calculator.
Cubic-EOS equations and property methods follow Smith, Van Ness, Abbott & Swihart,
Introduction to Chemical Engineering Thermodynamics, 8th ed. (2018).

Backend units: Pa, m^3/mol, J/mol, J/(mol K), K.
UI/API units: bar, L/mol, J/mol, J/(mol K), K.
"""
from __future__ import annotations

import math
from typing import Any

import numpy as np
from flask import Flask, jsonify, render_template, request

app = Flask(__name__)

R = 8.31446261815324  # J mol-1 K-1
T_REF_DEFAULT = 298.15  # K
P_REF_DEFAULT_BAR = 1.0  # bar

# Pure-component critical constants and acentric factors are taken from
# Table B.1, Appendix B, Smith et al. (8th ed., pp. 650–652).
# Ideal-gas Cp/R coefficients are from Table C.1, Appendix C (p. 656).
# Cp/R = A + B*T + C*T**2 + D*T**(-2), with T in K.
# The book table reports 10^3 B, 10^6 C and 10^-5 D; values below are converted.
# The Cp polynomial's stated validity interval begins at 298 K and ends at Tmax.
GASES: dict[str, dict[str, Any]] = {
    "methane": {
        "name": "Methane", "formula": "CH₄", "cas": "74-82-8", "family": "Alkane",
        "Tc": 190.6, "Pc_bar": 45.99, "omega": 0.012, "Zc": 0.286,
        "Vc_cm3": 98.6, "Tb": 111.4, "mw": 16.043,
        "cp_coeffs": {"A": 1.702, "B": 9.081e-3, "C": -2.164e-6, "D": 0.0},
        "cp_tmax": 1500.0, "cp298_R": 4.217,
        "accent": "#45cbe8", "theme": "blue-flame cyan",
        "visual_note": "Blue-flame palette inspiration; methane gas itself is colourless.",
    },
    "ethane": {
        "name": "Ethane", "formula": "C₂H₆", "cas": "74-84-0", "family": "Alkane",
        "Tc": 305.3, "Pc_bar": 48.72, "omega": 0.100, "Zc": 0.279,
        "Vc_cm3": 145.5, "Tb": 184.6, "mw": 30.070,
        "cp_coeffs": {"A": 1.131, "B": 19.225e-3, "C": -5.561e-6, "D": 0.0},
        "cp_tmax": 1500.0, "cp298_R": 6.369,
        "accent": "#8bd8ef", "theme": "glacier blue",
        "visual_note": "A cryogenic-blue palette; ethane gas itself is colourless.",
    },
    "propane": {
        "name": "Propane", "formula": "C₃H₈", "cas": "74-98-6", "family": "Alkane",
        "Tc": 369.8, "Pc_bar": 42.48, "omega": 0.152, "Zc": 0.276,
        "Vc_cm3": 200.0, "Tb": 231.1, "mw": 44.097,
        "cp_coeffs": {"A": 1.213, "B": 28.785e-3, "C": -8.824e-6, "D": 0.0},
        "cp_tmax": 1500.0, "cp298_R": 9.011,
        "accent": "#70b5ff", "theme": "torch blue",
        "visual_note": "Inspired by a propane burner flame; propane gas is colourless.",
    },
    "n-butane": {
        "name": "n-Butane", "formula": "C₄H₁₀", "cas": "106-97-8", "family": "Alkane",
        "Tc": 425.1, "Pc_bar": 37.96, "omega": 0.200, "Zc": 0.274,
        "Vc_cm3": 255.0, "Tb": 272.7, "mw": 58.123,
        "cp_coeffs": {"A": 1.935, "B": 36.915e-3, "C": -11.402e-6, "D": 0.0},
        "cp_tmax": 1500.0, "cp298_R": 11.928,
        "accent": "#ffad75", "theme": "warm flame amber",
        "visual_note": "Warm flame-inspired palette; n-butane gas itself is colourless.",
    },
    "n-pentane": {
        "name": "n-Pentane", "formula": "C₅H₁₂", "cas": "109-66-0", "family": "Alkane",
        "Tc": 469.7, "Pc_bar": 33.70, "omega": 0.252, "Zc": 0.270,
        "Vc_cm3": 313.0, "Tb": 309.2, "mw": 72.150,
        "cp_coeffs": {"A": 2.464, "B": 45.351e-3, "C": -14.111e-6, "D": 0.0},
        "cp_tmax": 1500.0, "cp298_R": 14.731,
        "accent": "#d59b75", "theme": "copper amber",
        "visual_note": "A warm hydrocarbon palette; n-pentane vapour is colourless.",
    },
    "nitrogen": {
        "name": "Nitrogen", "formula": "N₂", "cas": "7727-37-9", "family": "Atmospheric gas",
        "Tc": 126.2, "Pc_bar": 34.00, "omega": 0.038, "Zc": 0.289,
        "Vc_cm3": 89.2, "Tb": 77.3, "mw": 28.014,
        "cp_coeffs": {"A": 3.280, "B": 0.593e-3, "C": 0.0, "D": 0.040e5},
        "cp_tmax": 2000.0, "cp298_R": 3.502,
        "accent": "#8bbcff", "theme": "liquid-nitrogen ice blue",
        "visual_note": "Inspired by liquid nitrogen's pale blue appearance; nitrogen gas is colourless.",
    },
    "oxygen": {
        "name": "Oxygen", "formula": "O₂", "cas": "7782-44-7", "family": "Atmospheric gas",
        "Tc": 154.6, "Pc_bar": 50.43, "omega": 0.022, "Zc": 0.288,
        "Vc_cm3": 73.4, "Tb": 90.2, "mw": 31.999,
        "cp_coeffs": {"A": 3.639, "B": 0.506e-3, "C": 0.0, "D": -0.227e5},
        "cp_tmax": 2000.0, "cp298_R": 3.535,
        "accent": "#70d6f3", "theme": "liquid-oxygen blue",
        "visual_note": "Inspired by pale-blue liquid oxygen; gaseous oxygen is colourless.",
    },
    "carbon-dioxide": {
        "name": "Carbon dioxide", "formula": "CO₂", "cas": "124-38-9", "family": "Inorganic gas",
        "Tc": 304.2, "Pc_bar": 73.83, "omega": 0.224, "Zc": 0.274,
        "Vc_cm3": 94.0, "Tb": None, "mw": 44.010,
        "cp_coeffs": {"A": 5.457, "B": 1.045e-3, "C": 0.0, "D": -1.157e5},
        "cp_tmax": 2000.0, "cp298_R": 4.467,
        "accent": "#c4d6e6", "theme": "dry-ice frost",
        "visual_note": "Inspired by dry-ice frost; carbon dioxide gas is colourless.",
    },
    "ammonia": {
        "name": "Ammonia", "formula": "NH₃", "cas": "7664-41-7", "family": "Polar / refrigerant",
        "Tc": 405.7, "Pc_bar": 112.80, "omega": 0.253, "Zc": 0.242,
        "Vc_cm3": 72.5, "Tb": 239.7, "mw": 17.031,
        "cp_coeffs": {"A": 3.578, "B": 3.020e-3, "C": 0.0, "D": -0.186e5},
        "cp_tmax": 1800.0, "cp298_R": 4.269,
        "accent": "#586bd5", "theme": "deep refrigeration blue",
        "visual_note": "A deep industrial-refrigeration blue; ammonia gas is colourless.",
    },
    "water": {
        "name": "Water vapour", "formula": "H₂O", "cas": "7732-18-5", "family": "Polar / associating",
        "Tc": 647.1, "Pc_bar": 220.55, "omega": 0.345, "Zc": 0.229,
        "Vc_cm3": 55.9, "Tb": 373.2, "mw": 18.015,
        "cp_coeffs": {"A": 3.470, "B": 1.450e-3, "C": 0.0, "D": 0.121e5},
        "cp_tmax": 2000.0, "cp298_R": 4.038,
        "accent": "#3ca9ff", "theme": "ocean blue",
        "visual_note": "Ocean-blue palette inspiration; water vapour itself is colourless.",
    },
}

for _gas in GASES.values():
    _gas["Pc"] = _gas["Pc_bar"] * 1e5  # Pa
    _gas["cp"] = _gas["cp298_R"] * R  # J mol-1 K-1 at 298 K, for display only

# General cubic EOS form: P = RT/(V-b) - a(T)/(V^2 + u*b*V + w*b^2).
# Constants follow Smith et al., Table 3.1, p. 100.
EOS = {
    "vdw": {"name": "van der Waals", "omega_a": 27.0 / 64.0, "omega_b": 1.0 / 8.0,
            "u": 0.0, "w": 0.0, "alpha": "constant", "short": "vdW"},
    "rk": {"name": "Redlich–Kwong", "omega_a": 0.42748, "omega_b": 0.08664,
           "u": 1.0, "w": 0.0, "alpha": "rk", "short": "RK"},
    "srk": {"name": "Soave–Redlich–Kwong", "omega_a": 0.42748, "omega_b": 0.08664,
            "u": 1.0, "w": 0.0, "alpha": "soave", "short": "SRK"},
    "pr": {"name": "Peng–Robinson", "omega_a": 0.45724, "omega_b": 0.07780,
           "u": 2.0, "w": -1.0, "alpha": "soave", "short": "PR"},
}


def _finite_number(value: Any, field: str, *, positive: bool = False) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        raise ValueError(f"'{field}' must be a number.")
    if not math.isfinite(number):
        raise ValueError(f"'{field}' must be finite.")
    if positive and number <= 0:
        raise ValueError(f"'{field}' must be greater than zero.")
    return number


def cp_over_r(temperature: float, gas: dict[str, Any]) -> float:
    """Appendix C Table C.1 ideal-gas Cp/R correlation."""
    c = gas["cp_coeffs"]
    return c["A"] + c["B"] * temperature + c["C"] * temperature**2 + c["D"] / temperature**2


def cp_at_temperature(temperature: float, gas: dict[str, Any]) -> float:
    return R * cp_over_r(temperature, gas)


def ideal_gas_reference_properties(
    temperature: float, pressure: float, t_ref: float, p_ref: float,
    h_ref: float, s_ref: float, gas: dict[str, Any],
) -> dict[str, float]:
    """Integrate Cp(T) for ideal-gas H and S from the chosen reference state.

    Uses Eqs. (6.50) and (6.51) with Cp/R = A + BT + CT^2 + D/T^2.
    """
    c = gas["cp_coeffs"]
    a, b, cc, d = c["A"], c["B"], c["C"], c["D"]
    delta_h_over_r = (
        a * (temperature - t_ref)
        + 0.5 * b * (temperature**2 - t_ref**2)
        + (cc / 3.0) * (temperature**3 - t_ref**3)
        - d * (1.0 / temperature - 1.0 / t_ref)
    )
    delta_s_over_r = (
        a * math.log(temperature / t_ref)
        + b * (temperature - t_ref)
        + 0.5 * cc * (temperature**2 - t_ref**2)
        - 0.5 * d * (1.0 / temperature**2 - 1.0 / t_ref**2)
        - math.log(pressure / p_ref)
    )
    return {
        "H_ideal": h_ref + R * delta_h_over_r,
        "S_ideal": s_ref + R * delta_s_over_r,
        "Cp": cp_at_temperature(temperature, gas),
    }


def model_parameters(gas: dict[str, Any], eos_key: str, temperature: float) -> dict[str, float]:
    """Return a(T), da/dT, b and generalized cubic constants, all in SI."""
    if eos_key not in EOS:
        raise ValueError("Choose a supported cubic equation of state.")
    if temperature <= 0:
        raise ValueError("Temperature must be greater than zero kelvin.")

    spec = EOS[eos_key]
    tc, pc, omega = gas["Tc"], gas["Pc"], gas["omega"]
    ac = spec["omega_a"] * R**2 * tc**2 / pc
    b = spec["omega_b"] * R * tc / pc
    tr = temperature / tc

    if spec["alpha"] == "constant":
        alpha, dalpha_dt = 1.0, 0.0
    elif spec["alpha"] == "rk":
        alpha = tr ** -0.5
        dalpha_dt = -0.5 * alpha / temperature
    else:
        if eos_key == "srk":
            m = 0.480 + 1.574 * omega - 0.176 * omega**2
        else:
            m = 0.37464 + 1.54226 * omega - 0.26992 * omega**2
        factor = 1.0 + m * (1.0 - math.sqrt(tr))
        alpha = factor**2
        dalpha_dt = -m * factor / math.sqrt(temperature * tc)

    a_t = ac * alpha
    da_dt = ac * dalpha_dt
    return {"a": a_t, "da_dt": da_dt, "b": b, "u": spec["u"], "w": spec["w"],
            "ac": ac, "alpha": alpha}


def pressure_from_tv(temperature: float, volume: float, gas: dict[str, Any], eos_key: str) -> float:
    p = model_parameters(gas, eos_key, temperature)
    b, a_t, u, w = p["b"], p["a"], p["u"], p["w"]
    if volume <= b * (1.0 + 1e-10):
        raise ValueError(f"Molar volume must be greater than b = {b * 1e3:.6g} L/mol for this model.")
    denominator = volume**2 + u * b * volume + w * b**2
    if denominator <= 0:
        raise ValueError("The attractive-term denominator is non-positive at this volume.")
    pressure = R * temperature / (volume - b) - a_t / denominator
    if not math.isfinite(pressure) or pressure <= 0:
        raise ValueError("This T–V state gives non-positive pressure. Increase volume or temperature.")
    return pressure


def compressibility_roots(temperature: float, pressure: float, gas: dict[str, Any], eos_key: str) -> list[float]:
    """Return real, physically admissible roots Z > B of the generic cubic EOS."""
    if pressure <= 0:
        raise ValueError("Pressure must be greater than zero.")
    p = model_parameters(gas, eos_key, temperature)
    a_t, b, u, w = p["a"], p["b"], p["u"], p["w"]
    A = a_t * pressure / (R * temperature) ** 2
    B = b * pressure / (R * temperature)

    # Z^3 + c2*Z^2 + c1*Z + c0 = 0.
    c2 = (u - 1.0) * B - 1.0
    c1 = A + w * B**2 - u * B * (1.0 + B)
    c0 = -A * B - w * B**2 * (1.0 + B)
    roots = np.roots([1.0, c2, c1, c0])
    valid: list[float] = []
    for root in roots:
        if abs(float(np.imag(root))) <= 1e-8 * max(1.0, abs(float(np.real(root)))):
            z = float(np.real(root))
            if z > B + 1e-10 and math.isfinite(z) and not any(abs(z - old) < 1e-8 for old in valid):
                valid.append(z)
    valid.sort()
    if not valid:
        raise ValueError("No physically admissible real compressibility root was found for this state.")
    return valid


def departure_properties(
    temperature: float, pressure: float, volume: float,
    gas: dict[str, Any], eos_key: str,
) -> dict[str, float]:
    """Evaluate EOS residual H, residual S and pure-fluid fugacity coefficient.

    The formulation is the generic-cubic result in Eqs. (13.74)–(13.76),
    with ln(phi) = G^R/(RT) from Eq. (10.33) in Smith et al.
    """
    p = model_parameters(gas, eos_key, temperature)
    a_t, da_dt, b, u, w = p["a"], p["da_dt"], p["b"], p["u"], p["w"]
    z = pressure * volume / (R * temperature)
    A = a_t * pressure / (R * temperature) ** 2
    B = b * pressure / (R * temperature)

    discriminant = u**2 - 4.0 * w
    if discriminant < -1e-12:
        raise ValueError("This cubic model has an unsupported complex logarithm branch.")
    delta = math.sqrt(max(0.0, discriminant))

    # I = integral used in the book's generic-cubic residual expressions.
    # For epsilon != sigma, I = ln[(Z + sigma*B)/(Z + epsilon*B)]/(sigma-epsilon).
    # For van der Waals (epsilon = sigma = 0), the continuous limit is B/Z.
    if delta > 1e-10:
        numerator = 2.0 * z + B * (u + delta)
        denominator = 2.0 * z + B * (u - delta)
        if numerator <= 0 or denominator <= 0:
            raise ValueError("Logarithm domain error while evaluating EOS departure functions.")
        integral_i = math.log(numerator / denominator) / delta
    else:
        integral_i = B / z

    if z <= B:
        raise ValueError("Compressibility factor must be greater than B.")
    q = A / B  # a(T)/(bRT)
    ln_phi = z - 1.0 - math.log(z - B) - q * integral_i
    h_res = R * temperature * (z - 1.0) + ((temperature * da_dt - a_t) / b) * integral_i
    s_res = R * math.log(z - B) + (da_dt / b) * integral_i

    return {
        "Z": z, "A": A, "B": B, "ln_phi": ln_phi,
        "phi": math.exp(max(-700.0, min(700.0, ln_phi))),
        "h_res": h_res, "s_res": s_res,
        "a": a_t, "da_dt": da_dt, "b": b, "alpha": p["alpha"],
        "I": integral_i, "q": q,
    }


def _select_root(roots: list[float], temperature: float, pressure: float,
                 gas: dict[str, Any], eos_key: str, phase: str) -> tuple[float, str]:
    # A cubic may return three real roots below Tc. The middle root lies on
    # the mechanically unstable segment of the subcritical isotherm described
    # in Sec. 3.6; stable phase comparison is therefore between the outer roots.
    if len(roots) == 1:
        return roots[0], "Single physically admissible real root"
    if phase == "vapor":
        return roots[-1], "Vapour-like outer root (largest Z)"
    if phase == "liquid":
        return roots[0], "Liquid-like outer root (smallest Z)"

    stable_candidates = [roots[0], roots[-1]] if len(roots) >= 3 else roots
    candidates = []
    for root in stable_candidates:
        volume = root * R * temperature / pressure
        dep = departure_properties(temperature, pressure, volume, gas, eos_key)
        candidates.append((dep["ln_phi"], root))
    selected = min(candidates, key=lambda item: item[0])[1]
    return selected, "Lowest-fugacity stable outer root"


def solve_temperature(pressure: float, volume: float, gas: dict[str, Any], eos_key: str) -> float:
    """Find a positive T for specified P,V by bracketing and bisection."""
    lower, upper = 5.0, 5000.0
    samples = np.geomspace(lower, upper, 400)
    valid_points: list[tuple[float, float]] = []
    for temp in samples:
        try:
            residual = pressure_from_tv(float(temp), volume, gas, eos_key) - pressure
            valid_points.append((float(temp), residual))
        except ValueError:
            continue
    brackets: list[tuple[float, float]] = []
    for (t1, f1), (t2, f2) in zip(valid_points, valid_points[1:]):
        if abs(f1) < max(1e-5, pressure * 1e-10):
            return t1
        if f1 * f2 < 0:
            brackets.append((t1, t2))
    if not brackets:
        raise ValueError("Could not bracket a temperature for this P–V pair. Check the inputs and try another molar volume.")
    lo, hi = brackets[0]
    flo = pressure_from_tv(lo, volume, gas, eos_key) - pressure
    for _ in range(120):
        mid = 0.5 * (lo + hi)
        fm = pressure_from_tv(mid, volume, gas, eos_key) - pressure
        if abs(fm) <= max(1e-5, pressure * 1e-10):
            return mid
        if flo * fm <= 0:
            hi = mid
        else:
            lo, flo = mid, fm
    return 0.5 * (lo + hi)


def calculate(payload: dict[str, Any]) -> dict[str, Any]:
    gas_key = str(payload.get("gas", "methane"))
    eos_key = str(payload.get("eos", "pr"))
    mode = str(payload.get("mode", "volume"))
    phase = str(payload.get("phase", "stable"))
    if gas_key not in GASES:
        raise ValueError("Choose one of the ten gases in the database.")
    if eos_key not in EOS:
        raise ValueError("Choose a supported cubic equation of state.")
    if phase not in {"stable", "vapor", "liquid"}:
        raise ValueError("Choose a valid root-selection method.")
    gas = GASES[gas_key]

    t_ref = _finite_number(payload.get("t_ref", T_REF_DEFAULT), "Reference temperature", positive=True)
    p_ref_bar = _finite_number(payload.get("p_ref_bar", P_REF_DEFAULT_BAR), "Reference pressure", positive=True)
    h_ref = _finite_number(payload.get("h_ref", 0.0), "Reference enthalpy")
    s_ref = _finite_number(payload.get("s_ref", 0.0), "Reference entropy")
    if not (298.0 <= t_ref <= gas["cp_tmax"]):
        raise ValueError(f"Reference temperature must be from 298 K to {gas['cp_tmax']:.0f} K for the selected Appendix C Cp correlation.")
    p_ref = p_ref_bar * 1e5

    if mode == "volume":
        pressure = _finite_number(payload.get("p_bar"), "Pressure", positive=True) * 1e5
        temperature = _finite_number(payload.get("t_k"), "Temperature", positive=True)
        roots = compressibility_roots(temperature, pressure, gas, eos_key)
        z, phase_label = _select_root(roots, temperature, pressure, gas, eos_key, phase)
        volume = z * R * temperature / pressure
    elif mode == "pressure":
        temperature = _finite_number(payload.get("t_k"), "Temperature", positive=True)
        volume = _finite_number(payload.get("v_lmol"), "Molar volume", positive=True) * 1e-3
        pressure = pressure_from_tv(temperature, volume, gas, eos_key)
        roots = compressibility_roots(temperature, pressure, gas, eos_key)
        phase_label = "Specified T–V state"
    elif mode == "temperature":
        pressure = _finite_number(payload.get("p_bar"), "Pressure", positive=True) * 1e5
        volume = _finite_number(payload.get("v_lmol"), "Molar volume", positive=True) * 1e-3
        temperature = solve_temperature(pressure, volume, gas, eos_key)
        roots = compressibility_roots(temperature, pressure, gas, eos_key)
        phase_label = "Temperature solved from P–V"
    else:
        raise ValueError("Choose whether to solve for pressure, temperature, or volume.")

    if volume <= 0 or pressure <= 0:
        raise ValueError("Calculated state is outside the physical positive P–V region.")

    departure = departure_properties(temperature, pressure, volume, gas, eos_key)
    ideal = ideal_gas_reference_properties(temperature, pressure, t_ref, p_ref, h_ref, s_ref, gas)
    h_total = ideal["H_ideal"] + departure["h_res"]
    s_total = ideal["S_ideal"] + departure["s_res"]
    fugacity = departure["phi"] * pressure

    root_details = []
    for index, root in enumerate(roots):
        if len(roots) >= 3 and index == 0:
            root_kind = "liquid-like outer root"
        elif len(roots) >= 3 and index == len(roots) - 1:
            root_kind = "vapour-like outer root"
        elif len(roots) >= 3:
            root_kind = "intermediate mathematical root (not stable)"
        else:
            root_kind = "single admissible real root"
        root_details.append({
            "Z": root,
            "V_Lmol": root * R * temperature / pressure * 1e3,
            "kind": root_kind,
        })
    warnings = []
    if temperature < 298.0 or temperature > gas["cp_tmax"]:
        warnings.append(
            f"The state temperature is outside the stated 298 K–{gas['cp_tmax']:.0f} K range of the book's Appendix C Cp correlation; ideal-gas H/S contributions are extrapolated."
        )

    return {
        "gas": {"key": gas_key, "name": gas["name"], "formula": gas["formula"]},
        "eos": {"key": eos_key, "name": EOS[eos_key]["name"]},
        "mode": mode,
        "phase": phase_label,
        "inputs": {
            "T_K": temperature, "P_bar": pressure / 1e5, "V_Lmol": volume * 1e3,
            "T_ref_K": t_ref, "P_ref_bar": p_ref_bar,
            "H_ref_Jmol": h_ref, "S_ref_JmolK": s_ref,
        },
        "properties": {
            "Z": departure["Z"], "A": departure["A"], "B": departure["B"],
            "V_Lmol": volume * 1e3, "P_bar": pressure / 1e5, "T_K": temperature,
            "phi": departure["phi"], "ln_phi": departure["ln_phi"],
            "fugacity_bar": fugacity / 1e5,
            "H_Jmol": h_total, "H_ideal_Jmol": ideal["H_ideal"], "H_res_Jmol": departure["h_res"],
            "S_JmolK": s_total, "S_ideal_JmolK": ideal["S_ideal"], "S_res_JmolK": departure["s_res"],
            "Cp_JmolK": ideal["Cp"],
            "a": departure["a"], "da_dT": departure["da_dt"], "b": departure["b"],
            "alpha": departure["alpha"], "I": departure["I"], "q": departure["q"],
        },
        "roots": root_details,
        "gas_data": {
            "Tc_K": gas["Tc"], "Pc_bar": gas["Pc_bar"], "omega": gas["omega"],
            "Zc": gas["Zc"], "Vc_cm3mol": gas["Vc_cm3"], "Tb_K": gas["Tb"],
            "molecular_weight_gmol": gas["mw"], "Cp_298_JmolK": gas["cp"],
            "Cp_coeffs": gas["cp_coeffs"], "Cp_Tmax_K": gas["cp_tmax"],
        },
        "warnings": warnings,
        "notes": [
            "H = H_ideal + H_res and S = S_ideal + S_res; residual terms are evaluated from the selected cubic EOS.",
            "Ideal-gas Cp(T) is integrated using the Appendix C, Table C.1 polynomial rather than a constant-Cp approximation.",
            "H and S are relative to the user-selected ideal-gas reference state; H_ref and S_ref default to zero and are not standard formation properties.",
        ],
    }


@app.get("/")
def home():
    return render_template("index.html")


@app.get("/api/gases")
def gases_api():
    public = {}
    for key, gas in GASES.items():
        public[key] = {
            "name": gas["name"], "formula": gas["formula"], "cas": gas["cas"], "family": gas["family"],
            "Tc": gas["Tc"], "Pc_bar": gas["Pc_bar"], "omega": gas["omega"], "Zc": gas["Zc"],
            "Vc_cm3": gas["Vc_cm3"], "Tb": gas["Tb"], "mw": gas["mw"],
            "cp298": gas["cp"], "cp298_R": gas["cp298_R"], "cp_coeffs": gas["cp_coeffs"],
            "cp_tmax": gas["cp_tmax"], "accent": gas["accent"], "theme": gas["theme"],
            "visual_note": gas["visual_note"],
        }
    return jsonify({"gases": public, "equations": {k: v["name"] for k, v in EOS.items()}})


@app.post("/api/calculate")
def calculate_api():
    try:
        payload = request.get_json(force=True) or {}
        return jsonify(calculate(payload))
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 400
    except Exception as exc:
        app.logger.exception("Unexpected calculation error")
        return jsonify({"error": f"Calculation failed: {exc}"}), 500


@app.post("/api/curve")
def curve_api():
    try:
        payload = request.get_json(force=True) or {}
        gas_key = str(payload.get("gas", "methane"))
        eos_key = str(payload.get("eos", "pr"))
        if gas_key not in GASES or eos_key not in EOS:
            raise ValueError("Choose a valid gas and cubic equation of state.")
        gas = GASES[gas_key]
        temperature = _finite_number(payload.get("t_k"), "Temperature", positive=True)
        params = model_parameters(gas, eos_key, temperature)
        b = params["b"]
        v_min = b * 1.035
        v_max = max(14.0 * R * temperature / gas["Pc"], 22.0 * b)
        volumes = np.geomspace(v_min, v_max, 220)
        p_values, v_values = [], []
        for volume in volumes:
            try:
                pressure = pressure_from_tv(temperature, float(volume), gas, eos_key)
                if pressure > 0:
                    p_values.append(pressure / 1e5)
                    v_values.append(float(volume) * 1e3)
            except ValueError:
                continue
        if not p_values:
            raise ValueError("No positive-pressure section of the isotherm was found.")
        return jsonify({"T_K": temperature, "P_bar": p_values, "V_Lmol": v_values})
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 400
    except Exception as exc:
        app.logger.exception("Unexpected curve error")
        return jsonify({"error": f"Curve generation failed: {exc}"}), 500


if __name__ == "__main__":
    app.run(host="127.0.0.1", port=5000, debug=True)
