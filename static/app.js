const $ = (id) => document.getElementById(id);
let gasDatabase = {};
let lastResult = null;
let lastCurve = null;

const modelInfo = {
  pr: {
    code: 'PR', name: 'Peng–Robinson equation of state',
    equation: 'P = RT / (V − b) − a(T) / [V² + 2bV − b²]',
    note: 'A widely used cubic model for real-fluid PVT behaviour and vapour–liquid calculations.',
    alpha: 'α = [1 + m(1 − √Tr)]²; m = 0.37464 + 1.54226ω − 0.26992ω²',
  },
  srk: {
    code: 'SRK', name: 'Soave–Redlich–Kwong equation of state',
    equation: 'P = RT / (V − b) − a(T) / [V(V + b)]',
    note: 'A temperature-dependent attraction correction commonly used for vapour-phase and hydrocarbon calculations.',
    alpha: 'α = [1 + m(1 − √Tr)]²; m = 0.480 + 1.574ω − 0.176ω²',
  },
  rk: {
    code: 'RK', name: 'Redlich–Kwong equation of state',
    equation: 'P = RT / (V − b) − a(T) / [V(V + b)]',
    note: 'The original RK attraction term has α = Tr⁻¹ᐟ² and uses the same cubic denominator as SRK.',
    alpha: 'α = Tr⁻¹ᐟ²',
  },
  vdw: {
    code: 'vdW', name: 'van der Waals equation of state',
    equation: 'P = RT / (V − b) − a / V²',
    note: 'The foundational cubic EOS: b represents finite molecular size and a represents intermolecular attraction.',
    alpha: 'α = 1',
  },
};

function setAccent(hex) {
  document.documentElement.style.setProperty('--accent', hex);
  const clean = hex.replace('#', '');
  const r = parseInt(clean.slice(0, 2), 16);
  const g = parseInt(clean.slice(2, 4), 16);
  const b = parseInt(clean.slice(4, 6), 16);
  document.documentElement.style.setProperty('--accent-rgb', `${r}, ${g}, ${b}`);
}

function fmt(value, digits = 4) {
  if (value === null || value === undefined || !Number.isFinite(Number(value))) return '—';
  const n = Number(value);
  if (Math.abs(n) >= 1e5 || (Math.abs(n) > 0 && Math.abs(n) < 1e-4)) return n.toExponential(3);
  return n.toLocaleString('en-US', { maximumFractionDigits: digits, minimumFractionDigits: 0 });
}

function gasSelected() {
  return gasDatabase[$('gasSelect').value];
}

function cpFormulaText(coeffs) {
  const terms = [`${fmt(coeffs.A, 5)}`, `${fmt(coeffs.B, 7)}T`, `${fmt(coeffs.C, 9)}T²`, `${fmt(coeffs.D, 3)}/T²`];
  return `Cp/R = ${terms.join(' + ').replaceAll('+ -', '− ')}`;
}

function updateGasTheme() {
  const gas = gasSelected();
  if (!gas) return;
  setAccent(gas.accent);
  $('gasSymbol').textContent = gas.formula;
  $('gasName').textContent = gas.name;
  $('gasVisualNote').textContent = gas.visual_note;
  $('gasTc').innerHTML = `${fmt(gas.Tc, 1)} <i>K</i>`;
  $('gasPc').innerHTML = `${fmt(gas.Pc_bar, 2)} <i>bar</i>`;
  $('gasOmega').textContent = fmt(gas.omega, 3);
  $('gasTheme').textContent = gas.theme.toUpperCase();
  $('gasMeta').textContent = `${gas.family} · ${gas.formula} · ${fmt(gas.mw, 3)} g/mol`;
  $('gasMw').textContent = `${fmt(gas.mw, 3)} g/mol`;
  $('gasZc').textContent = fmt(gas.Zc, 3);
  $('gasVc').textContent = `${fmt(gas.Vc_cm3, 1)} cm³/mol`;
  $('gasTb').textContent = gas.Tb === null ? 'Not listed' : `${fmt(gas.Tb, 1)} K`;
  $('cpEquation').textContent = cpFormulaText(gas.cp_coeffs);
  $('cpRange').textContent = `Appendix C, Table C.1 · stated range: 298–${fmt(gas.cp_tmax, 0)} K · Cp/R at 298 K = ${fmt(gas.cp298_R, 3)}`;
  updateEquation();
}

function updateEquation() {
  const model = modelInfo[$('eosSelect').value];
  $('modelCode').textContent = model.code;
  $('equationName').textContent = model.name;
  $('equationDisplay').textContent = model.equation;
  $('modelNote').textContent = model.note;
  $('modelFoot').textContent = `${model.alpha}. Parameters a(T) and b are calculated from the selected gas's Tc, Pc and acentric factor ω where the model requires it.`;
}

function updateMode() {
  const mode = $('modeSelect').value;
  document.querySelectorAll('[data-input]').forEach((field) => {
    const key = field.dataset.input;
    const show = (mode === 'volume' && ['p', 't'].includes(key))
      || (mode === 'pressure' && ['t', 'v'].includes(key))
      || (mode === 'temperature' && ['p', 'v'].includes(key));
    field.style.display = show ? 'block' : 'none';
  });
  $('phaseRow').style.display = mode === 'volume' ? 'block' : 'none';
}

function payloadFromForm() {
  return {
    gas: $('gasSelect').value,
    eos: $('eosSelect').value,
    mode: $('modeSelect').value,
    phase: $('phaseSelect').value,
    p_bar: $('pInput').value,
    t_k: $('tInput').value,
    v_lmol: $('vInput').value,
    t_ref: $('tRefInput').value,
    p_ref_bar: $('pRefInput').value,
    h_ref: $('hRefInput').value,
    s_ref: $('sRefInput').value,
  };
}

function setError(message) {
  $('errorBox').textContent = message;
  $('errorBox').hidden = !message;
}

function updateResultUI(result) {
  const p = result.properties;
  lastResult = result;
  $('emptyState').hidden = true;
  $('resultsContent').hidden = false;
  $('thermoEmpty').hidden = true;
  $('thermoResults').hidden = false;
  $('copyBtn').disabled = false;
  $('phasePill').textContent = result.phase.toUpperCase();
  $('zValue').textContent = fmt(p.Z, 5);
  $('zCaption').textContent = p.Z > 1.015 ? 'Repulsive effects dominate' : (p.Z < 0.985 ? 'Attractive effects dominate' : 'Close to ideal-gas behaviour');
  $('vValue').textContent = `${fmt(p.V_Lmol, 5)} L/mol`;
  $('phiValue').textContent = fmt(p.phi, 5);
  $('phiCaption').textContent = `ln φ = ${fmt(p.ln_phi, 5)}`;
  $('fValue').textContent = `${fmt(p.fugacity_bar, 5)} bar`;
  $('stateSummary').textContent = `${fmt(p.T_K, 3)} K · ${fmt(p.P_bar, 4)} bar · ${result.eos.name}`;

  $('hValue').textContent = fmt(p.H_Jmol, 3);
  $('hIdeal').textContent = `${fmt(p.H_ideal_Jmol, 3)} J/mol`;
  $('hResidual').textContent = `${fmt(p.H_res_Jmol, 3)} J/mol`;
  $('sValue').textContent = fmt(p.S_JmolK, 4);
  $('sIdeal').textContent = `${fmt(p.S_ideal_JmolK, 4)} J/mol·K`;
  $('sResidual').textContent = `${fmt(p.S_res_JmolK, 4)} J/mol·K`;
  $('hBar').style.width = `${Math.min(100, Math.abs(p.H_res_Jmol) / (Math.abs(p.H_res_Jmol) + Math.abs(p.H_ideal_Jmol) + 1) * 100)}%`;
  $('sBar').style.width = `${Math.min(100, Math.abs(p.S_res_JmolK) / (Math.abs(p.S_res_JmolK) + Math.abs(p.S_ideal_JmolK) + 1) * 100)}%`;

  $('aValue').textContent = fmt(p.a, 5);
  $('bValue').textContent = fmt(p.b, 7);
  $('alphaValue').textContent = fmt(p.alpha, 6);
  $('daValue').textContent = fmt(p.da_dT, 7);
  $('chartTemp').textContent = `T = ${fmt(p.T_K, 3)} K`;
  $('rootSummary').textContent = result.roots.map((r, idx) => `Z${idx + 1} ${fmt(r.Z, 5)} · ${r.kind}`).join(' | ');
  if (result.warnings && result.warnings.length) {
    setError(result.warnings.join(' '));
  }
}

async function postJson(url, body) {
  const response = await fetch(url, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(body),
  });
  const data = await response.json();
  if (!response.ok) throw new Error(data.error || `Request failed (${response.status})`);
  return data;
}

async function runCalculation(event) {
  if (event) event.preventDefault();
  setError('');
  const btn = $('calculateBtn');
  btn.classList.add('loading');
  btn.querySelector('span').textContent = 'Solving equation…';
  try {
    const result = await postJson('/api/calculate', payloadFromForm());
    updateResultUI(result);
    await loadCurve(result.properties.T_K, result.properties.V_Lmol);
  } catch (error) {
    setError(error.message);
  } finally {
    btn.classList.remove('loading');
    btn.querySelector('span').textContent = 'Run calculation';
  }
}

async function loadCurve(temperature = Number($('tInput').value), pointVolume = null) {
  try {
    const data = await postJson('/api/curve', {
      gas: $('gasSelect').value,
      eos: $('eosSelect').value,
      t_k: temperature,
    });
    lastCurve = data;
    renderChart(data, pointVolume);
  } catch (error) {
    $('pvtChart').innerHTML = `<div class="chart-empty">${escapeHtml(error.message)}</div>`;
  }
}

function renderChart(data, pointVolume = null) {
  if (!window.Plotly) {
    $('pvtChart').innerHTML = '<div class="chart-empty">Plotly could not be loaded. Check your internet connection.</div>';
    return;
  }
  const accent = getComputedStyle(document.documentElement).getPropertyValue('--accent').trim();
  const traces = [{
    x: data.V_Lmol, y: data.P_bar, type: 'scatter', mode: 'lines',
    name: $('eosSelect').selectedOptions[0].text,
    line: { color: accent, width: 2.5, shape: 'spline', smoothing: 0.35 },
    hovertemplate: 'V = %{x:.5g} L/mol<br>P = %{y:.5g} bar<extra></extra>',
  }];
  if (pointVolume && lastResult) {
    traces.push({
      x: [lastResult.properties.V_Lmol], y: [lastResult.properties.P_bar],
      type: 'scatter', mode: 'markers', name: 'Calculated state',
      marker: { color: '#ffffff', size: 10, line: { color: accent, width: 2.5 } },
      hovertemplate: 'Calculated state<br>V = %{x:.5g} L/mol<br>P = %{y:.5g} bar<extra></extra>',
    });
  }
  const layout = {
    margin: { l: 55, r: 18, t: 12, b: 48 },
    paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
    font: { family: 'DM Mono, monospace', size: 9, color: '#8297b0' },
    xaxis: { title: { text: 'Molar volume, V (L/mol)', font: { size: 9, color: '#8297b0' } }, gridcolor: 'rgba(164,190,220,.09)', zerolinecolor: 'rgba(164,190,220,.15)', linecolor: 'rgba(164,190,220,.18)', type: 'log' },
    yaxis: { title: { text: 'Pressure, P (bar)', font: { size: 9, color: '#8297b0' } }, gridcolor: 'rgba(164,190,220,.09)', zerolinecolor: 'rgba(164,190,220,.15)', linecolor: 'rgba(164,190,220,.18)', rangemode: 'tozero' },
    showlegend: false,
    hoverlabel: { bgcolor: '#101b2b', bordercolor: accent, font: { color: '#e8f1fd', size: 9 } },
  };
  Plotly.react('pvtChart', traces, layout, { responsive: true, displayModeBar: false });
}

function escapeHtml(text) {
  return String(text).replace(/[&<>"']/g, (char) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#039;' }[char]));
}

async function copyResults() {
  if (!lastResult) return;
  const p = lastResult.properties;
  const text = [
    `PVT Lab — ${lastResult.gas.name} (${lastResult.gas.formula})`,
    `EOS: ${lastResult.eos.name}`,
    `T = ${p.T_K} K | P = ${p.P_bar} bar | V = ${p.V_Lmol} L/mol`,
    `Z = ${p.Z} | phi = ${p.phi} | fugacity = ${p.fugacity_bar} bar`,
    `H = ${p.H_Jmol} J/mol (ideal ${p.H_ideal_Jmol}; residual ${p.H_res_Jmol})`,
    `S = ${p.S_JmolK} J/mol K (ideal ${p.S_ideal_JmolK}; residual ${p.S_res_JmolK})`,
    `Cp(T) = ${p.Cp_JmolK} J/mol K`,
    `Root: ${lastResult.phase}`,
  ].join('\n');
  try {
    await navigator.clipboard.writeText(text);
    $('copyBtn').textContent = '✓';
    setTimeout(() => $('copyBtn').textContent = '⧉', 1400);
  } catch {
    setError('Clipboard permission was blocked by the browser.');
  }
}

async function init() {
  try {
    const response = await fetch('/api/gases');
    const data = await response.json();
    gasDatabase = data.gases;
    $('gasSelect').innerHTML = Object.entries(gasDatabase)
      .map(([key, gas]) => `<option value="${key}">${gas.name} · ${gas.formula}</option>`).join('');
    $('gasSelect').value = 'methane';
    updateGasTheme();
    updateMode();
  } catch (error) {
    setError(`Could not load the gas database: ${error.message}`);
  }
}

$('calcForm').addEventListener('submit', runCalculation);
$('gasSelect').addEventListener('change', updateGasTheme);
$('eosSelect').addEventListener('change', updateEquation);
$('modeSelect').addEventListener('change', updateMode);
$('copyBtn').addEventListener('click', copyResults);
$('advancedToggle').addEventListener('click', () => {
  const fields = $('advancedFields');
  fields.hidden = !fields.hidden;
  $('advancedToggle').querySelector('span').textContent = fields.hidden ? '＋' : '−';
});
$('refreshCurveBtn').addEventListener('click', () => {
  const t = lastResult ? lastResult.properties.T_K : Number($('tInput').value);
  loadCurve(t, lastResult ? lastResult.properties.V_Lmol : null);
});

init();
