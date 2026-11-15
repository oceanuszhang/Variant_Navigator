const $ = id => document.getElementById(id);
let reference = '';
let lastEstimate = null;

async function api(path, options) {
  const response = await fetch(path, options);
  const data = await response.json();
  if (!response.ok) throw new Error(data.error || 'Request failed');
  return data;
}
const pct = value => (100 * value).toFixed(1) + '%';
const label = value => value ? value[0].toUpperCase() + value.slice(1) : '—';

async function loadReference() {
  try {
    const data = await api('/api/reference');
    reference = data.sequence;
    $('sequence').value = reference;
  } catch (error) {
    $('message').textContent = error.message;
  }
}

function requestBody() {
  return {
    variant: $('variant').value.trim().toUpperCase(),
    sequence: $('sequence').value,
    model: $('model').value
  };
}

async function predict() {
  $('predict').disabled = true;
  $('message').textContent = 'Looking up saved ESM embedding and model output…';
  try {
    const data = await api('/api/predict', {
      method: 'POST', headers: {'Content-Type': 'application/json'},
      body: JSON.stringify(requestBody())
    });
    $('result-title').textContent = data.mutant + ' · ABL1';
    $('predicted').textContent = label(data.prediction);
    $('predicted').className = data.prediction;
    $('measured').textContent = label(data.measured_class);
    $('measured').className = data.measured_class;
    $('growth').textContent = (data.measured_netgrowth >= 0 ? '+' : '') + data.measured_netgrowth.toFixed(5);
    $('split').textContent = data.split + ' site';
    $('probabilities').replaceChildren();
    if (data.probabilities) {
      Object.entries(data.probabilities).forEach(([name, value]) => {
        const row = document.createElement('div');
        row.innerHTML = '<span>' + label(name) + '</span><div><i style="width:' + (100*value) + '%"></i></div><b>' + pct(value) + '</b>';
        $('probabilities').append(row);
      });
    }
    $('message').textContent = data.deployment_scope;
    lastEstimate = null;
    $('run-boltz').disabled = true;
  } catch (error) {
    $('message').textContent = error.message;
  } finally {
    $('predict').disabled = false;
  }
}

async function loadValidation() {
  try {
    const data = await api('/api/validation');
    const models = data.metrics.models;
    $('validation-count').textContent = data.metrics.split.validation_rows + ' variants · ' + data.metrics.split.validation_positions + ' sites';
    $('metrics').replaceChildren();
    [
      ['Majority baseline', models.majority_baseline],
      ['Ridge', models.ridge],
      ['XGBoost', models.xgboost]
    ].forEach(([name, values]) => {
      const card = document.createElement('div');
      card.innerHTML = '<span>' + name + '</span><strong>' + pct(values.balanced_accuracy) +
        '</strong><small>balanced accuracy · macro F1 ' + values.macro_f1.toFixed(3) + '</small>';
      $('metrics').append(card);
    });
    const matrix = models.xgboost.confusion_matrix;
    $('confusion').innerHTML =
      '<table><thead><tr><th>actual ↓ / predicted →</th><th>decrease</th><th>unchanged</th><th>increase</th></tr></thead><tbody>' +
      ['decrease','unchanged','increase'].map((name, i) =>
        '<tr><th>' + name + '</th>' + matrix[i].map(value => '<td>' + value + '</td>').join('') + '</tr>'
      ).join('') + '</tbody></table>';
  } catch (error) {
    $('validation-count').textContent = error.message;
  }
}

async function loadPtarmigan() {
  try {
    const data = await api('/api/ptarmigan/abl1');
    $('ptarmigan-note').textContent = data.scope;
    const hits = (data.hits.data || []).slice(0, 5);
    $('ptarmigan').innerHTML = '<table><thead><tr><th>Rank</th><th>Compound</th><th>Score</th></tr></thead><tbody>' +
      hits.map(hit => '<tr><td>' + hit.rank + '</td><td>' + hit.name + '</td><td>' +
        (hit.score ?? hit.prediction ?? '—') + '</td></tr>').join('') + '</tbody></table>';
  } catch (error) {
    $('ptarmigan-note').textContent = error.message;
  }
}

async function estimateBoltz() {
  $('estimate').disabled = true;
  $('boltz-message').textContent = 'Requesting a free estimate; no structure run is being submitted…';
  try {
    lastEstimate = await api('/api/boltz/estimate', {
      method: 'POST', headers: {'Content-Type': 'application/json'},
      body: JSON.stringify(requestBody())
    });
    $('boltz-message').textContent = 'Estimated cost: $' + lastEstimate.estimated_cost_usd.toFixed(2) +
      ' for ' + lastEstimate.num_samples + ' samples. No paid job submitted.';
    $('run-boltz').disabled = false;
  } catch (error) {
    $('boltz-message').textContent = error.message;
  } finally {
    $('estimate').disabled = false;
  }
}

async function runBoltz() {
  if (!lastEstimate) return;
  const amount = lastEstimate.estimated_cost_usd.toFixed(2);
  if (!window.confirm('Submit this Boltz job? It may charge up to $' + amount + ' in service credits.')) return;
  $('run-boltz').disabled = true;
  $('boltz-message').textContent = 'Boltz job running…';
  try {
    const body = {...requestBody(), confirm_paid_run: true, max_usd: Math.min(0.25, lastEstimate.estimated_cost_usd)};
    const data = await api('/api/boltz/run', {
      method: 'POST', headers: {'Content-Type': 'application/json'}, body: JSON.stringify(body)
    });
    $('boltz-message').textContent = 'Completed · confidence ' + data.confidence.toFixed(1) +
      ' · recorded cost $' + (data.cost_usd || 0).toFixed(2);
  } catch (error) {
    $('boltz-message').textContent = error.message;
    $('run-boltz').disabled = false;
  }
}

$('load').addEventListener('click', loadReference);
$('predict').addEventListener('click', predict);
$('estimate').addEventListener('click', estimateBoltz);
$('run-boltz').addEventListener('click', runBoltz);
$('variant').addEventListener('keydown', event => { if (event.key === 'Enter') predict(); });
document.querySelectorAll('[data-variant]').forEach(button => button.addEventListener('click', () => {
  $('variant').value = button.dataset.variant;
  predict();
}));
loadReference();
loadValidation();
loadPtarmigan();

