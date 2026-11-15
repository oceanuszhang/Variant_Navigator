const $ = id => document.getElementById(id);
const ns = 'http://www.w3.org/2000/svg';
const signed = n => Number.isFinite(n) ? (n >= 0 ? '+' : '') + n.toFixed(4) : '?';
let artifact = null;
let requestNumber = 0;
// Native disclosures provide a click/tap fallback without experimental APIs.
const interpretationText = [
  'Predicted corrected net growth under 1,200 nM imatinib. A saved held-out ESM estimate is shown when available; otherwise an unmeasured variant uses local RandomForest inference. Higher values mean more growth in this assay. Neither model output establishes clinical resistance, and ESM validation metrics do not apply to the runtime RandomForest.',
  'An experimental net-growth measurement from the BCR-ABL screen, not a prediction. Higher values mean more growth under this fixed drug exposure. Growth reflects both protein function and drug response. A dash means no measured value is available; replicate spread is not a confidence interval.'
];
document.querySelectorAll('.result-values > div').forEach((container, index) => {
  const help = document.createElement('details');
  help.className = 'result-help';
  const trigger = document.createElement('summary');
  trigger.textContent = 'i';
  trigger.setAttribute('aria-label', index === 0 ? 'Interpret model prediction' : 'Interpret measured response');
  const content = document.createElement('div');
  content.className = 'result-help-content';
  content.textContent = interpretationText[index];
  help.append(trigger, content);
  container.querySelector('.eyebrow').after(help);
  help.addEventListener('pointerenter', event => { if (event.pointerType !== 'touch') help.open = true; });
  help.addEventListener('pointerleave', () => { if (!help.contains(document.activeElement)) help.open = false; });
  trigger.addEventListener('focus', () => { help.open = true; });
  help.addEventListener('focusout', event => { if (!help.contains(event.relatedTarget)) help.open = false; });
  help.addEventListener('keydown', event => { if (event.key === 'Escape') { help.open = false; event.preventDefault(); event.stopPropagation(); } });
  // Focus/hover may have opened it before activation: activating must not hide it.
  trigger.addEventListener('click', event => { event.preventDefault(); help.open = true; });
});
async function api(path, options) {
  if (location.protocol === 'file:') throw new Error('Start python ui/server.py and open http://127.0.0.1:8000');
  const response = await fetch(path, options);
  const data = await response.json();
  if (!response.ok) throw new Error(data.error || 'Request failed');
  return data;
}
document.querySelectorAll('nav button').forEach(button => button.addEventListener('click', () => {
  document.querySelectorAll('.page').forEach(page => page.hidden = page.id !== button.dataset.view);
  document.querySelectorAll('nav button').forEach(b => {
    b.classList.toggle('active', b === button);
    if (b === button) b.setAttribute('aria-current', 'page'); else b.removeAttribute('aria-current');
  });
  $('breadcrumb').textContent = button.querySelector('span').textContent;
}));
async function explore() {
  const request = ++requestNumber;
  const variant = $('variant').value.trim().toUpperCase();
  $('explore').disabled = true;
  $('input-message').textContent = 'Running pipeline. First-time inference may load data and fit the local model.';
  try {
    const data = await api('/api/predict', {method:'POST', headers:{'Content-Type':'application/json'}, body:JSON.stringify({variant, sequence:$('aa-sequence').value})});
    if (request !== requestNumber) return;
    const result = data.result, held = data.held_out;
    const estimate = held ? held.predicted : result.predicted_netgrowth;
    $('result-variant').textContent = result.mutant;
    $('predicted').textContent = signed(estimate);
    $('measured').textContent = signed(result.measured_netgrowth);
    $('head-label').textContent = held ? `${held.head} ? saved held-out estimate` : result.kind === 'inference' ? 'RandomForest ? runtime inference ? unverified' : 'No model estimate for this variant';
    document.querySelector('.measured + .value-caption').textContent = result.kind === 'measurement' ? 'BCR-ABL DMS ? 1,200 nM imatinib' : 'No measured assay value';
    document.querySelector('.result-card .pill').textContent = held ? 'Measurement + held-out estimate' : result.kind;
    document.querySelectorAll('.sequence code').forEach((element, i) => {
      const b = document.createElement('b'); b.textContent = i ? data.sequence.alt : data.sequence.ref;
      if (i) b.className = 'alt';
      element.replaceChildren(document.createTextNode('? '+data.sequence.before+' '), b, document.createTextNode(' '+data.sequence.after+' ?'));
    });
    $('response-marker').hidden = !Number.isFinite(estimate);
    $('response-marker').style.left = Math.max(0, Math.min(100, (estimate+.03)/.1*100))+'%';
    $('variant').setAttribute('aria-invalid','false');
    $('input-message').textContent = `${result.source}. ${held ? held.provenance_limit : result.caveat || 'Measured assay lookup; not an inferred result.'}`;
  } catch (error) {
    if (request === requestNumber) {
      $('input-message').textContent = error.message + ' Previous results, if present, retain their original variant label.';
      $('variant').setAttribute('aria-invalid','true');
    }
  } finally { if (request === requestNumber) $('explore').disabled = false; }
}
$('explore').addEventListener('click', explore);
$('load-reference').addEventListener('click', async () => {
  try {
    const data = await api('/api/reference');
    $('aa-sequence').value = data.sequence;
    $('input-message').textContent = 'ABL1 study reference loaded. Enter a substitution and submit.';
  } catch (error) { $('input-message').textContent = error.message; }
});
$('variant').addEventListener('keydown', e => {if (e.key === 'Enter' && !e.isComposing && !$('explore').disabled) explore();});
document.querySelectorAll('[data-variant]').forEach(b => b.addEventListener('click', () => {$('variant').value = b.dataset.variant; if (!$('explore').disabled) explore();}));
function table() {
  const query = $('filter').value.trim().toUpperCase();
  $('rows').replaceChildren();
  const rows = (artifact?.rows || []).filter(r => r.mutant.toUpperCase().includes(query));
  rows.forEach(row => {
    const tr = document.createElement('tr');
    [row.mutant, signed(row.measured), signed(row.predicted), Math.abs(row.predicted-row.measured).toFixed(4), row.rep_low == null || row.rep_high == null ? 'Replicates unavailable' : row.in_range ? 'Within range' : 'Outside range'].forEach(value => {
      const td = document.createElement('td'); td.textContent = value; tr.append(td);
    });
    $('rows').append(tr);
  });
  if (!rows.length) { const tr=document.createElement('tr'),td=document.createElement('td');td.colSpan=5;td.textContent=artifact?'No matching variants':'No validation results loaded';tr.append(td);$('rows').append(tr); }
}
$('filter').addEventListener('input', table);
function plot() {
  $('scatter').replaceChildren();
  const rows = artifact?.rows || [];
  if (!rows.length) return;
  let low=Infinity, high=-Infinity;
  for (const r of rows) {low=Math.min(low,r.measured,r.predicted);high=Math.max(high,r.measured,r.predicted);}
  const pad=(high-low)*.06 || .01;low-=pad;high+=pad;
  const sx=v=>65+(v-low)/(high-low)*500, sy=v=>285-(v-low)/(high-low)*250;
  function add(tag,attrs,text) {const el=document.createElementNS(ns,tag);Object.entries(attrs).forEach(([k,v])=>el.setAttribute(k,v));if(text)el.textContent=text;$('scatter').append(el);}
  for(let i=0;i<=5;i++){const v=low+(high-low)*i/5;add('line',{x1:65,y1:sy(v),x2:565,y2:sy(v),stroke:'#dbe3d6'});add('text',{x:55,y:sy(v)+4,'text-anchor':'end',fill:'#526057','font-size':11},v.toFixed(3));add('text',{x:sx(v),y:306,'text-anchor':'middle',fill:'#526057','font-size':11},v.toFixed(3));}
  add('line',{x1:65,y1:285,x2:565,y2:35,stroke:'#536a57','stroke-dasharray':'5 5'});
  for(const r of rows){const highlight=r.mutant==='T315I';add('circle',{cx:sx(r.measured),cy:sy(r.predicted),r:highlight?5:2.5,fill:highlight?'#a94b30':'#14754b',opacity:highlight?1:.5});if(highlight)add('text',{x:sx(r.measured)+7,y:sy(r.predicted)-6,fill:'#82351f','font-size':12},r.mutant);}
  add('text',{x:315,y:333,'text-anchor':'middle',fill:'#526057','font-size':12},'Measured corrected net growth');
  add('text',{x:15,y:165,transform:'rotate(-90 15 165)','text-anchor':'middle',fill:'#526057','font-size':12},'Predicted corrected net growth');
}
$('download').addEventListener('click',()=>{
  if(!artifact)return;
  const rows=artifact.rows.filter(r=>r.mutant.toUpperCase().includes($('filter').value.trim().toUpperCase()));
  const keys=['mutant','position','measured','predicted','error','abs_error','rep_low','rep_high','in_range'];
  const csv=keys.join(',')+'\n'+rows.map(r=>keys.map(k=>JSON.stringify(r[k]??'')).join(',')).join('\n');
  const url=URL.createObjectURL(new Blob([csv],{type:'text/csv'}));const a=document.createElement('a');a.href=url;a.download='validation-predictions.csv';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
});
async function loadValidation(){
  try {
    artifact=await api('/api/validation');
    const m=artifact.metrics[artifact.head];
    const values=document.querySelectorAll('.metrics strong');
    values[0].textContent=m.mae.toFixed(4);values[1].textContent=m.spearman.toFixed(3);values[2].textContent=m.in_range_pct.toFixed(1)+'%';
    document.querySelector('#validation .notice').textContent=`${artifact.model} / ${artifact.head} ? ${artifact.rows.length} held-out variants ? seed ${artifact.seed}. ${artifact.provenance_limit}`;
    document.querySelector('.validation-grid .pill').textContent=artifact.head+' ? saved validation';
    document.querySelector('.split-card dd:nth-of-type(2)').textContent=artifact.seed;
    document.querySelector('#overview .notice').textContent=`Validation artifact SHA-256: ${artifact.artifact_sha256}. New unmeasured variants use the pipeline's local RandomForest, not this ESM head. Runtime validation is unverified.`;
    document.querySelector('#overview .pill').textContent='Local pipeline connected';
    $('download').disabled=false;
    table();plot();
  } catch(error){document.querySelector('#validation .notice').textContent=error.message;table();}
}
loadValidation();

// Two opposing sinusoidal backbones, joined by base-pair rungs.
// Draw each half-turn in depth order so the strands alternate at crossings.
const ribbon = document.querySelector('.dna-ribbon');
if (ribbon) {
  const helix = document.createElementNS(ns, 'svg');
  helix.setAttribute('viewBox', '0 0 720 64');
  helix.setAttribute('preserveAspectRatio', 'xMinYMid meet');
  helix.setAttribute('aria-hidden', 'true');
  helix.setAttribute('focusable', 'false');
  const center = 32, amplitude = 24, halfTurn = 72;
  const y = (x, strand) => center + strand * amplitude * Math.sin(x * Math.PI / halfTurn);
  const add = (tag, attrs) => {
    const element = document.createElementNS(ns, tag);
    Object.entries(attrs).forEach(([name, value]) => element.setAttribute(name, value));
    helix.append(element);
    return element;
  };
  // Rungs sit behind both strands; their length narrows toward each crossing.
  for (let x = 9; x < 720; x += 9) {
    add('line', {x1:x, y1:y(x, 1), x2:x, y2:y(x, -1),
      stroke:Math.floor(x / halfTurn) % 2 ? '#69b7be' : '#d9ad35',
      'stroke-width':2, opacity:0.75});
  }
  for (let turn = 0; turn < 10; turn++) {
    const start = turn * halfTurn;
    const front = turn % 2 ? -1 : 1;
    for (const strand of [-front, front]) {
      const points = [];
      for (let x = start; x <= start + halfTurn; x += 2) {
        points.push(`${x === start ? 'M' : 'L'}${x},${y(x, strand).toFixed(2)}`);
      }
      const d = points.join(' ');
      // A paper-colored halo makes the over/under crossing legible.
      add('path', {d, fill:'none', stroke:'#f7f7ee', 'stroke-width':7,
        'stroke-linecap':'round'});
      add('path', {d, fill:'none', stroke:strand === 1 ? '#168458' : '#e2b32d',
        'stroke-width':3.5, 'stroke-linecap':'round'});
    }
  }
  // Compose two complete helices: one across the top, one down the right.
  const composition = document.createElementNS(ns, 'svg');
  composition.setAttribute('viewBox', '0 0 720 160');
  composition.setAttribute('preserveAspectRatio', 'none');
  composition.setAttribute('aria-hidden', 'true');
  composition.setAttribute('focusable', 'false');
  const horizontal = document.createElementNS(ns, 'g');
  horizontal.setAttribute('transform', 'translate(0 8)');
  const vertical = document.createElementNS(ns, 'g');
  vertical.setAttribute('transform', 'translate(640 -260) rotate(90)');
  for (const element of helix.children) {
    horizontal.append(element.cloneNode(true));
    const copy = element.cloneNode(true);
    // Pink and blue distinguish the horizontal helix from green and gold.
    if (element.getAttribute('stroke') === '#168458') element.setAttribute('stroke', '#e65f88');
    if (element.getAttribute('stroke') === '#e2b32d') element.setAttribute('stroke', '#63b4d6');
    vertical.append(copy);
  }
  // Apply the horizontal palette after the independent vertical copy exists.
  horizontal.querySelectorAll('[stroke="#168458"]').forEach(e => e.setAttribute('stroke', '#e65f88'));
  horizontal.querySelectorAll('[stroke="#e2b32d"]').forEach(e => e.setAttribute('stroke', '#63b4d6'));
  composition.append(horizontal);
  ribbon.replaceChildren(composition);
  ribbon.closest('.page-top').append(ribbon);

  // Extend the same vertical helix through the header without moving the
  // horizontal strand out of its reserved space below the subtitle.
  const verticalExtension = document.createElementNS(ns, 'svg');
  verticalExtension.classList.add('dna-vertical-extension');
  verticalExtension.setAttribute('viewBox', '0 0 64 720');
  verticalExtension.setAttribute('preserveAspectRatio', 'xMidYMid slice');
  verticalExtension.setAttribute('aria-hidden', 'true');
  verticalExtension.setAttribute('focusable', 'false');
  const extensionGroup = document.createElementNS(ns, 'g');
  extensionGroup.setAttribute('transform', 'translate(64 0) rotate(90)');
  for (const element of helix.children) {
    const copy = element.cloneNode(true);
    if (copy.getAttribute('stroke') === '#e65f88') copy.setAttribute('stroke', '#168458');
    if (copy.getAttribute('stroke') === '#63b4d6') copy.setAttribute('stroke', '#e2b32d');
    extensionGroup.append(copy);
  }
  verticalExtension.append(extensionGroup);
  ribbon.closest('.page-top').append(verticalExtension);
}
