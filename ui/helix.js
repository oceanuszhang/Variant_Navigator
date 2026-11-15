// Decorative DNA: opposing sinusoidal backbones with base-pair rungs.
(() => {
  const host = document.querySelector('.dna-crossing');
  if (!host) return;
  const ns = 'http://www.w3.org/2000/svg';
  function element(tag, attrs) {
    const node = document.createElementNS(ns, tag);
    Object.entries(attrs).forEach(([key, value]) => node.setAttribute(key, value));
    return node;
  }
  function strands(colors) {
    const group = element('g', {});
    const y = (x, sign) => 32 + sign * 24 * Math.sin(x * Math.PI / 72);
    for (let x = 9; x < 1440; x += 9) {
      group.append(element('line', {x1:x, x2:x, y1:y(x,1), y2:y(x,-1),
        stroke:'#c4a852', 'stroke-width':2, opacity:.65}));
    }
    for (let turn = 0; turn < 20; turn++) {
      const front = turn % 2 ? -1 : 1;
      for (const sign of [-front, front]) {
        let d = '';
        for (let x = turn*72; x <= (turn+1)*72; x += 2) {
          d += `${x === turn*72 ? 'M' : 'L'}${x},${y(x,sign).toFixed(2)} `;
        }
        group.append(element('path', {d, fill:'none', stroke:'#f5f4ed', 'stroke-width':7, 'stroke-linecap':'round'}));
        group.append(element('path', {d, fill:'none', stroke:colors[sign === 1 ? 0 : 1], 'stroke-width':3.5, 'stroke-linecap':'round'}));
      }
    }
    return group;
  }
  const horizontal = element('svg', {viewBox:'0 0 1440 64', preserveAspectRatio:'xMinYMid slice', 'aria-hidden':'true', focusable:'false'});
  horizontal.append(strands(['#e65f88','#63b4d6']));
  host.append(horizontal);
  const vertical = element('svg', {viewBox:'0 0 64 1440', preserveAspectRatio:'xMidYMin slice', 'aria-hidden':'true', focusable:'false', class:'dna-upright'});
  const rotated = element('g', {transform:'translate(64 0) rotate(90)'});
  rotated.append(strands(['#168458','#e2b32d']));
  vertical.append(rotated);
  host.closest('.hero').append(vertical);
})();
