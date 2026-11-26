/* Local 2x2 operator and row-reduction experiments; no network access. */
(() => {
  const value = id => Number(document.getElementById(id)?.value);
  const write = (id, text) => {const node = document.getElementById(id); if (node) node.textContent = text;};
  const fmt = x => Math.abs(x) < 1e-12 ? '0' : Number(x.toFixed(6)).toString();
  function operator() {
    if (!document.getElementById('op-a')) return;
    const a=value('op-a'), b=value('op-b'), c=value('op-c'), d=value('op-d');
    if (![a,b,c,d].every(Number.isFinite) || [a,b,c,d].some(x=>Math.abs(x)>1e6)) {
      ['op-trace','op-det','op-eig','op-kind'].forEach(id=>write(id,'—'));
      write('op-error','Enter finite entries with absolute value at most 1e6.');
      return;
    }
    const tr=a+d, det=a*d-b*c, disc=tr*tr-4*det;
    write('op-trace',fmt(tr)); write('op-det',fmt(det));
    let kind;
    if (disc > 1e-12) {
      const r=Math.sqrt(disc);
      write('op-eig', `${fmt((tr+r)/2)} and ${fmt((tr-r)/2)}`);
      kind='Two distinct real eigenvalues: diagonalizable over R.';
    } else if (disc < -1e-12) {
      const r=Math.sqrt(-disc)/2;
      write('op-eig', `${fmt(tr/2)} ± ${fmt(r)}i`);
      kind='No real eigenvalue: not diagonalizable over R, diagonalizable over C (distinct complex eigenvalues).';
    } else {
      write('op-eig', `${fmt(tr/2)} (repeated)`);
      const scalar=Math.abs(b)<1e-12&&Math.abs(c)<1e-12&&Math.abs(a-d)<1e-12;
      kind=scalar ? 'Repeated eigenvalue and A = λI: diagonal already, minimal polynomial x − λ.' :
        'Repeated eigenvalue but A ≠ λI: not diagonalizable, minimal polynomial (x − λ)².';
    }
    write('op-kind', kind);
    const ch=[[a*a+b*c-tr*a+det, a*b+b*d-tr*b],[c*a+d*c-tr*c, c*b+d*d-tr*d+det]];
    const residual=Math.max(...ch.flat().map(Math.abs));
    write('op-error', `Cayley–Hamilton residual max|A² − (tr A)A + (det A)I| = ${residual.toExponential(2)}.`);
  }
  function rank() {
    if (!document.getElementById('rr-00')) return;
    const m=[0,1,2].map(i=>[0,1,2,3].map(j=>value(`rr-${i}${j}`)));
    if (!m.flat().every(Number.isFinite)) {write('rr-rank','—'); write('rr-pivots','—'); return;}
    let row=0; const pivots=[];
    for (let col=0; col<4 && row<3; col++) {
      let best=row;
      for (let i=row+1;i<3;i++) if (Math.abs(m[i][col])>Math.abs(m[best][col])) best=i;
      if (Math.abs(m[best][col])<1e-10) continue;
      [m[row],m[best]]=[m[best],m[row]];
      const p=m[row][col];
      m[row]=m[row].map(x=>x/p);
      for (let i=0;i<3;i++) if (i!==row) {const f=m[i][col]; m[i]=m[i].map((x,j)=>x-f*m[row][j]);}
      pivots.push(col+1); row++;
    }
    write('rr-rank', String(pivots.length));
    write('rr-pivots', pivots.length ? pivots.join(', ') : 'none');
    write('rr-null', String(4-pivots.length));
    write('rr-rref', m.map(r=>'[ '+r.map(fmt).join('  ')+' ]').join('\n'));
  }
  document.querySelectorAll('.lab-grid input').forEach(input=>input.addEventListener('input',()=>{operator();rank();}));
  operator(); rank();
})();
