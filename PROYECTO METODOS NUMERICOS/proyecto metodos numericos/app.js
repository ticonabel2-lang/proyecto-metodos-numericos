// ========================================================
// CONTROLADOR DE NAVEGACIÓN (TABS)
// ========================================================
function showSection(sectionId) {
    const contents = document.querySelectorAll('.tab-content');
    contents.forEach(content => {
        content.classList.add('hidden');
    });

    const buttons = document.querySelectorAll('.nav-btn');
    buttons.forEach(btn => {
        btn.classList.remove('active');
    });

    const targetSection = document.getElementById(sectionId);
    if (targetSection) {
        targetSection.classList.remove('hidden');
    }

    const targetBtn = Array.from(buttons).find(b => {
        const clickAttr = b.getAttribute('onclick');
        return clickAttr && clickAttr.includes(sectionId);
    });
    if (targetBtn) {
        targetBtn.classList.add('active');
    }
}

// ========================================================
// ESCENARIO A: RED MATRICIAL CON ANÁLISIS DE ESTABILIDAD
// ========================================================
function calcularEscenarioA() {
    const A = [
        [parseFloat(document.getElementById('eA-A00').value), parseFloat(document.getElementById('eA-A01').value), parseFloat(document.getElementById('eA-A02').value)],
        [parseFloat(document.getElementById('eA-A10').value), parseFloat(document.getElementById('eA-A11').value), parseFloat(document.getElementById('eA-A12').value)],
        [parseFloat(document.getElementById('eA-A20').value), parseFloat(document.getElementById('eA-A21').value), parseFloat(document.getElementById('eA-A22').value)]
    ];
    const b = [
        parseFloat(document.getElementById('eA-b0').value),
        parseFloat(document.getElementById('eA-b1').value),
        parseFloat(document.getElementById('eA-b2').value)
    ];
    const omega = parseFloat(document.getElementById('eA-omega').value) || 1.25;

    // LU Descomposición
    let L = [[1,0,0],[0,1,0],[0,0,1]], U = [[0,0,0],[0,0,0],[0,0,0]];
    for(let i=0; i<3; i++) {
        for(let k=i; k<3; k++) {
            let s = 0; for(let j=0; j<i; j++) s += L[i][j]*U[j][k]; U[i][k] = A[i][k] - s;
        }
        for(let k=i; k<3; k++) {
            if(i!==k) {
                let s = 0; for(let j=0; j<i; j++) s += L[k][j]*U[j][i]; L[k][i] = (A[k][i] - s) / U[i][i];
            }
        }
    }
    let y = [b[0], (b[1] - L[1][0]*b[0]), (b[2] - L[2][0]*b[0] - L[2][1]*(b[1] - L[1][0]*b[0]))];
    let xLU = [0,0,0];
    xLU[2] = y[2]/U[2][2]; xLU[1] = (y[1] - U[1][2]*xLU[2])/U[1][1]; xLU[0] = (y[0] - U[0][1]*xLU[1] - U[0][2]*xLU[2])/U[0][0];

    // Jacobi
    let xJ = [0,0,0], xOld = [0,0,0];
    for(let it=0; it<15; it++) {
        xJ[0] = (b[0] - A[0][1]*xOld[1] - A[0][2]*xOld[2])/A[0][0];
        xJ[1] = (b[1] - A[1][0]*xOld[0] - A[1][2]*xOld[2])/A[1][1];
        xJ[2] = (b[2] - A[2][0]*xOld[0] - A[2][1]*xOld[1])/A[2][2];
        xOld = [...xJ];
    }

    // Gauss-Seidel
    let xGS = [0,0,0];
    for(let it=0; it<15; it++) {
        xGS[0] = (b[0] - A[0][1]*xGS[1] - A[0][2]*xGS[2])/A[0][0];
        xGS[1] = (b[1] - A[1][0]*xGS[0] - A[1][2]*xGS[2])/A[1][1];
        xGS[2] = (b[2] - A[2][0]*xGS[0] - A[2][1]*xGS[1])/A[2][2];
    }

    // SOR
    let xSOR = [0,0,0];
    for(let it=0; it<15; it++) {
        let t0 = (b[0] - A[0][1]*xSOR[1] - A[0][2]*xSOR[2])/A[0][0]; xSOR[0] = (1-omega)*xSOR[0] + omega*t0;
        let t1 = (b[1] - A[1][0]*xSOR[0] - A[1][2]*xSOR[2])/A[1][1]; xSOR[1] = (1-omega)*xSOR[1] + omega*t1;
        let t2 = (b[2] - A[2][0]*xSOR[0] - A[2][1]*xSOR[1])/A[2][2]; xSOR[2] = (1-omega)*xSOR[2] + omega*t2;
    }

    // Gradiente Conjugado
    let xGC = [0,0,0], r = [b[0], b[1], b[2]], p = [...r];
    for(let it=0; it<3; it++) {
        let Ap = [A[0][0]*p[0]+A[0][1]*p[1]+A[0][2]*p[2], A[1][0]*p[0]+A[1][1]*p[1]+A[1][2]*p[2], A[2][0]*p[0]+A[2][1]*p[1]+A[2][2]*p[2]];
        let rDotR = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
        if (rDotR < 1e-10) break;
        let alpha = rDotR / (p[0]*Ap[0] + p[1]*Ap[1] + p[2]*Ap[2]);
        xGC[0] += alpha*p[0]; xGC[1] += alpha*p[1]; xGC[2] += alpha*p[2];
        let rNew = [r[0]-alpha*Ap[0], r[1]-alpha*Ap[1], r[2]-alpha*Ap[2]];
        let beta = (rNew[0]*rNew[0] + rNew[1]*rNew[1] + rNew[2]*rNew[2]) / rDotR;
        p = [rNew[0]+beta*p[0], rNew[1]+beta*p[1], rNew[2]+beta*p[2]]; r = [...rNew];
    }

    const fmt = arr => arr.map((v,i)=>`Zona ${i+1}: ${v.toFixed(1)} Tons`).join('<br>');
    document.getElementById('eA-out-lu').innerHTML = fmt(xLU);
    document.getElementById('eA-out-jacobi').innerHTML = fmt(xJ);
    document.getElementById('eA-out-gs').innerHTML = fmt(xGS);
    document.getElementById('eA-out-sor').innerHTML = fmt(xSOR);
    document.getElementById('eA-out-gc').innerHTML = fmt(xGC);

    let diag = "<strong>📊 ANÁLISIS AUTOMÁTICO DE ESTABILIDAD LINEAL:</strong><br>";
    let dom = true;
    for(let i=0; i<3; i++) {
        let sum = 0; for(let j=0; j<3; j++) { if(i!==j) sum += Math.abs(A[i][j]); }
        if(Math.abs(A[i][i]) <= sum) dom = false;
    }
    diag += `• Dominancia Diagonal: ${dom ? '<span style="color:#10b981">✅ CONVERGENCIA GARANTIZADA</span>' : '<span style="color:#f43f5e">⚠️ DIAGONAL NO DOMINANTE - Posible divergencia</span>'}.<br>`;
    diag += `• Recomendación: ${dom ? 'Los métodos iterativos convergerán a la solución exacta.' : 'Utilice LU o Gradiente Conjugado para evitar problemas de convergencia.'}`;
    document.getElementById('eA-analisis-automatico').innerHTML = diag;
    document.getElementById('res-escenarioA').classList.remove('hidden');
}

// ========================================================
// ESCENARIO B: VACIADO DINÁMICO DE CARBURANTES
// ========================================================
let chartBInstance = null;
function calcularEscenarioB() {
    const V0 = parseFloat(document.getElementById('eB-v0').value);
    const inj = parseFloat(document.getElementById('eB-inj').value);
    const cons = parseFloat(document.getElementById('eB-cons').value);

    function dVdt(t, V) { return inj - (cons * (1 + 0.04 * t)); }

    let h = 0.5, t = 0, rE = V0, rH = V0, rRK = V0;
    let tX = [], yE = [], yH = [], yRK = [], rows = "";
    let tCritico = null;

    while(t <= 8) {
        tX.push(t.toFixed(1) + "h"); yE.push(Math.max(0, rE)); yH.push(Math.max(0, rH)); yRK.push(Math.max(0, rRK));
        rows += `<tr><td>${t.toFixed(1)}</td><td>${rE.toFixed(0)} L</td><td>${rH.toFixed(0)} L</td><td>${rRK.toFixed(0)} L</td></tr>`;
        if (rRK <= 3000 && tCritico === null) tCritico = t;

        // Euler
        rE += h * dVdt(t, rE);
        
        // Heun
        let sP = dVdt(t, rH); let rP = rH + h * sP; rH += (h/2)*(sP + dVdt(t+h, rP));
        
        // RK4
        let k1 = dVdt(t, rRK); 
        let k2 = dVdt(t+h/2, rRK+(h*k1)/2); 
        let k3 = dVdt(t+h/2, rRK+(h*k2)/2); 
        let k4 = dVdt(t+h, rRK+h*k3);
        rRK += (h/6)*(k1 + 2*k2 + 2*k3 + k4);
        
        t += h;
    }
    document.getElementById('eB-tabla').innerHTML = rows;

    document.getElementById('eB-analisis-automatico').innerHTML = `
        <strong>🚨 ANÁLISIS DE RESILIENCIA EN RESERVAS:</strong><br>
        • Método RK4 (más preciso) establece colapso crítico (<3000L) a las: <strong>${tCritico !== null ? tCritico.toFixed(1) + ' horas' : 'Estable (No colapsa en 8h)'}</strong>.<br>
        • Error comparativo: Euler subestima el tiempo de colapso, Heun mejora pero RK4 es el recomendado.
    `;

    const ctx = document.getElementById('chart-eB').getContext('2d');
    if(chartBInstance) chartBInstance.destroy();
    chartBInstance = new Chart(ctx, {
        type: 'line',
        data: {
            labels: tX,
            datasets: [
                { label: 'Euler (Error O(h))', data: yE, borderColor: '#f43f5e', borderDash:[5,5], fill:false, tension:0.1 },
                { label: 'Heun (Error O(h²))', data: yH, borderColor: '#f59e0b', borderDash:[3,3], fill:false, tension:0.1 },
                { label: 'RK4 (Error O(h⁴))', data: yRK, borderColor: '#10b981', borderWidth:2.5, fill:false, tension:0.1 }
            ]
        },
        options: { responsive: true, maintainAspectRatio: false, scales: { x:{ticks:{color:'#64748b'}}, y:{ticks:{color:'#64748b'}} }, plugins: { legend: { labels: { color: '#94a3b8' } } } }
    });
    document.getElementById('res-escenarioB').classList.remove('hidden');
}

// ========================================================
// ESCENARIO C: INTERPOLACIÓN HISTÓRICA DE PRECIOS
// ========================================================
const cX = [1, 5, 10, 15, 20, 30]; 
const cY = [15.5, 16.2, 19.8, 24.5, 22.0, 31.0];
let chartCInstance = null;

function calcularSplinePaceño(xEval) {
    let n = cX.length - 1; let h = new Array(n);
    for (let i = 0; i < n; i++) h[i] = cX[i+1] - cX[i];
    let A_mat = []; for (let i = 0; i <= n; i++) A_mat.push(new Array(n + 1).fill(0));
    let B_vec = new Array(n + 1).fill(0); A_mat[0][0] = 1; A_mat[n][n] = 1;

    for (let i = 1; i < n; i++) {
        A_mat[i][i-1] = h[i-1]/6; A_mat[i][i] = (h[i-1]+h[i])/3; A_mat[i][i+1] = h[i]/6;
        B_vec[i] = ((cY[i+1]-cY[i])/h[i]) - ((cY[i]-cY[i-1])/h[i-1]);
    }
    let M = new Array(n+1).fill(0), cp = new Array(n+1).fill(0), dp = new Array(n+1).fill(0);
    cp[0] = A_mat[0][1]/A_mat[0][0]; dp[0] = B_vec[0]/A_mat[0][0];
    for (let i = 1; i <= n; i++) {
        let m = A_mat[i][i-1]; let den = A_mat[i][i] - m * (i<=n ? cp[i-1] : 0);
        if (i < n) cp[i] = A_mat[i][i+1]/den; dp[i] = (B_vec[i]-m*dp[i-1])/den;
    }
    M[n] = dp[n]; for (let i = n-1; i >= 0; i--) M[i] = dp[i] - cp[i]*M[i+1];

    for (let i = 0; i < n; i++) {
        if (xEval >= cX[i] && xEval <= cX[i+1]) {
            let d_izq = cX[i+1]-xEval, d_der = xEval-cX[i];
            return (M[i]*Math.pow(d_izq,3)/(6*h[i])) + (M[i+1]*Math.pow(d_der,3)/(6*h[i])) + 
                   (cY[i]-M[i]*Math.pow(h[i],2)/6)*(d_izq/h[i]) + (cY[i+1]-M[i+1]*Math.pow(h[i],2)/6)*(d_der/h[i]);
        }
    } return cY[n];
}

function calcularEscenarioC() {
    let tgt = parseFloat(document.getElementById('eC-eval').value);
    let n = cX.length;

    // Lagrange
    let resL = 0; for(let i=0; i<n; i++) { let t = cY[i]; for(let j=0; j<n; j++) { if(j!==i) t *= (tgt-cX[j])/(cX[i]-cX[j]); } resL += t; }
    
    // Newton (diferencias divididas)
    let dd = []; for(let i=0; i<n; i++) dd.push(new Array(n).fill(0)); for(let i=0; i<n; i++) dd[i][0] = cY[i];
    for(let j=1; j<n; j++) { for(let i=0; i<n-j; i++) dd[i][j] = (dd[i+1][j-1]-dd[i][j-1])/(cX[i+j]-cX[i]); }
    let resN = dd[0][0]; let p = 1; for(let i=1; i<n; i++) { p *= (tgt-cX[i-1]); resN += dd[0][i]*p; }
    
    // Splines
    let resS = calcularSplinePaceño(tgt);

    document.getElementById('eC-out-lagrange').innerText = resL.toFixed(2) + " Bs";
    document.getElementById('eC-out-newton').innerText = resN.toFixed(2) + " Bs";
    document.getElementById('eC-out-splines').innerText = resS.toFixed(2) + " Bs";
    
    // Interpretación añadida
    document.getElementById('eC-interpretacion-dinamica').innerHTML = `
        <strong>📈 ANÁLISIS DE CURVA DE PRECIOS:</strong><br>
        • La interpolación por <strong>splines cúbicos naturales</strong> proporciona la transición más suave entre los datos de mercado, minimizando oscilaciones artificiales (fenómeno de Runge).<br>
        • El precio interpolado para el día <strong>${tgt}</strong> es de <strong>${resS.toFixed(2)} Bs</strong>.<br>
        • Diferencia máxima entre métodos: <strong>${Math.max(Math.abs(resL-resS), Math.abs(resN-resS)).toFixed(2)} Bs</strong>.
    `;

    let lX = [], lS = [];
    for(let d=1; d<=30; d++) { lX.push("Día "+d); lS.push(calcularSplinePaceño(d)); }

    const ctx = document.getElementById('chart-eC').getContext('2d');
    if(chartCInstance) chartCInstance.destroy();
    chartCInstance = new Chart(ctx, {
        type: 'line',
        data: { labels: lX, datasets: [{ label: 'Spline Cúbico Natural', data: lS, borderColor: '#10b981', fill: false, tension: 0.2 }] },
        options: { responsive: true, maintainAspectRatio: false, scales:{x:{ticks:{color:'#64748b'}},y:{ticks:{color:'#64748b'}}} }
    });
    document.getElementById('res-escenarioC').classList.remove('hidden');
}

// ========================================================
// ESCENARIO D: GASTO INTEGRADO COMPUESTO
// ========================================================
function calcularEscenarioD() {
    let a = 1, b = 30;
    let nT = 29, hT = (b-a)/nT, sT = calcularSplinePaceño(a) + calcularSplinePaceño(b);
    for(let i=1; i<nT; i++) sT += 2 * calcularSplinePaceño(a + i*hT);
    let rT = (hT/2)*sT;

    let n13 = 28, h13 = (b-a)/n13, s13 = calcularSplinePaceño(a) + calcularSplinePaceño(b);
    for(let i=1; i<n13; i++) s13 += (i%2===0 ? 2 : 4) * calcularSplinePaceño(a + i*h13);
    let r13 = (h13/3)*s13;

    let n38 = 27, h38 = (b-a)/n38, s38 = calcularSplinePaceño(a) + calcularSplinePaceño(b);
    for(let i=1; i<n38; i++) s38 += (i%3===0 ? 2 : 3) * calcularSplinePaceño(a + i*h38);
    let r38 = ((3*h38)/8)*s38;

    let costoBase = 15.5 * 30; // Precio día 1 * 30 días
    
    document.getElementById('eD-trapecio').innerHTML = rT.toFixed(2) + " Bs";
    document.getElementById('eD-simpson13').innerHTML = r13.toFixed(2) + " Bs";
    document.getElementById('eD-simpson38').innerHTML = r38.toFixed(2) + " Bs";
    
    document.getElementById('eD-interpretacion-dinamica').innerHTML = `
        <strong>📊 EVALUACIÓN INTEGRAL DEL GASTO FAMILIAR:</strong><br>
        • Gasto real estimado (Simpson 3/8): <strong>${r38.toFixed(2)} Bs</strong><br>
        • Gasto sin inflación (precio constante): <strong>${costoBase.toFixed(2)} Bs</strong><br>
        • <strong style="color:#f43f5e;">Pérdida del poder adquisitivo: ${(r38 - costoBase).toFixed(2)} Bs (${((r38-costoBase)/costoBase*100).toFixed(1)}% más)</strong><br>
        • El cálculo por regla de Simpson de orden superior minimiza el error de truncamiento global.
    `;
    document.getElementById('res-escenarioD').classList.remove('hidden');
}

// ========================================================
// ESCENARIO E: ANÁLISIS DE PUNTOS DE INFLEXIÓN
// ========================================================
function f_eE(x) { return Math.pow(x,3) - 4.5*Math.pow(x,2) - 10*x + 15; }
function df_eE(x) { return 3*Math.pow(x,2) - 9*x - 10; }
let chartEInstance = null;

function calcularEscenarioE() {
    let x0 = parseFloat(document.getElementById('eE-x0').value);
    let x1 = parseFloat(document.getElementById('eE-x1').value);
    const mIt = 8, rV = 1.094568;

    // Bisección
    let a = x0, b = x1, hB = "";
    for(let i=1; i<=mIt; i++) {
        let c = (a+b)/2; hB += `<tr><td>${i}</td><td>${c.toFixed(6)}</td><td>${Math.abs(rV-c).toExponential(2)}</td></tr>`;
        if(f_eE(a)*f_eE(c)<0) b=c; else a=c;
    }
    document.getElementById('eE-tab-biseccion').innerHTML = hB;

    // Newton-Raphson
    let cN = x0, hN = "";
    for(let i=1; i<=mIt; i++) {
        let nxt = cN - f_eE(cN)/df_eE(cN); hN += `<tr><td>${i}</td><td>${nxt.toFixed(6)}</td><td>${Math.abs(rV-nxt).toExponential(2)}</td></tr>`;
        cN = nxt;
    }
    document.getElementById('eE-tab-newton').innerHTML = hN;

    // Secante
    let s0 = x0, s1 = x1, hS = "";
    for(let i=1; i<=mIt; i++) {
        let nxt = s1 - (f_eE(s1)*(s1-s0))/(f_eE(s1)-f_eE(s0)); hS += `<tr><td>${i}</td><td>${nxt.toFixed(6)}</td><td>${Math.abs(rV-nxt).toExponential(2)}</td></tr>`;
        s0 = s1; s1 = nxt;
    }
    document.getElementById('eE-tab-secante').innerHTML = hS;

    // Análisis de convergencia
    document.getElementById('eE-orden-convergencia').innerHTML = `
        <strong>📊 ANÁLISIS DE CONVERGENCIA:</strong><br>
        • <strong>Newton-Raphson</strong>: Convergencia cuadrática (más rápida, requiere derivada).<br>
        • <strong>Bisección</strong>: Convergencia lineal (robusta pero lenta).<br>
        • <strong>Secante</strong>: Convergencia superlineal (≈1.618 - relación áurea).<br>
        • La raíz real en el intervalo es aproximadamente: <strong>${rV.toFixed(6)}</strong>
    `;

    let gX = [], gY = []; for(let x=-1; x<=4; x+=0.1) { gX.push(x.toFixed(1)); gY.push(f_eE(x)); }
    const ctx = document.getElementById('chart-eE').getContext('2d');
    if(chartEInstance) chartEInstance.destroy();
    chartEInstance = new Chart(ctx, {
        type: 'line',
        data: { labels: gX, datasets: [{ label: 'f(x) = x³ - 4.5x² - 10x + 15', data: gY, borderColor: '#f43f5e', fill: false, tension: 0.1 }] },
        options: { responsive: true, maintainAspectRatio: false, plugins: { tooltip: { callbacks: { label: (ctx) => `f(x) = ${ctx.raw.toFixed(4)}` } } } }
    });
    document.getElementById('res-escenarioE').classList.remove('hidden');
}

// ========================================================
// ESCENARIO F: MAL CONDICIONAMIENTO MATRICIAL
// ========================================================
let chartFInstance = null;
function calcularEscenarioF() {
    const dB = parseFloat(document.getElementById('eF-demanda').value);
    const pR = parseFloat(document.getElementById('eF-rumor').value);

    const A = [[1.0, 1.0, 1.0], [1.0, 1.01, 1.0], [1.0, 1.0, 1.005]];
    const bN = [dB, dB * 1.005, dB * 1.002];
    const bR = bN.map(v => v * (1 + pR/100));

    function solve3x3(M, v) {
        let det = M[0][0]*(M[1][1]*M[2][2] - M[1][2]*M[2][1]) - M[0][1]*(M[1][0]*M[2][2] - M[1][2]*M[2][0]) + M[0][2]*(M[1][0]*M[2][1] - M[1][1]*M[2][0]);
        if(Math.abs(det) < 1e-7) return [NaN, NaN, NaN];
        let dx = v[0]*(M[1][1]*M[2][2] - M[1][2]*M[2][1]) - M[0][1]*(v[1]*M[2][2] - M[1][2]*v[2]) + M[0][2]*(v[1]*M[2][1] - M[1][1]*v[2]);
        let dy = M[0][0]*(v[1]*M[2][2] - M[1][2]*v[2]) - v[0]*(M[1][0]*M[2][2] - M[1][2]*M[2][0]) + M[0][2]*(M[1][0]*v[2] - v[1]*M[2][0]);
        let dz = M[0][0]*(M[1][1]*v[2] - v[1]*M[2][1]) - M[0][1]*(M[1][0]*v[2] - v[1]*M[2][0]) + v[0]*(M[1][0]*M[2][1] - M[1][1]*M[2][0]);
        return [dx/det, dy/det, dz/det];
    }

    const solN = solve3x3(A, bN); const solR = solve3x3(A, bR);
    const cambioOriente = ((solR[0]-solN[0])/solN[0]*100).toFixed(1);
    const cambioOccidente = ((solR[1]-solN[1])/solN[1]*100).toFixed(1);
    const cambioValles = ((solR[2]-solN[2])/solN[2]*100).toFixed(1);

    document.getElementById('eF-tab-normal').innerHTML = `Oriente: ${solN[0].toFixed(0)} T | Occidente: ${solN[1].toFixed(0)} T | Valles: ${solN[2].toFixed(0)} T`;
    document.getElementById('eF-tab-rumor').innerHTML = `Oriente: ${solR[0].toFixed(0)} T | Occidente: ${solR[1].toFixed(0)} T | Valles: ${solR[2].toFixed(0)} T`;
    
    document.getElementById('eF-interpretacion').innerHTML = `
        <strong>⚠️ ANÁLISIS DE SENSIBILIDAD POR RUMORES:</strong><br>
        • Una perturbación del <strong>${pR}% en la demanda</strong> genera cambios extremos:<br>
        &nbsp;&nbsp;→ Oriente: <strong>${cambioOriente}%</strong><br>
        &nbsp;&nbsp;→ Occidente: <strong>${cambioOccidente}%</strong><br>
        &nbsp;&nbsp;→ Valles: <strong>${cambioValles}%</strong><br>
        • <strong style="color:#f43f5e;">Conclusión: El sistema es MAL CONDICIONADO</strong>. Un rumor pequeño amplifica masivamente el desabastecimiento.
    `;

    const ctx = document.getElementById('chart-eF').getContext('2d');
    if(chartFInstance) chartFInstance.destroy();
    chartFInstance = new Chart(ctx, {
        type: 'bar',
        data: { labels: ['Oriente', 'Occidente', 'Valles'], datasets: [{ label: 'Demanda Normal', data: solN, backgroundColor: '#3b82f6' }, { label: 'Demanda con Rumor', data: solR, backgroundColor: '#f59e0b' }] },
        options: { responsive: true, maintainAspectRatio: false, plugins: { legend: { labels: { color: '#94a3b8' } } } }
    });
    document.getElementById('res-escenarioF').classList.remove('hidden');
}

// ========================================================
// ESCENARIO G: CINÉTICA DE MASAS (RK4 VECTORIAL)
// ========================================================
let chartGInstance = null;
function calcularEscenarioG() {
    const N0 = parseFloat(document.getElementById('eG-neut').value);
    const M0 = parseFloat(document.getElementById('eG-manif').value);
    const D0 = parseFloat(document.getElementById('eG-media').value);

    const alpha = 0.0008, beta = 0.035, gamma = 0.015;

    function dNdt(N, M, D) { return -alpha * N * M + gamma * D; }
    function dMdt(N, M, D) { return alpha * N * M - beta * M * D; }
    function dDdt(N, M, D) { return beta * M * D - gamma * D; }

    let diasX = [], sN = [], sM = [], sD = [], h = 0.2, t = 0, rows = "";
    let cN = N0, cM = M0, cD = D0;

    while(t <= 10) {
        if (Math.abs(t - Math.round(t)) < 1e-5) {
            let dE = Math.round(t); diasX.push(`Día ${dE}`); sN.push(cN); sM.push(cM); sD.push(cD);
            rows += `<tr><td>${dE}</td><td>${cN.toFixed(0)}</td><td style="color:var(--accent-rose); font-weight:bold;">${cM.toFixed(0)}</td><td style="color:var(--accent-cyan); font-weight:bold;">${cD.toFixed(0)}</td></tr>`;
        }

        let k1N = dNdt(cN, cM, cD), k1M = dMdt(cN, cM, cD), k1D = dDdt(cN, cM, cD);
        let k2N = dNdt(cN+h*k1N/2, cM+h*k1M/2, cD+h*k1D/2), 
            k2M = dMdt(cN+h*k1N/2, cM+h*k1M/2, cD+h*k1D/2), 
            k2D = dDdt(cN+h*k1N/2, cM+h*k1M/2, cD+h*k1D/2);
        let k3N = dNdt(cN+h*k2N/2, cM+h*k2M/2, cD+h*k2D/2), 
            k3M = dMdt(cN+h*k2N/2, cM+h*k2M/2, cD+h*k2D/2), 
            k3D = dDdt(cN+h*k2N/2, cM+h*k2M/2, cD+h*k2D/2);
        let k4N = dNdt(cN+h*k3N, cM+h*k3M, cD+h*k3D), 
            k4M = dMdt(cN+h*k3N, cM+h*k3M, cD+h*k3D), 
            k4D = dDdt(cN+h*k3N, cM+h*k3M, cD+h*k3D);

        cN += (h/6)*(k1N + 2*k2N + 2*k3N + k4N); 
        cM += (h/6)*(k1M + 2*k2M + 2*k3M + k4M); 
        cD += (h/6)*(k1D + 2*k2D + 2*k3D + k4D);
        
        if(cN<0) cN=0; if(cM<0) cM=0; if(cD<0) cD=0;
        t += h;
    }
    document.getElementById('eG-tabla').innerHTML = rows;

    let evalHtml = "<strong>📢 CONCLUSIONES CINÉTICAS DE LA DINÁMICA SOCIAL (RK4 VECTORIAL):</strong><br>";
    if(sM[sM.length-1] > M0 * 2) {
        evalHtml += `⚠️ El motor RK4 detecta una <strong>aceleración crítica</strong> en la movilización, alcanzando los <strong>${sM[sM.length-1].toFixed(0)} manifestantes activos</strong> (${((sM[sM.length-1]/M0-1)*100).toFixed(0)}% de aumento). El desabastecimiento actúa como un catalizador social directo.<br>`;
        evalHtml += `💡 <strong>Recomendación:</strong> Activar mecanismos de diálogo y mediación antes del día 5 para evitar la masificación del conflicto.`;
    } else {
        evalHtml += `✅ El sistema socioeconómico converge de forma asintótica hacia una resolución pacífica controlada. Los mediadores logran contener la escalada.`;
    }
    document.getElementById('eG-interpretacion').innerHTML = evalHtml;

    const ctx = document.getElementById('chart-eG').getContext('2d');
    if(chartGInstance) chartGInstance.destroy();
    chartGInstance = new Chart(ctx, {
        type: 'line',
        data: {
            labels: diasX,
            datasets: [
                { label: 'Neutrales N(t)', data: sN, borderColor: '#64748b', fill:false, pointRadius: 2, tension: 0.2 },
                { label: 'Manifestantes M(t)', data: sM, borderColor: '#f43f5e', borderWidth: 2.5, fill:false, pointRadius: 3, tension: 0.2 },
                { label: 'Mediadores D(t)', data: sD, borderColor: '#06b6d4', borderWidth: 1.5, fill:false, pointRadius: 2, tension: 0.2 }
            ]
        },
        options: { 
            responsive: true, 
            maintainAspectRatio: false,
            plugins: { legend: { labels: { color: '#94a3b8', font: { size: 11 } } } }, 
            scales: { 
                x: { grid: { color: 'rgba(255,255,255,0.03)' }, ticks: { color: '#64748b' } }, 
                y: { grid: { color: 'rgba(255,255,255,0.03)' }, ticks: { color: '#64748b' } } 
            } 
        }
    });
    document.getElementById('res-escenarioG').classList.remove('hidden');
}