// Gráficas interactivas del visor (Plotly). Tomado de web/graficas.js del visor de la compañera
// (github.com/anajul1a4/Hidrologia-Chicama-El-Tambo), sin su lógica de pestañas, tema y visor de imágenes,
// que aquí resuelve 22_visor_html.py. Datos en window.DATOS (los genera 22_visor_html.py).
// Uso: Graficas.dibujarEn(seccion) al mostrar una sección; Graficas.redibujar() al cambiar el tema.
(function () {
  "use strict";
  const D = window.DATOS;
  const MESES = ["Ene", "Feb", "Mar", "Abr", "May", "Jun", "Jul", "Ago", "Sep", "Oct", "Nov", "Dic"];
  const UNID = { P_L: "mm/mes", P_I: "mm/mes", R: "mm/mes", Q: "m³/s", T: "°C" };
  const NOMBRE = { P_L: "P<sub>L</sub>", P_I: "P<sub>I</sub>", R: "R", Q: "Q̄", T: "T" };
  const estado = { clim: "P_L", anom: "P_L", qc: "P_L", modelos: "PL", tend: "P_L", anios: false,
                   serie11: new Set(["P_L", "Q", "R", "T"]) };
  const dibujadas = new Set();

  // ------------------------------------------------------------ tema (tokens del CSS del visor)
  function tema() {
    const cs = getComputedStyle(document.documentElement);
    const v = (n) => cs.getPropertyValue(n).trim();
    return { card: v("--card"), ink: v("--ink"), ink2: v("--ink2"), muted: v("--muted"), line: v("--line"),
             axis: v("--axis"), pos: v("--pos"), neg: v("--neg"), excl: v("--excl"), vPL: v("--vPL"),
             vPI: v("--vPI"), vQ: v("--vQ"), vT: v("--vT"), vTE: v("--vTE"), apoyo: v("--apoyo") };
  }
  function colVar(t, v) {
    return { P_L: t.vPL, P_I: t.vPI, Q: t.vQ, R: t.vQ, T: t.vT, T_PISCOt: t.vT, T_ERA5L: t.vTE, PL: t.vPL }[v] || t.ink;
  }
  function rgba(hex, a) {
    const h = hex.replace("#", "");
    const n = parseInt(h.length === 3 ? h.split("").map((c) => c + c).join("") : h, 16);
    return `rgba(${(n >> 16) & 255},${(n >> 8) & 255},${n & 255},${a})`;
  }
  function eje(t, extra) {
    return Object.assign({ gridcolor: t.line, linecolor: t.axis, zerolinecolor: t.axis, showline: true,
      ticks: "outside", tickcolor: t.axis, ticklen: 4, tickfont: { color: t.muted }, automargin: true,
      title: { font: { color: t.ink2, size: 12 } } }, extra || {});
  }
  function base(t, extra) {
    return Object.assign({
      paper_bgcolor: t.card, plot_bgcolor: t.card,
      font: { family: 'system-ui,-apple-system,"Segoe UI",Roboto,sans-serif', size: 12, color: t.ink2 },
      margin: { l: 60, r: 18, t: 40, b: 40 }, hovermode: "x unified",
      hoverlabel: { bgcolor: t.card, bordercolor: t.line, font: { color: t.ink, size: 12 } },
      legend: { orientation: "h", x: 0, y: 1.02, yanchor: "bottom", font: { color: t.ink2 }, bgcolor: "rgba(0,0,0,0)" },
    }, extra);
  }
  const CONFIG = { responsive: true, displaylogo: false, scrollZoom: false,
    modeBarButtonsToRemove: ["lasso2d", "select2d", "autoScale2d"],
    toImageButtonOptions: { format: "png", scale: 2, filename: "el_tambo" } };
  function ejeTiempo(t, extra) {
    return eje(t, Object.assign({ type: "date", hoverformat: "%Y-%m", tickformat: "%Y" }, extra || {}));
  }
  function slider(t) {
    return { visible: true, thickness: 0.07, bgcolor: t.card, bordercolor: t.line, borderwidth: 1 };
  }
  function mesSiguiente(f) {
    let [y, m] = f.split("-").map(Number);
    m += 1; if (m > 12) { m = 1; y += 1; }
    return `${y}-${String(m).padStart(2, "0")}`;
  }
  function linea(x, y, nombre, color, extra) {
    return Object.assign({ type: "scatter", mode: "lines", x, y, name: nombre, line: { color, width: 2 },
      connectgaps: false }, extra || {});
  }

  // ------------------------------------------------------------ gráficas
  const G = {};

  // 1.1 · cuatro paneles alineados
  G.serie11 = function (t) {
    const x = D.fechas;
    const orden = ["P_L", "Q", "R", "T"].filter((v) => estado.serie11.has(v));
    const k = orden.length, gap = 0.05, h = (1 - 0.02 - gap * (k - 1)) / k;
    const dom = {}, ax = {};
    orden.forEach((v, i) => {
      const top = 1 - i * (h + gap);
      dom[v] = [top - h, top];
      ax[v] = i === 0 ? "y" : "y" + (i + 1);
    });
    const abajo = ax[orden[k - 1]];
    const data = [];
    const PANEL = {
      P_L: ["P<sub>L</sub> · precipitación acumulada mensual (CAMELS-PE)", "mm/mes"],
      Q: ["Q̄ · caudal medio mensual (SENAMHI, rellenado)", "m³/s"],
      R: ["R · escorrentía en lámina", "mm/mes"],
      T: ["T · temperatura media mensual", "°C"],
    };
    if (ax.P_L) data.push({ type: "bar", x, y: D.P_L, name: "P<sub>L</sub>", yaxis: ax.P_L, marker: { color: t.vPL },
      showlegend: false, hovertemplate: "%{y:.1f} mm/mes" });
    if (ax.Q) data.push(linea(x, D.Q, "Q̄", t.vQ, { yaxis: ax.Q, showlegend: false, line: { color: t.vQ, width: 1.6 },
      hovertemplate: "%{y:.2f} m³/s" }));
    if (ax.R) data.push(linea(x, D.R, "R", t.vQ, { yaxis: ax.R, showlegend: false, line: { color: t.vQ, width: 1.6 },
      hovertemplate: "%{y:.1f} mm/mes" }));
    if (ax.T) {
      data.push(linea(x, D.T_piscot, "T PISCOt", t.vT, { yaxis: ax.T, line: { color: t.vT, width: 1.6 },
        hovertemplate: "%{y:.2f} °C" }));
      data.push(linea(x, D.T_est, "T estimada (ERA5-Land + sesgo)", t.vT, { yaxis: ax.T,
        line: { color: t.vT, width: 1.8, dash: "dash" }, hovertemplate: "%{y:.2f} °C" }));
    }
    // sombreados los meses de caudal que sin relleno se excluían
    const shapes = [];
    ["Q", "R"].filter((v) => ax[v]).forEach((v) => {
      D.R_flag_obs.forEach((f, i) => {
        if (f === "excluido" || f === "sin_dato") {
          shapes.push({ type: "rect", xref: "x", yref: "paper", x0: x[i], x1: mesSiguiente(x[i]), y0: dom[v][0],
            y1: dom[v][1], fillcolor: t.excl, line: { width: 0 }, layer: "below" });
        }
      });
    });
    const annotations = orden.map((v) => ({ text: PANEL[v][0], xref: "paper", yref: "paper", x: 0, y: dom[v][1],
      xanchor: "left", yanchor: "bottom", showarrow: false, font: { color: t.ink, size: 12.5 } }));
    const layout = base(t, {
      height: 150 + 190 * k, margin: { l: 62, r: 18, t: 56, b: 30 }, bargap: 0.15, shapes, annotations,
      legend: { orientation: "h", x: 0, y: 1 + 34 / (150 + 190 * k), yanchor: "bottom", font: { color: t.ink2 } },
      xaxis: ejeTiempo(t, { anchor: abajo, rangeslider: slider(t) }),
    });
    orden.forEach((v) => {
      const nom = ax[v] === "y" ? "yaxis" : "yaxis" + ax[v].slice(1);
      layout[nom] = eje(t, { domain: dom[v], title: { text: PANEL[v][1] }, rangemode: v === "T" ? "normal" : "tozero" });
    });
    return { data, layout };
  };

  // 1.4 · disponibilidad (año × mes) y serie diaria de la variable elegida
  G.qc = function (t) {
    const v = estado.qc, dsp = D.disp[v], dd = D.diario;
    const anios = D.disp.anios;
    const nombreCond = { completo: "completo", escalado: "escalado (incompleto)", excluido: "excluido (> 10 % faltante)",
                         sin_dato: "sin dato" };
    const COND = ["completo", "escalado", "excluido", "sin_dato"];
    const COL = ["#256abf", "#9ec5f4", "#d03b3b", "#d9d9d9"];
    const diasMes = (y, m) => new Date(Date.UTC(y, m, 0)).getUTCDate();
    const z = [], etiqueta = [], hover = [];
    dsp.f.forEach((fila, i) => {
      z.push([]); etiqueta.push([]); hover.push([]);
      fila.forEach((f, j) => {
        const c = dsp.c[i][j];
        if (f == null || c == null) { z[i].push(null); etiqueta[i].push(""); hover[i].push("fuera del periodo común"); return; }
        const n = Math.round(f * diasMes(anios[i], j + 1));
        z[i].push(COND.indexOf(c));
        etiqueta[i].push(c === "completo" || c === "sin_dato" || v === "P_I" ? "" : String(n));
        hover[i].push(`${anios[i]}-${String(j + 1).padStart(2, "0")}<br>${n} días con dato (${Math.round(100 * f)} %)<br>` +
                      `condición: ${nombreCond[c] || c}`);
      });
    });
    const escala = [];
    COL.forEach((col, k) => { escala.push([k / 4, col], [(k + 1) / 4, col]); });
    const data = [{
      type: "heatmap", x: MESES, y: anios, z, customdata: hover, text: etiqueta, texttemplate: "%{text}",
      textfont: { size: 9, color: t.ink }, hovertemplate: "%{customdata}<extra></extra>",
      zmin: -0.5, zmax: 3.5, xgap: 1, ygap: 1, xaxis: "x", yaxis: "y", colorscale: escala,
      colorbar: { title: { text: "condición del mes<br>(antes del relleno)", font: { size: 11, color: t.ink2 } }, len: 0.52,
                  y: 1, yanchor: "top", thickness: 12, tickvals: [0, 1, 2, 3],
                  ticktext: ["completo", "escalado (≤ 10 % faltante)", "excluido (> 10 % faltante)", "sin dato"],
                  tickfont: { color: t.muted }, outlinewidth: 0 },
    }];
    const layout = base(t, {
      height: 640, hovermode: "closest", margin: { l: 56, r: 18, t: 30, b: 40 }, showlegend: v === "Q",
      legend: { orientation: "h", x: 0, y: 0.36, yanchor: "bottom", font: { color: t.ink2 } },
      xaxis: eje(t, { type: "category", domain: [0, 1], anchor: "y", showgrid: false, side: "top" }),
      yaxis: eje(t, { domain: [0.48, 1], autorange: "reversed", dtick: 2, showgrid: false, title: { text: "Año" } }),
      annotations: [],
    });
    if (v === "P_I") {
      layout.annotations.push({ text: `IMERG es un producto mensual: no tiene serie diaria. Los ${D.fechas.length} meses están completos.`,
        xref: "paper", yref: "paper", x: 0, y: 0.18, showarrow: false, xanchor: "left", font: { color: t.ink2, size: 13 } });
      return { data, layout };
    }
    const f0 = new Date(dd.ini + "T00:00:00Z");
    const xd = Array.from({ length: dd.n }, (_, i) => new Date(f0.getTime() + i * 864e5).toISOString().slice(0, 10));
    const u = { P_L: "mm/día", Q: "m³/s", T: "°C" }[v];
    const etq = { P_L: "P<sub>L</sub> diaria (PISCOp v2.1)", Q: "Caudal diario observado (vacíos sin rellenar)",
                  T: "T media diaria (PISCOt; termina en 2020-12)" }[v];
    data.push({ type: "scattergl", mode: "lines", x: xd, y: dd[v], name: v === "Q" ? "observado" : etq,
      xaxis: "x2", yaxis: "y2", connectgaps: false, line: { color: colVar(t, v), width: 1 },
      hovertemplate: `%{x}<br>%{y:.2f} ${u}<extra></extra>` });
    if (v === "Q") {
      data.push({ type: "scattergl", mode: "markers", x: xd, y: dd.Q_int, name: "días interpolados (relleno)",
        xaxis: "x2", yaxis: "y2", marker: { color: t.ink, size: 4 },
        hovertemplate: `%{x}<br>%{y:.2f} ${u} (interpolado)<extra></extra>` });
    }
    layout.xaxis2 = ejeTiempo(t, { domain: [0, 1], anchor: "y2", hoverformat: "%Y-%m-%d" });
    layout.yaxis2 = eje(t, { domain: [0, 0.34], title: { text: u }, rangemode: v === "T" ? "normal" : "tozero" });
    layout.annotations.push({ text: etq, xref: "paper", yref: "paper", x: 0, y: 0.345, xanchor: "left",
      yanchor: "bottom", showarrow: false, font: { color: t.ink, size: 12.5 } });
    return { data, layout };
  };

  // 1.1 · temperatura
  G.temp = function (t) {
    const x = D.fechas;
    // (b) ciclo anual medio 1998–2020 (meses con PISCOt y ERA5-Land a la vez)
    const ciclo = (serie) => {
      const suma = Array(12).fill(0), n = Array(12).fill(0);
      x.forEach((f, i) => {
        const [y, m] = f.split("-").map(Number);
        if (y >= 1998 && y <= 2020 && serie[i] != null && D.T_piscot[i] != null && D.T_era5[i] != null) {
          suma[m - 1] += serie[i]; n[m - 1] += 1;
        }
      });
      return suma.map((v, j) => (n[j] ? Math.round(100 * v / n[j]) / 100 : null));
    };
    const h2 = "%{x}: %{y:.2f} °C";
    return { data: [
      linea(x, D.T_piscot, "PISCOt v1.2 (CAMELS-PE)", t.vT, { legendgroup: "p", hovertemplate: "%{y:.2f} °C" }),
      linea(x, D.T_era5, "ERA5-Land sin corregir", t.vTE, { legendgroup: "e", line: { color: t.vTE, width: 1.6 },
        hovertemplate: "%{y:.2f} °C" }),
      linea(x, D.T_est, "Estimada: ERA5-Land + sesgo mensual", t.vT, { line: { color: t.vT, width: 2, dash: "dash" },
        hovertemplate: "%{y:.2f} °C" }),
      { type: "scatter", mode: "lines+markers", x: MESES, y: ciclo(D.T_piscot), xaxis: "x2", yaxis: "y2",
        name: "PISCOt", legendgroup: "p", showlegend: false, line: { color: t.vT, width: 2 },
        marker: { size: 7, color: t.vT, line: { color: t.card, width: 1.5 } }, hovertemplate: "PISCOt " + h2 + "<extra></extra>" },
      { type: "scatter", mode: "lines+markers", x: MESES, y: ciclo(D.T_era5), xaxis: "x2", yaxis: "y2",
        name: "ERA5-Land", legendgroup: "e", showlegend: false, line: { color: t.vTE, width: 2 },
        marker: { size: 7, color: t.vTE, line: { color: t.card, width: 1.5 } }, hovertemplate: "ERA5-Land " + h2 + "<extra></extra>" },
    ], layout: base(t, { height: 470, hovermode: "x",
      xaxis: ejeTiempo(t, { domain: [0, 0.7], rangeslider: slider(t) }),
      yaxis: eje(t, { title: { text: "Temperatura media mensual (°C)" } }),
      xaxis2: eje(t, { domain: [0.78, 1], anchor: "y2", type: "category" }),
      yaxis2: eje(t, { anchor: "x2", title: { text: "°C" } }),
      annotations: [
        { text: "(a) Serie mensual", xref: "paper", yref: "paper", x: 0, y: 1.0, xanchor: "left", yanchor: "bottom",
          showarrow: false, font: { color: t.ink, size: 12.5 }, yshift: 26 },
        { text: "(b) Ciclo anual medio 1998–2020", xref: "paper", yref: "paper", x: 0.78, y: 1.0, xanchor: "left",
          yanchor: "bottom", showarrow: false, font: { color: t.ink, size: 12.5 }, yshift: 26 }],
      margin: { l: 60, r: 18, t: 70, b: 40 } }) };
  };

  // 1.2 · fuentes de precipitación
  G.precip = function (t) {
    const x = D.fechas;
    const h = "%{y:.1f} mm/mes";
    return { data: [
      linea(x, D.P_L, "P<sub>L</sub> CAMELS-PE", t.vPL, { hovertemplate: h }),
      linea(x, D.P_I, "P<sub>I</sub> IMERG V07", t.vPI, { hovertemplate: h }),
      linea(x, D.chirps, "CHIRPS v2.0", t.apoyo, { visible: "legendonly", line: { color: t.apoyo, width: 1.6 }, hovertemplate: h }),
      linea(x, D.piscop, "PISCOp v3.0", t.muted, { visible: "legendonly", line: { color: t.muted, width: 1.6, dash: "dot" },
        hovertemplate: h }),
    ], layout: base(t, { height: 480, xaxis: ejeTiempo(t, { rangeslider: slider(t) }),
      yaxis: eje(t, { title: { text: "Precipitación (mm/mes)" }, rangemode: "tozero" }) }) };
  };

  // 1.2 · totales anuales
  G.anual = function (t) {
    const a = D.anual;
    const h = "%{y:.0f} mm/año";
    return { data: [
      { type: "bar", x: a.anios, y: a.P_L, name: "P<sub>L</sub> CAMELS-PE", marker: { color: t.vPL, cornerradius: 3 }, hovertemplate: h },
      { type: "bar", x: a.anios, y: a.P_I, name: "P<sub>I</sub> IMERG V07", marker: { color: t.vPI, cornerradius: 3 }, hovertemplate: h },
    ], layout: base(t, { height: 420, barmode: "group", bargap: 0.25, bargroupgap: 0.1,
      xaxis: eje(t, { dtick: 2, tickangle: 0 }),
      yaxis: eje(t, { title: { text: "Precipitación anual (mm/año)" }, rangemode: "tozero" }) }) };
  };

  // 1.5.a · climatología con bandas
  G.clim = function (t) {
    const v = estado.clim, c = D.clim[v], u = UNID[v], cv = colVar(t, v);
    const f = v === "T" ? ".2f" : (v === "Q" ? ".1f" : ".0f");
    const ht = (n) => `${n}: %{y:${f}} ${u}<extra></extra>`;
    const sin = { width: 0 };
    return { data: [
      { type: "scatter", mode: "lines", x: MESES, y: c.P90, line: sin, name: "P90", hovertemplate: ht("P90"), showlegend: false },
      { type: "scatter", mode: "lines", x: MESES, y: c.P10, line: sin, fill: "tonexty", fillcolor: rgba(cv, 0.14),
        name: "Percentiles 10–90", hovertemplate: ht("P10") },
      { type: "scatter", mode: "lines", x: MESES, y: c.Q3, line: sin, name: "Q3", hovertemplate: ht("Q3"), showlegend: false },
      { type: "scatter", mode: "lines", x: MESES, y: c.Q1, line: sin, fill: "tonexty", fillcolor: rgba(cv, 0.30),
        name: "Rango intercuartílico", hovertemplate: ht("Q1") },
      { type: "scatter", mode: "lines", x: MESES, y: c.mediana, name: "Mediana",
        line: { color: cv, width: 2, dash: "dash" }, hovertemplate: ht("Mediana") },
      { type: "scatter", mode: "lines+markers", x: MESES, y: c.media, name: "Media",
        line: { color: cv, width: 2.4 }, marker: { size: 8, color: cv, line: { color: t.card, width: 2 } },
        hovertemplate: ht("Media") },
    ], layout: base(t, { height: 460, xaxis: eje(t, { type: "category" }),
      yaxis: eje(t, { title: { text: `${NOMBRE[v]} (${u})` }, rangemode: v === "T" ? "normal" : "tozero" }) }) };
  };

  function colorSigno(t, v, y) {
    return y.map((a) => (a == null ? t.muted : ((a >= 0) !== (v === "T") ? t.pos : t.neg)));
  }

  // 1.5.b · anomalías mensuales
  G.anom = function (t) {
    const v = estado.anom, a = D.anom[v], u = UNID[v];
    const f = v === "T" ? ".2f" : ".1f";
    return { data: [
      { type: "bar", x: D.fechas, y: a.m, name: "Anomalía mensual", marker: { color: colorSigno(t, v, a.m) },
        hovertemplate: `%{y:+${f}} ${u}` },
      linea(D.fechas, a.mm12, "Media móvil de 12 meses", t.ink, { line: { color: t.ink, width: 1.8 },
        hovertemplate: `%{y:+${f}} ${u}` }),
    ], layout: base(t, { height: 440, bargap: 0, xaxis: ejeTiempo(t, { rangeslider: slider(t) }),
      yaxis: eje(t, { title: { text: `Anomalía de ${NOMBRE[v]} (${u})` }, zeroline: true }) }) };
  };

  // 1.5.b · anomalía anual
  G.anomanual = function (t) {
    const v = estado.anom, aa = D.anom_anual;
    const u = v === "T" ? "°C (media del año)" : "mm/año";
    const f = v === "T" ? ".2f" : ".0f";
    return { data: [
      { type: "bar", x: aa.anios, y: aa[v], name: "Anomalía anual", marker: { color: colorSigno(t, v, aa[v]), cornerradius: 3 },
        hovertemplate: `%{y:+${f}} ${v === "T" ? "°C" : "mm"}<extra></extra>` },
    ], layout: base(t, { height: 400, hovermode: "closest", showlegend: false, bargap: 0.2,
      xaxis: eje(t, { dtick: 2 }), yaxis: eje(t, { title: { text: `Anomalía de ${NOMBRE[v]} (${u})` }, zeroline: true }) }) };
  };

  // 2.3 · observado frente a modelado en los años de evaluación
  G.modelos = function (t) {
    const m = D.modelos, v = estado.modelos, d = m[v];
    const u = v === "PL" ? "mm/mes" : "m³/s";
    const estilos = [{ color: t.muted, width: 1.4, dash: "dot" }, { color: t.ink2, width: 1.4, dash: "dash" },
                     { color: colVar(t, v), width: 2.2 }];   // referencia, lineal simple, modelo elegido
    const data = [linea(m.fechas, d.obs, v === "PL" ? "Observado (P<sub>L</sub>)" : "Observado (Q̄)", t.ink,
      { line: { color: t.ink, width: 1.6 }, hovertemplate: `%{y:.1f} ${u}` })];
    d.mods.forEach((mo, i) => data.push(linea(m.fechas, mo.y, mo.nombre, estilos[i].color,
      { line: estilos[i], hovertemplate: `%{y:.1f} ${u}`, visible: i === 1 ? "legendonly" : true })));
    const sombra = (x0, x1) => ({ type: "rect", xref: "x", yref: "paper", x0, x1, y0: 0, y1: 1,
      fillcolor: rgba(t.muted, 0.10), line: { width: 0 }, layer: "below" });
    return { data, layout: base(t, { height: 470, shapes: [sombra("1998-01", "2007-01"), sombra("2016-01", "2025-01")],
      annotations: [["2002-06", "evaluación B"], ["2020-06", "evaluación A"], ["2011-06", "sin evaluar"]].map(([x, txt]) =>
        ({ x, y: 1, xref: "x", yref: "paper", yanchor: "bottom", text: txt, showarrow: false, font: { color: t.muted, size: 11 } })),
      xaxis: ejeTiempo(t, { rangeslider: slider(t) }),
      yaxis: eje(t, { title: { text: v === "PL" ? "P<sub>L</sub> (mm/mes)" : "Q̄ (m³/s)" }, rangemode: "tozero" }) }) };
  };

  // 3.4 · anomalías con la tendencia de Sen
  G.tend = function (t) {
    const v = estado.tend, d = D.tend[v];
    const esT = v.startsWith("T");
    const f = esT ? ".2f" : ".1f";
    return { data: [
      { type: "bar", x: d.fechas, y: d.a, name: "Anomalía mensual", marker: { color: colorSigno(t, esT ? "T" : v, d.a) },
        hovertemplate: `%{y:+${f}} ${d.unid}` },
      linea(d.fechas, d.mm12, "Media móvil de 12 meses", t.muted, { line: { color: t.ink2, width: 1.4 },
        hovertemplate: `%{y:+${f}} ${d.unid}` }),
      { type: "scatter", mode: "lines", x: d.sen_x, y: d.sen, name: d.txt, line: { color: t.ink, width: 2.6 },
        hoverinfo: "skip" },
    ], layout: base(t, { height: 450, bargap: 0, xaxis: ejeTiempo(t, { rangeslider: slider(t) }),
      yaxis: eje(t, { title: { text: `Anomalía (${d.unid})` }, zeroline: true }) }) };
  };

  // ------------------------------------------------------------ dibujo
  function dibujar(el) {
    const f = G[el.dataset.plot];
    if (!f || !window.Plotly) return;
    const g = f(tema());
    Plotly.react(el, g.data, g.layout, CONFIG);
    dibujadas.add(el);
  }

  document.querySelectorAll(".selector").forEach((sel) => {
    const para = sel.dataset.para;
    sel.addEventListener("click", (e) => {
      const b = e.target.closest("button");
      if (!b) return;
      if (sel.classList.contains("multi")) {                    // varios a la vez; siempre queda al menos uno
        const on = b.getAttribute("aria-pressed") === "true";
        if (on && estado[para].size === 1) return;
        b.setAttribute("aria-pressed", String(!on));
        if (on) estado[para].delete(b.dataset.v); else estado[para].add(b.dataset.v);
      } else {
        sel.querySelectorAll("button[data-v]").forEach((x) => x.setAttribute("aria-pressed", String(x === b)));
        estado[para] = b.dataset.v;
        if (para === "qc") {
          document.querySelectorAll(".qc-var").forEach((d) => { d.hidden = d.dataset.var !== b.dataset.v; });
        }
      }
      document.querySelectorAll(`.plot[data-plot="${para}"], .plot[data-plot="${para}anual"]`).forEach(dibujar);
    });
  });

  window.Graficas = {
    // dibuja (o ajusta al ancho) las gráficas de una sección recién mostrada
    dibujarEn(sec) {
      sec.querySelectorAll(".plot").forEach((el) => {
        if (dibujadas.has(el)) Plotly.Plots.resize(el); else dibujar(el);
      });
    },
    redibujar() { dibujadas.forEach(dibujar); },
  };
})();
