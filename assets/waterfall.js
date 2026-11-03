// Explorador interactivo de descomposición aditiva SHAP (Waterfall) por secuencia
(function () {
  const wfData = window.TARGETCNV_WF_DATA || [];
  const filterSelect = document.getElementById("wf-class-filter");
  const rowSelect = document.getElementById("wf-row-select");
  const prevBtn = document.getElementById("wf-prev-btn");
  const nextBtn = document.getElementById("wf-next-btn");
  const metaDiv = document.getElementById("wf-meta-info");

  if (!filterSelect || !rowSelect || !prevBtn || !nextBtn || !metaDiv) {
    return;
  }

  function renderWaterfall(rowNum) {
    const d = wfData.find(function (x) { return x.row === rowNum; }) || wfData[0];
    if (!d) {
      return;
    }

    const badgeCls = d.clase === "TP" ? "badge-tp" : "badge-fp";
    metaDiv.innerHTML =
      "<b>Fila #" + d.row + "</b> | Clase Real: <span class='" + badgeCls + "'>" + d.clase + "</span>" +
      " | <b>Locus:</b> <code>" + d.gene + "</code>" +
      " | <b>Ortogrupo:</b> <code>" + d.og + "</code>" +
      " | <b>ID:</b> <code>" + d.seq_id + "</code>" +
      " | <b>E[f(X)]:</b> " + d.base.toFixed(3) +
      " &rarr; <b>P(TP):</b> <b>" + (d.fx * 100).toFixed(2) + "%</b>";

    const ord = d.shaps
      .map(function (v, i) { return { i: i, abs: Math.abs(v) }; })
      .sort(function (a, b) { return a.abs - b.abs; })
      .map(function (o) { return o.i; });

    const yLabels = ord.map(function (i) { return d.features[i] + " = " + d.values[i]; });
    const xVals = ord.map(function (i) { return d.shaps[i]; });
    const txt = xVals.map(function (v) { return (v >= 0 ? "+" : "") + v.toFixed(4); });

    if (window.Plotly) {
      const trace = {
        type: "waterfall",
        orientation: "h",
        y: yLabels,
        x: xVals,
        text: txt,
        textposition: "outside",
        connector: { line: { color: "#7f8c8d", width: 1, dash: "dot" } },
        increasing: { marker: { color: "#27ae60" } },
        decreasing: { marker: { color: "#c0392b" } },
        base: d.base
      };

      const layout = {
        title: {
          text: "<b>Descomposición Aditiva SHAP: " + d.seq_id + " (" + d.gene + ")</b>",
          font: { size: 14 }
        },
        xaxis: {
          title: "Contribución sobre probabilidad P(TP)",
          range: [-0.1, 1.1]
        },
        margin: { l: 185, r: 45, t: 55, b: 50 }
      };

      Plotly.react("wf-plotly-container", [trace], layout, { responsive: true });
    }
  }

  function populateDropdown(filterValue) {
    const currentVal = rowSelect.value;
    rowSelect.innerHTML = "";

    const filtered = wfData.filter(function (d) {
      return filterValue === "ALL" || d.clase === filterValue;
    });

    filtered.forEach(function (d) {
      const opt = document.createElement("option");
      opt.value = d.row;
      opt.textContent =
        "Fila #" + d.row + " [" + d.clase + "] - " +
        d.gene + " | " + d.og + " (" + d.seq_id + ") | P(TP)=" +
        (d.fx * 100).toFixed(1) + "%";
      rowSelect.appendChild(opt);
    });

    if (filtered.some(function (d) { return String(d.row) === String(currentVal); })) {
      rowSelect.value = currentVal;
    } else if (filtered.length > 0) {
      rowSelect.value = filtered[0].row;
    }

    renderWaterfall(parseInt(rowSelect.value, 10));
  }

  filterSelect.addEventListener("change", function () {
    populateDropdown(filterSelect.value);
  });

  rowSelect.addEventListener("change", function () {
    renderWaterfall(parseInt(rowSelect.value, 10));
  });

  prevBtn.addEventListener("click", function () {
    if (rowSelect.selectedIndex > 0) {
      rowSelect.selectedIndex -= 1;
      renderWaterfall(parseInt(rowSelect.value, 10));
    }
  });

  nextBtn.addEventListener("click", function () {
    if (rowSelect.selectedIndex < rowSelect.options.length - 1) {
      rowSelect.selectedIndex += 1;
      renderWaterfall(parseInt(rowSelect.value, 10));
    }
  });

  setTimeout(function () {
    populateDropdown("ALL");
  }, 120);
})();

