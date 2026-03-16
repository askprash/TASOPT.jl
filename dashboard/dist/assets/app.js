(function () {
const COLORS = ["#71b0d4", "#df8f6a"];
const FILL_COLORS = ["rgba(95, 160, 199, 0.18)", "rgba(215, 123, 82, 0.18)"];
const MAP_FILL_COLORS = [
    [[0, "rgba(68, 119, 153, 0.18)"], [1, "rgba(126, 208, 154, 0.82)"]],
    [[0, "rgba(119, 78, 69, 0.18)"], [1, "rgba(215, 180, 106, 0.82)"]],
];
const TABS = [
    ["overview", "Overview"],
    ["mission", "Mission"],
    ["geometry", "Geometry"],
    ["weights", "Weights & Structures"],
    ["aero", "Aerodynamics"],
    ["engine", "Engine"],
];
const ENGINE_MAPS = [
    ["fan", "Fan"],
    ["lpc", "LPC"],
    ["hpc", "HPC"],
];
const state = {
    data: null,
    activeTab: "overview",
    pointIndex: 9,
    missionAxisMode: "range",
    geometryCompareMode: "full",
};
const root = document.getElementById("app");
document.addEventListener("DOMContentLoaded", bootstrap);
document.addEventListener("mathjax-ready", () => {
    typesetMath();
});
window.addEventListener("load", () => {
    typesetMath();
});
function bootstrap() {
    const embedded = readEmbeddedData();
    if (embedded) {
        setData(embedded);
        return;
    }
    const preloadUrl = readPreloadUrl();
    if (preloadUrl) {
        renderLoader(`Loading ${preloadUrl} ...`);
        loadDashboardUrl(preloadUrl);
        return;
    }
    renderLoader();
}
function readPreloadUrl() {
    try {
        return new URLSearchParams(window.location.search).get("input");
    }
    catch (_error) {
        return null;
    }
}
function readEmbeddedData() {
    const scriptNode = document.getElementById("taesopt-dashboard-data");
    if (scriptNode && scriptNode.textContent && scriptNode.textContent.trim()) {
        return JSON.parse(scriptNode.textContent);
    }
    if (window.__TAESOPT_DASHBOARD_DATA__) {
        return window.__TAESOPT_DASHBOARD_DATA__;
    }
    return null;
}
function setData(data) {
    state.data = trimHiddenMissionPoint(data);
    const defaultPoint = (data.meta && data.meta.default_point) || 1;
    const labels = (data.meta && data.meta.mission_short) || [];
    state.pointIndex = clamp(defaultPoint - 1, 0, Math.max(labels.length - 1, 0));
    renderDashboard();
}
function renderLoader(message) {
    root.innerHTML = `
    <div class="empty-state">
      <section class="empty-card">
        <h1>TASOPT Diagnostic Dashboard</h1>
        <p>Load a dashboard JSON export generated from one or two sized TASOPT aircraft. The browser only renders the snapshot; Julia stays responsible for reading .jld2 and exporting the diagnostic payload.</p>
        <div class="dropzone" id="dropzone">
          <strong>Drop a dashboard JSON file here</strong>
          <small>or use the picker below</small>
          <div style="margin-top:16px;">
            <button class="load-button" id="loader-button" type="button">Load Dashboard JSON</button>
            <input id="dashboard-file-input" type="file" accept=".json,application/json" hidden />
          </div>
        </div>
        ${message ? `<p style="margin-top:18px;color:#f0907a;">${escapeHtml(message)}</p>` : ""}
      </section>
    </div>
  `;
    wireFileInput();
    wireDropzone();
}
function renderDashboard() {
    const compare = isCompareMode();
    root.innerHTML = `
    <div class="dashboard-shell">
      <section class="toolbar top-nav">
        <nav class="tabs">
          ${TABS.map(([key, label]) => `<button class="tab-button ${state.activeTab === key ? "active" : ""}" data-tab="${key}" type="button">${label}</button>`).join("")}
        </nav>
        <div class="toolbar-side">
          <div class="toolbar-controls">
            <div class="segmented-control compact">
              <span class="control-label">Mission axis</span>
              <div class="segmented-buttons">
                <button class="segmented-button ${state.missionAxisMode === "range" ? "active" : ""}" data-axis-mode="range" type="button">Actual</button>
                <button class="segmented-button ${state.missionAxisMode === "normalized" ? "active" : ""}" data-axis-mode="normalized" type="button">Norm.</button>
              </div>
            </div>
            ${compare ? `
              <div class="segmented-control compact">
                <span class="control-label">Geometry</span>
                <div class="segmented-buttons">
                  <button class="segmented-button ${state.geometryCompareMode === "full" ? "active" : ""}" data-geometry-mode="full" type="button">Full</button>
                  <button class="segmented-button ${state.geometryCompareMode === "half" ? "active" : ""}" data-geometry-mode="half" type="button">Half</button>
                </div>
              </div>
            ` : ""}
          </div>
          <div class="point-control">
            <span>Mission point</span>
            <select id="point-select" class="point-select">
              ${state.data.meta.mission_short.map((shortLabel, index) => `
                <option value="${index}" ${index === state.pointIndex ? "selected" : ""}>
                  ${escapeHtml(`${shortLabel} / ${state.data.meta.mission_labels[index]}`)}
                </option>
              `).join("")}
            </select>
          </div>
        </div>
      </section>

      <header class="masthead">
        <div>
          <h1>${escapeHtml(state.data.meta.title || "TASOPT Diagnostic Dashboard")}</h1>
          <p>High-resolution SVG aircraft geometry with chart-focused mission, aero, weight, and engine diagnostics for a quick read on one design or a two-aircraft comparison.</p>
        </div>
        <div class="masthead-actions">
          <span class="pill"><strong>${compare ? "COMPARE" : "SINGLE"}</strong> ${state.data.aircraft.map((item) => escapeHtml(item.name)).join(" vs ")}</span>
          <button class="load-button" id="load-button" type="button">Load JSON</button>
          <input id="dashboard-file-input" type="file" accept=".json,application/json" hidden />
        </div>
      </header>

      ${renderHeroSection()}

      <main id="tab-content"></main>
    </div>
    <div class="tooltip" id="hover-tooltip"></div>
  `;
    wireFileInput();
    wireTabs();
    wireControls();
    wirePointSelect();
    renderTab();
    installSvgTooltips();
    typesetMath();
}
function renderHeroSection() {
    if (isCompareMode() && state.geometryCompareMode === "half") {
        return `
      <section class="hero-grid single">
        ${renderHalfCompareCard(state.data.aircraft[0], state.data.aircraft[1])}
      </section>
    `;
    }
    const compare = isCompareMode();
    return `
    <section class="hero-grid ${compare ? "compare" : "single"}">
      ${state.data.aircraft.map((aircraft, index) => renderAircraftCard(aircraft, index)).join("")}
    </section>
  `;
}
function renderAircraftCard(aircraft, index) {
    const accentClass = index === 1 ? "accent-b" : "";
    return `
    <article class="aircraft-card ${accentClass}">
      <div class="aircraft-head">
        <div>
          <h2 class="aircraft-title">${escapeHtml(aircraft.name)}</h2>
          <div class="aircraft-source">${aircraft.source ? `Source: ${escapeHtml(aircraft.source)}` : "Direct aircraft export"}</div>
        </div>
        <div class="pill"><strong>${index === 0 ? "A" : "B"}</strong> ${formatMetricValue(aircraft.summary.range_km, 0)} km range</div>
      </div>
      <div class="figure-wrap">
        ${renderAircraftSvg(aircraft, index)}
      </div>
      <div class="chip-row">
        ${renderMetricChips(aircraft)}
      </div>
    </article>
  `;
}
function renderHalfCompareCard(leftAircraft, rightAircraft) {
    return `
    <article class="aircraft-card compare-half-card">
      <div class="compare-half-head">
        <div class="compare-half-title">
          <span class="pill"><strong>A</strong> ${escapeHtml(leftAircraft.name)}</span>
          <div class="aircraft-source">${leftAircraft.source ? `Source: ${escapeHtml(leftAircraft.source)}` : "Direct aircraft export"}</div>
        </div>
        <div class="pill"><strong>HALF</strong> Shared centerline size comparison</div>
        <div class="compare-half-title right">
          <span class="pill"><strong>B</strong> ${escapeHtml(rightAircraft.name)}</span>
          <div class="aircraft-source">${rightAircraft.source ? `Source: ${escapeHtml(rightAircraft.source)}` : "Direct aircraft export"}</div>
        </div>
      </div>
      <div class="compare-half-body">
        <aside class="metric-rail left">
          ${renderMetricChips(leftAircraft)}
        </aside>
        <div class="figure-wrap compare-half-figure">
          ${renderHalfCompareSvg(leftAircraft, rightAircraft)}
        </div>
        <aside class="metric-rail right">
          ${renderMetricChips(rightAircraft)}
        </aside>
      </div>
    </article>
  `;
}
function renderMetricChips(aircraft) {
    const chips = [
        ["Span", `${formatMetricValue(aircraft.summary.wing_span_m, 1)} m`],
        ["MTOW", `${formatMetricValue(aircraft.summary.mtow_t, 1)} t`],
        ["Cruise L/D", formatMetricValue(aircraft.summary.cruise_lod, 2)],
        ["TSFC", `${formatMetricValue(aircraft.summary.cruise_tsfc_mg_ns, 3)} mg/N/s`],
    ];
    return chips
        .map(([label, value]) => `<div class="metric-chip"><span>${escapeHtml(label)}</span><strong>${escapeHtml(value)}</strong></div>`)
        .join("");
}
function renderAircraftSvg(aircraft, index) {
    const viewBox = aircraft.geometry.svg.view_box;
    return `
    <div class="figure-frame ${index === 1 ? "accent-b" : ""}">
      <svg viewBox="${viewBox.x} ${viewBox.y} ${viewBox.width} ${viewBox.height}" preserveAspectRatio="xMidYMid meet" aria-label="${escapeAttr(aircraft.name)} stick figure">
        <g transform="scale(1,-1)">
          <line class="figure-centerline" x1="${viewBox.x + 0.04 * viewBox.width}" y1="0" x2="${viewBox.x + 0.96 * viewBox.width}" y2="0"></line>
          ${renderSvgComponents(aircraft, index)}
        </g>
      </svg>
    </div>
  `;
}
function renderHalfCompareSvg(leftAircraft, rightAircraft) {
    const viewBox = rotatedSvgViewBox([leftAircraft, rightAircraft]);
    const leftClipWidth = Math.max(0 - viewBox.x, 0.001);
    const rightClipWidth = Math.max(viewBox.x + viewBox.width, 0.001);
    return `
    <div class="figure-frame half-compare">
      <svg viewBox="${viewBox.x} ${viewBox.y} ${viewBox.width} ${viewBox.height}" preserveAspectRatio="xMidYMid meet" aria-label="Half-view aircraft comparison">
        <defs>
          <clipPath id="compare-half-left-clip">
            <rect x="${viewBox.x}" y="${viewBox.y}" width="${leftClipWidth}" height="${viewBox.height}"></rect>
          </clipPath>
          <clipPath id="compare-half-right-clip">
            <rect x="0" y="${viewBox.y}" width="${rightClipWidth}" height="${viewBox.height}"></rect>
          </clipPath>
        </defs>
        <line class="figure-divider" x1="0" y1="${viewBox.y + 0.03 * viewBox.height}" x2="0" y2="${viewBox.y + 0.97 * viewBox.height}"></line>
        <g clip-path="url(#compare-half-left-clip)">
          ${renderSvgComponents(leftAircraft, 0, transformRotateToVertical)}
        </g>
        <g class="accent-b" clip-path="url(#compare-half-right-clip)">
          ${renderSvgComponents(rightAircraft, 1, transformRotateToVertical)}
        </g>
      </svg>
    </div>
  `;
}
function renderSvgComponents(aircraft, index, transformPoint = null) {
    return aircraft.geometry.svg.components
        .map((component) => {
        const points = component.points
            .map((point) => {
            const finalPoint = transformPoint ? transformPoint(point) : point;
            return `${finalPoint[0]},${finalPoint[1]}`;
        })
            .join(" ");
        const tip = escapeAttr(component.tooltip);
        return `<polygon class="aircraft-component kind-${component.kind}" points="${points}" data-tip="${tip}" data-accent="${index === 1 ? "b" : "a"}"></polygon>`;
    })
        .join("");
}
function rotatedSvgViewBox(aircraftList) {
    const allPoints = aircraftList.flatMap((aircraft) => aircraft.geometry.svg.components.flatMap((component) => component.points.map(transformRotateToVertical)));
    const xs = allPoints.map((point) => Number(point[0]));
    const ys = allPoints.map((point) => Number(point[1]));
    const minX = Math.min(...xs);
    const maxX = Math.max(...xs);
    const minY = Math.min(...ys);
    const maxY = Math.max(...ys);
    const xPad = Math.max(2.0, 0.10 * (maxX - minX));
    const yPad = Math.max(2.0, 0.06 * (maxY - minY));
    return {
        x: minX - xPad,
        y: minY - yPad,
        width: maxX - minX + 2 * xPad,
        height: maxY - minY + 2 * yPad,
    };
}
function transformRotateToVertical(point) {
    return [Number(point[1]), Number(point[0])];
}
function wireFileInput() {
    const fileInput = document.getElementById("dashboard-file-input");
    const loaderButton = document.getElementById("loader-button");
    const loadButton = document.getElementById("load-button");
    if (loaderButton) {
        loaderButton.addEventListener("click", () => fileInput && fileInput.click());
    }
    if (loadButton) {
        loadButton.addEventListener("click", () => fileInput && fileInput.click());
    }
    if (fileInput) {
        fileInput.addEventListener("change", (event) => {
            const target = event.target;
            const file = target.files && target.files[0];
            if (file) {
                readDashboardFile(file);
            }
        });
    }
}
function wireDropzone() {
    const dropzone = document.getElementById("dropzone");
    if (!dropzone) {
        return;
    }
    dropzone.addEventListener("dragover", (event) => {
        event.preventDefault();
        dropzone.classList.add("dragover");
    });
    dropzone.addEventListener("dragleave", () => {
        dropzone.classList.remove("dragover");
    });
    dropzone.addEventListener("drop", (event) => {
        event.preventDefault();
        dropzone.classList.remove("dragover");
        const file = event.dataTransfer && event.dataTransfer.files && event.dataTransfer.files[0];
        if (file) {
            readDashboardFile(file);
        }
    });
}
function readDashboardFile(file) {
    const reader = new FileReader();
    reader.onload = () => {
        try {
            const payload = JSON.parse(String(reader.result || "{}"));
            setData(payload);
        }
        catch (_error) {
            renderLoader(`Could not parse ${file.name}. Make sure this file came from export_dashboard_data().`);
        }
    };
    reader.readAsText(file);
}
function loadDashboardUrl(url) {
    fetch(url)
        .then((response) => {
        if (!response.ok) {
            throw new Error(`HTTP ${response.status}`);
        }
        return response.json();
    })
        .then((payload) => {
        setData(payload);
    })
        .catch((_error) => {
        renderLoader(`Could not load ${url}. Check that the file exists and is valid dashboard JSON.`);
    });
}
function wireTabs() {
    document.querySelectorAll("[data-tab]").forEach((node) => {
        node.addEventListener("click", () => {
            state.activeTab = node.getAttribute("data-tab");
            renderDashboard();
        });
    });
}
function wireControls() {
    document.querySelectorAll("[data-axis-mode]").forEach((node) => {
        node.addEventListener("click", () => {
            state.missionAxisMode = node.getAttribute("data-axis-mode");
            renderDashboard();
        });
    });
    document.querySelectorAll("[data-geometry-mode]").forEach((node) => {
        node.addEventListener("click", () => {
            state.geometryCompareMode = node.getAttribute("data-geometry-mode");
            renderDashboard();
        });
    });
}
function wirePointSelect() {
    const select = document.getElementById("point-select");
    if (!select) {
        return;
    }
    select.addEventListener("change", (event) => {
        state.pointIndex = Number(event.target.value);
        renderTab();
    });
}
function installSvgTooltips() {
    const tooltip = document.getElementById("hover-tooltip");
    if (!tooltip) {
        return;
    }
    document.querySelectorAll(".aircraft-component").forEach((node) => {
        node.addEventListener("mousemove", (event) => {
            tooltip.textContent = node.getAttribute("data-tip") || "";
            tooltip.style.opacity = "1";
            tooltip.style.left = `${event.clientX + 16}px`;
            tooltip.style.top = `${event.clientY + 16}px`;
        });
        node.addEventListener("mouseleave", () => {
            tooltip.style.opacity = "0";
        });
    });
}
function renderTab() {
    const container = document.getElementById("tab-content");
    if (!container || !state.data) {
        return;
    }
    if (state.activeTab === "overview") {
        container.innerHTML = renderOverviewTab();
        plotOverviewTab();
    }
    else if (state.activeTab === "mission") {
        container.innerHTML = renderMissionTab();
        plotMissionTab();
    }
    else if (state.activeTab === "geometry") {
        container.innerHTML = renderGeometryTab();
        plotGeometryTab();
    }
    else if (state.activeTab === "weights") {
        container.innerHTML = renderWeightsTab();
        plotWeightsTab();
    }
    else if (state.activeTab === "aero") {
        container.innerHTML = renderAeroTab();
        plotAeroTab();
    }
    else if (state.activeTab === "engine") {
        container.innerHTML = renderEngineTab();
        plotEngineTab();
    }
    typesetMath();
}
function renderOverviewTab() {
    return `
    <section class="tab-grid two-col">
      <article class="panel full-span">
        <div class="panel-head">
          <div class="eyebrow">Snapshot</div>
          <h2>Primary diagnostic metrics</h2>
          <div class="panel-note">Fast scalar comparison for the figures that usually drive the first design review conversation.</div>
        </div>
        <div class="panel-body summary-cards">${renderSummaryCards()}</div>
      </article>
      ${renderTablePanel("KPI table", "Full overview scalar readout for the loaded design set.", overviewSnapshotRows(), isCompareMode() ? "Comparison" : "Snapshot")}
      ${chartPanel("overview-efficiency", "Range efficiency over mission", "Folded BRE view: L/D divided by TSFC across the mission points.")}
      ${chartPanel("overview-drivers", "Cruise efficiency drivers", "Compact view of cruise L/D, 1 / TSFC, and inverse PFEI.")}
    </section>
  `;
}
function renderMissionTab() {
    const axisLabel = missionAxisLabel().toLowerCase();
    return `
    <section class="tab-grid three-col">
      ${renderTablePanel("Mission snapshot", `Selected-point mission state plus the current ${axisLabel} reference.`, missionSnapshotRows())}
      ${chartPanel("mission-alt", "Altitude and \\(\\gamma\\)", `Altitude against ${axisLabel}, with \\(\\gamma\\) overlaid on the secondary axis.`)}
      ${chartPanel("mission-mach", "Mach profile", `Mission Mach number against ${axisLabel}.`)}
      ${chartPanel("mission-lod", "\\(L/D\\)", `Computed \\(L/D\\) against ${axisLabel}.`)}
      ${chartPanel("mission-roc", "Rate of climb", `Useful for spotting climb and descent behavior shifts over ${axisLabel}.`)}
      ${chartPanel("mission-weight", "\\(W/W_{MTO}\\)", `Mission weight fraction positioned by ${axisLabel}.`)}
      ${chartPanel("mission-cd", "\\(C_D \\times 10^4\\)", "Total \\(C_D\\) scaled by \\(10^4\\) for readability.")}
    </section>
  `;
}
function renderGeometryTab() {
    const pointText = missionPointText(state.pointIndex);
    return `
    <section class="tab-grid three-col">
      ${renderTablePanel("Geometry snapshot", "Wing, fuselage, and propulsor dimensions for quick size comparison.", geometrySnapshotRows())}
      ${chartPanel("geometry-chord", "Wing chord distribution", "Planform chord values exported directly from the sized wing geometry.")}
      ${chartPanel("geometry-cl", "Spanwise section \\(C_l\\)", `Selected mission point: ${pointText}.`)}
      ${chartPanel("geometry-fuse", "Fuselage cross section", "Exported fuselage section geometry from the TASOPT layout, including multi-bubble shaping instead of a circular proxy.")}
      ${chartPanel("geometry-bars", "Geometry scalar comparison", "Span, area, AR, sweep, and fuselage dimensions in one scan.")}
    </section>
  `;
}
function renderWeightsTab() {
    const piePanels = state.data.aircraft
        .map((aircraft, index) => chartPanel(`weights-pie-${index}`, `${aircraft.name} MTOW breakdown`, "OEW, fuel, and payload as shares of MTOW."))
        .join("");
    return `
    <section class="tab-grid three-col">
      ${renderTablePanel("Weights snapshot", "Primary mass totals plus subsystem structure weights.", weightsSnapshotRows())}
      ${piePanels}
      ${chartPanel("weights-bar", "Primary weight totals", "MTOW, OEW, fuel, and payload side by side.")}
      ${chartPanel("weights-structure", "Structural breakdown", "Subsystem-level structure and installed system weights.", true)}
    </section>
  `;
}
function renderAeroTab() {
    const pointText = missionPointText(state.pointIndex);
    const hasTrefftz = state.data.aircraft.some((aircraft) => aircraft.aero.trefftz);
    const fractionPanels = state.data.aircraft
        .map((aircraft, index) => chartPanel(`aero-fractions-${index}`, `${aircraft.name} drag fractions`, "Component drag fractions stacked over the mission."))
        .join("");
    return `
    <section class="tab-grid three-col">
      ${renderTablePanel("Aerodynamic snapshot", `Cruise drag build-up plus selected-point aerodynamic state at ${pointText}.`, aeroSnapshotRows())}
      ${chartPanel("aero-drag", "Cruise drag areas", "\\(C_D S\\) at Cruise 1, shown as equivalent cm\\(^2\\).")}
      ${chartPanel("aero-lod", "Mission \\(L/D\\)", `Aerodynamic efficiency through the mission, positioned by ${missionAxisLabel().toLowerCase()}.`)}
      ${fractionPanels}
      ${hasTrefftz ? chartPanel("aero-circulation", "Trefftz circulation", `Selected mission point: ${pointText}.`) : ""}
      ${hasTrefftz ? chartPanel("aero-induced", "Trefftz induced drag density", `Selected mission point: ${pointText}; plotted as the exported \\(-\\Gamma v_n\\) distribution.`) : ""}
    </section>
  `;
}
function renderEngineTab() {
    const pointText = missionPointText(state.pointIndex);
    const mapPanels = ENGINE_MAPS
        .flatMap(([key, label]) => state.data.aircraft.map((aircraft, index) => chartPanel(`engine-map-${key}-${index}`, `${aircraft.name} ${label} map`, `Operating line over the normalized ${label} map. Selected mission point: ${pointText}.`, true)))
        .join("");
    return `
    <section class="tab-grid two-col">
      ${renderTablePanel("Engine snapshot", `Engine cycle metrics at ${pointText} plus mission-driving design values.`, engineSnapshotRows())}
      <div class="nested-grid two-col full-span">
        ${chartPanel("engine-temp", "Station temperatures", `Selected mission point: ${pointText}.`)}
        ${chartPanel("engine-pressure", "Station pressures", `Selected mission point: ${pointText}; pressures in kPa.`)}
        ${chartPanel("engine-tsfc", "\\(\\mathrm{TSFC}\\) and \\(T_{t4}\\)", `Mission trends for \\(\\mathrm{TSFC}\\) with \\(T_{t4}\\) overlaid on a secondary axis, positioned by ${missionAxisLabel().toLowerCase()}.`)}
        ${chartPanel("engine-ratios", "\\(\\mathrm{BPR}\\), \\(\\mathrm{OPR}\\), and \\(\\mathrm{FPR}\\)", "Mission trends for the main engine cycle ratios.")}
      </div>
      <div class="nested-grid two-col full-span">
        ${mapPanels}
      </div>
    </section>
  `;
}
function renderSummaryCards() {
    if (state.data.comparison) {
        return state.data.comparison.summary.slice(0, 8).map((row) => {
            const delta = deltaVisual(row);
            return `
        <div class="summary-card">
          <div class="label">${renderMathLabel(row.label)}</div>
          <div class="value">${formatMetricValue(row.right, row.digits)}</div>
          <div class="delta ${delta.className}">${delta.text}</div>
        </div>
      `;
        }).join("");
    }
    return state.data.aircraft[0].summary.kpis.slice(0, 8).map((row) => `
    <div class="summary-card">
      <div class="label">${renderMathLabel(row.label)}</div>
      <div class="value">${formatMetricValue(row.value, row.digits)}</div>
    </div>
  `).join("");
}
function renderTablePanel(title, note, rows, eyebrow = "Snapshot", fullSpan = false) {
    return `
    <article class="panel table-panel ${fullSpan ? "full-span" : ""}">
      <div class="panel-head">
        <div class="eyebrow">${escapeHtml(eyebrow)}</div>
        <h3>${renderRichText(title)}</h3>
        <div class="panel-note">${renderRichText(note)}</div>
      </div>
      <div class="panel-body">${renderMetricTable(rows)}</div>
    </article>
  `;
}
function renderMetricTable(rows) {
    if (!isCompareMode()) {
        const aircraft = state.data.aircraft[0];
        return `
      <table class="metric-table">
        <thead>
          <tr><th>Metric</th><th>${escapeHtml(aircraft.name)}</th></tr>
        </thead>
        <tbody>
          ${rows.map((row) => `
            <tr>
              <td>${renderMathLabel(row.label)}</td>
              <td class="number">${formatMetricValue(row.value, row.digits)}</td>
            </tr>
          `).join("")}
        </tbody>
      </table>
    `;
    }
    return `
    <table class="metric-table">
      <thead>
        <tr>
          <th>Metric</th>
          <th>${escapeHtml(state.data.aircraft[0].name)}</th>
          <th>${escapeHtml(state.data.aircraft[1].name)}</th>
          <th>Delta</th>
        </tr>
      </thead>
      <tbody>
        ${rows.map((row) => {
        const delta = deltaVisual(row);
        return `
            <tr>
              <td>${renderMathLabel(row.label)}</td>
              <td class="number">${formatMetricValue(row.left, row.digits)}</td>
              <td class="number">${formatMetricValue(row.right, row.digits)}</td>
              <td class="number delta-text ${delta.className}">${escapeHtml(delta.text)}</td>
            </tr>
          `;
    }).join("")}
      </tbody>
    </table>
  `;
}
function chartPanel(id, title, note, tall = false) {
    return `
    <article class="panel">
      <div class="panel-head">
        <h3>${renderRichText(title)}</h3>
        <div class="panel-note">${renderRichText(note)}</div>
      </div>
      <div class="panel-body">
        <div id="${id}" class="chart-frame ${tall ? "tall" : ""}"></div>
      </div>
    </article>
  `;
}
function plotOverviewTab() {
    missionPlot("overview-efficiency", "(L/D) / TSFC", (aircraft) => aircraft.mission.lod.map((lod, point) => {
        const tsfc = aircraft.engine.tsfc_mg_ns[point];
        return tsfc > 0 ? lod / tsfc : 0;
    }), "(L/D) / TSFC %{y:.3f}", {}, (_aircraft, index) => ({
        fill: "tozeroy",
        fillcolor: FILL_COLORS[index],
    }));
    const driverLabels = ["Cruise L/D", "1000 / TSFC", "10 / PFEI"];
    const driverTraces = state.data.aircraft.map((aircraft, index) => barTrace(driverLabels, [
        aircraft.summary.cruise_lod,
        aircraft.summary.cruise_tsfc_mg_ns > 0 ? 1000 / aircraft.summary.cruise_tsfc_mg_ns : 0,
        aircraft.summary.pfei_kj_kg_km > 0 ? 10 / aircraft.summary.pfei_kj_kg_km : 0,
    ], aircraft.name, index));
    plotChart("overview-drivers", driverTraces, layout({
        barmode: "group",
        yaxis: { title: "Scaled value" },
        margin: { l: 54, r: 20, t: 24, b: 60 },
    }));
}
function plotMissionTab() {
    const altTraces = state.data.aircraft.flatMap((aircraft, index) => ([
        missionTrace(aircraft, aircraft.mission.alt_km, `${aircraft.name} altitude`, index, {
            fill: "tozeroy",
            fillcolor: FILL_COLORS[index],
        }, "Altitude %{y:.2f} km"),
        missionTrace(aircraft, aircraft.mission.gamma_rad, `${aircraft.name} gamma`, index, {
            yaxis: "y2",
            mode: "lines",
            line: { color: COLORS[index], width: 1.6, dash: "dot" },
        }, "Gamma %{y:.4f} rad"),
    ]));
    plotChart("mission-alt", altTraces, layout({
        xaxis: missionAxisConfig(),
        yaxis: { title: "Altitude [km]" },
        yaxis2: { title: "$\\gamma$ [rad]", overlaying: "y", side: "right", showgrid: false },
        shapes: missionSelectionShapes(),
        hovermode: "closest",
    }));
    missionPlot("mission-mach", "Mach", (aircraft) => aircraft.mission.mach, "Mach %{y:.3f}");
    missionPlot("mission-lod", "$L/D$", (aircraft) => aircraft.mission.lod, "L/D %{y:.2f}");
    missionPlot("mission-roc", "ROC [m/s]", (aircraft) => aircraft.mission.roc_mps, "ROC %{y:.2f} m/s");
    missionPlot("mission-weight", "$W/W_{MTO}$", (aircraft) => aircraft.mission.weight_fraction, "W / WMTO %{y:.3f}");
    missionPlot("mission-cd", "$C_D \\times 10^4$", (aircraft) => aircraft.mission.cd_e4, "CD x 1e4 %{y:.1f}");
}
function plotGeometryTab() {
    const chordTraces = state.data.aircraft.map((aircraft, index) => lineTrace(aircraft.aero.wing.eta, aircraft.aero.wing.chord_m, aircraft.name, index, { mode: "lines+markers", marker: { size: 7, color: COLORS[index] } }));
    plotChart("geometry-chord", chordTraces, layout({
        xaxis: { title: "$\\eta = 2y/b$" },
        yaxis: { title: "Chord [m]" },
    }));
    const clTraces = state.data.aircraft.map((aircraft, index) => lineTrace(aircraft.aero.wing.cl_eta, aircraft.aero.wing.cl_by_point[state.pointIndex], aircraft.name, index, {
        mode: "lines+markers",
        marker: { size: 7, color: COLORS[index] },
        fill: "tozeroy",
        fillcolor: FILL_COLORS[index],
    }));
    plotChart("geometry-cl", clTraces, layout({
        xaxis: { title: "$\\eta = 2y/b$" },
        yaxis: { title: "Section $C_l$" },
    }));
    const fuseTraces = state.data.aircraft.map((aircraft, index) => ({
        x: aircraft.geometry.fuselage.cross_section.points.map((point) => point[0]),
        y: aircraft.geometry.fuselage.cross_section.points.map((point) => point[1]),
        name: aircraft.name,
        type: "scatter",
        mode: "lines",
        line: { color: COLORS[index], width: 2.2 },
        fill: "toself",
        fillcolor: FILL_COLORS[index],
        hovertemplate: "y %{x:.2f} m<br>z %{y:.2f} m<extra>" + aircraft.name + "</extra>",
    }));
    plotChart("geometry-fuse", fuseTraces, layout({
        xaxis: { title: "$y$ [m]", scaleanchor: "y" },
        yaxis: { title: "$z$ [m]" },
    }));
    const categories = ["Span [m]", "Area [m^2]", "AR", "Sweep [deg]", "Fuse L [m]", "Fuse W [m]"];
    const barTraces = state.data.aircraft.map((aircraft, index) => barTrace(categories, [
        aircraft.summary.wing_span_m,
        aircraft.summary.wing_area_m2,
        aircraft.summary.wing_ar,
        aircraft.summary.wing_sweep_deg,
        aircraft.geometry.fuselage.length_m,
        aircraft.geometry.fuselage.width_m,
    ], aircraft.name, index));
    plotChart("geometry-bars", barTraces, layout({
        barmode: "group",
        yaxis: { title: "Value" },
        margin: { l: 54, r: 20, t: 24, b: 72 },
    }));
}
function plotWeightsTab() {
    state.data.aircraft.forEach((aircraft, index) => {
        plotChart(`weights-pie-${index}`, [{
                labels: aircraft.weights.mtow_breakdown.labels,
                values: aircraft.weights.mtow_breakdown.values_t,
                type: "pie",
                hole: 0.48,
                marker: { colors: [COLORS[index], "#d7b46a", "#7ed09a"] },
                textinfo: "label+percent",
                hovertemplate: "%{label}: %{value:.2f} t<extra></extra>",
            }], {
            paper_bgcolor: "rgba(0,0,0,0)",
            plot_bgcolor: "rgba(0,0,0,0)",
            font: { family: "IBM Plex Sans, Avenir Next, Segoe UI, sans-serif", color: "#8fa5ae" },
            margin: { l: 4, r: 4, t: 8, b: 8 },
            showlegend: true,
            legend: { orientation: "h", y: -0.08, font: { size: 10 } },
            annotations: [{
                    showarrow: false,
                    x: 0.5,
                    y: 0.5,
                    text: `<b>${formatMetricValue(aircraft.weights.mtow_t, 1)}</b><br>MTOW`,
                    font: { color: "#e9efe8", size: 13 },
                }],
        });
    });
    const categories = ["MTOW", "OEW", "Fuel", "Payload"];
    const barTraces = state.data.aircraft.map((aircraft, index) => barTrace(categories, [aircraft.weights.mtow_t, aircraft.weights.oew_t, aircraft.weights.fuel_t, aircraft.weights.payload_t], aircraft.name, index));
    plotChart("weights-bar", barTraces, layout({
        barmode: "group",
        yaxis: { title: "Weight [t]" },
    }));
    const labels = state.data.aircraft[0].weights.structure_breakdown.labels;
    const structureTraces = state.data.aircraft.map((aircraft, index) => barTrace(labels, aircraft.weights.structure_breakdown.values_t, aircraft.name, index));
    plotChart("weights-structure", structureTraces, layout({
        barmode: "group",
        yaxis: { title: "Weight [t]" },
        margin: { l: 54, r: 20, t: 24, b: 78 },
    }));
}
function plotAeroTab() {
    const dragTraces = state.data.aircraft.map((aircraft, index) => barTrace(aircraft.aero.drag_keys, aircraft.aero.drag_area_cruise_cm2, aircraft.name, index));
    plotChart("aero-drag", dragTraces, layout({
        barmode: "group",
        yaxis: { title: "$C_D S$ [cm^2]" },
        margin: { l: 54, r: 20, t: 24, b: 62 },
    }));
    missionPlot("aero-lod", "$L/D$", (aircraft) => aircraft.aero.lod, "L/D %{y:.2f}");
    state.data.aircraft.forEach((aircraft, index) => {
        const traces = aircraft.aero.drag_keys.map((key, keyIndex) => ({
            x: missionPositions(aircraft),
            y: aircraft.aero.drag_fractions[keyIndex],
            customdata: missionCustomdata(aircraft),
            name: key,
            type: "scatter",
            mode: "lines",
            stackgroup: "one",
            fill: keyIndex === 0 ? "tozeroy" : "tonexty",
            line: {
                color: keyIndex % 2 === 0 ? COLORS[index] : "#d7b46a",
                width: 1.2,
                dash: keyIndex % 3 === 0 ? "solid" : "dot",
            },
            hovertemplate: "%{customdata[0]}<br>Range %{customdata[1]:.0f} km<br>R / Rdes %{customdata[2]:.3f}<br>Fraction %{y:.3f}<extra>" + key + "</extra>",
        }));
        plotChart(`aero-fractions-${index}`, traces, layout({
            xaxis: missionAxisConfig(),
            yaxis: { title: "Fraction of $C_D$" },
            legend: { orientation: "h", y: -0.22, font: { size: 10 } },
            shapes: missionSelectionShapes(index),
            hovermode: "closest",
            margin: { l: 54, r: 20, t: 24, b: 78 },
        }));
    });
    if (state.data.aircraft.some((aircraft) => aircraft.aero.trefftz)) {
        const circulationTraces = state.data.aircraft
            .map((aircraft, index) => {
            if (!aircraft.aero.trefftz) {
                return null;
            }
            return mirroredWingTrace(aircraft.aero.trefftz.wing_span_y_m, aircraft.aero.trefftz.circulation_by_point[state.pointIndex], aircraft.name, index, {
                mode: "lines+markers",
                marker: { size: 5, color: COLORS[index] },
            }, "Gamma %{y:.3f} m^2/s");
        })
            .filter(Boolean);
        plotChart("aero-circulation", circulationTraces, layout({
            xaxis: { title: "Spanwise location y [m]" },
            yaxis: { title: "Circulation $\\Gamma$ [m^2/s]" },
            hovermode: "closest",
        }));
        const inducedTraces = state.data.aircraft
            .map((aircraft, index) => {
            if (!aircraft.aero.trefftz) {
                return null;
            }
            return mirroredWingTrace(aircraft.aero.trefftz.wing_span_y_m, aircraft.aero.trefftz.drag_density_proxy[state.pointIndex], aircraft.name, index, {
                mode: "lines+markers",
                marker: { size: 5, color: COLORS[index], symbol: "square" },
            }, "-Gamma vn %{y:.3f}");
        })
            .filter(Boolean);
        plotChart("aero-induced", inducedTraces, layout({
            xaxis: { title: "Spanwise location y [m]" },
            yaxis: { title: "Induced drag density, $-\\Gamma v_n$" },
            hovermode: "closest",
        }));
    }
}
function plotEngineTab() {
    const pointIndex = state.pointIndex;
    const tempTraces = state.data.aircraft.map((aircraft, index) => ({
        x: aircraft.engine.temperature_stations,
        y: aircraft.engine.temperature_k[pointIndex],
        name: aircraft.name,
        type: "scatter",
        mode: "lines+markers",
        line: { color: COLORS[index], width: 2.4 },
        marker: { size: 7, color: COLORS[index] },
        hovertemplate: "%{y:.1f} K<extra>" + aircraft.name + "</extra>",
    }));
    plotChart("engine-temp", tempTraces, layout({
        yaxis: { title: "Temperature [K]" },
    }));
    const pressureTraces = state.data.aircraft.map((aircraft, index) => ({
        x: aircraft.engine.pressure_stations,
        y: aircraft.engine.pressure_kpa[pointIndex],
        name: aircraft.name,
        type: "scatter",
        mode: "lines+markers",
        line: { color: COLORS[index], width: 2.4 },
        marker: { size: 7, color: COLORS[index] },
        hovertemplate: "%{y:.1f} kPa<extra>" + aircraft.name + "</extra>",
    }));
    plotChart("engine-pressure", pressureTraces, layout({
        yaxis: { title: "Pressure [kPa]" },
    }));
    const tsfcTraces = state.data.aircraft.flatMap((aircraft, index) => ([
        missionTrace(aircraft, aircraft.engine.tsfc_mg_ns, `${aircraft.name} TSFC`, index, {}, "TSFC %{y:.3f} mg/N/s"),
        missionTrace(aircraft, aircraft.engine.tt4_k, `${aircraft.name} Tt4`, index, {
            yaxis: "y2",
            mode: "lines",
            line: { color: COLORS[index], width: 1.6, dash: "dot" },
        }, "Tt4 %{y:.1f} K"),
    ]));
    plotChart("engine-tsfc", tsfcTraces, layout({
        xaxis: missionAxisConfig(),
        yaxis: { title: "$\\mathrm{TSFC}$ [mg/N/s]" },
        yaxis2: { title: "$T_{t4}$ [K]", overlaying: "y", side: "right", showgrid: false },
        shapes: missionSelectionShapes(),
        hovermode: "closest",
    }));
    const ratioTraces = state.data.aircraft.flatMap((aircraft, index) => ([
        missionTrace(aircraft, aircraft.engine.bpr, `${aircraft.name} BPR`, index, {}, "BPR %{y:.3f}"),
        missionTrace(aircraft, aircraft.engine.opr.map((value) => value / 10), `${aircraft.name} OPR / 10`, index, {
            mode: "lines",
            line: { color: COLORS[index], width: 1.5, dash: "dot" },
        }, "OPR / 10 %{y:.2f}"),
        missionTrace(aircraft, aircraft.engine.fpr, `${aircraft.name} FPR`, index, {
            mode: "lines",
            line: { color: "#d7b46a", width: 1.4, dash: index === 0 ? "dash" : "dashdot" },
        }, "FPR %{y:.3f}"),
    ]));
    plotChart("engine-ratios", ratioTraces, layout({
        xaxis: missionAxisConfig(),
        yaxis: { title: "Cycle ratio" },
        legend: { orientation: "h", y: -0.22, font: { size: 10 } },
        shapes: missionSelectionShapes(),
        hovermode: "closest",
        margin: { l: 54, r: 20, t: 24, b: 80 },
    }));
    state.data.aircraft.forEach((aircraft, index) => {
        ENGINE_MAPS.forEach(([key]) => {
            plotCompressorMap(`engine-map-${key}-${index}`, aircraft, index, aircraft.engine.maps && aircraft.engine.maps[key]);
        });
    });
}
function plotCompressorMap(id, aircraft, index, mapPayload) {
    if (!mapPayload) {
        const node = document.getElementById(id);
        if (node) {
            node.innerHTML = `<div class="chart-warning">No ${escapeHtml(id)} data was exported for this aircraft.</div>`;
        }
        return;
    }
    const operatingCustomdata = missionCustomdata(aircraft).map((row, pointIndex) => [
        row[0],
        row[1],
        row[2],
        Number(mapPayload.operating_line.speed_frac[pointIndex] || 0),
        Number(mapPayload.operating_line.eff[pointIndex] || 0),
    ]);
    const traces = [{
            x: mapPayload.efficiency_grid.x,
            y: mapPayload.efficiency_grid.y,
            z: mapPayload.efficiency_grid.z,
            type: "contour",
            colorscale: MAP_FILL_COLORS[index],
            contours: {
                coloring: "fill",
                showlines: true,
                start: 0.6,
                end: 0.96,
                size: 0.02,
            },
            line: { width: 0.7, color: "rgba(8, 17, 24, 0.32)" },
            showscale: true,
            colorbar: {
                title: { text: "$\\eta_p$" },
                thickness: 10,
                len: 0.78,
                y: 0.5,
                tickfont: { size: 10, color: "#8fa5ae" },
            },
            hovertemplate: "mb / mbD %{x:.3f}<br>PR / piD %{y:.3f}<br>Eff %{z:.3f}<extra>" + aircraft.name + " map</extra>",
        }]
        .concat((mapPayload.speed_lines || []).map((line) => ({
        x: line.x,
        y: line.y,
        type: "scatter",
        mode: "lines",
        line: {
            color: Math.abs(Number(line.speed_frac) - 1.0) < 0.03 ? "rgba(233, 239, 232, 0.78)" : "rgba(143, 165, 174, 0.48)",
            width: Math.abs(Number(line.speed_frac) - 1.0) < 0.03 ? 2.0 : 1.0,
        },
        hovertemplate: "mb / mbD %{x:.3f}<br>PR / piD %{y:.3f}<br>N / ND " + formatMetricValue(line.speed_frac, 2) + "<extra>Speed line</extra>",
        showlegend: false,
    })))
        .concat([{
            x: mapPayload.operating_line.x,
            y: mapPayload.operating_line.y,
            customdata: operatingCustomdata,
            type: "scatter",
            mode: "markers",
            marker: {
                size: 7,
                color: COLORS[index],
                line: { color: "#081118", width: 1 },
            },
            name: `${aircraft.name} operating line`,
            hovertemplate: "%{customdata[0]}<br>Range %{customdata[1]:.0f} km<br>R / Rdes %{customdata[2]:.3f}<br>mb / mbD %{x:.3f}<br>PR / piD %{y:.3f}<br>N / ND %{customdata[3]:.3f}<br>Eff %{customdata[4]:.3f}<extra>" + aircraft.name + "</extra>",
            showlegend: false,
        }, {
            x: [mapPayload.operating_line.x[state.pointIndex]],
            y: [mapPayload.operating_line.y[state.pointIndex]],
            customdata: [operatingCustomdata[state.pointIndex]],
            type: "scatter",
            mode: "markers",
            marker: {
                size: 11,
                color: "#d7b46a",
                symbol: "diamond",
                line: { color: "#081118", width: 1.2 },
            },
            hovertemplate: "%{customdata[0]}<br>Range %{customdata[1]:.0f} km<br>R / Rdes %{customdata[2]:.3f}<br>mb / mbD %{x:.3f}<br>PR / piD %{y:.3f}<extra>Selected point</extra>",
            showlegend: false,
        }, {
            x: [1.0],
            y: [1.0],
            type: "scatter",
            mode: "markers",
            marker: {
                size: 10,
                color: "rgba(8, 17, 24, 0)",
                symbol: "star-diamond",
                line: { color: "rgba(233, 239, 232, 0.85)", width: 1.2 },
            },
            hovertemplate: "Design point<extra></extra>",
            showlegend: false,
        }]);
    plotChart(id, traces, layout({
        xaxis: { title: "$\\dot{m}_b/\\dot{m}_{bD}$", range: [0.45, 1.5] },
        yaxis: { title: "$PR/\\pi_D$", range: [0.5, 1.45] },
        shapes: [
            {
                type: "line",
                x0: 1,
                x1: 1,
                y0: 0.5,
                y1: 1.45,
                line: { color: "rgba(233, 239, 232, 0.18)", width: 1, dash: "dot" },
            },
            {
                type: "line",
                x0: 0.45,
                x1: 1.5,
                y0: 1,
                y1: 1,
                line: { color: "rgba(233, 239, 232, 0.18)", width: 1, dash: "dot" },
            },
        ],
        margin: { l: 58, r: 18, t: 18, b: 46 },
    }));
}
function missionAxisLabel() {
    return state.missionAxisMode === "normalized" ? "range normalized by design range" : "actual mission range";
}
function missionPosition(aircraft, pointIndex) {
    const rangeKm = Number(aircraft.mission.range_km && aircraft.mission.range_km[pointIndex]);
    if (!Number.isFinite(rangeKm)) {
        return Number(pointIndex);
    }
    if (state.missionAxisMode === "normalized") {
        const designRange = Number(aircraft.summary.range_km);
        return designRange > 0 ? rangeKm / designRange : 0;
    }
    return rangeKm;
}
function missionPositions(aircraft) {
    return state.data.meta.mission_short.map((_label, pointIndex) => missionPosition(aircraft, pointIndex));
}
function missionPointText(pointIndex) {
    return `${state.data.meta.mission_short[pointIndex]} / ${state.data.meta.mission_labels[pointIndex]}`;
}
function missionCustomdata(aircraft) {
    return state.data.meta.mission_short.map((_label, pointIndex) => [
        missionPointText(pointIndex),
        Number(aircraft.mission.range_km[pointIndex] || 0),
        normalizedMissionPosition(aircraft, pointIndex),
    ]);
}
function normalizedMissionPosition(aircraft, pointIndex) {
    const rangeKm = Number(aircraft.mission.range_km[pointIndex] || 0);
    const designRange = Number(aircraft.summary.range_km || 0);
    return designRange > 0 ? rangeKm / designRange : 0;
}
function missionAxisConfig() {
    return state.missionAxisMode === "normalized" ? {
        title: "$R/R_{des}$",
        tickformat: ".2f",
        automargin: true,
    } : {
        title: "$R$ [km]",
        tickformat: ".0f",
        automargin: true,
    };
}
function missionSelectionShapes(singleAircraftIndex = null) {
    const aircraftList = singleAircraftIndex === null
        ? state.data.aircraft.map((aircraft, index) => ({ aircraft, index }))
        : [{ aircraft: state.data.aircraft[singleAircraftIndex], index: singleAircraftIndex }];
    return aircraftList.map(({ aircraft, index }) => ({
        type: "line",
        x0: missionPosition(aircraft, state.pointIndex),
        x1: missionPosition(aircraft, state.pointIndex),
        yref: "paper",
        y0: 0,
        y1: 1,
        line: {
            color: index === 0 ? "rgba(113, 176, 212, 0.36)" : "rgba(223, 143, 106, 0.36)",
            width: 1.1,
            dash: "dot",
        },
    }));
}
function missionTrace(aircraft, y, name, index, extra, valueLine) {
    return Object.assign({
        x: missionPositions(aircraft),
        y,
        name,
        type: "scatter",
        mode: "lines+markers",
        line: { color: COLORS[index], width: 2.4 },
        marker: { size: 6, color: COLORS[index] },
        customdata: missionCustomdata(aircraft),
        hovertemplate: `%{customdata[0]}<br>Range %{customdata[1]:.0f} km<br>R / Rdes %{customdata[2]:.3f}<br>${valueLine}<extra>${name}</extra>`,
    }, extra || {});
}
function mirroredWingTrace(spanY, values, name, index, extra, valueLine) {
    const positiveSpan = spanY.slice();
    const mirroredSpan = positiveSpan.slice().reverse().map((value) => -value);
    const fullSpan = mirroredSpan.concat(positiveSpan);
    const fullValues = values.slice().reverse().concat(values);
    return Object.assign({
        x: fullSpan,
        y: fullValues,
        name,
        type: "scatter",
        mode: "lines",
        line: { color: COLORS[index], width: 2.2 },
        hovertemplate: `y %{x:.2f} m<br>${valueLine}<extra>${name}</extra>`,
    }, extra || {});
}
function missionPlot(id, ylabel, extractor, valueLine, customLayout, traceExtras) {
    const traces = state.data.aircraft.map((aircraft, index) => missionTrace(aircraft, extractor(aircraft), aircraft.name, index, traceExtras ? traceExtras(aircraft, index) : null, valueLine));
    const plotLayout = Object.assign({
        xaxis: missionAxisConfig(),
        yaxis: { title: ylabel },
        shapes: missionSelectionShapes(),
        hovermode: "closest",
    }, customLayout || {});
    plotChart(id, traces, layout(plotLayout));
}
function plotChart(id, traces, customLayout) {
    const node = document.getElementById(id);
    if (!node) {
        return;
    }
    if (!window.Plotly) {
        node.innerHTML = `<div class="chart-warning">Plotly.js is unavailable. Connect to the internet or vendor the pinned Plotly bundle locally if you need offline viewing.</div>`;
        return;
    }
    const plotLayout = normalizePlotlyMath(customLayout || layout());
    const plotTraces = Array.isArray(traces) ? traces.map((trace) => normalizePlotlyMath(trace)) : traces;
    window.Plotly.newPlot(id, plotTraces, plotLayout, {
        responsive: true,
        displaylogo: false,
        modeBarButtonsToRemove: ["lasso2d", "select2d", "autoScale2d"],
    });
}
function lineTrace(x, y, name, index, extra) {
    return Object.assign({
        x,
        y,
        name,
        type: "scatter",
        mode: "lines",
        line: { color: COLORS[index], width: 2.4 },
        hovertemplate: "%{y:.3f}<extra>" + name + "</extra>",
    }, extra || {});
}
function barTrace(x, y, name, index) {
    return {
        x,
        y,
        name,
        type: "bar",
        marker: {
            color: COLORS[index],
            opacity: 0.86,
            line: { color: COLORS[index], width: 1 },
        },
        hovertemplate: "%{y:.3f}<extra>" + name + "</extra>",
    };
}
function layout(overrides) {
    const base = {
        paper_bgcolor: "rgba(0,0,0,0)",
        plot_bgcolor: "#09141c",
        font: {
            family: "IBM Plex Sans, Avenir Next, Segoe UI, sans-serif",
            color: "#8fa5ae",
            size: 11,
        },
        margin: { l: 54, r: 20, t: 24, b: 42 },
        legend: {
            bgcolor: "rgba(0,0,0,0)",
            font: { size: 10 },
        },
        xaxis: {
            gridcolor: "rgba(214, 225, 232, 0.08)",
            linecolor: "rgba(214, 225, 232, 0.16)",
            zerolinecolor: "rgba(214, 225, 232, 0.12)",
            automargin: true,
        },
        yaxis: {
            gridcolor: "rgba(214, 225, 232, 0.08)",
            linecolor: "rgba(214, 225, 232, 0.16)",
            zerolinecolor: "rgba(214, 225, 232, 0.12)",
            automargin: true,
        },
        hovermode: "closest",
        hoverlabel: {
            bgcolor: "rgba(8, 17, 24, 0.96)",
            bordercolor: "rgba(215, 180, 106, 0.24)",
            font: { color: "#e9efe8", size: 11 },
        },
    };
    return Object.assign(base, overrides || {});
}
function overviewSnapshotRows() {
    if (state.data.comparison) {
        return state.data.comparison.summary;
    }
    return state.data.aircraft[0].summary.kpis.map((row) => ({
        label: row.label,
        value: row.value,
        digits: row.digits,
        lower_better: row.lower_better,
    }));
}
function missionSnapshotRows() {
    const pointText = missionPointText(state.pointIndex);
    return metricRows([
        { label: "\\(R_{des}\\) [km]", get: (aircraft) => aircraft.summary.range_km, digits: 0, lowerBetter: false },
        { label: `${pointText} \\(R\\) [km]`, get: (aircraft) => aircraft.mission.range_km[state.pointIndex], digits: 0, lowerBetter: false },
        { label: `${pointText} \\(R/R_{des}\\)`, get: (aircraft) => normalizedMissionPosition(aircraft, state.pointIndex), digits: 3, lowerBetter: false },
        { label: "\\(t\\) [hr]", get: (aircraft) => aircraft.mission.time_hr[state.pointIndex], digits: 2, lowerBetter: true },
        { label: "\\(h\\) [km]", get: (aircraft) => aircraft.mission.alt_km[state.pointIndex], digits: 2, lowerBetter: false },
        { label: "\\(M\\)", get: (aircraft) => aircraft.mission.mach[state.pointIndex], digits: 3, lowerBetter: false },
        { label: "\\(\\gamma\\) [deg]", get: (aircraft) => aircraft.mission.gamma_rad[state.pointIndex] * 180 / Math.PI, digits: 2, lowerBetter: false },
        { label: "\\(L/D\\)", get: (aircraft) => aircraft.mission.lod[state.pointIndex], digits: 2, lowerBetter: false },
        { label: "ROC [m/s]", get: (aircraft) => aircraft.mission.roc_mps[state.pointIndex], digits: 2, lowerBetter: false },
        { label: "\\(W/W_{MTO}\\)", get: (aircraft) => aircraft.mission.weight_fraction[state.pointIndex], digits: 3, lowerBetter: false },
        { label: "\\(C_D \\times 10^4\\)", get: (aircraft) => aircraft.mission.cd_e4[state.pointIndex], digits: 1, lowerBetter: true },
    ]);
}
function geometrySnapshotRows() {
    return metricRows([
        { label: "\\(b\\) [m]", get: (aircraft) => aircraft.summary.wing_span_m, digits: 2, lowerBetter: false },
        { label: "\\(S\\) [m^2]", get: (aircraft) => aircraft.summary.wing_area_m2, digits: 1, lowerBetter: false },
        { label: "\\(AR\\)", get: (aircraft) => aircraft.summary.wing_ar, digits: 2, lowerBetter: false },
        { label: "\\(\\Lambda\\) [deg]", get: (aircraft) => aircraft.summary.wing_sweep_deg, digits: 1, lowerBetter: false },
        { label: "\\(c_{root}\\) [m]", get: (aircraft) => aircraft.geometry.wing.root_chord_m, digits: 2, lowerBetter: false },
        { label: "\\(L_f\\) [m]", get: (aircraft) => aircraft.geometry.fuselage.length_m, digits: 2, lowerBetter: false },
        { label: "\\(W_f\\) [m]", get: (aircraft) => aircraft.geometry.fuselage.width_m, digits: 2, lowerBetter: false },
        { label: "\\(H_f\\) [m]", get: (aircraft) => aircraft.geometry.fuselage.height_m, digits: 2, lowerBetter: false },
        { label: "Cabin width [m]", get: (aircraft) => aircraft.geometry.fuselage.cabin_width_m, digits: 2, lowerBetter: false },
        { label: "\\(D_f\\) [m]", get: (aircraft) => aircraft.summary.fan_diameter_m, digits: 3, lowerBetter: false },
    ]);
}
function weightsSnapshotRows() {
    const structureLabels = state.data.aircraft[0].weights.structure_breakdown.labels;
    return metricRows([
        { label: "\\(W_{MTO}\\) [t]", get: (aircraft) => aircraft.weights.mtow_t, digits: 1, lowerBetter: true },
        { label: "\\(W_{OE}\\) [t]", get: (aircraft) => aircraft.weights.oew_t, digits: 1, lowerBetter: true },
        { label: "\\(W_f\\) [t]", get: (aircraft) => aircraft.weights.fuel_t, digits: 1, lowerBetter: true },
        { label: "\\(W_{pay}\\) [t]", get: (aircraft) => aircraft.weights.payload_t, digits: 1, lowerBetter: false },
    ].concat(structureLabels.map((label, index) => ({
        label: `${label} [t]`,
        get: (aircraft) => aircraft.weights.structure_breakdown.values_t[index],
        digits: 2,
        lowerBetter: true,
    }))));
}
function aeroSnapshotRows() {
    const pointText = missionPointText(state.pointIndex);
    const dragKeys = state.data.aircraft[0].aero.drag_keys;
    return metricRows([
        { label: `${pointText} \\(L/D\\)`, get: (aircraft) => aircraft.aero.lod[state.pointIndex], digits: 2, lowerBetter: false },
    ].concat(dragKeys.map((key, index) => ({
        label: `${dragKeyLabel(key)} fraction`,
        get: (aircraft) => aircraft.aero.drag_fractions[index][state.pointIndex],
        digits: 3,
        lowerBetter: index === 0,
    }))).concat(dragKeys.map((key, index) => ({
        label: `${dragKeyLabel(key)} @ CR1 [cm^2]`,
        get: (aircraft) => aircraft.aero.drag_area_cruise_cm2[index],
        digits: 1,
        lowerBetter: true,
    }))));
}
function engineSnapshotRows() {
    return metricRows([
        { label: "\\(F_{tot}\\) [kN]", get: (aircraft) => aircraft.engine.thrust_total_kN[state.pointIndex], digits: 1, lowerBetter: false },
        { label: "\\(F_e\\) [kN]", get: (aircraft) => aircraft.engine.thrust_per_engine_kN[state.pointIndex], digits: 1, lowerBetter: false },
        { label: "\\(\\mathrm{TSFC}\\) [mg/N/s]", get: (aircraft) => aircraft.engine.tsfc_mg_ns[state.pointIndex], digits: 3, lowerBetter: true },
        { label: "\\(\\mathrm{BPR}\\)", get: (aircraft) => aircraft.engine.bpr[state.pointIndex], digits: 2, lowerBetter: false },
        { label: "\\(\\mathrm{OPR}\\)", get: (aircraft) => aircraft.engine.opr[state.pointIndex], digits: 1, lowerBetter: false },
        { label: "Fan \\(\\pi_f\\)", get: (aircraft) => aircraft.engine.fan_pr[state.pointIndex], digits: 3, lowerBetter: false },
        { label: "LPC \\(\\pi_{LPC}\\)", get: (aircraft) => aircraft.engine.lpc_pr[state.pointIndex], digits: 3, lowerBetter: false },
        { label: "HPC \\(\\pi_{HPC}\\)", get: (aircraft) => aircraft.engine.hpc_pr[state.pointIndex], digits: 3, lowerBetter: false },
        { label: "Fan \\(\\eta_p\\)", get: (aircraft) => aircraft.engine.fan_eff[state.pointIndex], digits: 3, lowerBetter: false },
        { label: "LPC \\(\\eta_p\\)", get: (aircraft) => aircraft.engine.lpc_eff[state.pointIndex], digits: 3, lowerBetter: false },
        { label: "HPC \\(\\eta_p\\)", get: (aircraft) => aircraft.engine.hpc_eff[state.pointIndex], digits: 3, lowerBetter: false },
        { label: "\\(T_{t4}\\) [K]", get: (aircraft) => aircraft.engine.tt4_k[state.pointIndex], digits: 0, lowerBetter: true },
    ]);
}
function metricRows(defs) {
    if (!isCompareMode()) {
        return defs.map((definition) => ({
            label: definition.label,
            value: definition.get(state.data.aircraft[0]),
            digits: definition.digits,
            lower_better: definition.lowerBetter,
        }));
    }
    return defs.map((definition) => ({
        label: definition.label,
        left: definition.get(state.data.aircraft[0]),
        right: definition.get(state.data.aircraft[1]),
        digits: definition.digits,
        lower_better: definition.lowerBetter,
    }));
}
function isCompareMode() {
    return state.data && state.data.aircraft && state.data.aircraft.length > 1;
}
function deltaVisual(row) {
    if (!isFiniteMetric(row.left) || !isFiniteMetric(row.right)) {
        return { className: "flat", text: "n/a" };
    }
    const deltaAbs = row.right - row.left;
    const deltaPct = Math.abs(row.left) < Number.EPSILON ? null : 100 * deltaAbs / Math.abs(row.left);
    if (deltaPct === null || Math.abs(deltaPct) < 0.1) {
        return { className: "flat", text: "No material change" };
    }
    const improved = row.lower_better ? deltaPct < 0 : deltaPct > 0;
    return {
        className: improved ? "good" : "bad",
        text: `${deltaPct > 0 ? "+" : ""}${formatNumber(deltaPct, 1)}% vs baseline`,
    };
}
function dragKeyLabel(key) {
    return ({
        CDi: "\\(C_{D,i}\\)",
        CDfuse: "\\(C_{D,fuse}\\)",
        CDwing: "\\(C_{D,wing}\\)",
        CDhtail: "\\(C_{D,htail}\\)",
        CDvtail: "\\(C_{D,vtail}\\)",
        CDnace: "\\(C_{D,nace}\\)",
    }[key] || escapeHtml(String(key)));
}
function renderMathLabel(label) {
    return renderRichText(label);
}
function renderRichText(text) {
    return replaceInlineMath(String(text ?? ""), (content) => escapeHtml(texToReadableText(content)));
}
function isFiniteMetric(value) {
    return value !== null && value !== undefined && Number.isFinite(Number(value));
}
function formatMetricValue(value, digits) {
    return isFiniteMetric(value) ? formatNumber(Number(value), digits) : "n/a";
}
function formatNumber(value, digits) {
    return Number(value).toFixed(digits);
}
function clamp(value, min, max) {
    return Math.max(min, Math.min(max, value));
}
function trimHiddenMissionPoint(data) {
    const short = data && data.meta && Array.isArray(data.meta.mission_short) ? data.meta.mission_short : [];
    const labels = data && data.meta && Array.isArray(data.meta.mission_labels) ? data.meta.mission_labels : [];
    if (!short.length) {
        return data;
    }
    const lastShort = String(short[short.length - 1] || "").toUpperCase();
    const lastLabel = String(labels[labels.length - 1] || "").toLowerCase();
    if (lastShort !== "TST" && lastLabel !== "test") {
        return data;
    }
    data.meta.mission_short = short.slice(0, -1);
    data.meta.mission_labels = labels.slice(0, -1);
    data.meta.default_point = Math.min(Number(data.meta.default_point || 1), data.meta.mission_short.length);
    data.aircraft.forEach((aircraft) => {
        trimLastPoints(aircraft.mission, ["range_km", "time_hr", "gamma_rad", "alt_km", "mach", "lod", "roc_mps", "weight_fraction", "cd_e4"]);
        if (aircraft.aero) {
            aircraft.aero.lod = trimLast(aircraft.aero.lod);
            aircraft.aero.drag_fractions = (aircraft.aero.drag_fractions || []).map((series) => trimLast(series));
            if (aircraft.aero.wing) {
                aircraft.aero.wing.cl_by_point = trimLast(aircraft.aero.wing.cl_by_point);
            }
            if (aircraft.aero.trefftz) {
                trimLastPoints(aircraft.aero.trefftz, ["circulation_by_point", "wake_circulation_by_point", "induced_velocity_mps", "drag_density_proxy"]);
            }
        }
        if (aircraft.engine) {
            trimLastPoints(aircraft.engine, [
                "temperature_k", "pressure_kpa", "thrust_total_kN", "thrust_per_engine_kN",
                "tsfc_mg_ns", "bpr", "opr", "fpr", "fan_pr", "lpc_pr", "hpc_pr",
                "fan_eff", "lpc_eff", "hpc_eff", "fan_speed_frac", "lpc_speed_frac", "hpc_speed_frac",
                "fan_mass_flow_frac", "lpc_mass_flow_frac", "hpc_mass_flow_frac", "tt4_k",
            ]);
            if (aircraft.engine.maps) {
                Object.values(aircraft.engine.maps).forEach((mapPayload) => {
                    if (mapPayload && mapPayload.operating_line) {
                        trimLastPoints(mapPayload.operating_line, ["x", "y", "speed_frac", "eff"]);
                    }
                });
            }
        }
    });
    return data;
}
function trimLastPoints(target, keys) {
    if (!target) {
        return;
    }
    keys.forEach((key) => {
        if (Array.isArray(target[key])) {
            target[key] = trimLast(target[key]);
        }
    });
}
function trimLast(values) {
    return Array.isArray(values) ? values.slice(0, -1) : values;
}
function typesetMath() {
    return;
}
function escapeHtml(value) {
    return String(value)
        .replaceAll("&", "&amp;")
        .replaceAll("<", "&lt;")
        .replaceAll(">", "&gt;")
        .replaceAll("\"", "&quot;")
        .replaceAll("'", "&#39;");
}
function escapeAttr(value) {
    return escapeHtml(value).replaceAll("\n", "&#10;");
}
function replaceInlineMath(text, renderer) {
    const pattern = /\\\((.+?)\\\)|\$([^$]+?)\$/g;
    let result = "";
    let lastIndex = 0;
    let matched = false;
    let match;
    while ((match = pattern.exec(text)) !== null) {
        matched = true;
        if (match.index > lastIndex) {
            result += escapeHtml(text.slice(lastIndex, match.index));
        }
        result += renderer(match[1] || match[2] || "");
        lastIndex = pattern.lastIndex;
    }
    if (!matched) {
        return escapeHtml(text);
    }
    if (lastIndex < text.length) {
        result += escapeHtml(text.slice(lastIndex));
    }
    return result;
}
function texToReadableText(tex) {
    let value = String(tex ?? "").trim();
    let previous = "";
    while (value !== previous) {
        previous = value;
        value = value.replace(/\\frac\{([^{}]+)\}\{([^{}]+)\}/g, (_match, num, den) => `${texToReadableText(num)}/${texToReadableText(den)}`);
        value = value.replace(/\\mathrm\{([^{}]+)\}/g, "$1");
        value = value.replace(/\\text\{([^{}]+)\}/g, "$1");
        value = value.replace(/\\dot\{m\}/g, "ṁ");
        value = value.replace(/([A-Za-zΑ-Ωα-ωΓΛγληπṁ]+)_\{([^{}]+)\}/g, (_match, base, sub) => `${base}_${texToReadableText(sub)}`);
        value = value.replace(/([A-Za-zΑ-Ωα-ωΓΛγληπṁ]+)_([A-Za-z0-9,]+)/g, "$1_$2");
        value = value.replace(/([A-Za-zΑ-Ωα-ωΓΛγληπṁ0-9]+)\^\{([^{}]+)\}/g, (_match, base, sup) => `${base}^${texToReadableText(sup)}`);
        value = value.replace(/([A-Za-zΑ-Ωα-ωΓΛγληπṁ0-9]+)\^([A-Za-z0-9+-]+)/g, "$1^$2");
    }
    return value
        .replace(/\\gamma/g, "γ")
        .replace(/\\Gamma/g, "Γ")
        .replace(/\\Lambda/g, "Λ")
        .replace(/\\eta/g, "η")
        .replace(/\\pi/g, "π")
        .replace(/\\times/g, "×")
        .replace(/\\cdot/g, "·")
        .replace(/[{}]/g, "")
        .replace(/\s+/g, " ")
        .trim();
}
function plotlyText(text) {
    return replaceInlineMath(String(text ?? ""), (content) => texToReadableText(content));
}
function normalizePlotlyMath(value, key = "") {
    if (Array.isArray(value)) {
        return value.map((item) => normalizePlotlyMath(item, key));
    }
    if (value && typeof value === "object") {
        const normalized = {};
        Object.entries(value).forEach(([entryKey, entryValue]) => {
            normalized[entryKey] = normalizePlotlyMath(entryValue, entryKey);
        });
        return normalized;
    }
    if (typeof value === "string" && (key === "title" || key === "text")) {
        return plotlyText(value);
    }
    return value;
}
})();
