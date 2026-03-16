# TASOPT Diagnostic Dashboard

This dashboard is a static browser app for inspecting one or two saved TASOPT aircraft.
Julia is responsible for exporting the dashboard data. The browser only renders the
exported JSON snapshot.

## User workflow

### 1. Export dashboard JSON

From the TASOPT repo root:

```julia
using TASOPT

export_dashboard_data("ac1.jld2", "ac2.jld2"; filepath = "compare.json")
```

Single-aircraft export also works:

```julia
using TASOPT

export_dashboard_data("ac1.jld2"; filepath = "single.json")
```

You can also pass sized `aircraft` objects directly instead of file paths.

### 2. Open the dashboard

```julia
using TASOPT

TASOPT.start_dashboard()
```

This opens the shipped dashboard shell in your browser.

### 3. Load the exported JSON

In the browser:

1. Click `Load JSON`
2. Select the exported file such as `compare.json`

The dashboard will then render in:

- single-aircraft mode if the JSON contains one aircraft
- compare mode if the JSON contains two aircraft

## Notes

- The normal workflow does not require users to build the frontend.
- The shipped static app lives in `dashboard/dist/`.
- `export_dashboard_bundle(...)` still exists, but the default workflow is:
  1. `export_dashboard_data(...)`
  2. `TASOPT.start_dashboard()`
  3. choose the JSON in the browser

## Frontend development

If you are editing the dashboard UI:

```bash
cd dashboard
npm run dev
```

For a fresh shipped build:

```bash
cd dashboard
npm run build
```
