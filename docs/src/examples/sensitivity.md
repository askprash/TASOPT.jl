# Calculating Parameter Sensitivities

![SensitivityPlot](../assets/sensitivity_plot.svg)

Using the sensitivity module it is possible to calculate the gradients of TASOPT parameters over PFEI. The sensitivity is calculated using finite difference (central relative difference).

The function input takes parameters as a list of symbols. There are a variety of parameters that can be taken in shown in the example below. The sensitivities can also be plot in a bar chart.

```julia
using TASOPT
include(__TASOPTindices__)
# List of the parameters you want to update as symbols
params = [
    # parg index — single scalar
    :(ac.parg[igWfmax]),
    # para array — specific mission point and mission index
    :(ac.para[iaCL, ipcruise1, 1]),
    # Structural nested field (no array indices) — use dot-path notation
    :(ac.fuselage.layout.cross_section.radius)
]
# Engine operating parameters (e.g. Tt4, pif, pihc, epolf) live in typed EngineState:
#   ac.missions[1].points[ip].engine.Tt4
#   ac.missions[1].points[ip].engine.pif
#   ac.missions[1].points[ip].engine.design.epolf
# Set these via direct assignment before calling size_aircraft! rather than
# through the params list.
epsilon = 1e-5
default_model = load_default_model()
size_aircraft!(default_model)
TASOPT.get_sensitivity(params, model_state = default_model, eps = epsilon)
```

If you want the default model as the model state and epsilon as 1e-5 
you can also call the function directly with just the params:

```julia
sens = TASOPT.get_sensitivity(params)

#Plot Sensitivities
TASOPT.plot_sensitivities(sens)
```
