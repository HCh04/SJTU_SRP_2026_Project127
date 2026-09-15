# AISynbio Repository Instructions for Kimi Code

## Purpose of this file

This file tells Kimi Code how to work in `/home/ubuntu/AISynbio`.

It is a repository-level development guide, not the runtime system prompt for
the scientific agent. Runtime scientific instructions must live in a separate
file such as `prompts/aisynbio_system.md` and must be loaded explicitly by the
application.

Follow the instructions in this file for all work under this directory. Read
the relevant source files before editing them, keep changes narrowly scoped,
and verify behavior before reporting completion.

## Project mission

Build an AI agent for enzyme and metabolic-network optimization in synthetic
biology. The product should connect enzyme-level engineering to cell-level
metabolic modeling through a traceable Design-Build-Test-Learn (DBTL) loop.

This is not intended to be a general biology chatbot. It should coordinate
specialized databases, deterministic scientific tools, optimization methods,
and experimental results to produce ranked and testable engineering designs.

The language model may plan, select tools, reconcile evidence, explain tradeoffs,
and propose the next experiment. It must not invent biological measurements or
replace deterministic simulation with unsupported prose.

## Initial scope

Keep the first implementation deliberately narrow:

- Host: yeast, unless a task explicitly specifies another host.
- Input: target molecule plus optional host, feedstock, culture, laboratory,
  growth, yield, titer, productivity, and safety constraints.
- Student A: pathway discovery and enzyme-candidate pipeline.
- Student B: metabolic-network and cell-performance simulation pipeline.
- AISynbio: validates and connects both pipelines, coordinates Biomni, ranks
  interventions, records uncertainty, and produces an engineering report.

Do not claim that the prototype supports arbitrary hosts, complete de novo
enzyme design, automatic plasmid design, or autonomous wet-lab execution until
those capabilities are implemented and validated.

## Current repository state

The expected layout is:

```text
/home/ubuntu/AISynbio/
|-- AGENTS.md              # Instructions for Kimi Code
|-- agent.py               # Current orchestration prototype
|-- prompts/               # Runtime prompts; create when implementing them
|-- tests/                 # Automated tests; create when implementing them
|-- Biomni/                # Nested clone of upstream Biomni
|-- data/                  # Local Biomni data cache
|-- venv/                  # Python 3.12 virtual environment
`-- .env                   # Credentials and provider configuration
```

Important facts:

- The top-level AISynbio directory is not currently a Git repository.
- `Biomni/` is a nested Git clone, not a Git submodule.
- Biomni is installed editable from `Biomni/` as version `0.0.8`.
- The inspected Biomni revision is
  `400c1f366b96a35ca253e13c9b06c5076af41d65`.
- `agent.py` still contains placeholder Student A and Student B results.
- The local Biomni data lake is partial. Passing
  `expected_data_lake_files=[]` skips the data-lake check; it does not mean no
  data is required.
- Another collaborator may edit files through SSH. Never overwrite or revert
  changes that you did not make. Re-read a file immediately before patching it.

Do not initialize Git, replace the Biomni revision, download the complete data
lake, or install system packages unless the user asks or approves it.

## Architecture boundaries

Maintain these ownership boundaries:

1. `agent.py` or a future application module owns orchestration and dependency
   wiring.
2. `prompts/` owns runtime agent behavior and output contracts.
3. Student A owns pathway and enzyme prediction.
4. Student B owns metabolic simulation.
5. AISynbio adapters validate, normalize, and version Student A/B data.
6. Biomni supplies agent orchestration, tool retrieval, selected database
   wrappers, and selected scientific utilities.
7. Deterministic software performs numerical calculations. LLM text is never a
   substitute for a model run.

Prefer adapters around Student A, Student B, and Biomni. Do not put all business
logic into one large prompt or one `agent.py` file.

## Biomni integration

Upstream project:
[snap-stanford/Biomni](https://github.com/snap-stanford/Biomni)

Pinned source used during this review:
[revision 400c1f3](https://github.com/snap-stanford/Biomni/tree/400c1f366b96a35ca253e13c9b06c5076af41d65)

### Useful Biomni components

- [`biomni/agent/a1.py`](https://github.com/snap-stanford/Biomni/blob/400c1f366b96a35ca253e13c9b06c5076af41d65/biomni/agent/a1.py)
  provides the `A1` agent, LangGraph workflow, `add_tool()`, `add_data()`,
  `add_software()`, `go()`, `go_stream()`, and conversation export.
- [`biomni/tool/tool_registry.py`](https://github.com/snap-stanford/Biomni/blob/400c1f366b96a35ca253e13c9b06c5076af41d65/biomni/tool/tool_registry.py)
  stores and resolves tool schemas.
- [`biomni/model/retriever.py`](https://github.com/snap-stanford/Biomni/blob/400c1f366b96a35ca253e13c9b06c5076af41d65/biomni/model/retriever.py)
  selects tools and resources relevant to a task.
- [`biomni/tool/database.py`](https://github.com/snap-stanford/Biomni/blob/400c1f366b96a35ca253e13c9b06c5076af41d65/biomni/tool/database.py)
  includes wrappers for UniProt, AlphaFold, InterPro, PDB, KEGG, PubChem,
  ChEMBL, Reactome, and other databases.
- [`biomni/tool/systems_biology.py`](https://github.com/snap-stanford/Biomni/blob/400c1f366b96a35ca253e13c9b06c5076af41d65/biomni/tool/systems_biology.py)
  contains a baseline COBRA-based FBA helper.
- [`biomni/tool/synthetic_biology.py`](https://github.com/snap-stanford/Biomni/blob/400c1f366b96a35ca253e13c9b06c5076af41d65/biomni/tool/synthetic_biology.py)
  contains SBML construction and codon-optimization helpers.
- [`biomni/tool/literature.py`](https://github.com/snap-stanford/Biomni/blob/400c1f366b96a35ca253e13c9b06c5076af41d65/biomni/tool/literature.py)
  supports targeted PubMed, arXiv, web, PDF, and supplementary-material
  retrieval.
- [`biomni/config.py`](https://github.com/snap-stanford/Biomni/blob/400c1f366b96a35ca253e13c9b06c5076af41d65/biomni/config.py)
  defines provider, model, path, timeout, retrieval, and licensing settings.
- [`docs/mcp_integration.md`](https://github.com/snap-stanford/Biomni/blob/400c1f366b96a35ca253e13c9b06c5076af41d65/docs/mcp_integration.md)
  describes MCP integration. MCP may later expose Student A/B as services.

### Biomni limitations that must be preserved in code and reports

- A database wrapper is not a verified scientific conclusion. Validate IDs,
  species, assay conditions, units, evidence type, and source version.
- BRENDA, MetaCyc, SABIO-RK, ChEBI, eQuilibrator, and enzyme-constrained GEM
  support are not built-in capabilities in this checkout. Treat them as future
  custom tools or MCP integrations.
- `analyze_enzyme_kinetics_assay()` generates synthetic data from hard-coded
  parameters and random noise. Never report its output as measured, retrieved,
  or predicted enzyme kinetics.
- `simulate_metabolic_network_perturbation()` uses generic mass-action rates.
  Treat it as demonstration code, not a quantitative strain-design model.
- `perform_flux_balance_analysis()` is a basic FBA wrapper. It does not provide
  a validated full pipeline for FVA, pFBA, MOMA, ROOM, OptKnock, thermodynamic
  constraints, kinetic modeling, or enzyme-constrained modeling.
- Check the license and `commercial_mode` implications of datasets before
  including them in commercial output.

### Correct custom-tool pattern

Biomni `A1.add_tool()` inspects Python source to generate a schema. Custom tools
must therefore be ordinary top-level Python functions with descriptive
docstrings and JSON-serializable inputs and outputs.

Do not decorate AISynbio functions with LangChain `@tool`. That produces a
`StructuredTool` object that is incompatible with Biomni's source inspection.

Use this integration shape:

```python
def get_data_from_student_a(compound: str) -> dict:
    """Return validated pathway and enzyme candidates for a compound."""
    ...


def run_simulation_from_student_b(data_a: dict) -> dict:
    """Return a validated metabolic simulation for one design candidate."""
    ...


agent = BiomniAgent(
    path=data_path,
    llm=model_name,
    source=model_source,
    expected_data_lake_files=required_data_files,
)
agent.add_tool(get_data_from_student_a)
agent.add_tool(run_simulation_from_student_b)
result = agent.go(runtime_prompt)
```

Register tools before calling `go()`. Do not create an unused local `tools`
list and assume Biomni can access it.

## Data contracts

Do not pass loosely defined dictionaries between Student A, Student B, and the
agent once the adapters are implemented. Introduce typed models, preferably
with Pydantic, for at least:

- `EngineeringSpec`: target, host, feedstock, culture conditions, objectives,
  constraints, and available laboratory operations.
- `PathwayCandidate`: ordered reactions, metabolite IDs, enzyme candidates,
  cofactors, compartments, host compatibility, evidence, and uncertainty.
- `EnzymeEvidence`: sequence/accession, organism, reaction, kinetic value,
  units, temperature, pH, assay type, source, and confidence.
- `SimulationInput`: model version, medium, bounds, objective, interventions,
  enzyme constraints, and assumptions.
- `SimulationResult`: solver status, objective value, fluxes, yield, growth,
  byproducts, sensitivities, warnings, and artifact paths.
- `DesignPackage`: ranked intervention, supporting evidence, predicted benefit,
  uncertainty, risks, alternatives, and next experiment.

Every numeric field must have explicit units or be declared dimensionless.
Never silently convert IDs, units, compartments, organisms, or reaction
directions. Reject invalid or incomplete inputs with actionable errors.

## Scientific output contract

The runtime agent should return a structured design package, not only prose.
Every recommendation must distinguish:

- observed experimental data;
- database-retrieved data;
- deterministic model output;
- ML prediction;
- LLM inference or hypothesis.

For each result record the source, version, biological context, assumptions,
uncertainty, and artifact produced by the tool run. Never invent `kcat`, `Km`,
yield, titer, productivity, flux, citations, or experimental outcomes.

When data is missing, return a knowledge gap and propose the smallest useful
measurement or query. Do not fill the gap with a plausible-looking number.

Rank designs using explicit objectives and constraints. Prefer a Pareto set when
yield, growth, productivity, robustness, and experimental cost conflict.

## Runtime DBTL behavior

The runtime prompt, not this file, should instruct AISynbio to follow:

1. Design: assemble and rank pathway and enzyme candidates.
2. Build: produce a proposed construct manifest only when sequence and part data
   are validated.
3. Test: run deterministic models under recorded conditions.
4. Learn: compare predictions with observations, update parameters or beliefs,
   and propose the highest-information next experiment.

The first milestone only needs Design and in-silico Test orchestration with
mocked Student A/B adapters. Do not claim Build or Learn is complete until real
construct and experimental-data interfaces exist.

## Implementation priorities

Unless the user gives a different task, prioritize work in this order:

1. Move the runtime scientific prompt into `prompts/aisynbio_system.md` and load
   it explicitly.
2. Remove LangChain `@tool` wrappers from Student A/B adapters.
3. Register both functions using `agent.add_tool()`.
4. Remove hard-coded provider URLs and model settings from source. Read them
   from `.env` or typed configuration without printing secrets.
5. Add typed request and response models plus validation at both adapter
   boundaries.
6. Replace fixed biological values with clearly labelled mock fixtures in
   tests. Production code must fail clearly while a real adapter is absent.
7. Add a deterministic orchestration test that uses fake tool implementations
   and does not call an external LLM.
8. Add an opt-in live smoke test for one target after deterministic tests pass.
9. Add artifact and provenance storage for each DBTL run.
10. Only then expand database, modeling, UI, or wet-lab integrations.

## Coding rules

- Use the existing virtual environment: `/home/ubuntu/AISynbio/venv`.
- Prefer small modules with clear ownership over a growing monolithic
  `agent.py`.
- Use type hints for public functions and Pydantic models at external
  boundaries.
- Keep scientific computation deterministic where possible. Pin random seeds
  in tests and record solver/model versions.
- Do not edit upstream code inside `Biomni/` unless an adapter cannot solve the
  problem. If upstream modification is necessary, explain why and isolate it.
- Do not use `sudo`, install packages, download large datasets, change firewall
  rules, or start public services without explicit approval.
- Never read, print, log, commit, or expose values from `.env`.
- Do not hard-code credentials, API base URLs, user-specific absolute paths, or
  fabricated scientific results.
- Do not overwrite collaborator changes. Inspect timestamps and diffs before
  editing shared files.
- Avoid unrelated refactors while implementing a requested feature.
- Comments should explain scientific assumptions or non-obvious constraints,
  not narrate simple code.

## Verification

Before claiming a change works:

1. Run a syntax or import check for every changed Python module.
2. Run focused tests using `venv/bin/python -m pytest` when tests exist.
3. Run `venv/bin/pip check` after dependency changes.
4. Test invalid, incomplete, and tool-failure inputs as well as the happy path.
5. Confirm tests do not make paid API calls unless explicitly marked as live.
6. Confirm no credential values or generated data artifacts are included in
   source control or terminal output.
7. State exactly what was tested and what remains unverified.

## Definition of done

A task is complete only when:

- the requested behavior is implemented rather than described;
- Student A/B or Biomni boundaries are explicit and validated;
- deterministic tests cover the changed behavior;
- failures are returned clearly without invented fallback data;
- scientific provenance and units are preserved;
- no secrets or unrelated collaborator work were changed;
- documentation matches the behavior that actually exists.

## Communication

When reporting work, be concise and specific. Lead with blockers, scientific
risks, or failed validation. Then summarize changed files, tests run, and the
next unresolved integration boundary. Do not describe planned capabilities as
already implemented.
