# Swig 0.18.0

## Project directories

- Adds optional read/write project directories for supported secure-context browsers.
- Creates automatically numbered `runs/NNN-name/` directories instead of overwriting earlier analyses.
- Writes the AIRR table incrementally into the active run.
- Writes an inspectable `inputs/datasets.tsv` study manifest and a root layout README.
- Writes `state/latest.swig-session.json.gz` plus numbered, stage-named state checkpoints.
- Appends structured run and checkpoint events to `logs/events.jsonl`.
- Opening an existing project restores the latest session and its linked AIRR result without a second file-selection step.
- Portable session files remain available independently of project-directory support.

## Directory data loading

- Files or whole directory trees can be dragged onto the data panel; directory selection is also available as a button.
- Relative input paths are retained in the run manifest.
- For nested trees, the first directory below the selected root initializes `Donor / subject` for descendant files.
- A flat directory asks whether its files share one donor or represent separate donors.
- Inferred donor values remain editable before annotation and in the Results metadata editor.

## Terminology

- Replaces user-facing “upload” labels with “load”, “local file”, or “provided file” terminology. Internal compatibility keys used by earlier session files are unchanged.

## Verification

- Adds unit tests for flat and nested donor inference and for the project manifest/run/checkpoint/restore layout.
- Production TypeScript and static-site builds remain part of the release gate.
