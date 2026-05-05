# Cursor AI Agent Setup: Rust + Plonky3 Experiments

A step-by-step guide to configuring a Cursor AI agent context with `git-workflow-and-versioning` and `source-driven-development` skills from [addyosmani/agent-skills](https://github.com/addyosmani/agent-skills), targeting a Rust/Plonky3 experiment project.

---

## Step 1: Scaffold the Project

Create the Rust project and add Plonky3 as a git-sourced dependency:

```bash
cargo new plonky3-experiments
cd plonky3-experiments
```

Edit `Cargo.toml` to add the Plonky3 crates:

```toml
[package]
name = "plonky3-experiments"
version = "0.1.0"
edition = "2021"

[dependencies]
p3-field = "<check Cargo.toml>"
p3-matrix = "<check Cargo.toml>"
p3-goldilocks = "<check Cargo.toml>"
p3-fri = "<check Cargo.toml>"
p3-merkle-tree = "<check Cargo.toml>"

[dev-dependencies]
criterion = "0.5"
```

> **Why this matters for source-driven-development:** The agent reads `Cargo.toml` before writing any code to identify exact crate versions. With crates.io dependencies, the version numbers in `Cargo.toml` become the source of truth for which APIs and docs to verify before implementation.

---

## Step 2: Create Your `.cursor/rules` Directory

Cursor reads agent rules from `.cursor/rules/` at project root. Create the directory structure:

```
plonky3-experiments/
├── .cursor/
│   └── rules/
│       ├── git-workflow.mdc
│       ├── source-driven-dev.mdc
│       └── project-context.mdc
├── src/
│   └── lib.rs
├── Cargo.toml
└── .gitignore
```

### `git-workflow.mdc`

Paste the full content of [`skills/git-workflow-and-versioning/SKILL.md`](https://github.com/addyosmani/agent-skills/blob/main/skills/git-workflow-and-versioning/SKILL.md) into this file verbatim. This instructs the agent to:

- Work in short-lived feature branches (1–3 days max)
- Make atomic commits of ~100 lines per increment
- Use conventional commit message types (`feat`, `fix`, `refactor`, `test`, `docs`, `chore`)
- Emit a structured **CHANGES MADE / DIDN'T TOUCH / CONCERNS** summary after every modification
- Treat every passing test as a git save point

### `source-driven-dev.mdc`

Paste the full content of [`skills/source-driven-development/SKILL.md`](https://github.com/addyosmani/agent-skills/blob/main/skills/source-driven-development/SKILL.md) into this file verbatim. This instructs the agent to:

- Read `Cargo.toml` first and declare the detected stack explicitly
- Fetch official Plonky3 source/docs before writing any API call
- Cite every non-trivial decision with a full URL in a code comment
- Flag anything that cannot be verified as `UNVERIFIED` rather than guessing
- Surface conflicts between existing code and current documented patterns

### `project-context.mdc`

Write a short, project-specific context file:

```markdown
# Project: plonky3-experiments

## Stack
- Rust (check Cargo.toml for exact crate versions)
- Plonky3 crates from crates.io

## Goals
- Experiment with field arithmetic over the Goldilocks field
- Prototype FRI-based polynomial commitments
- Build small circuits and verify proofs end-to-end

## Module Conventions
- One module per concept: `fields/`, `fri/`, `circuits/`
- Every public `fn` has a doc comment with a usage example
- Tests live alongside source in `#[cfg(test)]` modules
- Benchmarks in `benches/` using Criterion

## Increment Size
Keep each agent increment to ~100 lines as per git-workflow skill.
Split into separate branches if a feature spans more than ~300 lines.
```

---

## Step 3: Initialise Git Properly

Set up git with a `.gitignore` and make the first save-point commit:

```bash
# .gitignore
cat <<'EOF' > .gitignore
/target
Cargo.lock
.env
*.pem
EOF

git init
git add Cargo.toml .gitignore .cursor/
git commit -m "chore: initialise plonky3-experiments with Cursor agent rules"
```

> **The save-point pattern:** Every commit is a stable state you can `git reset --hard HEAD` back to. This is the foundational safety net from the git-workflow skill — if the agent drifts or breaks the build, you lose at most one increment of work.

From this point, the recommended branching model is **trunk-based development**:

```
main ──●──────────●──────────●──────────●──  (always compiles + tests pass)
        ╲          ╱  ╲        ╱
         ●──●──●──╱    ●──●──╱
         feature/      feature/
         goldilocks-   poly-commit
         field-arith   (1-3 days)
```

---

## Step 4: Launching the Agent in Increment Mode

Open Cursor and activate **Agent mode** (`Cmd+I` / `Ctrl+I`). Send this structured opening prompt to prime the agent correctly:

```
Read all files in `.cursor/rules/` to understand your workflow constraints.
Then read `Cargo.toml` to detect exact crates.io versions and declare them.

First feature: implement a basic field arithmetic module over the
Goldilocks field using `p3-goldilocks`. Apply source-driven-development:
fetch the official docs or published crate source that matches the
versions in `Cargo.toml` before writing any code. Keep this increment
under 100 lines. After tests pass, emit a CHANGES MADE / DIDN'T TOUCH /
CONCERNS summary and propose a conventional commit message.
```

This single prompt activates both skills in concert:

| Prompt clause | Skill activated | Expected agent behaviour |
|---|---|---|
| "Read `.cursor/rules/`" | Both | Agent reads rules before acting |
| "Read `Cargo.toml`" | source-driven-dev | Declares detected stack explicitly |
| "Fetch Plonky3 source before writing" | source-driven-dev | No hallucinated APIs |
| "Keep under 100 lines" | git-workflow | Scoped, reviewable increment |
| "Emit CHANGES MADE summary" | git-workflow | Structured change audit |
| "Propose a commit message" | git-workflow | Closes the save-point loop |

---

## Step 5: The Increment Loop

Each feature session follows this repeating rhythm:

```
1. Agent reads Cargo.toml → declares stack (source-driven-dev)
       ↓
2. Agent fetches the crates.io-matching Plonky3 docs / relevant source for the feature
       ↓
3. Implements ≤100 lines with inline source-citation comments
       ↓
4. Runs: cargo test && cargo clippy
       ↓
5. Emits: CHANGES MADE / DIDN'T TOUCH / CONCERNS summary
       ↓
6. You review → approve → git commit -m "feat: ..."
       ↓
7. Next increment begins from clean HEAD  ──→ back to step 1
```

### Disciplined Prompt Patterns

Use these phrases consistently to keep the agent on-track across sessions:

- **`"Cite your sources before writing code"`** — triggers source-driven-dev behaviour
- **`"Keep this under 100 lines, split if needed"`** — enforces atomic increment size
- **`"Do not touch X, it is out of scope"`** — agent must list it in DIDN'T TOUCH
- **`"Run cargo test then propose a commit"`** — closes the save-point loop

### If the Agent Drifts

```bash
# Revert to the last known-good state instantly
git reset --hard HEAD

# Or inspect recent history to find a safe point
git log --oneline -10
git reset --hard <sha>
```

### Suggested Feature Sequence

| Increment | Module | Plonky3 crates |
|---|---|---|
| 1 | Goldilocks field arithmetic | `p3-field`, `p3-goldilocks` |
| 2 | Dense polynomial representation | `p3-field`, `p3-matrix` |
| 3 | FFT / NTT over the field | `p3-dft` |
| 4 | Merkle commitment to evaluations | `p3-merkle-tree` |
| 5 | FRI low-degree test | `p3-fri` |
| 6 | Simple constraint circuit | all of the above |

Each increment is a short-lived branch merged back to `main` once `cargo test` passes — consistent with the trunk-based model from the git-workflow skill.

---

*Skills sourced from [github.com/addyosmani/agent-skills](https://github.com/addyosmani/agent-skills). Project context updated for crates.io-based Plonky3 dependencies.*
