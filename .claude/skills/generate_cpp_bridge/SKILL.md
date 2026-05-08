---
name: generate_cpp_bridge
description: Wrap an existing MOM6 Fortran subroutine in a runtime-dispatched shim that selects between (a) the original Fortran code, (b) a binary capture mode that records inputs+outputs to disk for offline validation, and (c) a C++/AMReX bridge invoked through bind(C). Use when porting any MOM6 kernel to AMReX while keeping the Fortran caller unchanged and the Fortran truth available as a numerical reference. Mirrors the pattern established in TURBO-ESM/MOM6 PR #15.
user-invocable: true
argument-hint: <work-directory> <function-name>
---

# Generate C++ bridge for a MOM6 Fortran subroutine

This skill is the **execution checklist**. All templates, rationale,
type-mapping tables, and pitfalls live in [lessons.md](lessons.md) — read
that file once at the start of every run (Step 0 enforces this) and refer
back to its numbered sections from each step below. Do not reproduce
templates here.

## Help message

If `$ARGUMENTS` is empty, or equals `help`, or equals `--help`, or equals
`-h`, do NOT run any steps. Print the following help message verbatim
and stop:

```
Usage: /generate_cpp_bridge <work-directory> <function-name>

Wrap a MOM6 Fortran subroutine in a runtime-dispatched shim that selects
between the original Fortran code, a capture mode for offline validation,
and a C++/AMReX bridge. Mirrors TURBO-ESM/MOM6 PR #15.

Arguments:
  <work-directory>   Absolute path to an existing TURBO-ESM/MOM6 checkout
                     (must contain src/ and config_src/, and should be on
                     the dev/turbo-debug branch). Cloning is not performed.
  <function-name>    Name of the Fortran subroutine to wrap (case-insensitive
                     match against its declaration in the tree).

Example:
  /generate_cpp_bridge /glade/derecho/scratch/sunjian/MOM6 PPM_limit_pos
```

## Step 0 — validate inputs

`$0` = work-directory, `$1` = function-name. Run these checks **before any
other step** and stop on the first failure with a one-line, actionable
error. Do not retry, do not assume defaults, do not create anything.

1. **Argument count.** If `$0` empty OR `$1` empty → stop:
   `Error: missing arguments. Run "/generate_cpp_bridge --help" for usage.`
2. **Work directory is an existing MOM6 checkout.** If `$0` is not an
   existing directory → stop:
   `Error: work directory "<value>" does not exist.`
   The directory must already contain a TURBO-ESM/MOM6 checkout —
   cloning is not performed by this skill. Step 1 verifies the checkout
   identity and branch state.

The remaining validation (MOM6 layout, subroutine presence, `lessons.md`
present and loaded, plan confirmation) runs in Step 1. No wrapping work
below executes until both Step 0 and Step 1 validation pass.

## Settle these decisions (ask if not obvious from the tree)

1. **Bridge prefix** — default to whatever existing `_bridge) bind(C)`
   declarations in the tree use; otherwise `turbotmp_`.
2. **Env-var name** — default `<UPPERCASE_$1>_MODE`.
3. **AMReX gating** — default `#ifdef _TIM` around the `case (TIMH_runAMREX)`
   arm only; capture mode is never gated.

If the user already specified any of these, take their values as-is.

## Procedure

Each step is one action with a pointer to the lessons.md section that
holds the template or rationale.

### 1. Validate the existing checkout
   **Then validate the tree** before proceeding to Step 2 — stop on the
   first failure with a one-line, actionable error:
   - **Looks like MOM6.** If `$0/src` or `$0/config_src` is missing → stop:
     `Error: "<value>" is not a MOM6 source tree (missing src/ or config_src/).`
   - **Subroutine present.** Run `grep -irn "^[[:space:]]*subroutine[[:space:]]\+<function-name>\b" $0/src $0/config_src`.
     - 0 matches → stop: `Error: subroutine "<function-name>" not found under <work-directory>/{src,config_src}.`
     - >1 match → list candidates and ask the user which file to wrap.
   - **`lessons.md` present.** If `$0/.claude/skills/generate_cpp_bridge/lessons.md` is missing → stop:
     `Error: lessons.md not found at <work-directory>/.claude/skills/generate_cpp_bridge/lessons.md.`
   - **Read `lessons.md` in full.** Treat it as authoritative for naming,
     type mappings, and the dispatcher pattern. If any later step appears
     to conflict with lessons.md, prefer lessons.md and report the
     discrepancy before proceeding.
   - **Confirm the plan.** Print one paragraph naming: the resolved
     subroutine file path, the proposed env-var (`<UPPERCASE_$1>_MODE`),
     the proposed bridge symbol (`<prefix>_$1_bridge`; pick `<prefix>` by
     `grep -h "_bridge) bind(C)" $0/src $0/config_src` — fall back to
     `turbotmp_`), and the caller files found by
     `grep -irl "call[[:space:]]\+<function-name>(" $0/src $0/config_src`.
     Then proceed to Step 2.

### 2. Classify each dummy argument
   For every arg, record intent / kind / rank / optional. Map to bridge
   types using lessons.md §3.1–§3.2. If the original takes loop-bound
   integers, collapse them into a single `Box_t` (lessons.md §5).

### 3. Rename the original implementation
   Rename `$1` → `$1_fortran` in place; convert array dummies to
   `RealArray_t` / `IntArray_t`; rewrite loops as `do concurrent` over
   `bx%idxS`/`bx%idxE`. Template: lessons.md §12.

### 4. Add the `bind(C)` interface block
   At the top of the host module, declare `<prefix>_$1_bridge` per the
   template in lessons.md §3.3. Doc-comment every dummy with `!<`.

### 5. Write the shim subroutine `$1`
   Same public name and dummy list as the original (after array→container
   rewrite). Body is one `select case (mode)` over
   `getenv_mode("<ENV_VAR>", default=TIMH_runFORTRAN)` with three arms:
   `TIMH_capture`, `TIMH_runAMREX` (inside `#ifdef _TIM`), and `case
   default` falling through to `$1_fortran`. Template: lessons.md §2.
   Capture-record naming: lessons.md §13.

### 6. Update host-module `use` statements
   Add the imports listed in lessons.md §14, deduplicating with what the
   module already has.

### 7. Rewrite each caller
   Wrap raw arrays into `RealArray_t` containers, build the iteration
   `Box_t`, call the shim, copy results back, free. Recipe: lessons.md §4.
   For `intent(out)` args, skip the inbound `copy2Array`; for pure-input
   args, skip the outbound `copy2F`. Reuse adjacent-kernel containers if
   the caller already has them.

### 8. Relocate halo checks and CPU clocks
   Move any halo-sufficiency `MOM_error(FATAL,...)` from inside the
   kernel to the caller or the shim entry (lessons.md §15). Leave any
   pre-existing `cpu_clock_begin/end` at the caller, around the shim
   call — never inside the shim (lessons.md §16).

### 9. Verify
   Run the three-mode matrix in lessons.md §17. If only the Fortran
   shim + capture are being delivered, stop after CAPTURE verification
   and report that the C++ side of `<prefix>_$1_bridge` is the next
   deliverable.

### 10. Commit and push (optional — ask first, never create a branch)
   **Never create or switch branches.** All git work stays on whatever
   branch is currently checked out.

   **Ask the user first** before doing any git operations:
   > "Step 10 is optional. Should I stage the modified files, commit, and
   > push to origin on the current branch?"

   Only proceed if the user confirms. If they decline, skip this step
   entirely and report the commit message you would have used so they can
   run the git commands themselves.

   If the user confirms, stage every modified and newly created file and
   commit with a message that briefly summarizes the changes (the wrapped
   subroutine name, the new bridge symbol, and the affected caller files)
   and explicitly notes that the work is co-authored by Claude. Use a
   HEREDOC so the trailer is preserved verbatim:

   ```
   git add -A
   git commit -m "$(cat <<'EOF'
   <one-line summary of the wrapping work for $1>

   <optional 1–3 line body describing the bridge symbol, callers touched,
   and verification mode reached>

   Co-authored-by: Claude <noreply@anthropic.com>
   EOF
   )"
   git push
   ```

   If the push is rejected, stop and surface the conflict to the user
   rather than force-pushing.

## Hard rules

- Do not change the public name or argument list of `$1` (the shim must
  drop into existing call sites unchanged).
- Do not introduce a global mode variable; per-kernel env vars only.
- Do not gate capture mode behind `#ifdef _TIM`; only the AMReX arm is gated.
- Do not omit the `case default` arm.
- Do not pass default-kind `logical`/`integer` to a `bind(C)` routine —
  cast to `c_bool` / `c_int` (lessons.md §3.2).
- Do not call `c_loc(this%data(1))` on an unallocated container; allocate
  first.
- Do not re-implement `getenv_mode`, `already_recorded`, `mark_recorded`,
  or `io_recorder` — `use` them from `turbotmp_helperF`.
- Do not reuse a `kernel` string across kernels.

If something not covered here comes up, consult lessons.md §10
(recurring pitfalls) before improvising.

## Output to the user on success

Report:

1. The shim's env-var name and the three accepted values.
2. The `kernel` string used for capture filenames.
3. The bridge symbol name(s) the C++ side must implement.
4. The list of caller files that were rewritten.
5. What still needs to happen on the C++ side (if anything).
