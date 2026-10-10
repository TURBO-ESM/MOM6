---
name: convert_loops_to_box_iterators
version: "0.3.2"
description: Convert the scalar-range `do concurrent` loop headers in a MOM6 subroutine (or a whole call tree) to iterate over `Box_t` iteration boxes -- a small set of core boxes (one per staggering, built per tile by the tree root) plus derived boxes grown from them by at most a cell or two -- instead of scalar bounds such as `Isq:Ieq+1`. Every range is first reduced to core-box-plus-offset form and proven equal in both memory modes; loop bodies, index variables and their order never change; now-unused scalar bounds and box dummies are dropped. Use it after a tree has box/tile iteration (the tree root builds `bxH`/`bxQ`, etc. per tile and kernels read scalar bounds from them) and before Phase 3 bridging, so every loop maps one-to-one onto an AMReX `ParallelFor(box, ...)`. Distinct from convert_array_containers, which only changes argument types and leaves loop-bound expressions alone.
user-invocable: true
argument-hint: <work-directory> <function-name>[,<function-name>...] [--core <box>[,<box>...]] [--enable_git_commit] [--disable_git_commit]
---

# Convert scalar-range loops to box iterators

## Why this exists

A kernel that reads scalar bounds from its boxes once at the top
(`Isq = bxQ%idxS(1) ; …`, `array_container_lessons` §5) keeps its
original loop headers — `do concurrent (k=ksc:kec, j=Jsq:Jeq+1,
i=Isq:Ieq+1)` — which is the right low-risk choice while converting
argument types. But it hides the structure the AMReX port needs. There,
each loop is a `ParallelFor` over one `amrex::Box`, and the bridge
builds that box from a `Box_C`
(TIM `generate_amrex_code` lessons.md §3). A loop whose range is only a
set of scalar expressions has no box to hand over.

The ranges in a tree look numerous but are not. In CorAdCalc, 15
different spellings reduced to 11 distinct ranges, all of them one of
**5 core boxes** grown by at most one cell:

| Core | Staggering | Range |
|---|---|---|
| `bxH` | h, h | `[isc:iec, jsc:jec]` |
| `bxU` | B, h | `[IscB:IecB, jsc:jec]` |
| `bxV` | h, B | `[isc:iec, JscB:JecB]` |
| `bxQ` | B, B | `[IscB:IecB, JscB:JecB]` |
| `bxQs` | scheme-dependent vorticity range | `[Is_q:Ie_q, Js_q:Je_q]` |

This skill makes that structure explicit: each loop iterates over a
core box, or over a derived box built once from a core with
`grow`/`growLo`/`growHi`.

Natural place in a `convert_calltree` campaign: after the Stage 1
box/tile conversion, once the tree root builds per-tile boxes, and
before Phase 3. CorAdCalc ran it as steps 9c (U/V cores), 9d (grown
U/V) and 9e (Q, Qs, H), each its own bit-for-bit checkpoint. Splitting
by core box like that keeps each diff reviewable.

Assumes `array_container_lessons` (§5 `Box_t` API, memory-mode rule;
§9 #3 owned boxes; §9 #18 case-insensitive names) has already been
invoked this session.

## Scope

**In scope:** `do concurrent` headers whose every range (i, j and, when
present, k) is an expression in scalar bounds that trace back to a box
the subroutine has, or can be given.

**Out of scope, leave the header alone:**
- serial `do` loops (e.g. the OBC `do k=ksc,kec` / `do i=max(Is_q,
  OBC%segment(n)%HI%isd), …` loops): their ranges are clipped by
  non-box data;
- any range that depends on something other than a box plus a constant
  offset (segment bounds, a runtime `stencil` that is not already
  folded into a core box, a `min`/`max`);
- loop bodies, index variable names, their case and their order;
- `DO_LOCALITY(...)` and `!$omp` directives: carried over verbatim.

## The core technique

1. **Resolve every scalar bound to its box.** Use the extraction lines
   (`Isq = bxQ%idxS(1)`), the tree root's `%set` calls and the `G`
   fields behind them. Write each range as `core ± offset` per
   dimension, e.g. `j=Jsq:Jeq+1, I=Isq:Ieq` → `bxQ`, dim 2 `+1` at the
   end.

2. **Prove it in both memory modes.** `G%IscB`/`G%JscB` are `isc`/`jsc`
   with non-symmetric memory and `isc-1`/`jsc-1` with symmetric memory,
   but `IecB = iec` and `JecB = jec` always. So:
   - an **end** offset may cross staggerings. `Ieq` is `iec`, so the
     h-box range `is-1:Ieq` is `bxH%growLo(1,1)`;
   - a **start** offset may not. Never grow `bxH` to get `IscB`, or
     shrink `bxQ` to get `isc`, because that holds in one mode only.
     Use the core box of the right staggering, or add one (step 3).

3. **Add any missing core box in the tree root.** Allocate it once
   beside the others, `%set` it per tile from `G`, and free it after
   the tile loop (lessons §5 "Per-tile boxes"). Pass it to each kernel
   that needs it.
   - **Boxes come first** in every kernel's argument list, in the
     existing relative order `bxH, bxQ, bxU, bxV, bxQs` (the continuity
     `bxC`-first precedent). Update each call site.
   - Declare it with a `!<` doc comment giving its range, in the style
     of the neighbouring box dummies.

4. **Build each derived box once per subroutine, never per loop.**
   - Declare it as a local, after `! Local variables`, with a comment
     giving its range:
     `type(Box_t) :: bxQ_pij  ! bxQ grown by one point at the end in i and j,`.
   - Build it right after the scalar-bound extraction lines, with one
     `growBy` call per box: `lo`/`hi` give the amount to extend the start
     and the end in each dimension. No temporary box is needed.
     ```fortran
     bxQ_pij = bxQ%growBy(lo=[0,0,0], hi=[1,1,0])
     bxQ_pj  = bxQ%growHi(dim=2, n=1)
     bxUx    = bxU%growBy(lo=[0,1,0], hi=[1,0,0])   ! end in i, start in j
     ```
   - Free every derived box just before `end subroutine`. Check for an
     early `return` first: each one needs the same frees, or the box
     leaks (lessons §9 #3).
   - **Naming:** `<core>_<p|m|g><dims>`. `p` means grown at the end
     (`growHi`), `m` at the start (`growLo`), `g` on both sides
     (`grow`), e.g. `bxQ_pij`, `bxQs_pi`, `bxH_mij`, `bxQ_gij`. When the
     dimensions differ, use one group per dimension, e.g. `bxU_pi_mj`.
     (CorAdCalc's step 9d `bxUx`/`bxVx` predate this convention.) Use
     `n=2`… and a digit (`bxQ_p2i`) only when a range really needs it.

5. **Rewrite the header only.** Keep the index variables, their case
   and order, and whether a k index is present:
   ```fortran
   do concurrent (k=bxQ_pij%idxS(3):bxQ_pij%idxE(3), j=bxQ_pij%idxS(2):bxQ_pij%idxE(2), &
                  I=bxQ_pij%idxS(1):bxQ_pij%idxE(1)) DO_LOCALITY(local(hArea_q))
   ```
   If the second line would exceed 100 columns with the `DO_LOCALITY`
   tail, move the tail to a third line, indented 4 past `do`, after
   `) &`.
   - A **2-D loop** (no `k`) uses dims 1–2 of a 3-D box.
   - A **k-invariant loop before the tile loop** (e.g. CorAdCalc_TR's
     `Area_h`/`Area_q`) needs its box before the per-tile `%set`.
     Allocate and set that core box there, over `1:nz`, derive from it,
     and let the tile loop re-`set` it as before. Move its `safe_alloc`
     up rather than calling it twice.

6. **Drop what became unused.** Recompute usage over the whole
   subroutine body, ignoring the extraction and declaration lines:
   - remove each scalar bound nothing reads any more, from both its
     extraction line and its `integer ::` declaration. Bounds still
     read by out-of-scope loops stay;
   - remove each box dummy nothing references any more. A box used only
     in an automatic array's specification expression (e.g.
     `bxH%idxS(3):bxH%idxE(3)` for tile-shaped scratch) **is**
     referenced. Update the call site too.

## Steps

### 0. Validate inputs

`$0` = work-directory, `$1` = comma-separated function-name list.
`--core` optionally limits the run to loops whose core is one of the
named boxes (e.g. `--core bxU,bxV`), for staging as CorAdCalc did.

1. Help/empty-argument/flag handling identical to
   `convert_array_containers` Step 0 items 1–3.
2. Invoke `array_container_lessons` now if not already done this
   session.
3. Locate each subroutine (0 or >1 matches: stop). Every target must
   already have box iteration: box dummies, or a tile loop that builds
   them. A target with raw `G%isc`-style bounds and no box needs the
   Stage 1 box conversion first; stop and say so.

### 1. Inventory and reduce

List every `do concurrent` header in the targets with its subroutine,
line, raw range text and the range reduced to core ± offset (technique
steps 1–2). Report:
- the distinct reduced ranges and the core boxes they need;
- the count per subroutine;
- anything out of scope, and why.

**Show this table before editing** — it is the plan for the edit. If a
range reduces to no core box, or a new core box is needed that the
plan doc doesn't already record, ask.

### 2. Name-clash check

Check every new name (derived boxes, any new core box) against
every existing symbol in each target and in the module, ignoring case
(lessons §9 #18). On a clash, stop and ask.

### 3. Edit

In this order, one subroutine at a time:
1. add any core box to the tree root and to the kernel arguments
   (technique step 3);
2. add the derived-box declarations, builds and frees (step 4);
3. rewrite the headers (step 5);
4. drop unused bounds and dummies (step 6).

A small script is the reliable way to do this for more than a handful
of loops. Match each header against the exact range text from Step 1,
never a looser regex.

### 4. Verify

- **Canonicalised diff.** Join each new header back to one line,
  replace each `bx%idxS(d):bx%idxE(d)` and each old scalar range with
  a placeholder, and diff against the pre-edit file. The only
  remaining differences must be box declarations, builds, frees,
  dropped bounds, and argument-list and call-site changes. Zero body
  lines.
- **Range table.** Re-derive every converted loop's range from its box
  build and compare it with the Step 1 raw range, in **both** memory
  modes.
- **No leftovers.** No in-scope `do concurrent` still has a scalar
  range; every derived box is freed on every exit path.
- **Consistency.** Every call's argument count and order match its
  callee. Every dropped bound or dummy is gone from the argument list,
  the declaration, the extraction line and the call site. Grep each old
  name across the subroutine (repo CLAUDE.md "Verification
  discipline").
- **Declarations.** Every declaration precedes the first executable
  statement (new builds go after the extraction lines, never among the
  declarations).
- **Line length.** Every added line is ≤ 100 columns.
- **Build and test belong to the user** unless they say otherwise.
  Suggest double_gyre bit-for-bit under both infra layers, plus
  override runs that reach each converted loop: every scheme `select
  case` branch, the `KE_SCHEME` variants, and a scheme that changes the
  `bxQs` stencil. Then CI.

## Versioning marker

Every Fortran file this skill creates or modifies gets a `!!SKILLS: 0.3.2`
marker line — the shared version for this whole skill family. If
missing, add it right after the file's license/header block, before
`module`; if present, update it in place. Grep-able
(`grep -rn "!!SKILLS:"`), meant to be stripped later.

## Hard rules

- Never skip or duplicate the `!!SKILLS: 0.3.2` marker.
- Never change a loop body, an index variable's name or case, the index
  order, or a `DO_LOCALITY`/`!$omp` directive.
- Never derive a box whose **start** index crosses staggerings
  (`bxH` → `IscB`, `bxQ` → `isc`): it is bit-for-bit in one memory
  mode only. Build the right core box from `G` instead.
- Never build a derived box inside a loop, or per loop; build it once
  per subroutine, and free it on every exit path.
- Never convert a serial `do` loop or a range clipped by non-box data.
- Never leave an unused box dummy or scalar bound behind, and never
  drop one that an automatic array's bounds or an out-of-scope loop
  still reads.
- Box dummies always lead the argument list, in `bxH, bxQ, bxU, bxV,
  bxQs` relative order.
- Do not add comments about the conversion itself — same standing rule
  as the rest of this family.
- Do not attempt to install a compiler, or claim an unrun build passed.

## Commit gating

Same as the rest of the family: `--enable_git_commit`/
`--disable_git_commit` override; otherwise
`~/.claude/preferences.json`'s `git_commit_and_push` key. Branch:
`claude_<lowercased_first_function_name>_box_iterators`, or
`..._and_N_more_box_iterators` for a multi-name run.

## Output to the user on success

1. The Step 1 table: loops converted per subroutine, the core and
   derived boxes, and each derived box's build and range.
2. Out-of-scope loops left alone, and why.
3. Core boxes added, argument lists changed, and bounds and dummies
   dropped.
4. Verification results: canonicalised diff, range table in both
   memory modes, consistency, line length.
5. Build status, and the suggested override runs.
6. Whether committed, or the modified files for manual commit.
