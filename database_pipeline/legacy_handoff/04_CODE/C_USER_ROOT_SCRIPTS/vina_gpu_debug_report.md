# Vina-GPU parse_pdbqt.cpp:71 Assertion Crash — Debug Report

## Symptom
Vina-GPU crashes on **every docking attempt** during "Reading input ...":
```
Vina-GPU: ./lib/parse_pdbqt.cpp:71: std::string omit_whitespace(const std::string&, sz, sz): Assertion `i-1 < str.size()' failed.
Aborted (core dumped)
```

`--randomize_only` mode **works fine** (0.087s). Only docking mode crashes.

## Environment
- Server: A100-80GB, CUDA, Ubuntu
- Binary: `~/Vina-GPU/Vina-GPU`
- Compiled with: `-DBUILD_KERNEL_FROM_SOURCE -DBOOST_TIMER_ENABLE_DEPRECATED -DOPENCL_3_0 -DNVIDIA_PLATFORM`
- Boost: 1.85 (via conda env `docking`)
- Binary runs from `~/Vina-GPU/` directory (OpenCL kernel files at `./OpenCL/src/kernels/`)

## Source Code Analysis

### The crashing function (`parse_pdbqt.cpp:58-75`)
```cpp
std::string omit_whitespace(const std::string& str, sz i, sz j) {
    if(j < i-1) j = i-1; // i >= 1                    // line 60
    if(j < str.size()) j = str.size();                  // line 61
    
    while(i <= j && std::isspace(str[i-1])) ++i;        // line 64-65
    while(i <= j && std::isspace(str[j-1])) --j;        // line 68-69
    
    VINA_CHECK(i-1 < str.size());  // <--- LINE 71: CRASH HERE
    VINA_CHECK(j-i+1 < str.size());                     // line 72
    
    return str.substr(i-1, j-i+1);
}
```

### The only call site (`parse_pdbqt.cpp:109`)
```cpp
parsed_atom parse_pdbqt_atom_string(const std::string& str) {
    unsigned number = checked_convert_substring<unsigned>(str, 7, 11, "atom number");
    vec coords(checked_convert_substring<fl>(str, 31, 38, "coordinate"),
               checked_convert_substring<fl>(str, 39, 46, "coordinate"),
               checked_convert_substring<fl>(str, 47, 54, "coordinate"));
    fl charge = 0;
    if(!substring_is_blank(str, 69, 76))
        charge = checked_convert_substring<fl>(str, 69, 76, "charge");
    std::string name = omit_whitespace(str, 78, 79);  // <--- LINE 109: only call
    ...
}
```

### Key finding: `omit_whitespace` is called with `i=78, j=79`, requiring `str.size() >= 78`
For the assertion `VINA_CHECK(i-1 < str.size())` → `VINA_CHECK(77 < str.size())`, the string must be at least 78 characters long.

### Parser gating (how lines reach `parse_pdbqt_atom_string`)
All 3 call sites are gated by:
```cpp
else if(starts_with(str, "ATOM  ") || starts_with(str, "HETATM")) {
    parsed_atom a = parse_pdbqt_atom_string(str);
```

## What We've Checked

### 1. Receptor file (2gao.pdbqt) — CLEAN
- 3506 ATOM/HETATM lines, **ALL exactly 79 characters**
- No ROOT/BRANCH/ENDBRANCH/TORSDOF residue (grep confirmed blank)
- No empty lines, no non-ASCII characters
- Two duplicate `END` lines were present, **both removed** — crash persists

### 2. Ligand file (HMPD-CMPD-0175063.pdbqt) — NORMAL
- 77 lines total
- 36 ATOM lines, **ALL exactly 79 characters**
- Lines 1-16: REMARK (11 torsion descriptions + blank header) — properly handled by parser
- Line 17: ROOT
- Lines 18-22: ATOM (5 atoms in root)
- Line 23: ENDROOT
- Lines 24-76: BRANCH/ENDBRANCH (torsion tree)
- Line 77: TORSDOF 11

### 3. Receptor parser (`parse_pdbqt_rigid`, lines 258-283)
- Handles: empty, TER, WARNING, REMARK, ATOM/HETATM
- Does NOT handle: `END` keyword → would throw "Unknown or inappropriate tag", not assert
- Both END lines removed from receptor, crash persists

### 4. Ligand parser (`parse_pdbqt_root` → `parse_pdbqt_root_aux` → `parse_pdbqt_aux`)
- REMARK before ROOT: skipped (line 319)
- ROOT section: processes ATOM lines until ENDROOT
- BRANCH/ENDBRANCH/TORSDOF: all handled

### 5. Tested configurations
- Original config: crash
- Minimal config (no energy_range): crash
- Multiple receptors tested: all crash
- `--randomize_only`: WORKS

## The Puzzle

Every ATOM/HETATM line in both files is exactly **79 characters**, which should pass `VINA_CHECK(77 < 79)`. No line shorter than 78 chars exists in either file.

Yet the assertion fires. This means either:
1. **A non-ATOM line is somehow reaching `parse_pdbqt_atom_string`** — but all call sites have the `starts_with(str, "ATOM  ")` gate
2. **The compiled binary doesn't match the source** — maybe an older compiled version with a different bug
3. **There's a different code path** we haven't examined (e.g., flex residue parsing, grid map reading)

## Unchecked Leads

1. **Example files exist** at `~/Vina-GPU/input_file_example/2bm2_protein.pdbqt` and `2bm2_ligand.pdbqt` — do these dock successfully? (would prove whether binary itself is broken)

2. **`parse_pdbqt_atom_string` line 543 call site** — inside `parse_pdbqt_branch_aux`, also gated by `starts_with(str, "ATOM  ")`. Same protection.

3. **Recompile with `-g` and debug with GDB** (`apt install gdb`) to get exact stack trace

4. **Source/binary mismatch** — verify the binary was compiled from the current source files

## Recommended Next Steps (in order)
1. Test docking with the bundled `2bm2_protein.pdbqt` + `2bm2_ligand.pdbqt` example files
2. Install gdb and get a stack trace: `gdb -batch -ex run -ex bt ./Vina-GPU --config ...`
3. Confirm binary matches source: check compile timestamp vs source file modification times
4. If example files also crash, recompile clean and retry
