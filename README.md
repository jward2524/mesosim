# Simulation Input File Format

## 1. Example Input File

```text
systemsize NX NY NZ
temp T
struct FCC|BCC|SC
geometry (sheet N|cluster R|file path)
atomtype A B C ...
composition xA xB xC ...
dissolution true|false ...
nnlevels N
1nne eAA eAB ...
run time|iteration value
flavor KMC|MC
```

## 2. Conceptual Organization

### Geometry

| Command | Required | Description |
|----------|----------|-------------|
| `temp` | Yes | Simulation temperature in K. |
| `potential` | No | Constant or swept electric potential. Default: 0 |
| `nnlevels` | Yes | Number of nearest-neighbor shells. |
| `nne` | Yes | Nearest-neighbor energies for shell n (flattened upper-triangle). For a three-component system: AA AB AC BB BC CC |

### Output and Logging

| Command | Required | Description |
|----------|----------|-------------|
| `seed` | No | Selection of the random seed value. Default: default |
| `run` | Yes | Simulation end condition. |
| `flavor` | Yes | Simulation algorithm flavor. |

### Thermodynamics

| Command | Required | Description |
|----------|----------|-------------|
| `datalog` | No | Logging schedule configuration. |
| `output` | No | Output log filename. |
| `logtype` | No | Enable output formats. |
| `help` | No | Show documentation for commands. |

## 3. Detailed Command Reference

### `systemsize`

**Category:** Geometry  
**Required:** Yes

**Syntax**

```text
systemsize NX NY NZ
```

**Parameters**

- `NX` — required
- `NY` — required
- `NZ` — required

**Description**

Simulation box size in lattice units.

---

### `temp`

**Category:** Thermodynamics  
**Required:** Yes

**Syntax**

```text
temp T
```

**Parameters**

- `T` — required

**Description**

Simulation temperature in K.

---

### `seed`

**Category:** Run Control  
**Required:** No

**Syntax**

```text
seed random|default|N
```

**Parameters**

- `random|default|N` — required

**Description**

Selection of the random seed value. Default: default

---

### `potential`

**Category:** Thermodynamics  
**Required:** No

**Syntax**

```text
potential U0 [dUdt Umax]
```

**Parameters**

- `U0` — required
- `dUdt` — optional
- `Umax` — optional

**Description**

Constant or swept electric potential. Default: 0

---

### `datalog`

**Category:** Output and Logging  
**Required:** No

**Syntax**

```text
datalog (linear|ln|iteration) (interval a [b]|list ...)
```

**Parameters**

- `(linear|ln|iteration)` — required
- `(interval` — required
- `a` — required
- `b|list` — optional
- `...)` — required

**Description**

Logging schedule configuration.

---

### `struct`

**Category:** Geometry  
**Required:** Yes

**Syntax**

```text
struct FCC|BCC|SC
```

**Parameters**

- `FCC|BCC|SC` — required

**Description**

Crystal structure type.

---

### `output`

**Category:** Output and Logging  
**Required:** No

**Syntax**

```text
output path/outfile.out
```

**Parameters**

- `path/outfile.out` — required

**Description**

Output log filename.

---

### `geometry`

**Category:** Geometry  
**Required:** Yes

**Syntax**

```text
geometry (sheet N|cluster R|file path)
```

**Parameters**

- `(sheet` — required
- `N|cluster` — required
- `R|file` — required
- `path)` — required

**Description**

Initial geometry configuration.

---

### `atomtype`

**Category:** Geometry  
**Required:** Yes

**Syntax**

```text
atomtype A B [C ...]
```

**Parameters**

- `A` — required
- `B` — required
- `C` — optional
- `...` — optional

**Description**

Define atom types and their order.

---

### `composition`

**Category:** Geometry  
**Required:** Yes

**Syntax**

```text
composition xA xB [xC ...]
```

**Parameters**

- `xA` — required
- `xB` — required
- `xC` — optional
- `...` — optional

**Description**

Atomic composition fractions; order follows atomtype.

---

### `dissolution`

**Category:** Geometry  
**Required:** Yes

**Syntax**

```text
dissolution true|false ...
```

**Parameters**

- `true|false` — required
- `...` — required

**Description**

Dissolution flags per atom type; order follows atomtype.

---

### `nnlevels`

**Category:** Thermodynamics  
**Required:** Yes

**Syntax**

```text
nnlevels N
```

**Parameters**

- `N` — required

**Description**

Number of nearest-neighbor shells.

---

### `nne`

**Category:** Thermodynamics  
**Required:** Yes

**Syntax**

```text
1nne eAA eAB ...
```

**Parameters**

- `eAA` — required
- `eAB` — required
- `...` — required

**Description**

Nearest-neighbor energies for shell n (flattened upper-triangle). For a three-component system: AA AB AC BB BC CC

---

### `run`

**Category:** Run Control  
**Required:** Yes

**Syntax**

```text
run time|iteration value
```

**Parameters**

- `time|iteration` — required
- `value` — required

**Description**

Simulation end condition.

---

### `flavor`

**Category:** Run Control  
**Required:** Yes

**Syntax**

```text
flavor KMC|MC
```

**Parameters**

- `KMC|MC` — required

**Description**

Simulation algorithm flavor.

---

### `logtype`

**Category:** Output and Logging  
**Required:** No

**Syntax**

```text
logtype [iter] [csv] [xyz]
```

**Parameters**

- `iter` — optional
- `csv` — optional
- `xyz` — optional

**Description**

Enable output formats.

---

### `help`

**Category:** Output and Logging  
**Required:** No

**Syntax**

```text
help [command]
```

**Parameters**

- `command` — optional

**Description**

Show documentation for commands.

---

