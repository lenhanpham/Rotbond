<p align="center">
  <a href="https://github.com/lenhanpham/Rotbond">
    <picture>
      <img src="resources/rotbond-logo.svg" alt="Rotbond" style="width: 50%;">
    </picture>
  </a>
</p>
<p align="center">
  <a href="https://www.rust-lang.org/">
    <img src="https://img.shields.io/badge/Rust-1.70+-orange.svg" alt="Rust">
  </a>
  <a href="./LICENSE">
    <img src="https://img.shields.io/badge/License-MIT-green.svg" alt="License">
  </a>
  <a href="https://github.com/lenhanpham/Rotbond/actions/workflows/build.yml">
    <img src="https://github.com/lenhanpham/Rotbond/actions/workflows/build.yml/badge.svg" alt="Build Status">
  </a>
</p>

Rotbond is program for generating molecular conformers through systematic bond rotations. Rotbond takes an XYZ molecular structure file and rotation parameters, then systematically generates all possible conformers, outputting them as individual XYZ files and a trajectory file for further calculations and visualization.

## Features

### Core Capabilities

- **Flexible Rotation Specifications**
  
  - Step-based rotations - e.g., `e60` generates 0°, 60°, 120°, 180°, 240°, 300° (6 states); `e30` generates 12 states
  - Explicit angle lists - specify exact angles needed
  - Synchronous rotations - link multiple bonds to rotate together
  - Manual bond definitions - force or remove bonds as needed

- **Smart Output Generation**
  
  - Trajectory file with all conformers for easy visualization
  - Individual XYZ files with automatic smart padding
  - Progress reporting and statistics
 
- **Intelligent Validation**
  
  - Steric clash detection and filtering
  - Configurable bond detection parameters
  - Distance-based validation using covalent radii

- **Documentation**
  
  - Built-in help system with 8 comprehensive sections
  - Practical examples and troubleshooting guides
  - Complete reference materials

### Advanced Features

- Support for synchronous rotations (same or opposite direction)
- Manual bond forcing for special cases (metal complexes, coordinate bonds)
- False bond removal for accurate molecular graphs
- **Automatic template creation** - generates comprehensive .rp templates with examples
- **Flexible input formats** - accepts both `molecule` and `molecule.xyz` formats
- **Conformer generation limits** - prevent system crashes from large jobs
- **Interactive safety warnings** - alerts for memory-intensive generations
- **Auto-confirmation mode** - batch processing without manual intervention
- Comprehensive error messages with contextual help
- Modern Rust implementation for safety and performance

## Quick Start

### Installation

```bash
# Clone the repository
git clone https://github.com/lenhanpham/Rotbond.git
cd Rotbond

# Build in release mode for optimal performance
cargo build --release

# Or build in debug mode for development
cargo build
```

**A binary file is provided for windows users in Release**

Copy rotbond.exe to a directory and add this directory to your Windows *environment variable*. 

### Basic Usage

Export your the directory where the rotbond binary is located, and then use it.

```bash
export PATH=$PATH:your_dir/target/release

# Run with release build (accepts both formats)
rotbond <molecule_name>
rotbond <molecule_name>.xyz
```

### Example: Ethane Rotation

**Step 1: Prepare your stucture `ethane.xyz`**

```
8
Ethane molecule
C       -2.6972703886      0.5190638733     -0.0000003688                 
C       -1.1782073026      0.5190641624     -0.0000004138                 
H       -3.0738517143     -0.1361520096     -0.8135565479                 
H       -3.0738520178      1.5512320348     -0.1606558009                 
H       -3.0738517177      0.1421113795      0.9742112740                 
H       -0.8016259769      0.8960166481     -0.9742120609                 
H       -0.8016256734     -0.5131039980      0.1606550254                 
H       -0.8016259735      1.1742800518      0.8135557588         
```

**Step 2: Run Rotbond to create template**

```bash
rotbond ethane
```

**Output:**

```
Rotation parameters file 'ethane.rp' not found.
Creating template file with default parameters...
✓ Template file 'ethane.rp' created successfully!

Please edit 'ethane.rp' to specify your rotation parameters, then run rotbond again.
```

**Step 3: Edit the generated `ethane.rp` template**

The program creates a comprehensive template with examples. Simply uncomment and modify:

```
# =======================================================================
# ROTBOND ROTATION PARAMETERS FILE
# =======================================================================
# This is a template file with examples of all available features.
# Edit this file to specify your rotation parameters, then run rotbond again.

# Configuration parameters
bond_factor = 1.2
skip_factor = 0.7

# Add your rotation specifications (uncomment and modify):
1-2 e60        # Rotate C-C bond every 60°
```

**Step 4: Run Rotbond again**

```bash
rotbond ethane
```

**Step 5: Output Files**

- `ethane_traj.xyz` - Trajectory with all conformers
- `ethane_01.xyz`, `ethane_02.xyz`, etc. - Individual conformers

## Automatic Template Creation

Rotbond now automatically creates comprehensive template files to make getting started easier than ever!

### How It Works

1. **Run with just an XYZ file**: `rotbond molecule` or `rotbond molecule.xyz`
2. **Automatic detection**: If `molecule.rp` doesn't exist, Rotbond creates it automatically
3. **Comprehensive template**: Generated file includes examples of all features with detailed comments
4. **Edit and run**: Simply uncomment and modify the examples, then run Rotbond again

### Template Features

The automatically generated template includes:

- **Configuration parameters** with explanations and recommended values
- **Manual bond definitions** for special cases (metal complexes, etc.)
- **All rotation types** with examples:
  - Step-based rotations (`1-2 e60`)
  - Explicit angle lists (`1-2 0 60 120 180`)
  - Synchronous rotations (`3-4 syn 1`)
- **Real-world examples** for common molecules (ethane, butane, toluene, metal complexes)
- **Best practices** and tips for parameter selection
- **Detailed comments** explaining each feature

### Example Template Output

```bash
rotbond ethane

# Output:
# Rotation parameters file 'ethane.rp' not found.
# Creating template file with default parameters...
# ✓ Template file 'ethane.rp' created successfully!
# 
# Please edit 'ethane.rp' to specify your rotation parameters, then run rotbond again.
```

The generated template is ready to use - just uncomment the rotation specifications you need!

## Required Input Files

### 1. Molecular Structure File (`<name>.xyz`)

Standard XYZ format with atom coordinates:

```
n_atoms
comment_line
element x y z
element x y z
...
```

**Specifications:**

- Atom count: Integer on first line
- Comment: Any text on second line
- Coordinates: Element symbol followed by X, Y, Z in Angstroms
- Atom indices are 1-based in all input files
- Elements are case-insensitive
- Standard covalent radii used for bond detection

### 2. Rotation Parameters File (`<name>.rp`)

**Automatic Template Creation:** If this file doesn't exist, Rotbond will automatically create a comprehensive template with examples and documentation. Simply run `rotbond <molecule_name>` and edit the generated template.

Configuration and rotation specifications in **any order**:

```
# Configuration (optional - defaults provided)
bond_factor = 1.2    # Bond detection threshold multiplier
skip_factor = 0.7    # Validation threshold
maxgen = 500         # Maximum conformers to generate (or "max" for unlimited)
autoconfirm = false  # Skip interactive warnings for large jobs

# Manual bond definitions (optional)
8-12 bond            # Force bond between atoms 8 and 12
22-45 nobond         # Remove bond between atoms 22 and 45

# Rotation specifications
1-6 e60              # Step-based rotation
2-5 0 60 120 180     # Explicit angles
3-7 syn 1            # Synchronous with bond 1 (same direction)
4-8 syn -1           # Synchronous with bond 1 (opposite direction)
```

**Flexible Format Rules:**

- Comments start with `#` and continue to end of line
- Empty lines are ignored
- Format: `atom1-atom2` (use dash, not hyphen)
- All angles in degrees
- Angle range: -360° to +360°

**Example Flexible Files:**

```bash
# Everything mixed together
7-11 e30
bond_factor = 1.2
2-7 0 120
21-22 bond
skip_factor = 0.7

# Or in any other order
bond_factor = 1.3
2-5 0 90 180
7-11 e45
skip_factor = 0.8
```

All formats above work identically!

## Rotation Syntax Reference

### Step-Based Rotation

Generates evenly spaced angles from 0° to 360°-step:

```bash
atom1-atom2 e<step_angle>
```

**Examples:**

- `1-6 e60` → Generates: 0°, 60°, 120°, 180°, 240°, 300° (6 states)
- `2-5 e30` → Generates: 0°, 30°, 60°, 90°, 120°, 150°, 180°, 210°, 240°, 270°, 300°, 330° (12 states)
- `3-7 e90` → Generates: 0°, 90°, 180°, 270° (4 states)

**Note:** The number after 'e' can be any angle (e.g., e30, e45, e60, e90, e120).
The number of states is calculated as 360°/step. Smaller step angles generate more rotation states; larger step angles generate fewer.

**Use cases:**

- Systematic exploration of conformational space
- When you want evenly spaced angle coverage
- Default choice for most rotations

### Explicit Angle Lists

Specify exact angles to use:

```bash
atom1-atom2 angle1 angle2 angle3 ...
```

**Examples:**

- `1-6 0 60 120 180` → Uses exactly these 4 angles
- `2-5 30 45 60 90 120 150` → Uses 6 specific angles
- `3-7 -120 -60 0 60 120` → Includes negative angles

**Use cases:**

- Target specific conformations
- Reduce computational load
- Define symmetry-unique angles

### Synchronous Rotations

Link multiple bonds to rotate together:

```bash
# Same direction as reference
atom1-atom2 syn <reference_bond>

# Opposite direction to reference
atom1-atom2 syn -<reference_bond>
```

**Examples:**

- `1-2 e60` - Main independent bond
- `3-4 syn 1` - Bond 3 rotates with bond 1 (same angle)
- `5-6 syn -1` - Bond 5 rotates opposite to bond 1

**Use cases:**

- Methyl groups in aromatic compounds
- Symmetric molecular fragments
- Coupled conformational changes

**Important Notes:**

- Synchronous bonds reference bond numbers (1, 2, 3, etc.)
- Cannot reference themselves
- Cannot create circular references
- Use negative sign for opposite direction

## Configuration Parameters

### bond_factor

Controls bond detection sensitivity:

- **Default:** 1.0
- **Formula:** `threshold = (cov_radius₁ + cov_radius₂) × bond_factor + 0.45`
- **Lower values** (0.5-0.8): Stricter detection, fewer bonds
- **Higher values** (1.2-1.5): Looser detection, more bonds
- **Recommended ranges:**
  - Organic molecules: 1.0-1.2
  - Metal complexes: 1.2-1.5
  - Crystal structures: 0.8-1.0

### skip_factor

Controls steric clash validation:

- **Default:** 0.7
- **Formula:** `min_distance = (cov_radius₁ + cov_radius₂) × skip_factor`
- **Lower values** (0.4-0.6): More lenient, accept closer atoms
- **Higher values** (0.8-1.0): Stricter validation
- **Recommended ranges:**
  - Small rigid molecules: 0.8-1.0
  - Flexible molecules: 0.6-0.8
  - When in doubt: Start with 0.7

### maxgen

Controls the maximum number of conformers to generate:

- **Default:** 500
- **Values:** 
  - Integer (e.g., `maxgen = 300`): Limit to specific number
  - `max`, `maximum`, or `full`: Generate all theoretical conformers (no limit)
- **Safety:** Prevents system crashes from excessive memory usage
- **Interactive warning:** Displays warning for >500 theoretical conformers
- **Recommended ranges:**
  - Small molecules: 100-1000
  - Large molecules: 50-500
  - High-memory systems: 1000+
  - Unlimited generation: Use with caution

### autoconfirm

Controls automatic confirmation for large conformer generation jobs:

- **Default:** false
- **Values:** 
  - `true`, `yes`, `1`, `on`: Skip interactive warnings
  - `false`, `no`, `0`, `off`: Show interactive warnings (default)
- **Use cases:**
  - Batch processing without manual intervention
  - Automated workflows and scripts
  - When you're confident about system resources
- **Safety:** Only use when you understand the memory implications

## Command-Line Options

```bash
rotbond [OPTIONS] <molecule_name>
rotbond [OPTIONS] <molecule_name>.xyz
```

**Input Formats:**

- `rotbond molecule` - Uses `molecule.xyz` and `molecule.rp` (creates template if .rp missing)
- `rotbond molecule.xyz` - Same behavior, accepts .xyz extension

**Options:**

- `-h, --help` - Show basic help message
- `-v, --version` - Show version information
- `--help topics` - List all available help topics
- `--help <topic>` - Show specific help topic
- `--verbose` - Enable verbose output

**Help Topics:**

- `usage` - Command-line usage and options
- `features` - Feature overview
- `input` - Input file formats (XYZ, .rp)
- `output` - Output file formats
- `examples` - Practical usage examples
- `algorithms` - Algorithm details
- `troubleshoot` - Troubleshooting guide
- `reference` - Reference materials

## Output Files

### Trajectory File (`<name>_traj.xyz`)

Contains all conformers in a single file:

```
8
Conformer 1 of 56 - generated by Rotbond
C  0.000000  0.000000  0.000000
H  0.000000  0.000000  1.090000
...
8
Conformer 2 of 56 - generated by Rotbond
C  0.123456  0.234567  0.345678
H  0.234567  0.345678  0.456789
...
```

**Features:**

- Standard XYZ format
- Compatible with VMD, Chimera, PyMOL
- Each structure labeled with conformer number
- Easy for visualization and analysis

### Individual Conformer Files (`<name>_NN.xyz`)

One file per valid conformer with smart padding:

- 1-9 conformers: `molecule_1.xyz`, `molecule_2.xyz`, ..., `molecule_9.xyz`
- 10-99 conformers: `molecule_01.xyz`, `molecule_02.xyz`, ..., `molecule_99.xyz`
- 100-999 conformers: `molecule_001.xyz`, ..., `molecule_999.xyz`
- 1000+ conformers: `molecule_0001.xyz`, ..., `molecule_1234.xyz`

**Features:**

- Single XYZ structure per file
- 6-decimal precision
- Works with any molecular viewer
- Each file is self-contained

## Practical Examples

### Example 1: Simple Single Bond (Ethane)

**Input:** `ethane.rp`

```
bond_factor = 1.2
skip_factor = 0.7
1-2 e60
```

**Result:** 6 rotation states → 6 valid conformers (100% success rate)

### Example 2: Multiple Independent Bonds (Butane)

**Input:** `butane.rp`

```
bond_factor = 1.2
skip_factor = 0.7

# Rotate around C1-C2 bond
1-2 e120

# Rotate around C2-C3 bond
2-3 e60

# Rotate around C3-C4 bond
3-4 e180
```

**Result:** 3 bonds × [3, 6, 2] states = 36 total combinations

### Example 3: Synchronous Rotations (p-Xylene)

**Input:** `p_xylene.rp`

```
bond_factor = 1.2
skip_factor = 0.7

# Independent rotation
1-2 e60

# Methyl group 1 rotates with main bond (same direction)
3-4 syn 1

# Methyl group 2 rotates with main bond (same direction)
5-6 syn 1

# Methyl group 3 rotates with main bond (opposite direction)
7-8 syn -1
```

**Result:** Only 1 independent bond → 6 rotation states
All methyl groups rotate predictably with main bond

### Example 4: Complex Molecule with Manual Bonds

**Input:** `complex.rp`

```
bond_factor = 1.2
skip_factor = 0.7

# Force metal-ligand bonds
8-12 bond
8-15 bond
8-23 bond

# Remove false positive bond
22-45 nobond

# Rotations
1-5 e90
6-10 e120
12-18 syn 2
```

**Result:** Manual bonds ensure correct connectivity
~15-20 valid conformers

### Example 5: Large Molecule with Safety Limits

**Input:** `large_complex.rp`

```
bond_factor = 1.2
skip_factor = 0.7

# Safety configuration for large generation
maxgen = 300         # Limit to 300 conformers
autoconfirm = true   # Skip interactive warnings

# Force metal-ligand bonds
43-48 bond
48-47 bond

# Multiple rotations (would generate 1728 theoretical conformers)
1-5 e30    # 12 states
6-10 e30   # 12 states  
12-18 e30  # 12 states
20-25 330  # 1 state
```

**Result:** 12 × 12 × 12 × 1 = 1,728 theoretical conformers
Limited to 300 actual conformers, no interactive prompts

## Conformer Generation Limits & Safety

Rotbond includes a comprehensive safety system to prevent system crashes from excessive conformer generation while maintaining full flexibility for advanced users.

### Default Safety Behavior

- **Default limit:** 500 conformers maximum
- **Interactive warning:** Displays for >500 theoretical conformers
- **Dynamic limit adjustment:** Set custom limits during runtime
- **Memory protection:** Prevents system crashes from excessive memory usage
- **Flexible user control:** Multiple options for handling large jobs

### Example: Enhanced Interactive Warning

```bash
rotbond large_molecule

  WARNING: Large conformer generation detected!
   Theoretical conformers: 1728
   This may consume significant memory and processing time.
   Current limit: 500 conformers

   Options:
   - Press Enter or 'y' to continue with current limit
   - Enter a number (e.g., 300) to set a new limit
   - Enter 'max' for unlimited generation
   - Enter 'n' to cancel

   Your choice: 200
✓ Limit set to 200 conformers.
```

**Interactive Options:**

- **Continue with current limit**: Press Enter, 'y', or 'yes'
- **Set custom limit**: Enter any positive integer (e.g., 100, 300, 1000)
- **Unlimited generation**: Enter 'max', 'maximum', or 'unlimited'
- **Cancel generation**: Enter 'n' or 'no'

### Configuration Options

**Custom Limits:**

```bash
# In molecule.rp
maxgen = 300         # Limit to 300 conformers
autoconfirm = false  # Show warnings (default)
```

**Unlimited Generation:**

```bash
# In molecule.rp  
maxgen = max         # Generate all conformers (use with caution!)
autoconfirm = true   # Skip warnings for batch processing
```

**Batch Processing:**

```bash
# In molecule.rp
maxgen = 1000        # Higher limit for powerful systems
autoconfirm = true   # No interactive prompts
```

### Safety Guidelines

**Recommended Limits by System:**

- **8GB RAM:** maxgen = 200-500
- **16GB RAM:** maxgen = 500-1000  
- **32GB+ RAM:** maxgen = 1000+ or unlimited
- **Batch jobs:** Always use `autoconfirm = true`

**Memory Estimation:**

- ~100 bytes per atom per conformer
- 50 atoms × 500 conformers ≈ 2.5 MB
- 100 atoms × 1000 conformers ≈ 10 MB
- 200 atoms × 2000 conformers ≈ 40 MB

**Best Practices:**

1. Start with default limits (500)
2. Monitor system resources during generation
3. Use `maxgen = max` only when necessary
4. Enable `autoconfirm` for unattended runs
5. Consider using trajectory output for large jobs

### Theoretical Conformer Calculation

The number of theoretical conformers is calculated as the cartesian product of all independent rotation states:

```bash
# Example: 3 independent bonds
Bond 1: e30 → 12 states
Bond 2: e60 → 6 states  
Bond 3: e90 → 4 states
Total: 12 × 6 × 4 = 288 theoretical conformers
```

**Common Scenarios:**

- 2 bonds × e60 = 6 × 6 = 36 conformers (no warning)
- 3 bonds × e30 = 12 × 12 × 12 = 1,728 conformers (warning)
- 4 bonds × e30 = 12⁴ = 20,736 conformers (large warning)

## Built-in Help System

Rotbond includes a comprehensive built-in help system accessible via command-line:

```bash
# Basic help
rotbond --help

# List all topics
rotbond --help topics

# Specific topic
rotbond --help examples
rotbond --help input
rotbond --help troubleshoot
```

**Available Topics:**

1. **usage** - Command-line usage and options
2. **features** - Feature overview
3. **input** - Input file formats (XYZ, .rp)
4. **output** - Output file formats
5. **examples** - Practical usage examples
6. **algorithms** - Algorithm details
7. **troubleshoot** - Troubleshooting guide
8. **reference** - Reference materials

**Error Message Help:**
All error messages include contextual help suggestions:

```bash
ERROR: Input file 'molecule.xyz' not found

For help creating input files, use: --help input
```

## Algorithm Details

Rotbond implements 5 core algorithms:

### 1. Molecular Graph Construction

- Builds bond connectivity using covalent radii
- Applies bond_factor to threshold calculation
- Supports manual bond forcing/removal
- Complexity: O(n²) for distance matrix

### 2. Fragment Identification

- Uses Breadth-First Search (BFS)
- Identifies atoms belonging to each rotating fragment
- Determines rotation axis and pivot point
- Complexity: O(n + m) where m = bonds

### 3. 3D Coordinate Transformation

- Implements Rodrigues' Rotation Formula
- Transforms coordinates to align rotation axis
- Applies rotation around arbitrary 3D axis
- Maintains molecular integrity

### 4. Conformer Generation

- Cartesian product over independent bonds
- Handles synchronous bonds automatically
- Progress reporting for large jobs
- Complexity: O(C × n²) where C = conformers

### 5. Steric Clash Validation

- Distance-based validation
- Uses covalent radii and skip_factor
- Early exit on first violation
- Ensures chemically reasonable conformers

## Performance Considerations

### Computational Complexity

- **Bond Detection:** O(n²)
- **Rotation:** O(f) where f = fragment size
- **Validation:** O(n²) per conformer
- **Generation:** O(C × n²) total

### Memory Usage

- Approximately 100 bytes per atom per conformer
- 50 atoms × 100 conformers ≈ 0.5 MB
- 100 atoms × 1000 conformers ≈ 10 MB

### Optimization Tips

1. **Reduce angle states** - Use explicit angles instead of steps
2. **Use synchronous rotations** - Decrease independent bonds
3. **Adjust skip_factor** - Higher values filter early
4. **Consider trajectory output** - Lower memory usage
5. **Batch processing** - For very large molecules

### Performance Guidelines

- Typical: 100-1000 conformers/second
- Depends on molecule size and validation strictness
- Trajectory output faster than individual files
- SSD recommended for large jobs

## Troubleshooting

### Common Issues and Solutions

**1. "ERROR: Input file 'molecule.xyz' not found"**

- Verify file exists in current directory
- Check filename spelling
- Remember: `rotbond molecule` looks for `molecule.xyz`

**2. "Need to create rotation parameters file"**

- **Automatic solution**: Just run `rotbond molecule` - it will create a template automatically!
- Edit the generated `molecule.rp` template file
- Uncomment and modify the rotation examples
- Run `rotbond molecule` again

**3. "ERROR: No rotation bonds specified"**

- Check rotation syntax in .rp file
- Ensure manual bonds come before rotations
- Add at least one rotation specification

**4. "Generated 0 valid conformers!"**

- Lower skip_factor (try 0.5 or 0.6)
- Reduce rotation angles
- Adjust bond_factor
- Check bond detection

**5. "Too many conformers (over 100,000)"**

- Reduce number of independent bonds
- Use larger step sizes (e90 instead of e30)
- Use explicit angles
- Add constraints

**6. "Circular reference detected"**

- Check synchronous bond specifications
- Avoid: Bond A syn B, Bond B syn A
- Remove one reference
- Use explicit angles instead

**7. "System running out of memory"**

- Reduce `maxgen` value (try 100-300)
- Use trajectory output instead of individual files
- Close other applications
- Consider using explicit angles instead of steps

**8. "Generation stopped at limit before completion"**

- This is normal behavior when `maxgen` limit is reached
- Increase `maxgen` value if more conformers needed
- Use `maxgen = max` for unlimited generation (with caution)
- Check if theoretical count exceeds your limit

**9. "Interactive prompt blocking batch processing"**

- Add `autoconfirm = true` to .rp file
- Set appropriate `maxgen` limit for your system
- Use this for unattended/automated runs

**10. "Want to adjust limit during runtime"**

- Use the interactive prompt when it appears
- Enter a custom number (e.g., 200, 1000) 
- Enter 'max' for unlimited generation
- This overrides the .rp file setting for current run only

### Debugging Tips

1. **Use template creation** - Run `rotbond molecule` to auto-generate comprehensive .rp template
2. **Use --verbose flag** - Shows detailed processing
3. **Start simple** - Begin with one bond rotation
4. **Check intermediate output** - Examine trajectory file
5. **Adjust parameters** - Try bond_factor: 0.8-1.5, skip_factor: 0.5-1.0
6. **Visualize results** - Use molecular viewer to check conformers

### Getting Help

```bash
# Built-in help
rotbond --help
rotbond --help topics
rotbond --help examples
rotbond --help troubleshoot

# Error messages
# All errors include contextual help suggestions
```

## Reference Materials

### Supported Elements

All elements from H to Lr (atomic numbers 1-103) are supported.

**Most Common:**

- H, C, N, O - Organic molecules
- P, S - Organophosphorus, organosulfur
- F, Cl, Br, I - Halogenated compounds
- Si - Silanes, silicones
- B - Boron compounds
- Na, K, Mg, Ca - Salts, complexes
- Transition metals - Fe, Co, Ni, Cu, Zn, etc.

### Covalent Radii (Å)

Based on OpenBabel covalent radii values:

| Element    | Symbol | Radius | Notes              |
| ---------- | ------ | ------ | ------------------ |
| Hydrogen   | H      | 0.23   | Very common        |
| Carbon     | C      | 0.68   | Organic backbone   |
| Nitrogen   | N      | 0.68   | Common in organics |
| Oxygen     | O      | 0.68   | Very common        |
| Fluorine   | F      | 0.64   | Halogen            |
| Phosphorus | P      | 1.05   | Common nonmetal    |
| Sulfur     | S      | 1.02   | Common nonmetal    |
| Chlorine   | Cl     | 0.99   | Halogen            |

*[Complete table available via: `rotbond --help reference`]*

### Angle Conventions

- **Unit:** Degrees (°)
- **Range:** -360° to +360°
- **Positive rotation:** Clockwise when viewing from atom1 → atom2
- **Zero rotation:** No rotation (original structure)
- **Step-based:** Always includes 0°, generates 360°/step angles

## Best Practices

### Getting Started

1. **Use automatic templates** - Run `rotbond molecule` to generate comprehensive .rp template
2. **Start simple** - Begin with one bond rotation from the template examples
3. **Verify output** - Check conformers before adding complexity
4. **Test with known molecule** - Ethane is a good test case

### Parameter Selection

- **bond_factor:** Start with 1.0-1.2
- **skip_factor:** Start with 0.7
- **Step size:** e30-e60 for detail, e90-e120 for flexibility

### Rotation Strategy

- Use synchronous rotations for symmetry
- Explicit angles for control
- Manual bonds for special cases
- Consider chemical constraints

### Validation

- Visualize results with molecular viewer
- Check bond lengths reasonable
- Verify expected conformations
- Monitor skipped conformers

### Performance

- Estimate conformers before generation
- Use trajectory for analysis
- Process in batches if needed
- Save successful configurations

### Conformer Limits & Safety

- **Start conservative** - Use default 500 limit initially
- **Monitor resources** - Watch memory usage during generation
- **Calculate theoretical count** - Multiply rotation states before running
- **Use autoconfirm for batch jobs** - Set `autoconfirm = true` for scripts
- **Set appropriate limits** - Match `maxgen` to your system capabilities
- **Unlimited generation** - Only use `maxgen = max` when necessary
- **Memory planning** - ~100 bytes per atom per conformer

## License

This project is open source. See [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.

## Additional Resources

- **Built-in Help:** `rotbond --help`
- **Examples:** `rotbond --help examples`
- **Troubleshooting:** `rotbond --help troubleshoot`
- **Reference:** `rotbond --help reference`

---

**Rotbond** - Generating molecular conformers with precision and efficiency
