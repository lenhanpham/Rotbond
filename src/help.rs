//! # Help System and Documentation
//!
//! This module provides comprehensive help and documentation for the Rotbond
//! conformer generation tool. It includes detailed explanations of all features,
//! file formats, parameters, and usage examples.
//!
//! ## Help Topics
//!
//! - **Basic usage**: Command line interface and file requirements
//! - **File formats**: XYZ structure files and .rp parameter files
//! - **Rotation specifications**: Step-based and explicit angle definitions
//! - **Parameters**: Bond detection and validation thresholds
//! - **Examples**: Complete workflows for common use cases
//! - **Troubleshooting**: Common issues and solutions
//!
//! ## Interactive Help
//!
//! The help system supports topic-based queries and provides contextual
//! information based on user needs.

use crate::utils::print_separator;

/// Prints basic help information including usage, options, and quick start guide.
/// 
/// This is the main help function displayed when users request general help
/// or provide invalid command-line arguments. It provides essential information
/// to get started with Rotbond.
/// 
/// # Output Sections
/// 
/// - Command usage syntax
/// - Basic command-line options
/// - Quick start example
/// - References to detailed help topics
/// 
/// # Examples
/// 
/// ```rust
/// use rotbond::help::print_basic_help;
/// 
/// print_basic_help();
/// // Displays comprehensive basic help information
/// ```
pub fn print_basic_help() {
    print_header();
    println!("{}", CMD_USAGE);
    println!();
    println!("{}", BASIC_OPTIONS);
    println!();
    println!("{}", QUICK_START);
    println!();
    println!("For more detailed information, use:");
    println!("  --help topics     List all available help topics");
    println!("  --help examples   Show practical usage examples");
    println!("  --help reference  Show reference materials");
    println!("  \n");
}

/// Prints version information for the Rotbond tool with a decorative banner.
/// 
/// Displays a comprehensive program information banner including version,
/// developer information, repository link, and citation request.
/// Used when users request version details with --version flag.
/// 
/// # Examples
/// 
/// ```rust
/// use rotbond::help::print_version;
/// 
/// print_version();
/// // Displays a decorative banner with:
/// // - Program name and title
/// // - Version {CARGO_PKG_VERSION} and release date
/// // - Developer information
/// // - GitHub repository link
/// // - Citation request and format
/// ```
pub fn print_version() {
    // Program information banner
    println!("  ***********************************************************************");
    println!("                                ROTBOND                                  ");
    println!("  ***********************************************************************");
    println!("# -----------------------------------------------------------------------#");
    println!("# Version {}  Release date: 2025                                      #", env!("CARGO_PKG_VERSION"));
    println!("# Developer: Le Nhan Pham                                                #");
    println!("# https://github.com/lenhanpham/Rotbond                                  #");
    println!("# -----------------------------------------------------------------------#");
    println!();
    println!("                                                                          ");
    println!("                                                                          ");
    println!("Please cite this project if you use Rotbond for your research             ");
    println!("                                                                          ");
    println!("# -----------------------------------------------------------------------#");
    println!("# L.N Pham, \"Rotbond - A Systematic molecular conformer generator\"       #");
    println!("# https://github.com/lenhanpham/Rotbond                                  #");
    println!("# -----------------------------------------------------------------------#");
}

/// Lists all available help topics with brief descriptions.
/// 
/// Provides an overview of all help topics available through the --help system.
/// Each topic covers a specific aspect of Rotbond functionality.
/// 
/// # Available Topics
/// 
/// 1. **usage** - Command-line usage and options
/// 2. **features** - Feature overview and capabilities
/// 3. **input** - Input file formats (XYZ, .rp)
/// 4. **output** - Output file formats and naming
/// 5. **examples** - Practical usage examples
/// 6. **algorithms** - Algorithm details and complexity
/// 7. **troubleshoot** - Troubleshooting guide
/// 8. **reference** - Reference materials and tables
/// 
/// # Examples
/// 
/// ```rust
/// use rotbond::help::print_help_topics;
/// 
/// print_help_topics();
/// // Displays numbered list of all available help topics
/// ```
pub fn print_help_topics() {
    print_header();
    println!("Available Help Topics:\n");
    println!("1. usage         - Command-line usage and options");
    println!("2. features      - Feature overview");
    println!("3. input         - Input file formats (XYZ, .rp)");
    println!("4. output        - Output file formats");
    println!("5. examples      - Practical usage examples");
    println!("6. algorithms    - Algorithm details");
    println!("7. troubleshoot  - Troubleshooting guide");
    println!("8. reference     - Reference materials");
    println!();
    println!("Use: rotbond --help <topic>");
    println!("Example: rotbond --help examples");
}

/// Prints detailed help information for a specific topic.
/// 
/// Provides comprehensive documentation for the requested help topic.
/// Topics are matched case-insensitively and support multiple aliases.
/// 
/// # Arguments
/// 
/// * `topic` - The help topic to display (case-insensitive)
/// 
/// # Supported Topics
/// 
/// - **usage, help**: Command-line usage and basic options
/// - **features, feature**: Feature overview and capabilities
/// - **input, format**: Input file formats and syntax
/// - **output, files**: Output file formats and naming conventions
/// - **examples, example**: Practical usage examples and workflows
/// - **algorithms, algorithm, tech**: Algorithm details and performance
/// - **troubleshoot, troubleshooting, errors**: Common issues and solutions
/// - **reference, ref**: Reference materials and lookup tables
/// 
/// # Examples
/// 
/// ```rust
/// use rotbond::help::print_help_topic;
/// 
/// print_help_topic("examples");
/// // Displays detailed examples and usage patterns
/// 
/// print_help_topic("input");
/// // Shows input file format documentation
/// ```
/// 
/// # Error Handling
/// 
/// If an unknown topic is requested, displays an error message and
/// shows the list of available topics.
pub fn print_help_topic(topic: &str) {
    let topic = topic.to_lowercase();

    match topic.as_str() {
        "usage" | "help" => {
            print_header();
            println!("{}", CMD_USAGE);
            println!();
            println!("{}", BASIC_OPTIONS);
            println!();
            println!("{}", CMD_EXAMPLES);
        }
        "features" | "feature" => {
            print_header();
            println!("FEATURE OVERVIEW\n");
            println!("{}", FEATURE_OVERVIEW);
        }
        "input" | "format" => {
            print_header();
            println!("INPUT FILE FORMATS\n");
            println!("{}", INPUT_FORMAT);
        }
        "output" | "files" => {
            print_header();
            println!("OUTPUT FILE FORMATS\n");
            println!("{}", OUTPUT_FORMAT);
        }
        "examples" | "example" => {
            print_header();
            println!("PRACTICAL EXAMPLES\n");
            println!("{}", EXAMPLES);
        }
        "algorithms" | "algorithm" | "tech" => {
            print_header();
            println!("ALGORITHM DETAILS\n");
            println!("{}", ALGORITHMS);
        }
        "troubleshoot" | "troubleshooting" | "errors" => {
            print_header();
            println!("TROUBLESHOOTING GUIDE\n");
            println!("{}", TROUBLESHOOTING);
        }
        "reference" | "ref" => {
            print_header();
            println!("REFERENCE MATERIALS\n");
            println!("{}", REFERENCE);
        }
        _ => {
            eprintln!("Unknown help topic: '{}'", topic);
            println!();
            print_help_topics();
        }
    }
}

/// Prints all help information in a comprehensive format.
/// 
/// Displays every help topic in sequence, separated by visual dividers.
/// This creates a complete reference manual that can be saved or printed.
/// Currently unused but available for future comprehensive help modes.
/// 
/// # Output Structure
/// 
/// 1. Basic help and usage
/// 2. Feature overview
/// 3. Input file formats
/// 4. Output file formats
/// 5. Practical examples
/// 6. Algorithm details
/// 7. Troubleshooting guide
/// 8. Reference materials
/// 
/// # Examples
/// 
/// ```rust
/// use rotbond::help::print_all_help;
/// 
/// print_all_help();
/// // Displays complete documentation (very long output)
/// ```
#[allow(dead_code)]
pub fn print_all_help() {
    print_basic_help();
    println!();
    print_separator();
    println!();
    print_help_topic("features");
    println!();
    print_separator();
    println!();
    print_help_topic("input");
    println!();
    print_separator();
    println!();
    print_help_topic("output");
    println!();
    print_separator();
    println!();
    print_help_topic("examples");
    println!();
    print_separator();
    println!();
    print_help_topic("algorithms");
    println!();
    print_separator();
    println!();
    print_help_topic("troubleshoot");
    println!();
    print_separator();
    println!();
    print_help_topic("reference");
}

/// Prints a formatted header for help sections.
/// 
/// Creates a consistent header format with separators and title information
/// for all help topic displays. Used internally by help functions.
fn print_header() {
    println!("  ***********************************************************************");
    println!("                                ROTBOND");
    println!("  ***********************************************************************");
    println!("# -----------------------------------------------------------------------#");
    println!("# Version {}  Release date: 2025                                      #", env!("CARGO_PKG_VERSION"));
    println!("# Developer: Le Nhan Pham                                                #");
    println!("# https://github.com/lenhanpham/Rotbond                                  #");
    println!("# -----------------------------------------------------------------------#");
    println!();
    println!();
    println!("Please cite this project if you use Rotbond for your research");
    println!();
    println!("# -----------------------------------------------------------------------#");
    println!("# L.N Pham, \"Rotbond - A Systematic molecular conformer generator\"       #");
    println!("# https://github.com/lenhanpham/Rotbond                                  #");
    println!("# -----------------------------------------------------------------------#");
    println!();
}

// =======================================================================
// HELP CONTENT CONSTANTS
// =======================================================================

pub const CMD_USAGE: &str = r#"
USAGE:
    rotbond [OPTIONS] <molecule_name>
    rotbond [OPTIONS] <molecule_name>.xyz

DESCRIPTION:
    Generate molecular conformers through bond rotations.

    Reads molecular structure from <molecule_name>.xyz and rotation
    parameters from <molecule_name>.rp, then generates all possible
    conformers as individual XYZ files and a trajectory file.

    If the .rp file doesn't exist, a comprehensive template file
    will be created automatically with examples and documentation.

REQUIRED INPUT FILES:
    <molecule_name>.xyz    - Molecular structure (XYZ format)
    <molecule_name>.rp     - Rotation parameters (auto-created if missing)

OUTPUT FILES:
    <molecule_name>_traj.xyz        - All conformers (trajectory)
    <molecule_name>_NN.xyz          - Individual conformers (smart padding)
"#;

pub const BASIC_OPTIONS: &str = r#"
OPTIONS:
    -h, --help              Show this help message
    -v, --version           Show version information
    --help topics           List all available help topics
    --help <topic>          Show specific help topic
                            (usage, features, input, output, examples,
                             algorithms, troubleshoot, reference)
    --verbose               Enable verbose output

ARGUMENTS:
    <molecule_name>         Base name for input/output files
                            (e.g., 'ethane' uses ethane.xyz and ethane.rp)
                            Also accepts 'molecule.xyz' format
"#;

pub const QUICK_START: &str = r#"
QUICK START:

1. Create XYZ file (ethane.xyz):
       8
       Ethane molecule
       C  0.000000  0.000000  0.000000
       H  0.000000  0.000000  1.090000
       ...

2. Run rotbond to create template:
       rotbond ethane
       (Creates ethane.rp template automatically)

3. Edit rotation file (ethane.rp):
       bond_factor = 1.2
       skip_factor = 0.7
       1-2 e60

4. Run rotbond again:
       rotbond ethane

5. Output:
       ethane_traj.xyz       - All conformers
       ethane_01.xyz ...     - Individual files
"#;

pub const CMD_EXAMPLES: &str = r#"
EXAMPLES:

    rotbond molecule_a           # Basic usage
    rotbond --help               # Show help
    rotbond --help examples      # Show examples
    rotbond --version            # Show version
    rotbond --verbose ethane     # Verbose output
"#;

pub const FEATURE_OVERVIEW: &str = r#"
FEATURE OVERVIEW

1. STEP-BASED ROTATIONS
   Format: atom1-atom2 e<step_angle>
   Example: 1-6 e60
   Generates: 0°, 60°, 120°, 180°, 240°, 300° (6 states)

2. EXPLICIT ANGLE ROTATIONS
   Format: atom1-atom2 angle1 angle2 angle3 ...
   Example: 2-5 0 60 120 180
   Uses: Exactly these 4 angles

3. SYNCHRONOUS ROTATIONS
   Format: atom1-atom2 syn <reference_bond>
   Example: 3-7 syn 1
   Effect: Bond 3 rotates exactly when bond 1 rotates

   Opposite direction:
   Format: atom1-atom2 syn -<reference_bond>
   Example: 4-8 syn -1
   Effect: Bond 4 rotates in opposite direction to bond 1

4. MANUAL BOND DEFINITIONS
   Force a bond:
   Format: atom1-atom2 bond
   Example: 8-12 bond
   Use: Force bonds that distance detection misses

   Remove a bond:
   Format: atom1-atom2 nobond
   Example: 22-45 nobond
   Use: Remove false positive bonds

5. CONFIGURATION PARAMETERS
   bond_factor = <value>
   - Multiplier for covalent radius threshold (default: 1.0)
   - Higher = detect more bonds, Lower = detect fewer

   skip_factor = <value>
   - Minimum distance validation (default: 0.7)
   - Filters conformers with steric clashes
   - Formula: d[atom1][atom2] >= (cov_radius1 + cov_radius2) × skip_factor
"#;

pub const INPUT_FORMAT: &str = r#"
INPUT FILE FORMATS

1. XYZ MOLECULAR STRUCTURE FILE (<molecule_name>.xyz)

   Format:
       n_atoms
       comment_line
       element x y z
       element x y z
       ...

   Example:
       8
       Ethane molecule - optimized geometry
       C       -2.6972703886      0.5190638733     -0.0000003688                 
       C       -1.1782073026      0.5190641624     -0.0000004138                 
       H       -3.0738517143     -0.1361520096     -0.8135565479                 
       H       -3.0738520178      1.5512320348     -0.1606558009                 
       H       -3.0738517177      0.1421113795      0.9742112740                 
       H       -0.8016259769      0.8960166481     -0.9742120609                 
       H       -0.8016256734     -0.5131039980      0.1606550254                 
       H       -0.8016259735      1.1742800518      0.8135557588 

   Notes:
   - Atom indices are 1-based in input files
   - Coordinates in Angstroms
   - Elements case-insensitive (C, c, carbon all valid)
   - Comment line can be any text

2. ROTATION PARAMETERS FILE (<molecule_name>.rp)

   Configuration Section (optional):
       bond_factor = 1.2    # Bond detection threshold
       skip_factor = 0.7    # Validation threshold

   Manual Bond Definitions (optional):
       # Force a bond
       8-12 bond

       # Remove a bond
       22-45 nobond

   Rotation Specifications:
       # Step-based rotation
       1-6 e60

       # Explicit angles
       2-5 0 60 120 180

       # Synchronous rotation (same direction)
       3-7 syn 1

       # Synchronous rotation (opposite direction)
       4-8 syn -1

   Flexible Format - Parameters can be in ANY order:
   - Comments start with '#'
   - Empty lines ignored
   - atom1-atom2 format (dash, not hyphen)
   - All parameters can be mixed - NO ordering restrictions!
   - All angles in degrees
   - Angles range: -360° to +360°

   Example Flexible File (all formats work the same):
       7-11 e30
       bond_factor = 1.2
       2-5 0 60 120
       skip_factor = 0.7
       8-12 bond
       22-45 nobond

   File Example:
       bond_factor = 1.2
       skip_factor = 0.7

       # Manual bond definitions
       8-12 bond
       15-23 bond
       22-45 nobond

       # Rotation specifications
       1-6 e60
       13-15 e120
       17-25 0 60
       24-26 syn 2
       25-78 e-60
       17-23 syn -3
"#;

pub const OUTPUT_FORMAT: &str = r#"
OUTPUT FILE FORMATS

1. TRAJECTORY FILE (<molecule_name>_traj.xyz)

   Contains all conformers in a single file.
   Format: Multiple XYZ structures concatenated

   Example:
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
       (pattern repeats for all conformers)

   Notes:
   - Standard XYZ format
   - Compatible with VMD, Chimera, PyMOL
   - Each structure has n_atoms + comment + coordinates
   - Comment includes conformer number for identification

2. INDIVIDUAL CONFORMER FILES (<molecule_name>_NN.xyz)

   One file per valid conformer.

   Naming Convention (Smart Padding):
       1-9 conformers:   molecule_1.xyz, molecule_2.xyz, ..., molecule_9.xyz
       10-99 conformers: molecule_01.xyz, molecule_02.xyz, ..., molecule_99.xyz
       100-999 conformers: molecule_001.xyz, ..., molecule_999.xyz
       1000-9999 conformers: molecule_0001.xyz, ..., molecule_9999.xyz
       10000+ conformers: molecule_00001.xyz, ...

   Example (56 conformers):
       molecule_a_traj.xyz     - Trajectory (all 56 conformers)
       molecule_a_01.xyz       - Conformer 1
       molecule_a_02.xyz       - Conformer 2
       ...
       molecule_a_56.xyz       - Conformer 56

   File Format:
       n_atoms
       Generated by Rotbond
       element x y z
       element x y z
       ...

   Notes:
   - 6-decimal precision for coordinates
   - Standard XYZ format
   - Each file is independent and complete
   - Files work with any molecular viewer

3. INTERACTIVE WARNINGS (for large jobs)

   When theoretical conformers > 500:

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

4. OUTPUT SUMMARY

   Console Output Example:
   ***********************************************************************
                                ROTBOND
   ***********************************************************************
   # -----------------------------------------------------------------------#
   # Version {CARGO_PKG_VERSION}  Release date: 2025                      #
   # Developer: Le Nhan Pham                                                #
   # https://github.com/lenhanpham/Rotbond                                  #
   # -----------------------------------------------------------------------#



   Please cite this project if you use Rotbond for your research

   # -----------------------------------------------------------------------#
   # L.N Pham, "Rotbond - A Systematic molecular conformer generator"       #
   # https://github.com/lenhanpham/Rotbond                                  #
   # -----------------------------------------------------------------------#

   Configuration:
     Bond factor:   1.20
     Skip factor:   0.70

   Rotation specifications:
     Bond 1: atoms 1-2 (6 rotation states)
              Angles: [0.0, 60.0, 120.0, 180.0, 240.0, 300.0]

   Generating conformers...

   =======================================================================
   Generation Summary:
     Theoretical combinations: 6
     Rejected (steric clashes): 0
     Valid conformers generated: 6
     Success rate: 100.0%
   =======================================================================

   Conformer generation completed successfully!

   =======================================================================
                           ROTBOND FOUND A BOOK :-)
   =======================================================================

   Title:  Pride and Prejudice
   Author: Jane Austen

   Summary:
      Elizabeth Bennet is Austen's most liberated and unambiguously
      appealing heroine, and Pride and Prejudice has remained over most
      of the past two centuries Austen's most popular novel...

   =======================================================================
                               Happy reading :-)
   =======================================================================

   Statistics:
   - Theoretical combinations: Total possible conformer combinations
   - Rejected (steric clashes): Conformers filtered due to atomic overlaps
   - Valid conformers generated: Conformers that passed validation
   - Success rate: Percentage of valid conformers (valid/theoretical × 100%)
   
   Additional Features:
   - Interactive conformer limit warnings for large generation jobs
   - Random book recommendations after successful completion
   - Professional banner with version and citation information
   - Consistent separator lines (71 characters) matching banner style
"#;

pub const EXAMPLES: &str = r#"
PRACTICAL EXAMPLES

1. SIMPLE SINGLE BOND ROTATION (Ethane)

   Input Files:
   -----------
   8
   Ethane molecule - optimized geometry
   C       -2.6972703886      0.5190638733     -0.0000003688                 
   C       -1.1782073026      0.5190641624     -0.0000004138                 
   H       -3.0738517143     -0.1361520096     -0.8135565479                 
   H       -3.0738520178      1.5512320348     -0.1606558009                 
   H       -3.0738517177      0.1421113795      0.9742112740                 
   H       -0.8016259769      0.8960166481     -0.9742120609                 
   H       -0.8016256734     -0.5131039980      0.1606550254                 
   H       -0.8016259735      1.1742800518      0.8135557588 

   ethane.rp:
       bond_factor = 1.2
       skip_factor = 0.7
       1-6 e60

   Command:
   --------
   rotbond ethane

   Output:
   -------
   ethane_traj.xyz          - Trajectory with all conformers
   ethane_1.xyz ...         - Individual conformer files

   Result:
   -------
   7 rotation states → 6 valid conformers (1 skipped due to steric clash)

---

2. MULTIPLE INDEPENDENT BONDS (Butane)

   Input Files:
   -----------
   butane.xyz: (25 atoms)

   butane.rp:
       bond_factor = 1.2
       skip_factor = 0.7

       # Rotate around C1-C2 bond
       1-2 e120

       # Rotate around C2-C3 bond
       2-3 e60

       # Rotate around C3-C4 bond
       3-4 e180

   Command:
   --------
   rotbond butane

   Result:
   -------
   3 bonds × [3, 6, 2] states = 36 total combinations
   ~30-35 valid conformers (after validation)

---

3. SYNCHRONOUS ROTATIONS (p-Xylene)

   Input Files:
   -----------
   p_xylene.xyz: (20 atoms)

   p_xylene.rp:
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

   Command:
   --------
   rotbond p_xylene

   Result:
   -------
   Only 1 independent bond → 7 rotation states
   Synchronous bonds follow automatically
   All methyl groups rotate predictably

---

4. COMPLEX MOLECULE WITH MANUAL BONDS

   Input Files:
   -----------
   complex.xyz: (45 atoms, includes metal center)

   complex.rp:
       bond_factor = 1.2
       skip_factor = 0.7

       # Force metal-ligand bonds (distance-based detection fails)
       8-12 bond
       8-15 bond
       8-23 bond

       # Remove false positive bond (steric proximity, not bonded)
       22-45 nobond

       # Rotations
       1-5 e90
       6-10 e120
       12-18 syn 2

   Command:
   --------
   rotbond complex

   Result:
   -------
   Manual bonds ensure correct connectivity
   False bonds removed to prevent incorrect rotations
   ~15-20 valid conformers

---

EXAMPLE FILE PATTERNS:

Small molecules (1-9 conformers):
    molecule_1.xyz, molecule_2.xyz, ..., molecule_9.xyz

Medium molecules (10-99 conformers):
    molecule_01.xyz, molecule_02.xyz, ..., molecule_56.xyz

Large molecules (100-999 conformers):
    molecule_001.xyz, molecule_002.xyz, ..., molecule_216.xyz

Very large (10000+ conformers):
    molecule_00001.xyz, molecule_00002.xyz, ..., molecule_12345.xyz

---

TIPS FOR EFFECTIVE USAGE:

1. Start Simple
   - Begin with one bond rotation to test
   - Verify output before adding more bonds

2. Choose Appropriate Step Sizes
   - Small molecules: e30, e45, e60
   - Large molecules: e60, e90, e120
   - Very flexible: e30 for fine detail

3. Use Explicit Angles for Control
   - Specify exact angles you need
   - Avoid unnecessary conformers

4. Synchronous Rotations for Symmetry
   - Methyl groups rotate together
   - Symmetric fragments
   - Maintains molecular symmetry

5. Manual Bonds for Special Cases
   - Metal complexes
   - Coordinate bonds
   - Aromatic systems
   - Remove false positives

6. Adjust Validation
   - Strict: skip_factor = 0.8-1.0
   - Relaxed: skip_factor = 0.5-0.7
   - Consider molecular flexibility

7. Estimate Before Generation
   - Count independent bonds
   - Multiply angle states
   - Check for combinatorial explosion
   - 5 bonds × 7 states = 16,807 combinations!

8. Use Trajectory Files
   - Easier for visualization
   - Better for analysis
   - Compatible with MD software
   - Use for initial screening

9. Validate Results
   - Check bond lengths reasonable
   - Verify angles make chemical sense
   - Look for unexpected conformations
   - Use chemical knowledge

10. Performance Optimization
    - Reduce unnecessary angle states
    - Use skip_factor to filter early
    - Consider parallel processing
    - Save intermediate results
"#;

pub const ALGORITHMS: &str = r#"
ALGORITHM DETAILS

1. MOLECULAR GRAPH CONSTRUCTION

   Purpose: Determine which atoms are bonded

   Process:
   1. Calculate distance matrix: d[i][j] = √[(xᵢ-xⱼ)² + (yᵢ-yⱼ)² + (zᵢ-zⱼ)²]
   2. Initialize bonds based on covalent radii:
      if d[i][j] < (cov_radius(i) + cov_radius(j)) × bond_factor
         then atoms i and j are bonded
   3. Apply manual bond definitions:
      - Add forced bonds (even if distance > threshold)
      - Remove forbidden bonds (even if distance < threshold)
   4. Build adjacency list representation

   Covalent Radii (Å):
      H: 0.23, C: 0.68, N: 0.68, O: 0.68, F: 0.64
      P: 1.05, S: 1.02, Cl: 0.99, Br: 1.21, I: 1.40
      (Full table in Reference section)

   Bond Factor Effect:
      - bond_factor = 1.0: Standard threshold
      - bond_factor > 1.0: Looser detection (more bonds)
      - bond_factor < 1.0: Stricter detection (fewer bonds)

---

2. FRAGMENT IDENTIFICATION

   Purpose: Identify atoms belonging to each rotating fragment

   Input: Bond between atoms a and b

   Algorithm (Breadth-First Search):
   1. Initialize fragment1 (atoms bonded to a, excluding b)
   2. Initialize fragment2 (atoms bonded to b, excluding a)
   3. For each atom in fragment1:
      - Add all bonded neighbors (except a and b)
      - Continue until no new atoms added
   4. For each atom in fragment2:
      - Add all bonded neighbors (except a and b)
      - Continue until no new atoms added

   Result: Two disjoint atom sets (fragments)

   Rotation Convention:
      - Fragment = atoms attached to atom2
      - Atom1 serves as rotation pivot
      - Bond defines rotation axis

---

3. 3D COORDINATE TRANSFORMATION

   Purpose: Rotate fragment around bond axis

   Method: Rodrigues' Rotation Formula

   Steps:
   1. Translation: Move bond to origin
      T = translation_matrix(-p₁)

   2. Axis Alignment: Rotate bond to align with Z-axis
      k = normalized(p₂ - p₁)  (bond direction vector)
      Calculate Rₓ, Rᵧ, Rₓ⁻¹, Rᵧ⁻¹

   3. Rotation: Apply rotation around Z-axis by angle θ
      R_z(θ) = rotation_matrix_z(θ)

   4. Restoration: Apply inverse transformations
      M = T⁻¹ × Rₓ⁻¹ × Rᵧ⁻¹ × R_z(θ) × Rᵧ × Rₓ × T

   5. Apply transformation to all fragment atoms

   Mathematical Formula:
      v_rot = v × cos(θ) + (k × v) × sin(θ) + k × (k · v) × (1 - cos(θ))
      where:
        v = point to rotate
        k = unit vector along rotation axis
        θ = rotation angle

---

4. CONFORMER GENERATION

   Purpose: Generate all possible conformer combinations

   Method: Cartesian Product

   Steps:
   1. Separate independent bonds from synchronous bonds
   2. Calculate total combinations:
      total = ∏(angle_sets[i].length) for independent bonds i
   3. Generate combinations recursively:
      For each independent bond:
        For each angle:
          Apply rotation
          Recurse to next bond
   4. Apply synchronous bond rotations:
      Use reference bond's angle × direction
   5. Validate conformer (steric clash check)
   6. Save if valid

   Complexity:
      - Time: O(N × M × A)
        N = number of independent bonds
        M = average angle states per bond
        A = validation cost (O(n²) for n atoms)
      - Space: O(C) for storing conformers
        C = number of valid conformers

---

5. STERIC CLASH VALIDATION

   Purpose: Filter out unrealistic conformers

   Method: Distance-based validation

   Algorithm:
   For each atom pair (i, j):
      distance = √[(xᵢ-xⱼ)² + (yᵢ-yⱼ)² + (zᵢ-zⱼ)²]
      min_distance = (cov_radius(i) + cov_radius(j)) × skip_factor
      if distance < min_distance:
         return INVALID (steric clash)

   If all pairs pass:
      return VALID

   Skip Factor Effect:
      - skip_factor = 1.0: Standard minimum distance
      - skip_factor > 1.0: Stricter validation (fewer conformers)
      - skip_factor < 0.5: Very relaxed (may include clashes)

   Trade-offs:
      - Lower skip_factor: More conformers, faster generation
      - Higher skip_factor: Fewer conformers, slower generation
      - Recommended: 0.6-0.8 for most molecules

---

PERFORMANCE CONSIDERATIONS:

1. Combinatorial Explosion
   - Each independent bond multiplies conformers
   - 5 bonds × 7 angles = 16,807 conformers
   - 6 bonds × 7 angles = 117,649 conformers!
   - Consider reducing angle states

2. Validation Cost
   - O(n²) distance checks per conformer
   - n = number of atoms
   - 50 atoms = 1,225 distance checks per conformer
   - Optimize by:
     * Early exit on first violation
     * Spatial indexing for large molecules
     * Parallel processing

3. Memory Usage
   - Store all valid conformers in memory
   - For large sets, stream to disk
   - Consider trajectory-only output

4. I/O Performance
   - Individual files: Many small writes
   - Trajectory file: Single large write
   - Use buffered I/O

5. Optimization Strategies
   - Reduce unnecessary angle states
   - Use explicit angles instead of steps
   - Filter early with higher skip_factor
   - Generate in batches
   - Parallel generation (independent processes)

ALGORITHM COMPLEXITY SUMMARY:

Operation                Time Complexity      Space Complexity
-----------------       -----------------   -----------------
Bond Detection          O(n²)               O(n²)
Fragment ID             O(n + m)            O(n)
Rotation               O(f)                O(1)
Validation             O(n²)               O(1)
Generation             O(C × n²)           O(C × n)

n = number of atoms
m = number of bonds
f = fragment size
C = number of conformers
"#;

pub const TROUBLESHOOTING: &str = r#"
TROUBLESHOOTING GUIDE

COMMON ERRORS AND SOLUTIONS

1. "ERROR: Input file 'molecule.xyz' not found"

   Cause: Missing XYZ input file

   Solution:
   - Verify file exists: ls molecule.xyz
   - Check filename spelling
   - Ensure file is in current directory
   - Remember: rotbond molecule expects molecule.xyz (not molecule_input.xyz)

---

2. "ERROR: Rotation parameters file 'molecule.rp' not found"

   Cause: Missing .rp rotation parameters file

   Solution:
   - Create rotation parameters file
   - Use exact molecule name (case-sensitive)
   - File must be: <molecule_name>.rp
   - Minimum file can be empty, but should have bond_factor/skip_factor

---

3. "ERROR: Invalid atom index"

   Cause: Rotation specification references non-existent atom

   Solution:
   - Check atom indices (1-based in input)
   - Verify XYZ file has correct number of atoms
   - Count atoms carefully (including all H atoms)
   - Example: "1-6" rotates bond between atom 1 and atom 6

---

4. "ERROR: Invalid bond specification"

   Cause: Bond format incorrect

   Solution:
   - Use dash (-), not hyphen (–) or em-dash (—)
   - Format: atom1-atom2 (e.g., 1-6)
   - No spaces around dash
   - Both atom indices required

---

5. "ERROR: No rotation bonds specified"

   Cause: Rotation file has no valid rotation specifications

   Solution:
   - Check rotation syntax
   - Ensure manual bonds come before rotations
   - Add at least one rotation line:
     * Step-based: 1-6 e60
     * Explicit: 1-6 0 60 120
     * Synchronous: 3-7 syn 1

---

6. "ERROR: Invalid step angle"

   Cause: Step angle format incorrect

   Solution:
   - Format: e<number> (e.g., e60, e90, e120)
   - No space between 'e' and number
   - Angle must be numeric (60, not sixty)
   - Range: 0-360 degrees

---

7. "Generated 0 valid conformers!"

   Cause: All conformers rejected by validation

   Solution:
   - Lower skip_factor (try 0.5 or 0.6)
   - Reduce rotation angles
   - Check for impossible rotations
   - Verify bond detection (adjust bond_factor)
   - Use explicit angles instead of steps

---

8. "Too many conformers (over 100,000)"

   Cause: Combinatorial explosion

   Solution:
   - Reduce number of independent bonds
   - Use larger step sizes (e90 instead of e30)
   - Use explicit angles
   - Add constraints or filters

---

9. "ERROR: Bond references itself"

   Cause: Synchronous rotation self-reference

   Solution:
   - Check synchronous bond specification
   - Cannot use: 3-7 syn 3
   - Must reference different bond
   - Example: 3-7 syn 1 (bond 3 syncs with bond 1)

---

10. "ERROR: Circular reference detected"

    Cause: Synchronous bonds in a loop

    Solution:
    - Check all syn specifications
    - Avoid: Bond A syn B, Bond B syn A
    - Fix by removing one reference
    - Use explicit angles instead

---

11. "Conformer has steric clashes"

    Cause: Atoms too close together

    Solution:
    - Lower skip_factor (e.g., 0.5)
    - Reduce rotation angles
    - Check skip_factor value (default 0.7)
    - May be chemically valid (check manually)
    - Use --verbose to see which conformers skipped

---

12. "Memory error during generation"

    Cause: Too many conformers for available memory

    Solution:
    - Reduce conformer count
    - Use trajectory-only output
    - Generate in batches
    - Increase system memory
    - Use 64-bit system

---

13. "Cannot create output file"

    Cause: Permission or disk space issue

    Solution:
    - Check write permissions in directory
    - Verify sufficient disk space
    - Close file if open in another program
    - Check filename validity

---

14. "DistMatrix function error"

    Cause: Fragment identification failed

    Solution:
    - Check bond connectivity
    - Adjust bond_factor
    - Verify manual bonds
    - Ensure molecule is connected

---

15. "Invalid angle value"

    Cause: Non-numeric angle in explicit list

    Solution:
    - Check explicit angle list
    - All values must be numeric
    - Format: 0 60 120 180
    - Range: -360 to +360 degrees

---

16. "Too many rotation specifications"

    Cause: Exceeded maximum bonds (usually 26, letters a-z)

    Solution:
    - Combine rotations
    - Use synchronous bonds
    - Reduce number of rotating bonds
    - Process in multiple runs

---

PERFORMANCE ISSUES:

1. Very Slow Generation

   Causes & Solutions:
   - Too many combinations: Reduce independent bonds
   - Large molecule: Increase skip_factor
   - Many validations: Lower skip_factor
   - Use explicit angles (fewer than steps)
   - Consider parallel processing

---

2. High Memory Usage

   Causes & Solutions:
   - Many conformers: Stream to disk
   - Trajectory-only output
   - Batch processing
   - 64-bit system recommended

---

3. I/O Bottleneck

   Causes & Solutions:
   - Many individual files: Use trajectory
   - Slow disk: Use faster storage
   - Buffered I/O: Already implemented
   - SSD recommended for large jobs

---

DEBUGGING TIPS:

1. Use --verbose Flag
   - Shows detailed processing information
   - Progress updates
   - Which conformers skipped
   - Fragment identification

2. Start Simple
   - Begin with one bond rotation
   - Verify output before adding complexity
   - Test with known molecule

3. Check Intermediate Output
   - Examine trajectory file
   - Verify conformers visually
   - Check bond lengths
   - Validate angles

4. Adjust Parameters
   - bond_factor: Try 0.8-1.5 range
   - skip_factor: Try 0.5-1.0 range
   - Step size: Larger for flexibility

5. Use Chemical Knowledge
   - Check if rotations are chemically reasonable
   - Verify expected conformations present
   - Look for symmetry

---

GETTING HELP:

1. Built-in Help System:
   rotbond --help                    # Basic help
   rotbond --help topics            # List topics
   rotbond --help examples          # Show examples
   rotbond --help reference         # Reference materials

2. Documentation:
   README.md                        # Complete guide
   Rotbond_plan.md                  # Specification

3. Common File Locations:
   Input:  ./molecule.xyz, ./molecule.rp
   Output: ./molecule_traj.xyz, ./molecule_NN.xyz

4. Log Files:
   Use --verbose for detailed console output
   No log files (console-only)

---

BEST PRACTICES:

1. Start Small
   - One bond rotation first
   - Verify output
   - Add complexity gradually

2. Validate Results
   - Check bond lengths reasonable
   - Verify angles make sense
   - Visualize with molecular viewer

3. Choose Parameters Carefully
   - bond_factor: Start with 1.0-1.2
   - skip_factor: Start with 0.7
   - Adjust based on results

4. Use Synchronous Rotations
   - Reduce independent bonds
   - Maintain molecular symmetry
   - More realistic conformers

5. Monitor Progress
   - Use --verbose
   - Watch for skipped conformers
   - Adjust parameters as needed

6. Backup Input Files
   - Preserve original structures
   - Keep different parameter sets
   - Document successful configurations

7. Test on Known Cases
   - Ethane rotation (standard test)
   - Butane (multiple bonds)
   - Compare with literature

8. Optimize for Your Molecule
   - Flexible: Smaller step sizes
   - Rigid: Larger step sizes
   - Consider chemical constraints

9. Use Trajectory for Analysis
   - Easier visualization
   - Better for statistics
   - Compatible with MD software

10. Document Your Workflow
    - Record successful parameters
    - Note issues and solutions
    - Share with colleagues
"#;

pub const REFERENCE: &str = r#"
REFERENCE MATERIALS

1. COVALENT RADII TABLE (Å)

   Based on OpenBabel covalent radii values
   with 0.4 Å tolerance for bond detection

   Element    Symbol    Radius    Notes
   --------   ------    ------    -----
   Hydrogen  H         0.23      Very common
   Helium    He        1.50      Noble gas
   Lithium   Li        0.68      Alkali metal
   Beryllium Be        0.35      Alkaline earth
   Boron     B         0.83      Metalloid
   Carbon    C         0.68      Organic backbone
   Nitrogen  N         0.68      Common in organics
   Oxygen    O         0.68      Very common
   Fluorine  F         0.64      Halogen
   Neon      Ne        1.50      Noble gas
   Sodium    Na        0.97      Alkali metal
   Magnesium Mg        1.10      Alkaline earth
   Aluminum  Al        1.35      Metallic
   Silicon   Si        1.20      Metalloid
   Phosphorus P        1.05      Common nonmetal
   Sulfur    S         1.02      Common nonmetal
   Chlorine  Cl        0.99      Halogen
   Argon     Ar        1.51      Noble gas
   Potassium K         1.33      Alkali metal
   Calcium   Ca        0.99      Alkaline earth
   ... (extends to all elements)
   [Full table includes all 118 elements]

   Default for unknown elements: 1.50 Å

   Usage in bond detection:
      threshold = (radius₁ + radius₂) × bond_factor

---

2. SUPPORTED ELEMENTS

   All elements from H to Lr (atomic numbers 1-103)
   Full periodic table support

   Most Common in Rotbond:
   - H, C, N, O: Organic molecules
   - P, S: Organophosphorus, organosulfur
   - F, Cl, Br, I: Halogenated compounds
   - Si: Silanes, silicones
   - B: Boron compounds
   - Na, K, Mg, Ca: Salts, complexes
   - Transition metals: Fe, Co, Ni, Cu, Zn, etc.
   - Heavy elements: Br, I, Pd, Pt, Au, Hg, etc.

   Element Format:
   - Case-insensitive: C, c, CARBON all valid
   - Standard symbols: Use periodic table notation
   - One or two characters: H, He, Li, etc.

---

3. ANGLE CONVENTIONS

   Unit: Degrees (°)
   Range: -360° to +360°

   Positive Rotation:
   - Clockwise when viewing bond from atom1 → atom2
   - Right-handed coordinate system
   - Standard chemistry convention

   Negative Rotation:
   - Counterclockwise
   - Opposite of positive

   Zero Rotation:
   - No rotation (original structure)
   - Included by default in step rotations

   Step-Based Generation:
   - Generates: 0°, step, 2×step, ..., 360°-step
   - Example (e60): 0°, 60°, 120°, 180°, 240°, 300° (6 states)
   - Always includes 0°, number of states = 360°/step

   Explicit Angles:
   - User specifies exact angles
   - No automatic generation
   - Duplicates automatically removed
   - Sorted for consistency

---

4. COORDINATE SYSTEM

   Standard: Cartesian XYZ
   Units: Angstroms (Å)
   Origin: As specified in input file

   Rotation Convention:
   - Axis: Bond direction (atom1 → atom2)
   - Pivot: atom1 position
   - Direction: From atom1 to atom2

   Example:
      Bond 1-6: atom1=1, atom2=6
      Rotation axis: Vector from atom1 to atom2
      Pivot point: Position of atom1
      Positive angle: Clockwise when looking from atom1 to atom6

---

5. BOND DETECTION PARAMETERS

   Default bond_factor = 1.0
   Formula: threshold = (radius₁ + radius₂) × bond_factor

   Effects:
   bond_factor < 1.0: Stricter detection (fewer bonds)
      0.5: Very strict, only obvious bonds
      0.8: Moderately strict
      1.0: Standard OpenBabel threshold
      1.2: Looser, more bonds detected
      1.5: Very loose, may include false positives
      2.0: Maximum, detect all proximities

   Recommendations:
   - Organic molecules: 1.0-1.2
   - Metal complexes: 1.2-1.5
   - Crystal structures: 0.8-1.0
   - Flexible molecules: 1.0
   - When in doubt: Start with 1.0

---

6. VALIDATION PARAMETERS

   Default skip_factor = 0.7
   Formula: min_distance = (radius₁ + radius₂) × skip_factor

   Effects:
   skip_factor < 0.7: More lenient (accept closer atoms)
      0.4: Very lenient, may include clashes
      0.5-0.6: Moderately lenient
      0.7: Standard validation
      0.8-0.9: Strict validation
      1.0: Very strict (rarely used)
      >1.0: Excessively strict (reject valid conformers)

   Recommendations:
   - Small rigid molecules: 0.8-1.0
   - Flexible molecules: 0.6-0.8
   - When in doubt: Start with 0.7
   - If too many rejected: Lower to 0.5-0.6
   - If too many accepted: Raise to 0.8

---

7. FILE FORMAT SPECIFICATIONS

   XYZ Format (Input):
      Line 1: Integer (atom count)
      Line 2: Comment string
      Lines 3-n+2: element x y z (floating point)
      Units: Ångstroms
      Delimiter: Whitespace
      Precision: 6+ decimal places recommended

   .rp Format (Input):
      Line-based text format
      Encoding: UTF-8
      Line endings: Unix (LF) or Windows (CRLF)
      Comments: '#' to end of line
      Delimiter: Whitespace
      Case-sensitive keywords

   Trajectory XYZ (Output):
      Multiple XYZ structures concatenated
      No blank lines between structures
      Each structure: n_atoms + comment + coordinates
      Total size: n_atoms × (n_conformers + 1) lines

   Individual XYZ (Output):
      Single XYZ structure per file
      Same format as input
      Filenames with zero-padding

---

8. ROTATION SYNTAX REFERENCE

   Step-based:
      Format: atom1-atom2 e<number>
      Examples:
        1-6 e60      → 6 states (0°, 60°, 120°, 180°, 240°, 300°)
        2-5 e30      → 12 states (0°, 30°, 60°, 90°, 120°, 150°, 180°, 210°, 240°, 270°, 300°, 330°)
        3-7 e90      → 4 states (0°, 90°, 180°, 270°)
        4-8 e120     → 3 states (0°, 120°, 240°)

   Explicit:
      Format: atom1-atom2 angle1 angle2 angle3 ...
      Examples:
        1-6 0 60 120 180              → 4 states
        2-5 0 90 180 270              → 4 states
        3-7 -120 -60 0 60 120         → 5 states
        4-8 30 45 60 90 120 150 180   → 7 states

   Synchronous (same direction):
      Format: atom1-atom2 syn <reference>
      Examples:
        3-7 syn 1       → Bond 3 syncs with bond 1
        5-9 syn 2       → Bond 5 syncs with bond 2

   Synchronous (opposite direction):
      Format: atom1-atom2 syn -<reference>
      Examples:
        4-8 syn -1      → Bond 4 syncs with bond 1 (opposite)
        6-10 syn -2     → Bond 6 syncs with bond 2 (opposite)

---

9. PERFORMANCE GUIDELINES

   Memory Usage:
      ~100 bytes per atom per conformer
      Example: 50 atoms × 100 conformers = 0.5 MB
      100 atoms × 1000 conformers = 10 MB
      200 atoms × 10000 conformers = 200 MB

   Processing Speed:
      Typical: 100-1000 conformers/second
      Depends on: molecule size, validation strictness
      Use skip_factor to control speed

   Optimization Tips:
      - Reduce angle states
      - Use explicit angles
      - Increase skip_factor
      - Use trajectory output
      - Process in batches

---

10. CHEMISTRY BACKGROUND

    Conformers:
    - Different spatial arrangements of same molecule
    - Connected by bond rotations
    - Same connectivity, different geometry
    - Important for:
      * Energy calculations
      * Spectral properties
      * Reaction pathways
      * Molecular recognition

    Steric Effects:
    - Repulsion between electron clouds
    - Prevent atoms from overlapping
    - skip_factor validates physical reasonableness
    - Based on van der Waals radii

    Bond Rotation:
    - Most common conformational change
    - Around single (σ) bonds
    - Low energy barrier
    - Temperature-dependent populations

    Synchronous Rotations:
    - Physically linked motions
    - Methyl groups
    - Symmetric fragments
    - Coupled bond rotations

    Validation Importance:
    - Prevents impossible conformations
    - Filters high-energy structures
    - Ensures chemical reasonableness
    - Saves computational resources
"#;
