#![deny(missing_docs)]
#![warn(clippy::all)]
//! # Rotbond - Molecular Conformer Generation Tool
//!
//! Rotbond is a high-performance molecular conformer generation tool that systematically
//! explores conformational space by rotating around rotatable bonds. It uses dihedral
//! angle-based rotation algorithms to generate chemically meaningful conformers.
//!
//! ## Features
//!
//! - **Automatic template creation** - generates comprehensive .rp files with examples
//! - **Flexible input formats** - accepts both `molecule` and `molecule.xyz` formats
//! - **Automatic bond detection** using OpenBabel-style covalent radii
//! - **Dihedral-based rotation** with Rodrigues rotation formula for mathematical precision
//! - **Steric clash detection** to filter invalid conformers
//! - **Multiple output formats** including individual XYZ files and trajectory files
//! - **Flexible rotation specifications** supporting step-based (e30, e60) and explicit angles
//! - **Synchronous rotations** for coordinated bond movements
//! - **Interactive safety warnings** for large conformer generation jobs
//! - **Professional XYZ formatting** with proper column alignment
//! - **High performance** with optimized algorithms and memory management
//!
//! ## Performance Considerations
//!
//! - **Complexity**: O(n^k) where n = average angles per bond, k = number of rotatable bonds
//! - **Memory usage**: Scales with the number of valid conformers generated
//! - **Recommended limits**: < 10 rotatable bonds for exhaustive search
//! - **Large molecules**: Use larger step sizes (e120, e180) to reduce search space
//!
//! ## Best Practices
//!
//! - Use automatic template creation - run `rotbond molecule` to generate comprehensive .rp template
//! - Start with `bond_factor = 1.0` and `skip_factor = 0.7`
//! - Use `e60` or `e120` for initial exploration, `e30` for detailed studies
//! - Check success rates - low rates may indicate parameter adjustment needed
//! - For flexible molecules, consider focusing on key rotatable bonds only
//!
//! ## Usage
//!
//! ```bash
//! rotbond molecule_name
//! rotbond molecule_name.xyz
//! ```
//!
//! This expects `molecule_name.xyz` (structure) and `molecule_name.rp` (rotation parameters).
//! If the .rp file doesn't exist, a comprehensive template will be created automatically
//! with examples and documentation. Simply edit the template and run again.
//!
//! ## Quick Start Example (New Simplified Workflow)
//!
//! 1. **Create structure file** (`butane.xyz`):
//! ```text
//!14
//!Butane molecule
//!C    -4.390308655000     -0.632987590700      0.043677475300
//!H    -5.330057936000     -0.527510732400     -0.503444569700
//!H    -4.389058916600      0.077648432000      0.873491606400
//!H    -4.345674372500     -1.643033486700      0.457086407500
//!C    -3.188670330900     -0.384837924800     -0.886583891400
//!C    -3.226558849100      1.036845318100     -1.485516216000
//!H    -2.262792259600     -0.522845611400     -0.322588193200
//!H    -3.198108077300     -1.124349892400     -1.691208307400
//!H    -3.217121102700      1.776357285700     -0.680891800000
//!C    -2.024920525000      1.284994984100     -2.415777582700
//!H    -4.152436920400      1.174853004700     -2.049511914300
//!H    -1.085171244000      1.179518125700     -1.868655537800
//!H    -2.026170263400      0.574358961400     -3.245591713900
//!H    -2.069554807400      2.295040880000     -2.829186515000
//! ```
//!
//! 2. **Run rotbond to create template**:
//! ```bash
//! rotbond butane
//! ```
//! Output: "✓ Template file 'butane.rp' created successfully!"
//!
//! 3. **Edit the generated template** (`butane.rp`):
//! ```text
//! bond_factor = 1.0
//! skip_factor = 0.7
//! 1-5 e60  # Rotate C-C bond every 60°
//! ```
//!
//! 4. **Run conformer generation**:
//! ```bash
//! rotbond butane
//! ```
//!
//! 5. **Output files**:
//! - `butane_1.xyz`, `butane_2.xyz`, ..., `butane_6.xyz` (individual conformers)
//! - `butane_traj.xyz` (all conformers in trajectory format)

mod algorithms;
mod booklist;
mod errors;
mod generate;
mod help;
mod io;
mod molecule;
mod rotation;
mod utils;

#[cfg(test)]
mod tests;

use booklist::print_random_book;
use generate::generate_conformers;
use help::{print_basic_help, print_help_topic, print_help_topics, print_version};
use molecule::{Bond, Molecule, RotationSpec};
use rotation::{generate_angle_sets, parse_rotation_file, validate_synchronous_references};
use std::env;
use std::io::{Write, stdin};
use std::process;
use utils::{
    print_configuration, print_rotation_summary, print_scanning_summary, print_summary, print_usage,
};

/// Creates a template rotation parameters (.rp) file with comprehensive examples and documentation.
///
/// Generates a well-documented template file that includes:
/// - Configuration parameters with explanations
/// - Examples of all rotation specification types
/// - Manual bond definition examples
/// - Comments explaining each feature
///
/// # Arguments
///
/// * `filename` - Path where the template .rp file should be created
///
/// # Returns
///
/// Result indicating success or failure of file creation
///
/// # Examples
///
/// ```rust
/// create_template_rp_file("molecule.rp")?;
/// ```
fn create_template_rp_file(filename: &str) -> std::io::Result<()> {
    use std::fs::File;
    use std::io::Write;

    let template_content = r#"# =======================================================================
# ROTBOND PARAMETERS FILE
# =======================================================================
# This is a template file with examples of all available features.
# Edit this file to specify your parameters, then run rotbond again.
#
# IMPORTANT: Choose either ROTATION or SCANNING mode - they cannot be mixed!
#
# For detailed help: rotbond --help input
# =======================================================================

# -----------------------------------------------------------------------
# CONFIGURATION PARAMETERS (Optional)
# -----------------------------------------------------------------------
# These parameters apply to both rotation and scanning modes

bond_factor = 1.0
# - Multiplier for covalent radius threshold (default: 1.0)
# - Higher values detect more bonds, lower values detect fewer
# - Recommended: 1.0-1.2 for organic molecules, 1.2-1.5 for metal complexes

skip_factor = 0.7
# - Minimum distance validation threshold (default: 0.7)
# - Filters conformers with steric clashes
# - Lower values = more lenient, higher values = stricter
# - Recommended: 0.6-0.8 for most molecules

maxgen = 500
# - Maximum number of conformers to generate (default: 500)
# - Use 'max' for unlimited generation
# - Higher values may consume significant memory

autoconfirm = false
# - Automatically confirm large generation jobs (default: false)
# - Set to 'true' to skip interactive warnings for large jobs

# -----------------------------------------------------------------------
# MANUAL BOND DEFINITIONS (Optional)
# -----------------------------------------------------------------------
# Force or remove bonds that automatic detection might miss or incorrectly identify
# These work with both rotation and scanning modes

# Force a bond (useful for metal-ligand bonds, weak interactions)
# 8-12 bond
# 15-23 bond

# Remove a false positive bond (atoms close but not bonded)
# 22-45 nobond

# =======================================================================
# ROTATION MODE SPECIFICATIONS
# =======================================================================
# Use these for traditional dihedral angle rotation conformer generation
# Uncomment and modify the examples below for your molecule

# STEP-BASED ROTATIONS (generates evenly spaced angles)
# Format: atom1-atom2 e<step_angle>
# Examples:
# 1-6 e60        # 6 states: 0°, 60°, 120°, 180°, 240°, 300°
# 2-5 e30        # 12 states: 0°, 30°, 60°, 90°, 120°, 150°, 180°, 210°, 240°, 270°, 300°, 330°
# 3-7 e90        # 4 states: 0°, 90°, 180°, 270°
# 4-8 e120       # 3 states: 0°, 120°, 240°

# EXPLICIT ANGLE ROTATIONS (specify exact angles)
# Format: atom1-atom2 angle1 angle2 angle3 ...
# Examples:
# 1-6 0 60 120 180                    # 4 specific angles
# 2-5 -120 -60 0 60 120               # 5 angles including negatives
# 3-7 30 45 60 90 120 150 180         # 7 custom angles

# SYNCHRONOUS ROTATIONS (coordinate with other bonds)
# Format: atom1-atom2 syn <reference_bond_number>
# Examples:
# 3-7 syn 1      # Bond 3 rotates with bond 1 (same direction)
# 5-9 syn 2      # Bond 5 rotates with bond 2 (same direction)
# 4-8 syn -1     # Bond 4 rotates opposite to bond 1
# 6-10 syn -2    # Bond 6 rotates opposite to bond 2

# =======================================================================
# SCANNING MODE SPECIFICATIONS
# =======================================================================
# Use these for bond length scanning conformer generation
# Format: atom1-atom2 scan steps step_size
# - steps: number of scanning steps (positive integer)
# - step_size: increment in Angstroms (positive = stretch, negative = compress)

# SINGLE BOND SCANNING
# Examples:
# 1-2 scan 10 0.1     # Stretch C-C bond in 10 steps of 0.1 Å each
# 3-4 scan 5 -0.05    # Compress bond in 5 steps of 0.05 Å each
# 5-6 s 15 0.2        # Alternative syntax: 's' instead of 'scan'

# MULTI-DIMENSIONAL SCANNING (multiple bonds simultaneously)
# Examples:
# 1-2 scan 8 0.1      # First bond: 8 steps, +0.1 Å
# 3-4 scan 6 -0.05    # Second bond: 6 steps, -0.05 Å
# Total combinations: 8 × 6 = 48 conformers

# -----------------------------------------------------------------------
# EXAMPLE CONFIGURATIONS FOR COMMON MOLECULES
# -----------------------------------------------------------------------

# ROTATION MODE EXAMPLES:

# ETHANE (simple single bond rotation)
# 1-2 e60

# BUTANE (multiple independent bonds)
# 1-2 e120       # C1-C2 bond
# 2-3 e60        # C2-C3 bond
# 3-4 e180       # C3-C4 bond

# TOLUENE (aromatic + methyl rotation)
# 7-8 e60        # Methyl group rotation

# P-XYLENE (coordinated methyl rotations)
# 1-7 e60        # First methyl (independent)
# 8-9 syn 1      # Second methyl (same as first)
# 10-11 syn -1   # Third methyl (opposite to first)

# METAL COMPLEX (with manual bonds)
# 1-25 bond      # Force metal-ligand bond
# 1-30 bond      # Force another metal-ligand bond
# 15-20 nobond   # Remove false positive
# 2-3 e90        # Ligand rotation
# 4-5 syn 1      # Coordinated ligand rotation

# SCANNING MODE EXAMPLES:

# ETHANE (C-C bond length scanning)
# 1-2 scan 10 0.1    # Stretch C-C bond from 1.5 to 2.5 Å

# HYDROGEN MOLECULE (H-H bond compression and stretching)
# 1-2 scan 20 0.05   # Fine scanning of H-H bond

# WATER MOLECULE (O-H bond scanning)
# 1-2 scan 8 0.1     # First O-H bond
# 1-3 scan 8 0.1     # Second O-H bond (2D scanning)
# Total: 8 × 8 = 64 conformers

# METAL COMPLEX (with forced bonds and scanning)
# 1-25 bond          # Force metal-ligand bond
# 1-30 bond          # Force another metal-ligand bond
# 15-20 nobond       # Remove false positive
# 1-25 scan 6 0.2    # Scan first metal-ligand bond
# 1-30 scan 6 -0.1   # Scan second metal-ligand bond (compress)

# -----------------------------------------------------------------------
# NOTES
# -----------------------------------------------------------------------
# - Atom indices are 1-based (first atom is 1, not 0)
# - Use dashes (-) not hyphens or em-dashes in bond specifications
# - Rotation angles are in degrees, range: -360° to +360°
# - Scanning step sizes are in Angstroms
# - Comments start with '#' and continue to end of line
# - Empty lines are ignored
# - Parameters can be in any order
# - CANNOT mix rotation and scanning in the same file

# -----------------------------------------------------------------------
# QUICK START
# -----------------------------------------------------------------------

# FOR ROTATION MODE:
# 1. Identify rotatable bonds in your molecule
# 2. Uncomment rotation examples above
# 3. Start with e60 or e90 for initial exploration
# 4. Run: rotbond <molecule_name>
# 5. Check success rate and adjust parameters if needed

# FOR SCANNING MODE:
# 1. Identify bonds to scan in your molecule
# 2. Uncomment scanning examples above
# 3. Start with 5-10 steps and small step sizes (0.05-0.1 Å)
# 4. Run: rotbond <molecule_name>
# 5. Check bond lengths are chemically reasonable

# DELETE THIS LINE AND ADD YOUR SPECIFICATIONS BELOW:
# (The program requires at least one specification to run)
"#;

    let mut file = File::create(filename)?;
    file.write_all(template_content.as_bytes())?;
    Ok(())
}

/// Calculates the theoretical number of conformers from angle sets.
///
/// Computes the cartesian product of all angle sets to determine
/// the total number of conformer combinations that would be generated.
///
/// # Arguments
///
/// * `angle_sets` - Vector of angle vectors for each rotatable bond
///
/// # Returns
///
/// Total theoretical conformer count
fn calculate_theoretical_conformers(angle_sets: &[Vec<f64>]) -> usize {
    let mut total = 1;
    for angles in angle_sets {
        if !angles.is_empty() {
            total *= angles.len();
        }
    }
    total
}

/// Calculates the theoretical number of conformers from scanning specifications.
///
/// Computes the cartesian product of all scanning steps to determine
/// the total number of conformer combinations that would be generated.
///
/// # Arguments
///
/// * `rotation_specs` - Vector of rotation specifications (containing scanning specs)
///
/// # Returns
///
/// Total theoretical conformer count for scanning mode
fn calculate_theoretical_scanning_conformers(rotation_specs: &[molecule::RotationSpec]) -> usize {
    let mut total = 1;
    for spec in rotation_specs {
        if let molecule::RotationSpec::Scanning { steps, .. } = spec {
            total *= steps;
        }
    }
    total
}

/// Main entry point for the Rotbond conformer generation tool.
///
/// Parses command line arguments, loads molecular structure and rotation parameters,
/// then generates all possible conformers using dihedral angle rotation.
/// Supports both `molecule` and `molecule.xyz` input formats and automatically
/// creates comprehensive .rp template files when they don't exist.
///
/// # Process Flow
///
/// 1. Parse command line arguments for molecule name and options (supports .xyz extension)
/// 2. Check for XYZ structure file existence
/// 3. Create comprehensive .rp template if rotation parameter file doesn't exist
/// 4. Load rotation parameters and detect rotatable bonds using covalent radii
/// 5. Generate conformers by systematic rotation with safety limits
/// 6. Apply steric clash filtering and interactive warnings for large jobs
/// 7. Write output files with proper formatting and display book recommendations
///
/// # Exit Codes
///
/// - `0`: Success - all conformers generated successfully
/// - `1`: Error - invalid arguments, file not found, or generation failure
fn main() {
    // Parse command line arguments
    let args: Vec<String> = env::args().collect();

    // Handle help and version flags before any processing
    if args.len() == 2 {
        let arg = &args[1];

        // Version flag
        if arg == "--version" || arg == "-v" || arg == "version" {
            print_version();
            process::exit(0);
        }

        // Help flag
        if arg == "--help" || arg == "-h" || arg == "help" {
            print_basic_help();
            process::exit(0);
        }

        // List help topics
        if arg == "--help" && args.len() == 2 {
            print_basic_help();
            process::exit(0);
        }
    }

    // Handle topic-specific help
    if args.len() == 3 {
        let arg2 = &args[2];

        if arg2 == "topics" || arg2 == "topic" {
            print_help_topics();
            process::exit(0);
        }

        // Topic-specific help
        if &args[1] == "--help" || &args[1] == "-h" {
            print_help_topic(arg2);
            process::exit(0);
        }
    }

    // Handle --help <topic> format
    if args.len() == 3 && (&args[1] == "--help" || &args[1] == "-h") {
        print_help_topic(&args[2]);
        process::exit(0);
    }

    // Parse standard arguments (molecule name)
    if args.len() != 2 {
        print_usage();
        eprintln!("\nFor help, use: --help or --help topics");
        process::exit(1);
    }

    // Display program banner
    print_version();
    println!();

    // Handle both "molecule" and "molecule.xyz" formats
    let base_name = if args[1].ends_with(".xyz") {
        &args[1][..args[1].len() - 4] // Remove .xyz extension
    } else {
        &args[1]
    };

    let xyz_filename = format!("{}.xyz", base_name);
    let rp_filename = format!("{}.rp", base_name);

    // Check if XYZ file exists
    if !std::path::Path::new(&xyz_filename).exists() {
        eprintln!("ERROR: Input file '{}' not found", xyz_filename);
        eprintln!("\nFor help creating input files, use: --help input");
        process::exit(1);
    }

    // Check if RP file exists, create template if not
    if !std::path::Path::new(&rp_filename).exists() {
        println!("Rotation parameters file '{}' not found.", rp_filename);
        println!("Creating template file with default parameters...");

        match create_template_rp_file(&rp_filename) {
            Ok(()) => {
                println!("✓ Template file '{}' created successfully!", rp_filename);
                println!(
                    "\nPlease edit '{}' to specify your rotation parameters, then run rotbond again.",
                    rp_filename
                );
                println!("\nFor help with rotation parameters, use: --help input");
                process::exit(0);
            }
            Err(e) => {
                eprintln!(
                    "ERROR: Could not create template file '{}': {}",
                    rp_filename, e
                );
                eprintln!("\nFor help creating rotation files, use: --help input");
                process::exit(1);
            }
        }
    }

    // Read XYZ file
    // Reading molecule...
    let (atoms, _comment) = match io::read_xyz_file(&xyz_filename) {
        Ok(result) => result,
        Err(e) => {
            eprintln!("ERROR reading XYZ file: {}", e);
            eprintln!("\nFor XYZ format help, use: --help input");
            process::exit(1);
        }
    };
    // Molecule loaded

    // Parse rotation/scanning parameters file
    // Parsing parameters...
    let (
        bond_factor,
        skip_factor,
        forced_bonds,
        forbidden_bonds,
        rotation_specs,
        conformer_config,
        operation_mode,
    ) = match parse_rotation_file(&rp_filename) {
        Ok(result) => result,
        Err(e) => {
            eprintln!("ERROR parsing parameters: {}", e);
            // Check if error message mentions scanning to provide appropriate help
            if e.to_string().to_lowercase().contains("scan") {
                eprintln!("\nFor scanning file format help, use: --help input");
            } else {
                eprintln!("\nFor rotation file format help, use: --help input");
            }
            process::exit(1);
        }
    };

    // Use bond_factor as specified in file or default (1.0)

    // Rotation specifications parsed

    print_configuration(bond_factor, skip_factor);

    // Validate specifications based on operation mode
    match operation_mode {
        molecule::OperationMode::Rotation => {
            // Validate synchronous references for rotation mode
            if let Err(e) = validate_synchronous_references(&rotation_specs) {
                eprintln!("ERROR validating synchronous references: {}", e);
                eprintln!("\nFor synchronous rotation help, use: --help features");
                process::exit(1);
            }
        }
        molecule::OperationMode::Scanning => {
            // Validate scanning specifications with integration to existing validation systems
            for (i, spec) in rotation_specs.iter().enumerate() {
                if let molecule::RotationSpec::Scanning {
                    atom1: _,
                    atom2: _,
                    steps,
                    step_size,
                } = spec
                {
                    if *steps == 0 {
                        eprintln!("ERROR: Scanning bond {} has zero steps", i + 1);
                        eprintln!("\nFor scanning syntax help, use: --help input");
                        process::exit(1);
                    }
                    if *step_size == 0.0 {
                        eprintln!("ERROR: Scanning bond {} has zero step size", i + 1);
                        eprintln!("\nFor scanning syntax help, use: --help input");
                        process::exit(1);
                    }

                    // Enhanced validation with bond_factor and skip_factor integration
                    // Warn about very large step sizes that might create unrealistic bond lengths
                    if step_size.abs() > 2.0 {
                        println!(
                            "WARNING: Scanning bond {} has large step size ({:.2} Å)",
                            i + 1,
                            step_size
                        );
                        println!("         This may create unrealistic bond lengths.");
                        println!(
                            "         Consider the bond_factor ({:.2}) and skip_factor ({:.2}) settings.",
                            bond_factor, skip_factor
                        );
                    }

                    // Validate scanning limits with maxgen parameter integration
                    let total_change = *steps as f64 * step_size.abs();
                    if total_change > 4.0 {
                        println!(
                            "WARNING: Scanning bond {} would change length by {:.2} Å total",
                            i + 1,
                            total_change
                        );
                        println!(
                            "         Large changes may violate steric constraints (skip_factor = {:.2})",
                            skip_factor
                        );
                        println!(
                            "         Consider reducing steps or step_size for better results."
                        );
                    }
                }
            }

            // Validate scanning parameters for critical errors only (overflow, invalid params)
            // Note: Conformer limit checking is handled later with interactive prompts
            if let Err(e) = errors::validate_scanning_specs(&rotation_specs) {
                eprintln!("ERROR: {}", errors::format_user_error(&e));
                process::exit(1);
            }

            // Validate forced and forbidden bonds compatibility with scanning
            let scanning_bond_indices: Vec<usize> = rotation_specs
                .iter()
                .enumerate()
                .filter_map(|(i, spec)| {
                    if matches!(spec, molecule::RotationSpec::Scanning { .. }) {
                        Some(i)
                    } else {
                        None
                    }
                })
                .collect();

            // Create a temporary molecule to validate forced/forbidden bonds
            let temp_molecule = molecule::Molecule {
                atoms: atoms.clone(),
                bonds: Vec::new(), // We don't need bonds for this validation
                fragments: Vec::new(),
                bond_factor,
                skip_factor,
                forced_bonds: forced_bonds.clone(),
                forbidden_bonds: forbidden_bonds.clone(),
                operation_mode: operation_mode.clone(),
            };

            if let Err(e) =
                errors::validate_forced_forbidden_bonds(&temp_molecule, &scanning_bond_indices)
            {
                eprintln!("ERROR: {}", errors::format_user_error(&e));
                process::exit(1);
            }
        }
    }

    // Generate angle sets (only for rotation mode)
    let angle_sets = match operation_mode {
        molecule::OperationMode::Rotation => {
            // Generating angle sets...
            match generate_angle_sets(rotation_specs.clone()) {
                Ok(result) => result,
                Err(e) => {
                    eprintln!("ERROR generating angle sets: {}", e);
                    eprintln!("\nFor rotation syntax help, use: --help features");
                    process::exit(1);
                }
            }
        }
        molecule::OperationMode::Scanning => {
            // Scanning mode doesn't use angle sets - create empty vector
            Vec::new()
        }
    };

    // Calculate theoretical conformer count and check limits
    let theoretical_count = match operation_mode {
        molecule::OperationMode::Rotation => calculate_theoretical_conformers(&angle_sets),
        molecule::OperationMode::Scanning => {
            calculate_theoretical_scanning_conformers(&rotation_specs)
        }
    };

    // Check if we need to warn about large generation jobs
    let mut final_max_conformers = conformer_config.max_conformers;

    if theoretical_count > 500 && !conformer_config.auto_confirm {
        match operation_mode {
            molecule::OperationMode::Rotation => {
                println!("\n   WARNING: Large rotation job detected!");
                println!("   Theoretical conformers: {}", theoretical_count);
                println!("   This may consume significant memory and processing time.");
            }
            molecule::OperationMode::Scanning => {
                println!("\n   WARNING: Large scanning job detected!");
                println!("   Theoretical conformers: {}", theoretical_count);
                println!("   This may consume significant memory and processing time.");
                println!("   Consider reducing the number of scanning steps or bonds.");
            }
        }

        if let Some(max_limit) = conformer_config.max_conformers {
            println!("   Current limit: {} conformers", max_limit);
        } else {
            println!(
                "   Current limit: unlimited (all {} conformers)",
                theoretical_count
            );
        }

        println!("\n   Options:");
        println!("   - Press Enter or 'y' to continue with current limit");
        println!("   - Enter a number (e.g., 300) to set a new limit");
        println!("   - Enter 'max' for unlimited generation");
        println!("   - Enter 'n' to cancel");
        print!("\n   Your choice: ");
        std::io::stdout().flush().unwrap();

        let mut input = String::new();
        match stdin().read_line(&mut input) {
            Ok(_) => {
                let response = input.trim();

                if response.is_empty()
                    || response.to_lowercase() == "y"
                    || response.to_lowercase() == "yes"
                {
                    // Continue with current limit
                    println!("✓ Continuing with current limit.");
                } else if response.to_lowercase() == "n" || response.to_lowercase() == "no" {
                    println!("Generation cancelled by user.");
                    process::exit(0);
                } else if response.to_lowercase() == "max"
                    || response.to_lowercase() == "maximum"
                    || response.to_lowercase() == "unlimited"
                {
                    final_max_conformers = None;
                    println!(
                        "✓ Limit set to unlimited - all {} conformers will be generated.",
                        theoretical_count
                    );
                } else if let Ok(new_limit) = response.parse::<usize>() {
                    if new_limit == 0 {
                        println!("ERROR: Limit must be greater than 0");
                        process::exit(1);
                    }
                    final_max_conformers = Some(new_limit);
                    println!("✓ Limit set to {} conformers.", new_limit);
                } else {
                    println!(
                        "ERROR: Invalid input '{}'. Please enter a number, 'max', 'y', or 'n'.",
                        response
                    );
                    process::exit(1);
                }
            }
            Err(_) => {
                eprintln!("ERROR: Failed to read user input");
                process::exit(1);
            }
        }
    } else if theoretical_count > 500 {
        match operation_mode {
            molecule::OperationMode::Rotation => {
                println!(
                    "INFO: Large rotation job ({}), auto-confirmed by configuration.",
                    theoretical_count
                );
            }
            molecule::OperationMode::Scanning => {
                println!(
                    "INFO: Large scanning job ({}), auto-confirmed by configuration.",
                    theoretical_count
                );
            }
        }
        if let Some(max_limit) = conformer_config.max_conformers {
            println!(
                "      Generation will be limited to {} conformers.",
                max_limit
            );
        }
    }

    // Create molecule
    let mut molecule = Molecule::new(atoms);
    molecule.bond_factor = bond_factor;
    molecule.skip_factor = skip_factor;
    molecule.forced_bonds = forced_bonds;
    molecule.forbidden_bonds = forbidden_bonds;
    molecule.operation_mode = operation_mode.clone();

    // Create bonds from specifications
    match operation_mode {
        molecule::OperationMode::Rotation => {
            // Creating rotation bonds...
            for (i, spec) in rotation_specs.iter().enumerate() {
                match spec {
                    RotationSpec::Step { step: _ } | RotationSpec::Explicit(_) => {
                        // These will be handled by angle_sets
                        let angles = if i < angle_sets.len() {
                            angle_sets[i].clone()
                        } else {
                            Vec::new()
                        };
                        molecule.add_bond(Bond {
                            atom1: 0,
                            atom2: 0,
                            angles,
                            is_synchronous: false,
                            reference_bond: None,
                            direction: 1.0,
                        });
                    }
                    RotationSpec::Synchronous {
                        reference,
                        direction,
                    } => {
                        molecule.add_bond(Bond {
                            atom1: 0,
                            atom2: 0,
                            angles: Vec::new(), // Synchronous bonds use reference bond's angles
                            is_synchronous: true,
                            reference_bond: Some(*reference),
                            direction: *direction,
                        });
                    }
                    RotationSpec::Scanning { .. } => {
                        // This shouldn't happen in rotation mode due to mutual exclusion
                        eprintln!("ERROR: Scanning specification found in rotation mode");
                        process::exit(1);
                    }
                }
            }
        }
        molecule::OperationMode::Scanning => {
            // Scanning bonds are created in generate.rs during conformer generation
            // No need to create bonds here for scanning mode
        }
    }

    // Re-parse rotation file to get atom indices (we need them for bonds)
    // This is a bit of a hack - we should refactor to get atom indices during parsing
    let atom_pairs = extract_atom_pairs_from_rp(&rp_filename).unwrap_or_else(|e| {
        eprintln!("ERROR: {}", e);
        process::exit(1);
    });

    // Update bonds with atom indices
    for (i, (atom1, atom2)) in atom_pairs.iter().enumerate() {
        if i < molecule.bonds.len() {
            molecule.bonds[i].atom1 = *atom1 - 1; // Convert to 0-indexed
            molecule.bonds[i].atom2 = *atom2 - 1;
        }
    }

    // Print summary based on operation mode
    match operation_mode {
        molecule::OperationMode::Rotation => {
            print_rotation_summary(&molecule.bonds, &angle_sets);
        }
        molecule::OperationMode::Scanning => {
            print_scanning_summary(&molecule.bonds, &rotation_specs);
        }
    }

    // Validate bonds (only for rotation mode - scanning bonds are created in generate.rs)
    if molecule.bonds.is_empty() && matches!(operation_mode, molecule::OperationMode::Rotation) {
        eprintln!("ERROR: No rotation bonds specified");
        eprintln!("\nFor rotation syntax help, use: --help features");
        process::exit(1);
    }

    for (i, bond) in molecule.bonds.iter().enumerate() {
        if bond.atom1 >= molecule.atoms.len() || bond.atom2 >= molecule.atoms.len() {
            eprintln!(
                "ERROR: Bond {} references non-existent atoms ({}, {})",
                i + 1,
                bond.atom1 + 1,
                bond.atom2 + 1
            );
            eprintln!("\nFor atom indexing help, use: --help input");
            process::exit(1);
        }
    }

    // Generate conformers
    match operation_mode {
        molecule::OperationMode::Rotation => {
            println!("Generating conformers through bond rotation...");
        }
        molecule::OperationMode::Scanning => {
            println!("Generating conformers through bond scanning...");
        }
    }
    println!();

    let (total_combinations, valid_conformers, min_length_violations) = match generate_conformers(
        &mut molecule,
        &angle_sets,
        &rotation_specs,
        base_name,
        final_max_conformers,
    ) {
        Ok((total, valid, min_violations)) => (total, valid, min_violations),
        Err(e) => {
            eprintln!("ERROR generating conformers: {}", e);
            eprintln!("\nFor troubleshooting help, use: --help troubleshoot");
            process::exit(1);
        }
    };

    // Print summary
    print_summary(total_combinations, valid_conformers, min_length_violations);

    println!("\nConformer generation completed successfully!");

    // Display a random book recommendation
    print_random_book();
}

/// Extracts atom pairs from rotation parameter file.
///
/// This helper function parses the .rp file to extract the atom indices
/// for each rotatable bond specification. It handles various rotation
/// formats including step-based, explicit angles, and synchronous rotations.
///
/// # Arguments
///
/// * `filename` - Path to the .rp parameter file
///
/// # Returns
///
/// A vector of tuples containing atom index pairs (1-based indexing as in file).
///
/// # Errors
///
/// Returns an error if:
/// - File cannot be read
/// - Atom indices cannot be parsed as numbers
/// - Invalid bond specification format
///
/// # Examples
///
/// ```rust
/// let pairs = extract_atom_pairs_from_rp("molecule.rp")?;
/// for (atom1, atom2) in pairs {
///     println!("Bond between atoms {} and {}", atom1, atom2);
/// }
/// ```
fn extract_atom_pairs_from_rp(
    filename: &str,
) -> Result<Vec<(usize, usize)>, Box<dyn std::error::Error>> {
    let contents = std::fs::read_to_string(filename)?;
    let mut atom_pairs = Vec::new();

    for line in contents.lines() {
        let line = line.trim();

        // Skip empty lines and comments
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let parts: Vec<&str> = line.split_whitespace().collect();

        if parts.len() < 2 {
            continue;
        }

        // Skip configuration parameters
        if line.starts_with("bond_factor") || line.starts_with("skip_factor") {
            continue;
        }

        // Parse bond specification (atom1-atom2 format)
        let bond_spec = parts[0];
        let atom_pair: Vec<&str> = bond_spec.split('-').collect();

        if atom_pair.len() != 2 {
            continue;
        }

        // Check if this is a valid rotation specification or manual bond
        // Valid lines are:
        // - atom1-atom2 eXX (step-based)
        // - atom1-atom2 X Y Z (explicit angles)
        // - atom1-atom2 syn N (synchronous)
        // - atom1-atom2 bond (manual bond)
        // - atom1-atom2 nobond (manual bond removal)

        let is_rotation_or_bond = if parts[1] == "syn" || parts[1] == "bond" || parts[1] == "nobond"
        {
            // Synchronous or manual bond - valid
            parts.len() >= 3
        } else if parts[1].starts_with('e') {
            // Step-based rotation - valid
            true
        } else {
            // Explicit angles or manual bond - check if all remaining parts are numbers
            parts[1..].iter().all(|p| p.parse::<f64>().is_ok())
        };

        if is_rotation_or_bond {
            let atom1 = atom_pair[0].parse::<usize>()?;
            let atom2 = atom_pair[1].parse::<usize>()?;
            atom_pairs.push((atom1, atom2));
        }
    }

    Ok(atom_pairs)
}
