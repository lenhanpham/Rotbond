//! # Utility Functions for Output Formatting and Display
//!
//! This module provides utility functions for formatting and displaying
//! information about molecular conformer generation processes.
//!
//! ## Key Functions
//!
//! - **Configuration display**: Shows bond detection and validation parameters
//! - **Rotation specification display**: Lists all rotatable bonds and their angles
//! - **Summary statistics**: Reports generation success rates and conformer counts
//! - **Usage information**: Provides help and usage instructions
//!
//! ## Output Format
//!
//! All output is designed to be clear and professional, suitable for both
//! interactive use and log file analysis.

/// Prints a separator line for output formatting.
/// 
/// Creates a visual separator using 80 equal signs to improve readability
/// of console output by clearly delineating different sections.
/// 
/// # Examples
/// 
/// ```rust
/// use rotbond::utils::print_separator;
/// 
/// print_separator();
/// // Output: =======================================================================
/// ```
pub fn print_separator() {
    println!("{}", "=".repeat(71));
}

/// Prints the program header with version information.
/// 
/// Displays the program name, version, and decorative separators in the
/// official Rotbond styling format. This provides a consistent header
/// that matches the main program output.
/// 
/// # Examples
/// 
/// ```rust
/// use rotbond::utils::print_header;
/// 
/// print_header();
/// // Output:
/// // ***********************************************************************
/// //                                 ROTBOND
/// // ***********************************************************************
/// // # -----------------------------------------------------------------------#
/// // # Version {CARGO_PKG_VERSION}  Release date: 2025                        #
/// // # Developer: Le Nhan Pham                                                #
/// // # https://github.com/lenhanpham/Rotbond                                  #
/// // # -----------------------------------------------------------------------#
/// ```
#[allow(dead_code)]
pub fn print_header() {
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
    println!();
    println!("Please cite this project if you use Rotbond for your research");
    println!();
    println!();
    println!("# -----------------------------------------------------------------------#");
    println!("# L.N Pham, \"Rotbond - A Systematic molecular conformer generator\"       #");
    println!("# https://github.com/lenhanpham/Rotbond                                  #");
    println!("# -----------------------------------------------------------------------#");
    println!();
}

/// Prints comprehensive usage information and command-line help.
/// 
/// Displays detailed information about command-line options, required files,
/// output files, and usage examples. This is the main help function shown
/// when users provide invalid arguments or request help.
/// 
/// # Output Sections
/// 
/// - **Options**: Command-line flags and their descriptions
/// - **Arguments**: Required positional arguments
/// - **Required Files**: Input file format descriptions
/// - **Output Files**: Generated file descriptions
/// - **Examples**: Common usage patterns
/// 
/// # Examples
/// 
/// ```rust
/// use rotbond::utils::print_usage;
/// 
/// print_usage();
/// // Displays comprehensive help information
/// ```
pub fn print_usage() {
    println!("Usage: rotbond [OPTIONS] <molecule_name>");
    println!();
    println!("OPTIONS:");
    println!("  -h, --help              Show this help message");
    println!("  -v, --version           Show version information");
    println!("  --help topics           List all available help topics");
    println!("  --help <topic>          Show specific help topic");
    println!();
    println!("ARGUMENTS:");
    println!("  <molecule_name>         Base name for input/output files");
    println!("                          (e.g., 'ethane' uses ethane.xyz and ethane.rp)");
    println!();
    println!("REQUIRED FILES:");
    println!("  <molecule_name>.xyz    - Input molecule structure (XYZ format)");
    println!("  <molecule_name>.rp     - Rotation parameters file");
    println!();
    println!("OUTPUT FILES:");
    println!("  <molecule_name>_traj.xyz   - All conformers in one trajectory file");
    println!("  <molecule_name>_NN.xyz     - Individual conformer files (smart padding)");
    println!();
    println!("EXAMPLES:");
    println!("  rotbond molecule_a         # Basic usage");
    println!("  rotbond --help           # Show help");
    println!("  rotbond --help examples  # Show examples");
    println!("  rotbond --version        # Show version");
    println!();
    println!("For more information, use: --help topics");
}

/// Prints a summary of rotation specifications in a clear, formatted display.
/// 
/// Shows all rotatable bonds and their associated rotation angles,
/// providing a complete overview of the conformational search space.
/// Handles both independent and synchronous bond rotations.
/// 
/// # Arguments
/// 
/// * `bonds` - Slice of bond specifications
/// * `angle_sets` - Corresponding angle sets for each bond
/// 
/// # Output Format
/// 
/// For each bond, displays:
/// - Bond number and atom indices (1-based for user clarity)
/// - Number of rotation states
/// - List of rotation angles
/// - Synchronous bond relationships if applicable
/// 
/// # Examples
/// 
/// ```rust
/// use rotbond::utils::print_rotation_summary;
/// 
/// print_rotation_summary(&bonds, &angle_sets);
/// // Output:
/// // Rotation specifications:
/// //   Bond 1: atoms 1-2 (6 rotation states)
/// //            Angles: [0.0, 60.0, 120.0, 180.0, 240.0, 300.0]
/// //   Bond 2: atoms 3-4 (synchronous with bond 1)
/// ```
pub fn print_rotation_summary(bonds: &[crate::molecule::Bond], angle_sets: &[Vec<f64>]) {
    println!("Rotation specifications:");
    for (i, (bond, angles)) in bonds.iter().zip(angle_sets.iter()).enumerate() {
        let atom1_idx = bond.atom1 + 1; // Convert to 1-indexed for display
        let atom2_idx = bond.atom2 + 1;

        if bond.is_synchronous {
            println!("  Bond {}: atoms {}-{} (synchronous with bond {})",
                     i + 1, atom1_idx, atom2_idx, bond.reference_bond.unwrap_or(0));
        } else if angles.is_empty() {
            println!("  Bond {}: atoms {}-{} (no rotations)",
                     i + 1, atom1_idx, atom2_idx);
        } else {
            let num_states = angles.len();
            println!("  Bond {}: atoms {}-{} ({} rotation states)",
                     i + 1, atom1_idx, atom2_idx, num_states);
            println!("           Angles: {:?}", angles);
        }
    }
    println!();
}

/// Prints the current configuration parameters in a formatted display.
/// 
/// Shows the bond detection and steric clash parameters being used
/// for the current conformer generation run.
/// 
/// # Arguments
/// 
/// * `bond_factor` - Bond detection threshold multiplier
/// * `skip_factor` - Steric clash detection threshold
/// 
/// # Examples
/// 
/// ```rust
/// use rotbond::utils::print_configuration;
/// 
/// print_configuration(1.0, 0.7);
/// // Output:
/// // Configuration:
/// //   Bond factor:   1.00
/// //   Skip factor:   0.70
/// ```
pub fn print_configuration(bond_factor: f64, skip_factor: f64) {
    println!("Configuration:");
    println!("  Bond factor:   {:.2}", bond_factor);
    println!("  Skip factor:   {:.2}", skip_factor);
    println!();
}

/// Prints a comprehensive summary of the conformer generation results.
/// 
/// Displays statistics including total combinations attempted, number of
/// valid conformers generated, rejected conformers, and success rate percentage.
/// Uses visual separators for clear presentation.
/// 
/// # Arguments
/// 
/// * `total_theoretical` - Total number of conformer combinations attempted
/// * `valid_conformers` - Number of conformers that passed validation
/// 
/// # Examples
/// 
/// ```rust
/// use rotbond::utils::print_summary;
/// 
/// print_summary(24, 21);
/// // Output:
/// // =======================================================================
/// // Generation Summary:
/// //   Theoretical combinations: 24
/// //   Rejected (steric clashes): 3
/// //   Valid conformers generated: 21
/// //   Success rate: 87.5%
/// // =======================================================================
/// ```
pub fn print_summary(total_theoretical: usize, valid_conformers: usize) {
    print_separator();
    println!("Generation Summary:");
    println!("  Theoretical combinations: {}", total_theoretical);
    
    if valid_conformers < total_theoretical {
        let rejected = total_theoretical - valid_conformers;
        println!("  Rejected (steric clashes): {}", rejected);
        println!("  Valid conformers generated: {}", valid_conformers);
        println!("  Success rate: {:.1}%",
                 (valid_conformers as f64 / total_theoretical as f64) * 100.0);
    } else {
        println!("  Valid conformers generated: {}", valid_conformers);
        println!("  Success rate: 100.0%");
    }
    print_separator();
}

/// Prints an error message and exits the program with status code 1.
/// 
/// This is a utility function for fatal error handling that ensures
/// consistent error message formatting and proper program termination.
/// 
/// # Arguments
/// 
/// * `msg` - The error message to display
/// 
/// # Examples
/// 
/// ```rust
/// use rotbond::utils::print_error_and_exit;
/// 
/// print_error_and_exit("File not found");
/// // Output: ERROR: File not found
/// // Program exits with code 1
/// ```
#[allow(dead_code)]
pub fn print_error_and_exit(msg: &str) -> ! {
    eprintln!("ERROR: {}", msg);
    std::process::exit(1);
}
