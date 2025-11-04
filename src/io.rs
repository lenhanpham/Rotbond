//! # Input/Output Operations for Molecular Files
//!
//! This module handles reading and writing molecular structure files,
//! with a focus on the XYZ format commonly used in computational chemistry.
//!
//! ## Supported Formats
//!
//! - **XYZ Format**: Standard molecular coordinate format
//!   - First line: number of atoms
//!   - Second line: comment/title
//!   - Remaining lines: element symbol and x, y, z coordinates
//!
//! ## Output Features
//!
//! - **Professional formatting**: Properly aligned columns for readability
//! - **Individual files**: Separate XYZ file for each conformer
//! - **Trajectory files**: All conformers in a single multi-structure file
//! - **Smart numbering**: Automatic padding for consistent file ordering

use crate::molecule::Atom;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

/// Loads a molecular structure from an XYZ format file.
///
/// Parses the standard XYZ format where the first line contains the number
/// of atoms, the second line is a comment, and subsequent lines contain
/// element symbols and coordinates.
///
/// # Arguments
///
/// * `filename` - Path to the XYZ file to load
///
/// # Returns
///
/// A tuple containing:
/// - Vector of [`Atom`] structures representing the molecule
/// - Comment string from the second line of the file
///
/// # Errors
///
/// Returns an error if:
/// - File cannot be opened or read
/// - File format is invalid (wrong number of atoms, malformed lines)
/// - Coordinates cannot be parsed as floating-point numbers
///
/// # Examples
///
/// ```rust
/// use rotbond::io::read_xyz_file;
///
/// let (atoms, comment) = read_xyz_file("molecule.xyz")?;
/// println!("Loaded {} atoms: {}", atoms.len(), comment);
/// ```
///
/// # XYZ Format
///
/// ```text
/// 3
/// Water molecule
/// O    0.000000    0.000000    0.000000
/// H    0.757000    0.586000    0.000000
/// H   -0.757000    0.586000    0.000000
/// ```
pub fn read_xyz_file(filename: &str) -> Result<(Vec<Atom>, String), Box<dyn std::error::Error>> {
    let file = File::open(filename).map_err(|e| format!("Failed to open {}: {}", filename, e))?;

    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    // Read number of atoms
    let nat_line = lines
        .next()
        .ok_or("Empty file")?
        .map_err(|e| format!("Failed to read atom count: {}", e))?;

    let nat = nat_line
        .trim()
        .parse::<usize>()
        .map_err(|e| format!("Invalid atom count: {}", e))?;

    // Read comment line
    let comment = lines
        .next()
        .ok_or("Missing comment line")?
        .map_err(|e| format!("Failed to read comment: {}", e))?;

    // Read atom coordinates
    let mut atoms = Vec::with_capacity(nat);
    for i in 0..nat {
        let line = lines
            .next()
            .ok_or(format!("Expected {} atoms, found fewer", nat))?
            .map_err(|e| format!("Failed to read atom {}: {}", i + 1, e))?;

        let parts: Vec<&str> = line.split_whitespace().collect();

        if parts.len() < 4 {
            return Err(format!(
                "Invalid atom format on line {}: expected 4 columns (element x y z), got {}",
                i + 1,
                parts.len()
            )
            .into());
        }

        let element = parts[0].to_string();
        let x = parts[1]
            .parse::<f64>()
            .map_err(|e| format!("Invalid x coordinate on line {}: {}", i + 1, e))?;
        let y = parts[2]
            .parse::<f64>()
            .map_err(|e| format!("Invalid y coordinate on line {}: {}", i + 1, e))?;
        let z = parts[3]
            .parse::<f64>()
            .map_err(|e| format!("Invalid z coordinate on line {}: {}", i + 1, e))?;

        atoms.push(Atom::new(element, x, y, z));
    }

    Ok((atoms, comment))
}

/// Writes a single conformer to an individual XYZ format file.
///
/// Creates a properly formatted XYZ file with professional column alignment.
/// The output uses 6 decimal places for coordinates and ensures consistent
/// spacing for easy reading and processing by other tools.
///
/// # Arguments
///
/// * `filename` - Path where the XYZ file should be written
/// * `atoms` - Slice of atoms representing the conformer
///
/// # Returns
///
/// `Ok(())` on success, or an error if the file cannot be created or written.
///
/// # Errors
///
/// Returns an error if:
/// - File cannot be created (permissions, disk space, etc.)
/// - Write operations fail
///
/// # Format
///
/// The output format uses professional alignment:
/// - Element symbols: left-aligned, 2 characters wide
/// - Coordinates: right-aligned, 15 characters wide, 6 decimal places
///
/// # Examples
///
/// ```rust
/// use rotbond::io::write_individual_xyz_file;
///
/// write_individual_xyz_file("conformer_1.xyz", &atoms)?;
/// ```
pub fn write_individual_xyz_file(
    filename: &str,
    atoms: &[Atom],
) -> Result<(), Box<dyn std::error::Error>> {
    let mut file =
        File::create(filename).map_err(|e| format!("Failed to create {}: {}", filename, e))?;

    // Write header
    writeln!(file, "{}", atoms.len())?;
    writeln!(file, "Generated by Rotbond")?;

    // Write atoms with proper column alignment
    for atom in atoms {
        writeln!(
            file,
            "{:<2} {:>15.6} {:>15.6} {:>15.6}",
            atom.element, atom.x, atom.y, atom.z
        )?;
    }

    Ok(())
}

/// Writes all conformers to a single trajectory XYZ file.
///
/// Creates a multi-structure XYZ file containing all valid conformers,
/// which is useful for visualization and analysis of conformational ensembles.
/// Each conformer is separated by its own header with conformer numbering.
///
/// # Arguments
///
/// * `filename` - Path where the trajectory file should be written
/// * `_base_name` - Base name for the molecule (unused, for future features)
/// * `conformers` - Slice of conformer atom vectors to write
/// * `total_conformers` - Total number of conformers for progress tracking
///
/// # Returns
///
/// `Ok(())` on success, or an error if file operations fail.
///
/// # Format
///
/// Each conformer in the trajectory file has the format:
/// ```text
/// 3
/// Conformer 1 of 5 - generated by Rotbond
/// O    0.000000    0.000000    0.000000
/// H    0.757000    0.586000    0.000000
/// H   -0.757000    0.586000    0.000000
/// 3
/// Conformer 2 of 5 - generated by Rotbond
/// ...
/// ```
///
/// # Examples
///
/// ```rust
/// use rotbond::io::write_trajectory_file;
///
/// write_trajectory_file("molecule_traj.xyz", "molecule", &conformers, 24)?;
/// ```
pub fn write_trajectory_file(
    filename: &str,
    _base_name: &str,
    conformers: &[Vec<Atom>],
    total_conformers: usize,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut file = File::create(filename)
        .map_err(|e| format!("Failed to create trajectory file {}: {}", filename, e))?;

    for (idx, atoms) in conformers.iter().enumerate() {
        let conformer_num = idx + 1;

        // Write header
        writeln!(file, "{}", atoms.len())?;
        writeln!(
            file,
            "Conformer {} of {} - generated by Rotbond",
            conformer_num, total_conformers
        )?;

        // Write atoms with proper column alignment
        for atom in atoms {
            writeln!(
                file,
                "{:<2} {:>15.6} {:>15.6} {:>15.6}",
                atom.element, atom.x, atom.y, atom.z
            )?;
        }
    }

    Ok(())
}

/// Calculates the number of padding digits needed for conformer filenames.
///
/// Determines how many digits are needed to represent the total number of
/// conformers, ensuring consistent filename sorting in directory listings.
///
/// # Arguments
///
/// * `total_conformers` - Total number of conformers to be generated
///
/// # Returns
///
/// Number of digits needed for zero-padding.
///
/// # Examples
///
/// ```rust
/// use rotbond::io::calculate_padding_digit_count;
///
/// assert_eq!(calculate_padding_digit_count(5), 1);    // 1-5
/// assert_eq!(calculate_padding_digit_count(50), 2);   // 01-50
/// assert_eq!(calculate_padding_digit_count(500), 3);  // 001-500
/// ```
pub fn calculate_padding_digit_count(total_conformers: usize) -> usize {
    if total_conformers < 10 {
        1
    } else if total_conformers < 100 {
        2
    } else if total_conformers < 1000 {
        3
    } else if total_conformers < 10000 {
        4
    } else {
        total_conformers.to_string().len()
    }
}

/// Formats a conformer filename with smart padding for consistent ordering.
///
/// Generates filenames with appropriate zero-padding so that files sort
/// correctly in directory listings (e.g., "mol_01.xyz", "mol_02.xyz", etc.).
///
/// # Arguments
///
/// * `base_name` - Base name for the file (e.g., "butane")
/// * `index` - Current conformer number (1-based)
/// * `total` - Total number of conformers for padding calculation
///
/// # Returns
///
/// A formatted filename string with appropriate padding.
///
/// # Examples
///
/// ```rust
/// use rotbond::io::format_conformer_filename;
///
/// let filename = format_conformer_filename("butane", 5, 100);
/// assert_eq!(filename, "butane_005.xyz");
///
/// let filename = format_conformer_filename("mbutane", 3, 10);
/// assert_eq!(filename, "butane_03.xyz");
/// ```
pub fn format_conformer_filename(base_name: &str, index: usize, total: usize) -> String {
    let padding = calculate_padding_digit_count(total);
    let index_str = format!("{:0width$}", index, width = padding);
    format!("{}_{}.xyz", base_name, index_str)
}
