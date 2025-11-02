//! # Rotation Parameter Parsing and Management
//!
//! This module handles parsing of rotation parameter files (.rp format)
//! and manages rotation specifications for conformer generation.
//!
//! ## File Format
//!
//! The .rp file format is flexible and supports parameters in any order:
//!
//! ```text
//! bond_factor = 1.0     # Bond detection threshold
//! skip_factor = 0.7     # Steric clash threshold
//! 17-19 bond           # Force a bond between atoms 17 and 19
//! 22-45 nobond         # Remove a bond between atoms 22 and 45
//! 2-7 e120             # Step rotation: 0°, 120°, 240°
//! 1-6 0 60 120 180     # Explicit angles
//! 7-11 e30             # Step rotation every 30°: 0°, 30°, 60°, ..., 330°
//! 24-26 syn 2          # Synchronous with bond 2, same direction
//! 25-27 syn -2         # Synchronous with bond 2, opposite direction
//! ```
//!
//! ## Rotation Types
//!
//! - **Step rotations**: `e30`, `e60`, `e90`, `e120`, `e180` generate evenly spaced angles
//! - **Explicit angles**: List specific angles to use for rotation
//! - **Synchronous rotations**: Follow another bond's rotation with optional direction reversal
//! - **Forced bonds**: Override automatic bond detection for specific atom pairs
//! - **Forbidden bonds**: Remove automatically detected bonds
//!
//! ## Parameter Defaults
//!
//! - `bond_factor`: 1.0 (standard OpenBabel threshold)
//! - `skip_factor`: 0.7 (moderate steric clash detection)

use crate::molecule::RotationSpec;

/// Parses rotation parameters from a .rp format file.
/// 
/// Supports flexible parameter ordering and multiple rotation specification formats.
/// The parser is robust and handles comments, empty lines, and various spacing.
/// 
/// # Arguments
/// 
/// * `filename` - Path to the .rp parameter file
/// 
/// # Returns
/// 
/// A tuple containing:
/// - `bond_factor`: Bond detection threshold multiplier
/// - `skip_factor`: Steric clash detection threshold
/// - `forced_bonds`: List of user-forced bonds (atom index pairs, 1-based)
/// - `forbidden_bonds`: List of user-forbidden bonds (atom index pairs, 1-based)
/// - `rotation_specs`: List of rotation specifications
/// 
/// # Errors
/// 
/// Returns an error if:
/// - File cannot be read
/// - Invalid syntax in parameter file
/// - Atom indices are invalid
/// - Angle values cannot be parsed
/// - Step angles are zero or exceed ±360°
/// 
/// # Examples
/// 
/// ```rust
/// use rotbond::rotation::parse_rotation_file;
/// 
/// let (bond_factor, skip_factor, forced, forbidden, rotations) = 
///     parse_rotation_file("molecule.rp")?;
/// println!("Found {} rotation specifications", rotations.len());
/// ```
/// 
/// # File Format Examples
/// 
/// ```text
/// # Configuration parameters
/// bond_factor = 1.0
/// skip_factor = 0.7
/// 
/// # Bond overrides
/// 17-19 bond           # Force this bond
/// 22-45 nobond         # Remove this bond
/// 
/// # Rotation specifications
/// 1-2 e60              # Every 60°: 0°, 60°, 120°, 180°, 240°, 300°
/// 3-4 0 90 180         # Explicit angles
/// 5-6 syn 1            # Synchronous with bond 1, same direction
/// 7-8 syn -1           # Synchronous with bond 1, opposite direction
/// ```
/// Configuration for conformer generation limits and behavior.
#[derive(Debug, Clone)]
pub struct ConformerConfig {
    /// Maximum number of conformers to generate (None = unlimited)
    pub max_conformers: Option<usize>,
    /// Whether to automatically confirm large generation jobs without prompting
    pub auto_confirm: bool,
}

impl Default for ConformerConfig {
    fn default() -> Self {
        ConformerConfig {
            max_conformers: Some(500), // Default limit of 500 conformers
            auto_confirm: false,
        }
    }
}

pub fn parse_rotation_file(filename: &str) -> Result<(f64, f64, Vec<(usize, usize)>, Vec<(usize, usize)>, Vec<RotationSpec>, ConformerConfig), Box<dyn std::error::Error>> {
    let contents = std::fs::read_to_string(filename)?;

    let mut bond_factor = 1.0;
    let mut skip_factor = 0.7;
    let mut forced_bonds = Vec::new();
    let mut forbidden_bonds = Vec::new();
    let mut rotation_specs = Vec::new();
    let mut conformer_config = ConformerConfig::default();

    for (line_num, line) in contents.lines().enumerate() {
        let line = line.trim();

        // Skip empty lines and comments
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        // Check for configuration parameters
        if line.starts_with("bond_factor") {
            if let Some(value) = line.split('=').nth(1) {
                bond_factor = value.trim().parse::<f64>()
                    .map_err(|e| format!("Line {}: Invalid bond_factor: {}", line_num + 1, e))?;
            }
            continue;
        }

        if line.starts_with("skip_factor") {
            if let Some(value) = line.split('=').nth(1) {
                skip_factor = value.trim().parse::<f64>()
                    .map_err(|e| format!("Line {}: Invalid skip_factor: {}", line_num + 1, e))?;
            }
            continue;
        }

        if line.starts_with("maxgen") {
            if let Some(value) = line.split('=').nth(1) {
                let value = value.trim().to_lowercase();
                if value == "max" || value == "maximum" || value == "full" {
                    conformer_config.max_conformers = None; // Unlimited
                } else {
                    let max_val = value.parse::<usize>()
                        .map_err(|e| format!("Line {}: Invalid maxgen value '{}': {}", line_num + 1, value, e))?;
                    conformer_config.max_conformers = Some(max_val);
                }
            }
            continue;
        }

        if line.starts_with("autoconfirm") {
            if let Some(value) = line.split('=').nth(1) {
                let value = value.trim().to_lowercase();
                conformer_config.auto_confirm = match value.as_str() {
                    "true" | "yes" | "1" | "on" => true,
                    "false" | "no" | "0" | "off" => false,
                    _ => return Err(format!("Line {}: Invalid autoconfirm value '{}' (use true/false)", line_num + 1, value).into()),
                };
            }
            continue;
        }

        // Parse bond definitions or rotation specifications
        let parts: Vec<&str> = line.split_whitespace().collect();

        if parts.len() < 2 {
            continue;
        }

        // Parse bond specification (atom1-atom2 format)
        let bond_spec = parts[0];
        let atom_pair: Vec<&str> = bond_spec.split('-').collect();

        if atom_pair.len() != 2 {
            return Err(format!("Line {}: Invalid bond specification '{}'", line_num + 1, bond_spec).into());
        }

        let atom1 = atom_pair[0].parse::<usize>()
            .map_err(|e| format!("Line {}: Invalid atom index '{}': {}", line_num + 1, atom_pair[0], e))?;
        let atom2 = atom_pair[1].parse::<usize>()
            .map_err(|e| format!("Line {}: Invalid atom index '{}': {}", line_num + 1, atom_pair[1], e))?;

        // Check if this is a manual bond definition
        if parts.len() >= 2 && (parts[1] == "bond" || parts[1] == "nobond") {
            if parts[1] == "bond" {
                forced_bonds.push((atom1, atom2));
            } else {
                forbidden_bonds.push((atom1, atom2));
            }
            continue;
        }

        // Parse rotation specification
        if parts[1] == "syn" {
            // Synchronous rotation
            if parts.len() != 3 {
                return Err(format!("Line {}: Invalid synchronous rotation specification (expected: atom1-atom2 syn <reference>)", line_num + 1).into());
            }

            let ref_str = parts[2];
            let direction;
            let reference;

            // Check for direction sign (+ or -)
            if ref_str.starts_with('-') {
                direction = -1.0;
                let ref_num = &ref_str[1..];
                reference = ref_num.parse::<usize>()
                    .map_err(|e| format!("Line {}: Invalid reference bond index: {}", line_num + 1, e))?;
            } else {
                direction = 1.0;
                reference = ref_str.parse::<usize>()
                    .map_err(|e| format!("Line {}: Invalid reference bond index: {}", line_num + 1, e))?;
            }

            rotation_specs.push(RotationSpec::Synchronous {
                reference,
                direction,
            });
        } else if parts[1].starts_with('e') {
            // Step-based rotation (e.g., e20, e-60)
            let step_str = &parts[1][1..];
            let step = step_str.parse::<f64>()
                .map_err(|e| format!("Line {}: Invalid step angle '{}': {}", line_num + 1, step_str, e))?;

            if step == 0.0 || step.abs() > 360.0 {
                return Err(format!("Line {}: Step angle must be non-zero and between -360 and 360 degrees", line_num + 1).into());
            }

            rotation_specs.push(RotationSpec::Step { step });
        } else {
            // Explicit angles
            let mut angles = Vec::new();

            for i in 1..parts.len() {
                let angle = parts[i].parse::<f64>()
                    .map_err(|e| format!("Line {}: Invalid angle '{}': {}", line_num + 1, parts[i], e))?;

                if angle < -360.0 || angle > 360.0 {
                    return Err(format!("Line {}: Angle must be between -360 and 360 degrees", line_num + 1).into());
                }

                angles.push(angle);
            }

            if angles.is_empty() {
                return Err(format!("Line {}: No angles specified", line_num + 1).into());
            }

            rotation_specs.push(RotationSpec::Explicit(angles));
        }
    }

    Ok((bond_factor, skip_factor, forced_bonds, forbidden_bonds, rotation_specs, conformer_config))
}

/// Generates angle sets from rotation specifications.
/// 
/// Converts rotation specifications into concrete lists of angles for each bond.
/// Handles step-based rotations, explicit angles, and synchronous rotations.
/// 
/// # Arguments
/// 
/// * `specs` - Vector of rotation specifications to process
/// 
/// # Returns
/// 
/// A vector of vectors, where each inner vector contains the rotation angles
/// for one bond. Synchronous bonds have empty angle vectors since they follow
/// their reference bonds.
/// 
/// # Algorithm
/// 
/// - **Step rotations**: Generates evenly spaced angles from 0° to 360°-step
/// - **Explicit angles**: Uses provided angles, sorted and deduplicated
/// - **Synchronous rotations**: Creates empty angle set (uses reference bond's angles)
/// 
/// # Examples
/// 
/// ```rust
/// use rotbond::rotation::{generate_angle_sets, RotationSpec};
/// 
/// let specs = vec![
///     RotationSpec::Step { step: 60.0 },
///     RotationSpec::Explicit(vec![0.0, 90.0, 180.0]),
/// ];
/// let angle_sets = generate_angle_sets(specs)?;
/// // angle_sets[0] = [0.0, 60.0, 120.0, 180.0, 240.0, 300.0]
/// // angle_sets[1] = [0.0, 90.0, 180.0]
/// ```
/// 
/// # Step Rotation Details
/// 
/// For step rotations, the number of angles is calculated as `360° / |step|`.
/// The direction of rotation is preserved:
/// - Positive step: clockwise rotation (0°, +step, +2*step, ...)
/// - Negative step: counterclockwise rotation (0°, -step, -2*step, ...)
pub fn generate_angle_sets(specs: Vec<RotationSpec>) -> Result<Vec<Vec<f64>>, Box<dyn std::error::Error>> {
    let mut angle_sets = Vec::new();

    for spec in specs {
        match spec {
            RotationSpec::Step { step } => {
                let mut angles = Vec::new();
                
                // Always use absolute step size for angle calculation
                let abs_step = step.abs();
                
                // Calculate number of rotations: 360 / |step|
                let num_rotations = (360.0 / abs_step).round() as i32;
                
                // Generate sequential angles: 0, step, 2*step, 3*step, ...
                // Apply direction by using the sign of the original step
                for i in 0..num_rotations {
                    let base_angle = i as f64 * abs_step;
                    
                    // Apply direction: positive step = clockwise, negative step = counterclockwise
                    let angle = if step > 0.0 {
                        // Clockwise: use positive angles as-is
                        base_angle
                    } else {
                        // Counterclockwise: negate the angle
                        -base_angle
                    };
                    
                    angles.push(angle);
                }
                
                // Keep angles in sequential order (do NOT sort or normalize)
                // This preserves the trajectory order for users
                angles.dedup();
                angle_sets.push(angles);
            }
            RotationSpec::Explicit(angles) => {
                // Remove duplicates and sort
                let mut unique_angles = angles;
                unique_angles.sort_by(|a, b| a.partial_cmp(b).unwrap());
                unique_angles.dedup();
                angle_sets.push(unique_angles);
            }
            RotationSpec::Synchronous { .. } => {
                // Synchronous bonds don't have their own angles
                // They'll use the reference bond's angles
                angle_sets.push(Vec::new());
            }
        }
    }

    Ok(angle_sets)
}

/// Validates synchronous rotation references for correctness.
/// 
/// Ensures that all synchronous rotation specifications reference valid bonds
/// and do not create invalid reference patterns like self-references.
/// 
/// # Arguments
/// 
/// * `specs` - Slice of rotation specifications to validate
/// 
/// # Returns
/// 
/// `Ok(())` if all references are valid, or an error describing the problem.
/// 
/// # Validation Rules
/// 
/// - Reference bond indices must be within valid range (1 to number of bonds)
/// - Bonds cannot reference themselves
/// - Reference bonds must exist in the specification list
/// 
/// # Errors
/// 
/// Returns an error if:
/// - A bond references a non-existent bond index
/// - A bond references itself (circular reference)
/// - Reference indices are out of bounds
/// 
/// # Examples
/// 
/// ```rust
/// use rotbond::rotation::{validate_synchronous_references, RotationSpec};
/// 
/// let specs = vec![
///     RotationSpec::Step { step: 60.0 },
///     RotationSpec::Synchronous { reference: 1, direction: 1.0 },
/// ];
/// validate_synchronous_references(&specs)?; // OK
/// 
/// let bad_specs = vec![
///     RotationSpec::Synchronous { reference: 1, direction: 1.0 }, // Self-reference
/// ];
/// // validate_synchronous_references(&bad_specs) returns error
/// ```
/// 
/// # Future Enhancements
/// 
/// Currently performs basic validation. Future versions may include:
/// - Circular reference chain detection
/// - Reference dependency graph analysis
/// - Validation of reference bond types
pub fn validate_synchronous_references(specs: &[RotationSpec]) -> Result<(), Box<dyn std::error::Error>> {
    let mut reference_counts = vec![0; specs.len()];

    for (i, spec) in specs.iter().enumerate() {
        if let RotationSpec::Synchronous { reference, .. } = spec {
            let ref_idx = *reference - 1; // Convert to 0-indexed

            if ref_idx >= specs.len() {
                return Err(format!("Bond {} references non-existent bond {}", i + 1, reference).into());
            }

            if ref_idx == i {
                return Err(format!("Bond {} references itself", i + 1).into());
            }

            // Check for circular references would require a graph traversal
            // For now, just count references
            reference_counts[ref_idx] += 1;
        }
    }

    Ok(())
}
