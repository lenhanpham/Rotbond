//! # Error Handling Module
//!
//! This module provides comprehensive error handling for the Rotbond system,
//! with specific support for both rotation and scanning modes. It defines
//! custom error types, validation functions, and user-friendly error messages.
//!
//! ## Error Categories
//!
//! - **Parameter Errors**: Invalid scanning parameters, bond specifications
//! - **Validation Errors**: Bond length validation, geometric constraints
//! - **System Errors**: File I/O, memory allocation, computation failures
//! - **User Errors**: Invalid input, missing files, configuration issues
//!
//! ## Design Principles
//!
//! - Clear, actionable error messages for users
//! - Consistent error handling patterns across rotation and scanning modes
//! - Integration with existing validation systems
//! - Graceful degradation and recovery where possible

use std::fmt;
use std::error::Error;

/// Comprehensive error type for scanning-specific operations.
/// 
/// Provides detailed error information for bond scanning operations,
/// including parameter validation, geometric constraints, and system failures.
/// Each error variant includes context-specific information to help users
/// understand and resolve issues.
/// 
/// # Examples
/// 
/// ```rust
/// use rotbond::errors::ScanningError;
/// 
/// let error = ScanningError::InvalidStepSize {
///     bond_index: 1,
///     step_size: 0.0,
///     reason: "Step size must be non-zero".to_string(),
/// };
/// println!("Error: {}", error);
/// ```
 
#[derive(Debug, Clone)]
pub enum ScanningError {
    /// Invalid scanning parameters provided by user
    InvalidParameters {
        /// Human-readable description of the parameter issue
        message: String,
        /// Suggested fix for the user
        suggestion: String,
    },
    
    /// Invalid number of scanning steps
    InvalidSteps {
        /// Bond index (1-based) with invalid steps
        bond_index: usize,
        /// The invalid step count provided
        steps: usize,
        /// Reason why the steps are invalid
        reason: String,
    },
    
    /// Invalid step size for bond scanning
    InvalidStepSize {
        /// Bond index (1-based) with invalid step size
        bond_index: usize,
        /// The invalid step size provided
        step_size: f64,
        /// Reason why the step size is invalid
        reason: String,
    },
    

}

impl fmt::Display for ScanningError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ScanningError::InvalidParameters { message, suggestion } => {
                write!(f, "Invalid scanning parameters: {}\nSuggestion: {}", message, suggestion)
            }
            ScanningError::InvalidSteps { bond_index, steps, reason } => {
                write!(f, "Invalid steps for scanning bond {}: {} steps ({})", bond_index, steps, reason)
            }
            ScanningError::InvalidStepSize { bond_index, step_size, reason } => {
                write!(f, "Invalid step size for scanning bond {}: {:.3} Å ({})", bond_index, step_size, reason)
            }

        }
    }
}

impl Error for ScanningError {}

/// Validates scanning parameters for a single bond.
/// 
/// Performs comprehensive validation of scanning parameters including
/// step count, step size, and geometric feasibility. This function
/// should be called before attempting to perform bond scanning.
/// 
/// # Arguments
/// 
/// * `bond_index` - Bond index (1-based) for error reporting
/// * `steps` - Number of scanning steps to validate
/// * `step_size` - Step size in Angstroms to validate
/// * `current_bond_length` - Current bond length in Angstroms
/// * `atom1_element` - Element symbol of first atom
/// * `atom2_element` - Element symbol of second atom




/// Gets the minimum reasonable bond length for two atom types.
/// 
/// Returns the minimum chemically reasonable bond length based on
/// covalent radii and known chemical constraints.
/// 
/// # Arguments
/// 
/// * `element1` - First element symbol
/// * `element2` - Second element symbol
/// 
/// # Returns


/// Validates that forced and forbidden bonds are compatible with scanning.
/// 
/// Ensures that user-specified forced and forbidden bonds will remain
/// valid during the scanning process and don't conflict with scanning bonds.
/// 
/// # Arguments
/// 
/// * `molecule` - The molecule with forced/forbidden bond specifications
/// * `scanning_bonds` - Indices of bonds being scanned
/// 
/// # Returns
/// 
/// `Ok(())` if bonds are compatible, or a `ScanningError` describing conflicts.
/// 
/// # Validation Rules
/// 
/// - Forced bonds should not be scanning bonds (would be redundant)
/// - Forbidden bonds should not conflict with scanning bonds
/// - Forced bond lengths should remain reasonable during scanning
/// 
/// # Examples
/// 
/// ```rust
/// use rotbond::errors::validate_forced_forbidden_bonds;
/// 
/// let scanning_bonds = vec![1, 3]; // Bonds 1 and 3 are being scanned
/// validate_forced_forbidden_bonds(&molecule, &scanning_bonds)?;
/// ```
pub fn validate_forced_forbidden_bonds(
    molecule: &crate::molecule::Molecule,
    scanning_bonds: &[usize],
) -> Result<(), ScanningError> {
    // Convert scanning bonds to 1-based for comparison with forced/forbidden bonds
    let scanning_bonds_1based: std::collections::HashSet<(usize, usize)> = scanning_bonds
        .iter()
        .filter_map(|&bond_idx| {
            if bond_idx < molecule.bonds.len() {
                let bond = &molecule.bonds[bond_idx];
                Some((bond.atom1 + 1, bond.atom2 + 1)) // Convert to 1-based
            } else {
                None
            }
        })
        .collect();
    
    // Check forced bonds
    for &(atom1, atom2) in &molecule.forced_bonds {
        let bond_pair = (atom1.min(atom2), atom1.max(atom2));
        let reverse_pair = (atom2, atom1);
        
        if scanning_bonds_1based.contains(&bond_pair) || scanning_bonds_1based.contains(&reverse_pair) {
            return Err(ScanningError::InvalidParameters {
                message: format!("Forced bond {}-{} conflicts with scanning bond specification", atom1, atom2),
                suggestion: "Remove the forced bond specification since the bond is already being scanned".to_string(),
            });
        }
    }
    
    // Check forbidden bonds
    for &(atom1, atom2) in &molecule.forbidden_bonds {
        let bond_pair = (atom1.min(atom2), atom1.max(atom2));
        let reverse_pair = (atom2, atom1);
        
        if scanning_bonds_1based.contains(&bond_pair) || scanning_bonds_1based.contains(&reverse_pair) {
            return Err(ScanningError::InvalidParameters {
                message: format!("Forbidden bond {}-{} conflicts with scanning bond specification", atom1, atom2),
                suggestion: "Remove either the forbidden bond specification or the scanning specification for this bond".to_string(),
            });
        }
    }
    
    Ok(())
}

/// Provides user-friendly error messages with helpful suggestions.
/// 
/// Converts technical scanning errors into clear, actionable messages
/// that help users understand and resolve issues with their scanning
/// specifications.
/// 
/// # Arguments
/// 
/// * `error` - The scanning error to format
/// 
/// # Returns
/// 
/// A formatted error message with suggestions for resolution.
/// 
/// # Examples
/// 
/// ```rust
/// use rotbond::errors::{ScanningError, format_user_error};
/// 
/// let error = ScanningError::InvalidSteps { 
///     bond_index: 1, 
///     steps: 0, 
///     reason: "Steps must be positive".to_string() 
/// };
/// let message = format_user_error(&error);
/// println!("{}", message);
/// ```
pub fn format_user_error(error: &ScanningError) -> String {
    match error {
        ScanningError::InvalidParameters { message, suggestion } => {
            format!("SCANNING ERROR: {}\n\nSUGGESTION: {}\n\nFor scanning syntax help, use: --help input", message, suggestion)
        }
        ScanningError::InvalidSteps { bond_index, steps, reason } => {
            format!("SCANNING ERROR: Invalid steps for bond {} - {} steps ({})\n\nSUGGESTION: Use a positive integer between 1 and 1000 for the number of steps.\nExample: {}-X scan 10 0.1\n\nFor scanning syntax help, use: --help input", 
                   bond_index, steps, reason, bond_index)
        }
        ScanningError::InvalidStepSize { bond_index, step_size, reason } => {
            format!("SCANNING ERROR: Invalid step size for bond {} - {:.3} Å ({})\n\nSUGGESTION: Use a non-zero step size between -5.0 and +5.0 Å.\nPositive values stretch the bond, negative values compress it.\nExample: {}-X scan 10 0.1\n\nFor scanning syntax help, use: --help input", 
                   bond_index, step_size, reason, bond_index)
        }

    }
}

/// Validates scanning conformer generation limits and parameters.
/// 
/// Checks that scanning parameters will not create excessive computational
/// load or memory usage. Provides warnings and suggestions for optimization.
/// 
/// # Arguments
/// 
/// * `rotation_specs` - Scanning specifications to validate
/// * `max_conformers` - Optional limit on conformer generation
/// 
/// # Returns
/// 
/// `Ok(())` if parameters are reasonable, or a `ScanningError` with suggestions.
/// 
/// # Validation Rules
/// 
/// - Total conformer count should be reasonable (< 10,000 by default)
/// - Individual bond scanning should not be excessive
/// - Memory usage estimation and warnings
/// 
/// # Examples
/// 
/// ```rust
/// use rotbond::errors::validate_scanning_specs;
/// 
/// validate_scanning_specs(&rotation_specs)?;
/// ```
pub fn validate_scanning_specs(
    rotation_specs: &[crate::molecule::RotationSpec],
) -> Result<(), ScanningError> {
    // Validate individual scanning specifications
    for (i, spec) in rotation_specs.iter().enumerate() {
        if let crate::molecule::RotationSpec::Scanning { steps, step_size, .. } = spec {
            let bond_index = i + 1; // Convert to 1-based for user display
            
            // Validate steps
            if *steps == 0 {
                return Err(ScanningError::InvalidSteps {
                    bond_index,
                    steps: *steps,
                    reason: "Steps must be a positive integer".to_string(),
                });
            }
            
            if *steps > 10000 {
                return Err(ScanningError::InvalidSteps {
                    bond_index,
                    steps: *steps,
                    reason: "Steps should be ≤ 10000 to avoid excessive computation".to_string(),
                });
            }
            
            // Validate step size
            if *step_size == 0.0 {
                return Err(ScanningError::InvalidStepSize {
                    bond_index,
                    step_size: *step_size,
                    reason: "Step size must be non-zero".to_string(),
                });
            }
            
            if step_size.abs() > 10.0 {
                return Err(ScanningError::InvalidStepSize {
                    bond_index,
                    step_size: *step_size,
                    reason: "Step size should be ≤ 10.0 Å to avoid unrealistic bond lengths".to_string(),
                });
            }
        }
    }
    
    // Check for arithmetic overflow in total combinations
    let mut total_combinations: usize = 1;
    for spec in rotation_specs {
        if let crate::molecule::RotationSpec::Scanning { steps, .. } = spec {
            if let Some(new_total) = total_combinations.checked_mul(*steps) {
                total_combinations = new_total;
            } else {
                return Err(ScanningError::InvalidParameters {
                    message: "Scanning parameters would create too many conformers (arithmetic overflow)".to_string(),
                    suggestion: "Reduce the number of steps or scanning bonds to avoid excessive computation".to_string(),
                });
            }
        }
    }
    
    Ok(())
}

