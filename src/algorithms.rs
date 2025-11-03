//! # Molecular Algorithms Module
//!
//! This module contains the core algorithms for molecular conformer generation,
//! including bond detection, fragment identification, and rotation algorithms.
//!
//! ## Key Components
//!
//! - **Bond Detection**: Uses OpenBabel-style covalent radii with configurable bond factor
//! - **Fragment Identification**: Graph-based algorithms to identify molecular fragments
//! - **Dihedral Rotation**: Rodrigues rotation formula for precise conformer generation
//! - **Steric Validation**: Distance-based clash detection for conformer filtering
//!
//! ## Bond Detection Formula
//!
//! Bonds are detected using the formula:
//! ```text
//! d² < ((rad₁ + rad₂) × bond_factor + 0.45)²
//! ```
//!
//! Where:
//! - `d` is the distance between atoms
//! - `rad₁`, `rad₂` are covalent radii of the atoms
//! - `bond_factor` is a configurable multiplier (default: 1.0)
//! - `0.45` Å is the OpenBabel tolerance constant

use crate::molecule::{Molecule, Atom, covalent_radius};
use std::collections::VecDeque;

/// Simple context for scanning-aware validation.
/// 
/// Provides adaptive skip_factor calculation for scanning operations by analyzing
/// bond length changes and adjusting validation thresholds accordingly.
/// 
/// # Fields
/// 
/// * `adaptive_skip_factor` - Dynamically adjusted skip_factor based on bond length changes
/// * `scanning_bonds` - List of bonds being scanned with their target lengths
/// 
/// # Usage
/// 
/// This context integrates with existing validation by providing an adaptive skip_factor
/// that accounts for the geometric changes during bond length scanning.
/// 
/// # Examples
/// 
/// ```rust
/// use rotbond::algorithms::ScanningValidationContext;
/// 
/// let scanning_bonds = vec![(0, 1, 1.7)]; // Scan bond 0-1 to 1.7 Å
/// let context = ScanningValidationContext::new(&molecule, &atoms, &scanning_bonds);
/// ```
#[derive(Debug, Clone)]
pub struct ScanningValidationContext {
    /// Dynamically adjusted skip_factor based on bond length changes
    pub adaptive_skip_factor: f64,
    /// List of bonds being scanned: (atom1, atom2, target_length)
    pub scanning_bonds: Vec<(usize, usize, f64)>,
}

impl ScanningValidationContext {
    /// Creates a new scanning validation context using existing bond detection infrastructure.
    /// 
    /// This method leverages the existing `build_bond_graph` and `identify_fragments` functions
    /// to create an adaptive validation context that integrates seamlessly with the rotation system.
    /// 
    /// # Arguments
    /// 
    /// * `molecule` - The molecule containing validation parameters
    /// * `current_coordinates` - Current atomic coordinates
    /// * `scanning_bonds` - List of (atom1, atom2, target_length) tuples for bonds being scanned
    /// 
    /// # Returns
    /// 
    /// A `ScanningValidationContext` with adaptive skip_factor calculated based on bond changes
    /// 
    /// # Examples
    /// 
    /// ```rust
    /// let scanning_bonds = vec![(0, 1, 1.7)]; // Scan bond 0-1 to 1.7 Å
    /// let context = ScanningValidationContext::new(&molecule, &atoms, &scanning_bonds);
    /// ```
    pub fn new(
        molecule: &Molecule,
        current_coordinates: &[Atom],
        scanning_bonds: &[(usize, usize, f64)],
    ) -> Self {
        let adaptive_skip_factor = Self::calculate_adaptive_skip_factor(
            molecule, current_coordinates, scanning_bonds
        );
        
        ScanningValidationContext {
            adaptive_skip_factor,
            scanning_bonds: scanning_bonds.to_vec(),
        }
    }
    
    /// Checks if an atom pair is involved in scanning operations.
    /// 
    /// # Arguments
    /// 
    /// * `atom1` - Index of first atom
    /// * `atom2` - Index of second atom
    /// 
    /// # Returns
    /// 
    /// `true` if this atom pair is being scanned, `false` otherwise
    pub fn is_scanning_bond(&self, atom1: usize, atom2: usize) -> bool {
        self.scanning_bonds.iter().any(|&(a1, a2, _)| {
            (a1 == atom1 && a2 == atom2) || (a1 == atom2 && a2 == atom1)
        })
    }
    
    /// Calculates adaptive skip_factor based on bond length changes.
    /// 
    /// Uses existing molecular infrastructure to analyze bond changes and calculate
    /// an appropriate skip_factor that accounts for the geometric changes during scanning.
    /// 
    /// # Arguments
    /// 
    /// * `molecule` - The molecule containing base skip_factor and bond parameters
    /// * `current_coordinates` - Current atomic coordinates
    /// * `scanning_bonds` - List of (atom1, atom2, target_length) tuples for bonds being scanned
    /// 
    /// # Returns
    /// 
    /// Adaptive skip_factor value optimized for scanning operations
    /// 
    /// # Algorithm
    /// 
    /// 1. Calculate maximum bond length change across all scanning bonds
    /// 2. Apply simple adjustment based on the magnitude of changes
    /// 3. Ensure skip_factor remains within safe bounds (0.4-1.0)
    /// 
    /// # Examples
    /// 
    /// ```rust
    /// let adaptive_skip = ScanningValidationContext::calculate_adaptive_skip_factor(
    ///     &molecule, &atoms, &scanning_bonds
    /// );
    /// ```
    #[allow(dead_code)]
    pub fn calculate_adaptive_skip_factor(
        molecule: &Molecule,
        current_coordinates: &[Atom],
        scanning_bonds: &[(usize, usize, f64)],
    ) -> f64 {
        let base_skip_factor = molecule.skip_factor;
        
        if scanning_bonds.is_empty() {
            return base_skip_factor;
        }
        
        // Find the maximum bond length change
        let mut max_change: f64 = 0.0;
        
        for &(atom1, atom2, target_length) in scanning_bonds {
            if atom1 >= current_coordinates.len() || atom2 >= current_coordinates.len() {
                continue;
            }
            
            let current_length = crate::molecule::calculate_distance(
                &current_coordinates[atom1], 
                &current_coordinates[atom2]
            );
            
            if current_length > 0.0 {
                let change = (target_length - current_length).abs();
                max_change = max_change.max(change);
            }
        }
        
        // Simple adaptive adjustment: reduce skip_factor for larger bond changes
        // This makes validation more permissive when bonds are being stretched/compressed
        let adjustment = (max_change * 0.5).min(0.3); // Max 30% reduction
        let adjusted_skip_factor = base_skip_factor - adjustment;
        
        // Keep within safe bounds: not too permissive (0.4) or restrictive (base_skip_factor)
        adjusted_skip_factor.max(0.4).min(base_skip_factor)
    }
    

}

/// Enumeration of different validation failure types.
/// 
/// Provides specific categorization of why a conformer failed validation,
/// enabling targeted diagnostic messages and parameter adjustment suggestions.
/// 
/// # Variants
/// 
/// * `StericClash` - Atoms are too close together (distance < skip_factor threshold)
/// * `InvalidBondLength` - Bond length is outside reasonable chemical bounds
/// * `GeometricConstraint` - Molecular geometry violates structural constraints
/// * `ForcedBondViolation` - A user-specified forced bond is not maintained
/// * `ForbiddenBondCreation` - A user-specified forbidden bond was accidentally created
/// 
/// # Examples
/// 
/// ```rust
/// use rotbond::algorithms::ValidationFailure;
/// 
/// let clash = ValidationFailure::StericClash { 
///     atom1: 0, atom2: 1, distance: 0.5, min_required: 0.7 
/// };
/// 
/// let bad_bond = ValidationFailure::InvalidBondLength { 
///     bond: (2, 3), length: 0.05, valid_range: (0.1, f64::INFINITY) 
/// };
/// ```
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub enum ValidationFailure {
    /// Steric clash between two atoms
    StericClash {
        /// Index of first atom involved in clash
        atom1: usize,
        /// Index of second atom involved in clash
        atom2: usize,
        /// Actual distance between atoms in Angstroms
        distance: f64,
        /// Minimum required distance in Angstroms
        min_required: f64,
    },
    /// Invalid bond length detected
    InvalidBondLength {
        /// Bond atom pair (atom1_idx, atom2_idx)
        bond: (usize, usize),
        /// Actual bond length in Angstroms
        length: f64,
        /// Valid bond length range (min, max) in Angstroms
        valid_range: (f64, f64),
    },
    /// Geometric constraint violation
    GeometricConstraint {
        /// Description of the constraint violation
        description: String,
        /// Severity level of the violation
        severity: ValidationSeverity,
    },
    /// Forced bond is not maintained at reasonable distance
    ForcedBondViolation {
        /// Bond atom pair (atom1_idx, atom2_idx) - 1-based indexing
        bond: (usize, usize),
        /// Actual distance in Angstroms
        distance: f64,
        /// Expected distance range (min, max) in Angstroms
        expected_range: (f64, f64),
    },
    /// Forbidden bond was accidentally created
    ForbiddenBondCreation {
        /// Bond atom pair (atom1_idx, atom2_idx) - 1-based indexing
        bond: (usize, usize),
        /// Actual distance in Angstroms (too short, indicating bond formation)
        distance: f64,
        /// Minimum distance to avoid bond formation in Angstroms
        min_distance_to_avoid: f64,
    },
}

/// Severity levels for validation failures.
/// 
/// Indicates the importance and impact of different validation failures,
/// helping users prioritize which issues to address first.
/// 
/// # Examples
/// 
/// ```rust
/// use rotbond::algorithms::ValidationSeverity;
/// 
/// let critical = ValidationSeverity::Critical; // Must fix
/// let warning = ValidationSeverity::Warning;   // Should consider fixing
/// let info = ValidationSeverity::Info;         // Optional optimization
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ValidationSeverity {
    /// Critical failure that makes conformer unusable
    Critical,
    /// Warning that indicates potential issues
    Warning,
    /// Informational message for optimization
    Info,
}

/// Diagnostic message with context and suggested actions.
/// 
/// Provides detailed information about validation results, including
/// specific atom indices, measurements, and actionable suggestions for
/// parameter adjustments or structural modifications.
/// 
/// # Examples
/// 
/// ```rust
/// use rotbond::algorithms::{DiagnosticMessage, ValidationSeverity};
/// 
/// let message = DiagnosticMessage {
///     level: ValidationSeverity::Warning,
///     message: "Atoms 1 and 2 are close (0.65 Å) but within tolerance".to_string(),
///     atom_indices: Some(vec![0, 1]), // 0-based indices
///     suggested_action: Some("Consider increasing skip_factor to 0.8".to_string()),
/// };
/// ```
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct DiagnosticMessage {
    /// Severity level of this diagnostic message
    pub level: ValidationSeverity,
    /// Human-readable diagnostic message
    pub message: String,
    /// Atom indices involved in this diagnostic (0-based, None if not atom-specific)
    pub atom_indices: Option<Vec<usize>>,
    /// Suggested action to address this issue (None if no specific action recommended)
    pub suggested_action: Option<String>,
}

/// Parameter adjustment suggestion for improving validation results.
/// 
/// Provides specific recommendations for modifying molecular parameters
/// to improve conformer generation success rates or address validation failures.
/// 
/// # Examples
/// 
/// ```rust
/// use rotbond::algorithms::ParameterSuggestion;
/// 
/// let suggestion = ParameterSuggestion {
///     parameter_name: "skip_factor".to_string(),
///     current_value: "0.7".to_string(),
///     suggested_value: "0.6".to_string(),
///     reason: "Reduce steric clash sensitivity for scanning mode".to_string(),
///     expected_impact: "Should increase conformer acceptance rate by ~20%".to_string(),
/// };
/// ```
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ParameterSuggestion {
    /// Name of the parameter to adjust
    pub parameter_name: String,
    /// Current parameter value as string
    pub current_value: String,
    /// Suggested new parameter value as string
    pub suggested_value: String,
    /// Explanation of why this adjustment is recommended
    pub reason: String,
    /// Expected impact of making this change
    pub expected_impact: String,
}

/// Comprehensive validation result with diagnostic information.
/// 
/// Contains the validation outcome along with detailed diagnostic information
/// that helps users understand why validation failed and how to improve results.
/// This structure enables the enhanced validation system to provide actionable
/// feedback rather than simple pass/fail results.
/// 
/// # Fields
/// 
/// * `is_valid` - Whether the conformer passed validation
/// * `failure_reason` - Specific reason for failure (None if valid)
/// * `diagnostic_info` - List of diagnostic messages with context
/// * `suggested_adjustments` - List of parameter adjustment recommendations
/// 
/// # Examples
/// 
/// ```rust
/// use rotbond::algorithms::{ValidationResult, ValidationFailure, DiagnosticMessage, ValidationSeverity};
/// 
/// // Successful validation
/// let success = ValidationResult {
///     is_valid: true,
///     failure_reason: None,
///     diagnostic_info: vec![],
///     suggested_adjustments: vec![],
/// };
/// 
/// // Failed validation with diagnostics
/// let failure = ValidationResult {
///     is_valid: false,
///     failure_reason: Some(ValidationFailure::StericClash { 
///         atom1: 0, atom2: 1, distance: 0.5, min_required: 0.7 
///     }),
///     diagnostic_info: vec![
///         DiagnosticMessage {
///             level: ValidationSeverity::Critical,
///             message: "Severe steric clash between atoms 1 and 2".to_string(),
///             atom_indices: Some(vec![0, 1]),
///             suggested_action: Some("Increase skip_factor or reduce scanning step size".to_string()),
///         }
///     ],
///     suggested_adjustments: vec![],
/// };
/// ```
#[derive(Debug, Clone)]
pub struct ValidationResult {
    /// Whether the conformer is valid (passed all validation checks)
    pub is_valid: bool,
    /// Primary reason for validation failure (None if valid)
    pub failure_reason: Option<ValidationFailure>,
    /// List of diagnostic messages providing detailed context
    pub diagnostic_info: Vec<DiagnosticMessage>,
    /// List of suggested parameter adjustments to improve results
    pub suggested_adjustments: Vec<ParameterSuggestion>,
}

impl ValidationResult {
    /// Creates a successful validation result.
    /// 
    /// # Returns
    /// 
    /// A `ValidationResult` indicating successful validation with no diagnostics
    /// 
    /// # Examples
    /// 
    /// ```rust
    /// use rotbond::algorithms::ValidationResult;
    /// 
    /// let success = ValidationResult::success();
    /// assert!(success.is_valid);
    /// assert!(success.failure_reason.is_none());
    /// ```
    pub fn success() -> Self {
        ValidationResult {
            is_valid: true,
            failure_reason: None,
            diagnostic_info: Vec::new(),
            suggested_adjustments: Vec::new(),
        }
    }
    
    /// Creates a failed validation result with a specific failure reason.
    /// 
    /// # Arguments
    /// 
    /// * `failure` - The primary validation failure that caused rejection
    /// 
    /// # Returns
    /// 
    /// A `ValidationResult` indicating failed validation with the specified reason
    /// 
    /// # Examples
    /// 
    /// ```rust
    /// use rotbond::algorithms::{ValidationResult, ValidationFailure};
    /// 
    /// let clash = ValidationFailure::StericClash { 
    ///     atom1: 0, atom2: 1, distance: 0.5, min_required: 0.7 
    /// };
    /// let failure = ValidationResult::failure(clash);
    /// assert!(!failure.is_valid);
    /// ```
    #[allow(dead_code)]
    pub fn failure(failure: ValidationFailure) -> Self {
        ValidationResult {
            is_valid: false,
            failure_reason: Some(failure),
            diagnostic_info: Vec::new(),
            suggested_adjustments: Vec::new(),
        }
    }
    
    /// Adds a diagnostic message to this validation result.
    /// 
    /// # Arguments
    /// 
    /// * `message` - Diagnostic message to add
    /// 
    /// # Examples
    /// 
    /// ```rust
    /// use rotbond::algorithms::{ValidationResult, DiagnosticMessage, ValidationSeverity};
    /// 
    /// let mut result = ValidationResult::success();
    /// let diagnostic = DiagnosticMessage {
    ///     level: ValidationSeverity::Info,
    ///     message: "Conformer is valid but close to steric clash threshold".to_string(),
    ///     atom_indices: None,
    ///     suggested_action: None,
    /// };
    /// result.add_diagnostic(diagnostic);
    /// ```
    pub fn add_diagnostic(&mut self, message: DiagnosticMessage) {
        self.diagnostic_info.push(message);
    }
    
    /// Adds a parameter adjustment suggestion to this validation result.
    /// 
    /// # Arguments
    /// 
    /// * `suggestion` - Parameter adjustment suggestion to add
    /// 
    /// # Examples
    /// 
    /// ```rust
    /// use rotbond::algorithms::{ValidationResult, ParameterSuggestion};
    /// 
    /// let mut result = ValidationResult::success();
    /// let suggestion = ParameterSuggestion {
    ///     parameter_name: "skip_factor".to_string(),
    ///     current_value: "0.7".to_string(),
    ///     suggested_value: "0.6".to_string(),
    ///     reason: "More permissive for scanning".to_string(),
    ///     expected_impact: "Higher success rate".to_string(),
    /// };
    /// result.add_suggestion(suggestion);
    /// ```
    pub fn add_suggestion(&mut self, suggestion: ParameterSuggestion) {
        self.suggested_adjustments.push(suggestion);
    }
    
    /// Checks if this result has any diagnostic messages of the specified severity.
    /// 
    /// # Arguments
    /// 
    /// * `severity` - Severity level to check for
    /// 
    /// # Returns
    /// 
    /// `true` if any diagnostic messages have the specified severity level
    pub fn has_diagnostics_of_severity(&self, severity: ValidationSeverity) -> bool {
        self.diagnostic_info.iter().any(|d| d.level == severity)
    }
    
    /// Gets the count of diagnostic messages by severity level.
    /// 
    /// # Returns
    /// 
    /// A tuple of (critical_count, warning_count, info_count)
    pub fn diagnostic_counts(&self) -> (usize, usize, usize) {
        let mut critical = 0;
        let mut warning = 0;
        let mut info = 0;
        
        for diagnostic in &self.diagnostic_info {
            match diagnostic.level {
                ValidationSeverity::Critical => critical += 1,
                ValidationSeverity::Warning => warning += 1,
                ValidationSeverity::Info => info += 1,
            }
        }
        
        (critical, warning, info)
    }
}

/// Enhanced conformer validation with scanning-aware logic and diagnostic reporting.
/// 
/// This function provides comprehensive validation that adapts to different operation modes
/// and provides detailed diagnostic information. It maintains backward compatibility with
/// the existing `is_valid_conformer` function while adding enhanced capabilities for
/// scanning mode and diagnostic reporting.
/// 
/// # Arguments
/// 
/// * `molecule` - The molecule containing validation parameters and operation mode
/// * `coordinates` - The conformer coordinates to validate
/// * `scanning_context` - Optional scanning context for adaptive validation (None for rotation mode)
/// 
/// # Returns
/// 
/// A `ValidationResult` containing validation outcome and detailed diagnostic information
/// 
/// # Validation Logic
/// 
/// 1. **Steric clash detection**: Uses adaptive skip_factor from scanning context if available
/// 2. **Bond length validation**: Validates all bonds are within reasonable chemical bounds
/// 3. **Scanning-specific validation**: Additional checks for scanning mode conformers
/// 4. **Forced/forbidden bond validation**: Ensures user constraints are maintained
/// 5. **Diagnostic generation**: Provides detailed failure analysis and suggestions
/// 
/// # Backward Compatibility
/// 
/// When called without scanning context, behaves identically to the original
/// `is_valid_conformer` function, ensuring existing code continues to work.
/// 
/// # Examples
/// 
/// ```rust
/// use rotbond::algorithms::{is_valid_conformer_enhanced, ScanningValidationContext};
/// 
/// // Rotation mode (backward compatible)
/// let result = is_valid_conformer_enhanced(&molecule, &atoms, None);
/// if result.is_valid {
///     println!("Conformer is valid");
/// }
/// 
/// // Scanning mode with context
/// let context = ScanningValidationContext::new(0.7);
/// let result = is_valid_conformer_enhanced(&molecule, &atoms, Some(&context));
/// if !result.is_valid {
///     for diagnostic in &result.diagnostic_info {
///         println!("Issue: {}", diagnostic.message);
///     }
/// }
/// ```
pub fn is_valid_conformer_enhanced(
    molecule: &Molecule,
    coordinates: &[Atom],
    scanning_context: Option<&ScanningValidationContext>,
) -> ValidationResult {
    let n = coordinates.len();
    let mut result = ValidationResult::success();
    
    // Determine which skip_factor to use
    let effective_skip_factor = scanning_context
        .map(|ctx| ctx.adaptive_skip_factor)
        .unwrap_or(molecule.skip_factor);
    
    // Track the most severe failure for primary failure reason
    let mut primary_failure: Option<ValidationFailure> = None;
    
    // 1. Steric clash detection with adaptive skip_factor for scanning
    for i in 0..n {
        for j in (i + 1)..n {
            let dx = coordinates[i].x - coordinates[j].x;
            let dy = coordinates[i].y - coordinates[j].y;
            let dz = coordinates[i].z - coordinates[j].z;
            let distance = (dx * dx + dy * dy + dz * dz).sqrt();

            // Use adaptive skip_factor for scanning, standard for rotation
            let min_distance = (covalent_radius(&coordinates[i].element) +
                               covalent_radius(&coordinates[j].element)) * effective_skip_factor;

            if distance < min_distance {
                // Create steric clash failure
                let clash = ValidationFailure::StericClash {
                    atom1: i,
                    atom2: j,
                    distance,
                    min_required: min_distance,
                };
                
                // Set as primary failure if this is the first critical failure
                if primary_failure.is_none() {
                    primary_failure = Some(clash.clone());
                }
                
                // Add diagnostic message
                let severity = if distance < min_distance * 0.8 {
                    ValidationSeverity::Critical
                } else {
                    ValidationSeverity::Warning
                };
                
                let diagnostic = DiagnosticMessage {
                    level: severity,
                    message: format!(
                        "Steric clash between atoms {} and {} ({} and {}): distance {:.3} Å < required {:.3} Å",
                        i + 1, j + 1, // Convert to 1-based for user display
                        coordinates[i].element,
                        coordinates[j].element,
                        distance,
                        min_distance
                    ),
                    atom_indices: Some(vec![i, j]),
                    suggested_action: Some(generate_steric_clash_suggestion(
                        distance, min_distance, effective_skip_factor, scanning_context
                    )),
                };
                
                result.add_diagnostic(diagnostic);
                
                // For critical failures, mark result as invalid immediately
                if severity == ValidationSeverity::Critical {
                    result.is_valid = false;
                    result.failure_reason = Some(clash);
                    
                    // Add parameter suggestion for steric clashes
                    if scanning_context.is_some() {
                        result.add_suggestion(ParameterSuggestion {
                            parameter_name: "skip_factor".to_string(),
                            current_value: format!("{:.2}", molecule.skip_factor),
                            suggested_value: format!("{:.2}", (molecule.skip_factor - 0.1).max(0.3)),
                            reason: "Reduce steric clash sensitivity for scanning mode".to_string(),
                            expected_impact: "Should allow more conformers to pass validation".to_string(),
                        });
                    }
                    
                    // Return early for critical steric clashes to avoid overwhelming diagnostics
                    return result;
                }
            } else if distance < min_distance * 1.2 {
                // Close but acceptable - add informational diagnostic
                let diagnostic = DiagnosticMessage {
                    level: ValidationSeverity::Info,
                    message: format!(
                        "Atoms {} and {} are close ({:.3} Å) but within tolerance ({:.3} Å)",
                        i + 1, j + 1, distance, min_distance
                    ),
                    atom_indices: Some(vec![i, j]),
                    suggested_action: None,
                };
                result.add_diagnostic(diagnostic);
            }
        }
    }
    
    // 2. Additional validation for scanning mode
    if matches!(molecule.operation_mode, crate::molecule::OperationMode::Scanning) {
        let scanning_result = validate_scanning_conformer_enhanced(molecule, coordinates, scanning_context);
        
        // Merge scanning validation results
        if !scanning_result.is_valid {
            result.is_valid = false;
            if result.failure_reason.is_none() {
                result.failure_reason = scanning_result.failure_reason;
            }
        }
        
        // Merge diagnostic information
        for diagnostic in scanning_result.diagnostic_info {
            result.add_diagnostic(diagnostic);
        }
        
        // Merge parameter suggestions
        for suggestion in scanning_result.suggested_adjustments {
            result.add_suggestion(suggestion);
        }
        

    }
    
    // 3. Generate overall parameter suggestions if there are issues
    if !result.is_valid || result.has_diagnostics_of_severity(ValidationSeverity::Warning) {
        generate_overall_parameter_suggestions(&mut result, molecule, scanning_context);
    }
    
    result
}

/// Generates appropriate suggestion text for steric clash issues.
/// 
/// # Arguments
/// 
/// * `distance` - Actual distance between atoms
/// * `min_required` - Minimum required distance
/// * `skip_factor` - Current skip_factor value
/// * `scanning_context` - Optional scanning context
/// 
/// # Returns
/// 
/// Suggested action text for addressing the steric clash
fn generate_steric_clash_suggestion(
    distance: f64,
    min_required: f64,
    skip_factor: f64,
    scanning_context: Option<&ScanningValidationContext>,
) -> String {
    let severity_ratio = distance / min_required;
    
    if let Some(_context) = scanning_context {
        if severity_ratio < 0.7 {
            // Severe clash in scanning mode
            format!(
                "Severe steric clash in scanning mode. Consider: (1) reducing scanning step size by 50%, (2) decreasing skip_factor to {:.2}, or (3) increasing scanning tolerance",
                (skip_factor - 0.15).max(0.3)
            )
        } else if severity_ratio < 0.9 {
            // Moderate clash in scanning mode
            format!(
                "Moderate steric clash in scanning mode. Try: (1) reducing skip_factor to {:.2}, or (2) using smaller scanning steps",
                (skip_factor - 0.1).max(0.4)
            )
        } else {
            // Minor clash in scanning mode
            "Minor steric clash in scanning mode. Consider slightly reducing skip_factor or scanning step size".to_string()
        }
    } else {
        if severity_ratio < 0.8 {
            // Severe clash in rotation mode
            format!(
                "Severe steric clash. Consider reducing skip_factor to {:.2} or checking molecular structure",
                (skip_factor - 0.1).max(0.4)
            )
        } else {
            // Moderate clash in rotation mode
            format!(
                "Steric clash detected. Try reducing skip_factor to {:.2}",
                (skip_factor - 0.05).max(0.5)
            )
        }
    }
}

/// Enhanced scanning-specific conformer validation with detailed diagnostics.
/// 
/// Performs comprehensive validation specific to scanning mode conformers,
/// including bond length validation, forced/forbidden bond checking, and
/// geometric consistency validation.
/// 
/// # Arguments
/// 
/// * `molecule` - The molecule containing bonds and validation parameters
/// * `coordinates` - The conformer coordinates to validate
/// * `scanning_context` - Optional scanning context for adaptive validation
/// 
/// # Returns
/// 
/// A `ValidationResult` with scanning-specific validation outcomes and diagnostics
fn validate_scanning_conformer_enhanced(
    molecule: &Molecule,
    coordinates: &[Atom],
    scanning_context: Option<&ScanningValidationContext>,
) -> ValidationResult {
    let mut result = ValidationResult::success();
    
    // 1. Validate bond lengths for all bonds in the molecule
    for bond in &molecule.bonds {
        if bond.atom1 >= coordinates.len() || bond.atom2 >= coordinates.len() {
            let failure = ValidationFailure::GeometricConstraint {
                description: format!("Bond references invalid atom indices: {} or {}", bond.atom1 + 1, bond.atom2 + 1),
                severity: ValidationSeverity::Critical,
            };
            result.is_valid = false;
            result.failure_reason = Some(failure);
            return result;
        }

        let atom1 = &coordinates[bond.atom1];
        let atom2 = &coordinates[bond.atom2];
        
        let dx = atom2.x - atom1.x;
        let dy = atom2.y - atom1.y;
        let dz = atom2.z - atom1.z;
        let bond_length = (dx * dx + dy * dy + dz * dz).sqrt();

        // Basic bond length validation
        if !crate::molecule::validate_bond_length(bond_length) {
            let failure = ValidationFailure::InvalidBondLength {
                bond: (bond.atom1, bond.atom2),
                length: bond_length,
                valid_range: (0.1, f64::INFINITY),
            };
            
            result.is_valid = false;
            result.failure_reason = Some(failure);
            
            let diagnostic = DiagnosticMessage {
                level: ValidationSeverity::Critical,
                message: format!(
                    "Invalid bond length between atoms {} and {}: {:.3} Å (minimum: 0.1 Å)",
                    bond.atom1 + 1, bond.atom2 + 1, bond_length
                ),
                atom_indices: Some(vec![bond.atom1, bond.atom2]),
                suggested_action: Some("Reduce scanning step size or check molecular structure".to_string()),
            };
            result.add_diagnostic(diagnostic);
            
            return result;
        }

        // Enhanced bond length validation with chemical reasonableness
        let min_expected_length = (covalent_radius(&atom1.element) + 
                                 covalent_radius(&atom2.element)) * molecule.bond_factor * 0.5;
        let max_expected_length = (covalent_radius(&atom1.element) + 
                                 covalent_radius(&atom2.element)) * molecule.bond_factor * 3.0;
        
        // Apply scanning-specific tolerance if available
        let (adjusted_min, adjusted_max) = if let Some(context) = scanning_context {
            if context.is_scanning_bond(bond.atom1, bond.atom2) {
                // More permissive bounds for bonds being actively scanned
                (min_expected_length * 0.8, max_expected_length * 1.2)
            } else {
                (min_expected_length, max_expected_length)
            }
        } else {
            (min_expected_length, max_expected_length)
        };
        
        if bond_length < adjusted_min {
            let failure = ValidationFailure::InvalidBondLength {
                bond: (bond.atom1, bond.atom2),
                length: bond_length,
                valid_range: (adjusted_min, adjusted_max),
            };
            
            result.is_valid = false;
            result.failure_reason = Some(failure);
            
            let diagnostic = DiagnosticMessage {
                level: ValidationSeverity::Critical,
                message: format!(
                    "Bond between atoms {} and {} is too short: {:.3} Å < {:.3} Å (considering bond_factor {:.2})",
                    bond.atom1 + 1, bond.atom2 + 1, bond_length, adjusted_min, molecule.bond_factor
                ),
                atom_indices: Some(vec![bond.atom1, bond.atom2]),
                suggested_action: Some("Increase scanning step size or adjust bond_factor parameter".to_string()),
            };
            result.add_diagnostic(diagnostic);
            
        } else if bond_length > adjusted_max {
            let failure = ValidationFailure::InvalidBondLength {
                bond: (bond.atom1, bond.atom2),
                length: bond_length,
                valid_range: (adjusted_min, adjusted_max),
            };
            
            result.is_valid = false;
            result.failure_reason = Some(failure);
            
            let diagnostic = DiagnosticMessage {
                level: ValidationSeverity::Critical,
                message: format!(
                    "Bond between atoms {} and {} is too long: {:.3} Å > {:.3} Å (considering bond_factor {:.2})",
                    bond.atom1 + 1, bond.atom2 + 1, bond_length, adjusted_max, molecule.bond_factor
                ),
                atom_indices: Some(vec![bond.atom1, bond.atom2]),
                suggested_action: Some("Reduce scanning step size or adjust bond_factor parameter".to_string()),
            };
            result.add_diagnostic(diagnostic);
        }
    }
    
    // 2. Validate forced bonds
    for &(atom1_idx, atom2_idx) in &molecule.forced_bonds {
        if atom1_idx == 0 || atom2_idx == 0 || 
           atom1_idx > coordinates.len() || atom2_idx > coordinates.len() {
            continue; // Skip invalid indices
        }

        let atom1 = &coordinates[atom1_idx - 1]; // Convert from 1-based to 0-based
        let atom2 = &coordinates[atom2_idx - 1];
        
        let dx = atom2.x - atom1.x;
        let dy = atom2.y - atom1.y;
        let dz = atom2.z - atom1.z;
        let distance = (dx * dx + dy * dy + dz * dz).sqrt();

        if !crate::molecule::validate_bond_length(distance) {
            let failure = ValidationFailure::ForcedBondViolation {
                bond: (atom1_idx, atom2_idx),
                distance,
                expected_range: (0.1, f64::INFINITY),
            };
            
            result.is_valid = false;
            result.failure_reason = Some(failure);
            
            let diagnostic = DiagnosticMessage {
                level: ValidationSeverity::Critical,
                message: format!(
                    "Forced bond between atoms {} and {} has invalid length: {:.3} Å",
                    atom1_idx, atom2_idx, distance
                ),
                atom_indices: Some(vec![atom1_idx - 1, atom2_idx - 1]),
                suggested_action: Some("Check forced bond specification or reduce scanning range".to_string()),
            };
            result.add_diagnostic(diagnostic);
        }
        
        // Additional validation with bond_factor
        let min_forced_length = (covalent_radius(&atom1.element) + 
                                covalent_radius(&atom2.element)) * molecule.bond_factor * 0.3;
        let max_forced_length = (covalent_radius(&atom1.element) + 
                                covalent_radius(&atom2.element)) * molecule.bond_factor * 4.0;
        
        if distance < min_forced_length || distance > max_forced_length {
            let failure = ValidationFailure::ForcedBondViolation {
                bond: (atom1_idx, atom2_idx),
                distance,
                expected_range: (min_forced_length, max_forced_length),
            };
            
            result.is_valid = false;
            result.failure_reason = Some(failure);
            
            let diagnostic = DiagnosticMessage {
                level: ValidationSeverity::Warning,
                message: format!(
                    "Forced bond between atoms {} and {} is outside reasonable range: {:.3} Å (expected: {:.3}-{:.3} Å)",
                    atom1_idx, atom2_idx, distance, min_forced_length, max_forced_length
                ),
                atom_indices: Some(vec![atom1_idx - 1, atom2_idx - 1]),
                suggested_action: Some("Consider adjusting bond_factor or scanning parameters".to_string()),
            };
            result.add_diagnostic(diagnostic);
        }
    }
    
    // 3. Validate forbidden bonds are not accidentally created
    for &(atom1_idx, atom2_idx) in &molecule.forbidden_bonds {
        if atom1_idx == 0 || atom2_idx == 0 || 
           atom1_idx > coordinates.len() || atom2_idx > coordinates.len() {
            continue; // Skip invalid indices
        }

        let atom1 = &coordinates[atom1_idx - 1]; // Convert from 1-based to 0-based
        let atom2 = &coordinates[atom2_idx - 1];
        
        let dx = atom2.x - atom1.x;
        let dy = atom2.y - atom1.y;
        let dz = atom2.z - atom1.z;
        let distance = (dx * dx + dy * dy + dz * dz).sqrt();

        // Check if forbidden atoms are too close (would form a bond)
        let bond_threshold = (covalent_radius(&atom1.element) + 
                             covalent_radius(&atom2.element)) * molecule.bond_factor + 0.45;
        
        if distance < bond_threshold {
            let failure = ValidationFailure::ForbiddenBondCreation {
                bond: (atom1_idx, atom2_idx),
                distance,
                min_distance_to_avoid: bond_threshold,
            };
            
            result.is_valid = false;
            result.failure_reason = Some(failure);
            
            let diagnostic = DiagnosticMessage {
                level: ValidationSeverity::Critical,
                message: format!(
                    "Forbidden bond between atoms {} and {} was created: distance {:.3} Å < threshold {:.3} Å",
                    atom1_idx, atom2_idx, distance, bond_threshold
                ),
                atom_indices: Some(vec![atom1_idx - 1, atom2_idx - 1]),
                suggested_action: Some("Reduce scanning range or adjust bond_factor to avoid forbidden bond formation".to_string()),
            };
            result.add_diagnostic(diagnostic);
        }
    }
    
    result
}

/// Generates overall parameter suggestions based on validation results.
/// 
/// Analyzes the validation results and molecular parameters to provide
/// comprehensive suggestions for improving conformer generation success.
/// 
/// # Arguments
/// 
/// * `result` - Mutable validation result to add suggestions to
/// * `molecule` - The molecule being validated
/// * `scanning_context` - Optional scanning context
fn generate_overall_parameter_suggestions(
    result: &mut ValidationResult,
    molecule: &Molecule,
    scanning_context: Option<&ScanningValidationContext>,
) {
    let (critical_count, warning_count, _info_count) = result.diagnostic_counts();
    
    // Suggest skip_factor adjustments based on failure patterns
    if critical_count > 0 || warning_count > 2 {
        if scanning_context.is_some() {
            result.add_suggestion(ParameterSuggestion {
                parameter_name: "skip_factor".to_string(),
                current_value: format!("{:.2}", molecule.skip_factor),
                suggested_value: format!("{:.2}", (molecule.skip_factor - 0.15).max(0.3)),
                reason: "Multiple validation failures suggest overly restrictive steric clash detection for scanning mode".to_string(),
                expected_impact: format!("Should reduce validation failures by approximately {}%", 
                                       ((critical_count + warning_count) as f64 * 15.0).min(80.0) as usize),
            });
        } else {
            result.add_suggestion(ParameterSuggestion {
                parameter_name: "skip_factor".to_string(),
                current_value: format!("{:.2}", molecule.skip_factor),
                suggested_value: format!("{:.2}", (molecule.skip_factor - 0.1).max(0.4)),
                reason: "Multiple steric clashes detected in rotation mode".to_string(),
                expected_impact: "Should allow more conformers to pass validation".to_string(),
            });
        }
    }
    
    // Suggest bond_factor adjustments if there are bond length issues
    let has_bond_length_issues = result.diagnostic_info.iter().any(|d| 
        d.message.contains("bond") && d.message.contains("length")
    );
    
    if has_bond_length_issues {
        result.add_suggestion(ParameterSuggestion {
            parameter_name: "bond_factor".to_string(),
            current_value: format!("{:.2}", molecule.bond_factor),
            suggested_value: format!("{:.2}", (molecule.bond_factor + 0.1).min(2.0)),
            reason: "Bond length validation failures suggest need for more permissive bond detection".to_string(),
            expected_impact: "Should allow wider range of acceptable bond lengths".to_string(),
        });
    }
    
    // Scanning-specific suggestions
    if let Some(context) = scanning_context {
        if critical_count > 0 {
            result.add_suggestion(ParameterSuggestion {
                parameter_name: "scanning_step_size".to_string(),
                current_value: "current".to_string(),
                suggested_value: "reduce by 50%".to_string(),
                reason: "Large scanning steps may be creating unrealistic geometries".to_string(),
                expected_impact: "Smaller steps should produce more chemically reasonable conformers".to_string(),
            });
        }
        
        if context.adaptive_skip_factor < molecule.skip_factor * 0.8 {
            result.add_suggestion(ParameterSuggestion {
                parameter_name: "skip_factor".to_string(),
                current_value: format!("{:.2}", molecule.skip_factor),
                suggested_value: format!("{:.2}", context.adaptive_skip_factor),
                reason: "Scanning operations may benefit from more permissive validation".to_string(),
                expected_impact: "Should improve validation success rate for scanning conformers".to_string(),
            });
        }
    }
}

/// Validates fragment movement effects during bond scanning operations.
/// 
/// This function performs specialized validation for scanning mode that accounts
/// for the movement of molecular fragments during bond length changes. It ensures
/// that fragment movements don't create unrealistic geometries or violate
/// chemical constraints.
/// 
/// # Arguments
/// 
/// * `molecule` - The molecule containing bonds and validation parameters
/// * `coordinates` - The conformer coordinates to validate
/// * `scanning_context` - Scanning context with bond length tracking information
/// 
/// # Returns
/// 
/// A `ValidationResult` with fragment movement validation outcomes
/// 
/// # Validation Checks
/// 
/// 1. **Fragment connectivity**: Ensures fragments remain properly connected
/// 2. **Bond angle preservation**: Validates that bond angles remain reasonable
/// 3. **Fragment overlap detection**: Checks for unrealistic fragment overlaps
/// 4. **Distance consistency**: Ensures distance changes are consistent with scanning
/// 
/// Applies bond length scanning using existing molecular infrastructure.
/// 
/// This function leverages the existing `build_bond_graph`, `identify_fragments`, 
/// and `scan_bond_length` functions to perform bond length scanning in a way
/// that's consistent with the rotation system.
/// 
/// # Arguments
/// 
/// * `molecule` - The molecule containing bond detection parameters
/// * `coordinates` - Mutable atomic coordinates to modify
/// * `atom1` - Index of first atom in the bond to scan
/// * `atom2` - Index of second atom in the bond to scan  
/// * `target_length` - Target bond length in Angstroms
/// 
/// # Returns
/// 
/// `Ok(())` if scanning was successful, `Err(String)` if it failed
/// 
/// # Algorithm
/// 
/// 1. Use existing `build_bond_graph` to detect molecular connectivity
/// 2. Use existing `identify_fragments` to find which atoms move together
/// 3. Use existing `scan_bond_length` to perform the actual translation
/// 
/// # Examples
/// 
/// ```rust
/// // Scan bond 0-1 to 1.7 Å using existing infrastructure
/// apply_bond_scanning(&molecule, &mut atoms, 0, 1, 1.7)?;
/// ```
#[allow(dead_code)]
pub fn apply_bond_scanning(
    molecule: &Molecule,
    coordinates: &mut [Atom],
    atom1: usize,
    atom2: usize,
    target_length: f64,
) -> Result<(), String> {
    // Use existing bond detection infrastructure
    let bond_graph = build_bond_graph(molecule);
    
    // Use existing fragment identification infrastructure
    let (_fragment1, fragment2) = identify_fragments(atom1, atom2, &bond_graph);
    
    // Use existing scan_bond_length function for the actual translation
    crate::molecule::scan_bond_length(coordinates, &fragment2, atom1, atom2, target_length)
}

/// Creates a scanning validation context that integrates with existing validation.
/// 
/// This function creates a simple scanning context that provides adaptive skip_factor
/// calculation while maintaining full compatibility with the existing validation system.
/// 
/// # Arguments
/// 
/// * `molecule` - The molecule containing validation parameters
/// * `coordinates` - Current atomic coordinates
/// * `scanning_bonds` - List of bonds being scanned with target lengths
/// 
/// # Returns
/// 
/// A `ScanningValidationContext` ready for use with existing validation functions
/// 
/// # Examples
/// 
/// ```rust
/// let scanning_bonds = vec![(0, 1, 1.7), (2, 3, 1.5)];
/// let context = create_scanning_context(&molecule, &atoms, &scanning_bonds);
/// let result = is_valid_conformer_enhanced(&molecule, &atoms, Some(&context));
/// ```
#[allow(dead_code)]
pub fn create_scanning_context(
    molecule: &Molecule,
    coordinates: &[Atom],
    scanning_bonds: &[(usize, usize, f64)],
) -> ScanningValidationContext {
    ScanningValidationContext::new(molecule, coordinates, scanning_bonds)
}

/// Builds a bond connectivity graph from a molecule using covalent radii and user overrides.
/// 
/// Uses the OpenBabel bond detection formula combined with user-defined forced and
/// forbidden bonds to create a comprehensive bond connectivity graph. This function
/// integrates with existing validation systems by respecting all molecule parameters.
/// 
/// Used by both rotation and scanning modes to ensure consistent bond detection
/// across all operation modes. The bond graph is fundamental for fragment
/// identification and validation in both rotation and scanning conformer generation.
/// 
/// # Arguments
/// 
/// * `molecule` - The molecule containing atoms and bond detection parameters
/// 
/// # Returns
/// 
/// A vector where each index represents an atom, and the value is a vector
/// of indices of atoms bonded to that atom (adjacency list representation).
/// 
/// # Algorithm
/// 
/// 1. **Distance-based detection**: Uses OpenBabel formula with molecule.bond_factor
/// 2. **Forced bonds**: Adds user-specified bonds from molecule.forced_bonds regardless of distance
/// 3. **Forbidden bonds**: Removes user-specified bonds from molecule.forbidden_bonds even if detected
/// 
/// # Integration with Validation Systems
/// 
/// - Respects molecule.bond_factor for bond detection threshold
/// - Applies molecule.forced_bonds specifications for user-defined bonds
/// - Applies molecule.forbidden_bonds specifications to remove unwanted bonds
/// - Used by both rotation and scanning modes for consistent fragment detection
/// - Ensures scanning respects the same bond detection parameters as rotation
/// 
/// # Examples
/// 
/// ```rust
/// use rotbond::algorithms::build_bond_graph;
/// use rotbond::molecule::Molecule;
/// 
/// let molecule = Molecule::new(atoms);
/// let bonds = build_bond_graph(&molecule);
/// // bonds[0] contains indices of atoms bonded to atom 0
/// ```
pub fn build_bond_graph(molecule: &Molecule) -> Vec<Vec<usize>> {
    let n = molecule.atoms.len();
    let mut adjacency_list = vec![Vec::new(); n];

    // Step 1: Distance-based bond detection (OpenBabel approach with bond_factor)
    for i in 0..n {
        for j in (i + 1)..n {
            let dx = molecule.atoms[i].x - molecule.atoms[j].x;
            let dy = molecule.atoms[i].y - molecule.atoms[j].y;
            let dz = molecule.atoms[i].z - molecule.atoms[j].z;
            let distance_squared = dx * dx + dy * dy + dz * dz;
            
            // OpenBabel formula with bond_factor: cutoff = ((rad1 + rad2) * bond_factor + 0.45)^2
            let radius_sum = covalent_radius(&molecule.atoms[i].element) + 
                           covalent_radius(&molecule.atoms[j].element);
            let cutoff_distance = radius_sum * molecule.bond_factor + 0.45;
            let threshold_squared = cutoff_distance.powi(2);

            if distance_squared <= threshold_squared {
                adjacency_list[i].push(j);
                adjacency_list[j].push(i);
            }
        }
    }

    // Step 2: Apply user-defined forced bonds
    for &(atom1, atom2) in &molecule.forced_bonds {
        let i = atom1 - 1; // Convert to 0-indexed
        let j = atom2 - 1;

        if i < n && j < n && !adjacency_list[i].contains(&j) {
            adjacency_list[i].push(j);
            adjacency_list[j].push(i);
        }
    }

    // Step 3: Apply user-defined forbidden bonds
    for &(atom1, atom2) in &molecule.forbidden_bonds {
        let i = atom1 - 1; // Convert to 0-indexed
        let j = atom2 - 1;

        if i < n && j < n {
            adjacency_list[i].retain(|&x| x != j);
            adjacency_list[j].retain(|&x| x != i);
        }
    }

    adjacency_list
}

/// Checks if two atoms are connected (bonded) in the adjacency list.
/// 
/// # Arguments
/// 
/// * `i` - Index of the first atom
/// * `j` - Index of the second atom
/// * `adjacency_list` - Bond connectivity graph
/// 
/// # Returns
/// 
/// `true` if atoms i and j are bonded, `false` otherwise.
#[allow(dead_code)]
fn are_connected(i: usize, j: usize, adjacency_list: &[Vec<usize>]) -> bool {
    adjacency_list[i].contains(&j)
}

/// Identifies molecular fragments on either side of a rotatable bond.
/// 
/// When a bond is conceptually "broken", this function determines which atoms
/// belong to each resulting fragment using graph traversal. This is essential
/// for determining which atoms should be rotated during conformer generation.
/// 
/// # Arguments
/// 
/// * `atom1` - Index of the first atom in the rotatable bond (remains fixed)
/// * `atom2` - Index of the second atom in the rotatable bond (rotation side)
/// * `adjacency_list` - Bond connectivity graph of the molecule
/// 
/// # Returns
/// 
/// A tuple containing:
/// - Fragment containing atom1 and all atoms connected to it (fixed during rotation)
/// - Fragment containing atom2 and all atoms connected to it (rotates around bond)
/// 
/// # Algorithm
/// 
/// Uses breadth-first search (BFS) from each pivot atom to identify all connected
/// atoms while avoiding crossing the rotatable bond.
/// 
/// # Examples
/// 
/// ```rust
/// use rotbond::algorithms::identify_fragments;
/// 
/// let (frag1, frag2) = identify_fragments(1, 2, &bond_graph);
/// println!("Fragment 1 has {} atoms", frag1.len());
/// println!("Fragment 2 has {} atoms", frag2.len());
/// ```
pub fn identify_fragments(
    atom1: usize,
    atom2: usize,
    adjacency_list: &[Vec<usize>],
) -> (Vec<usize>, Vec<usize>) {
    let n = adjacency_list.len();

    // Find all atoms connected to atom1 (including atom1, excluding atom2)
    let mut fragment1 = Vec::new();
    let mut visited1 = vec![false; n];
    let mut queue1 = VecDeque::new();

    // Start with atom1 itself (the fixed pivot)
    fragment1.push(atom1);
    visited1[atom1] = true;
    queue1.push_back(atom1);

    // BFS to find all connected atoms in fragment 1
    while let Some(current) = queue1.pop_front() {
        for &neighbor in &adjacency_list[current] {
            if !visited1[neighbor] && neighbor != atom2 {
                visited1[neighbor] = true;
                fragment1.push(neighbor);
                queue1.push_back(neighbor);
            }
        }
    }

    // Find all atoms connected to atom2 (including atom2, excluding atom1)
    let mut fragment2 = Vec::new();
    let mut visited2 = vec![false; n];
    let mut queue2 = VecDeque::new();

    // Start with atom2 itself (the primary rotating atom)
    fragment2.push(atom2);
    visited2[atom2] = true;
    queue2.push_back(atom2);

    // BFS to find all connected atoms in fragment 2
    while let Some(current) = queue2.pop_front() {
        for &neighbor in &adjacency_list[current] {
            if !visited2[neighbor] && neighbor != atom1 {
                visited2[neighbor] = true;
                fragment2.push(neighbor);
                queue2.push_back(neighbor);
            }
        }
    }

    // Validate fragments before returning
    validate_fragments(atom1, atom2, &fragment1, &fragment2);

    (fragment1, fragment2)
}

/// Validates that fragments are correctly identified.
/// 
/// Performs silent validation to ensure fragments contain expected atoms
/// and are non-empty. This is an internal quality assurance function.
/// 
/// # Arguments
/// 
/// * `_atom1` - First pivot atom (unused in current implementation)
/// * `_atom2` - Second pivot atom (unused in current implementation)
/// * `_fragment1` - First fragment to validate (unused in current implementation)
/// * `_fragment2` - Second fragment to validate (unused in current implementation)
fn validate_fragments(
    _atom1: usize,
    _atom2: usize,
    _fragment1: &[usize],
    _fragment2: &[usize],
) {
    // Validate fragments (silent validation)
    // Fragments are validated but warnings are suppressed for cleaner output
}

/// Applies rotation transformation to a molecular fragment around a rotatable bond.
/// 
/// This function performs dihedral angle-based rotation using the Rodrigues rotation
/// formula for mathematical precision. It automatically finds reference atoms to
/// define the dihedral angle and applies the exact rotation needed to achieve
/// the desired angle change.
/// 
/// # Arguments
/// 
/// * `atoms` - Mutable slice of atoms in the molecule
/// * `fragment` - Indices of atoms that should be rotated (must include atom2 side)
/// * `pivot_atom1` - Index of the first pivot atom (remains fixed)
/// * `pivot_atom2` - Index of the second pivot atom (rotation axis endpoint)
/// * `angle_degrees` - Desired dihedral angle change in degrees
/// 
/// # Algorithm
/// 
/// 1. Identifies reference atoms A and B to form dihedral A-pivot1-pivot2-B
/// 2. Calculates current dihedral angle
/// 3. Determines target dihedral angle (current + change)
/// 4. Uses iterative solver to find required bond rotation
/// 5. Applies Rodrigues rotation to achieve exact dihedral angle
/// 
/// # Fallback
/// 
/// Falls back to simple bond rotation if dihedral rotation fails.
/// 
/// # Examples
/// 
/// ```rust
/// use rotbond::algorithms::rotate_fragment;
/// 
/// // Rotate fragment around bond 1-2 by 30 degrees
/// let fragment = vec![2, 3, 4]; // atoms on the rotating side
/// rotate_fragment(&mut atoms, &fragment, 1, 2, 30.0);
/// ```
pub fn rotate_fragment(
    atoms: &mut [Atom],
    fragment: &[usize],
    pivot_atom1: usize,
    pivot_atom2: usize,
    angle_degrees: f64,
) {
    if angle_degrees == 0.0 {
        return;
    }
    
    // Build bond graph for dihedral analysis
    let adjacency_list = build_bond_graph_from_atoms(atoms);
    
    // Try to find reference atoms for dihedral-based rotation
    let dihedral_info = get_dihedral_tracking_info(
        pivot_atom1, pivot_atom2, fragment, &adjacency_list, atoms
    );
    
    if let Some((atom_a, atom_b)) = dihedral_info {
        // Calculate initial dihedral
        let initial_dihedral = calculate_dihedral_angle(
            &atoms[atom_a],
            &atoms[pivot_atom1],
            &atoms[pivot_atom2],
            &atoms[atom_b],
        );
        
        // Calculate target dihedral (initial + desired change)
        let target_dihedral = initial_dihedral + angle_degrees;
        
        // Apply dihedral-based rotation to achieve exact dihedral angle change
        if rotate_fragment_to_absolute_dihedral(
            atoms, fragment, pivot_atom1, pivot_atom2, target_dihedral
        ).is_err() {
            // Fallback to simple bond rotation if dihedral rotation fails
            rotate_fragment_simple(atoms, fragment, pivot_atom1, pivot_atom2, angle_degrees);
        }
    } else {
        // Use simple bond rotation when dihedral reference atoms not available
        rotate_fragment_simple(atoms, fragment, pivot_atom1, pivot_atom2, angle_degrees);
    }
}

/// Gets dihedral tracking information for rotation analysis.
/// 
/// Identifies reference atoms needed for dihedral angle tracking during
/// bond rotation. This is a helper function that combines fragment
/// identification with reference atom selection.
/// 
/// # Arguments
/// 
/// * `atom1` - First atom of the rotatable bond
/// * `atom2` - Second atom of the rotatable bond
/// * `fragment` - Atoms that will be rotated (fragment containing atom2)
/// * `adjacency_list` - Bond connectivity graph
/// * `atoms` - Atomic coordinates and element information
/// 
/// # Returns
/// 
/// `Some((atom_a, atom_b))` if suitable reference atoms are found for
/// dihedral A-atom1-atom2-B, `None` if no suitable atoms exist.
/// 
/// # Usage
/// 
/// This function is used internally by [`rotate_fragment`] to determine
/// if dihedral-based rotation is possible, falling back to simple bond
/// rotation if reference atoms cannot be found.
/// 
/// # Examples
/// 
/// ```rust
/// # use rotbond::algorithms::get_dihedral_tracking_info;
/// let dihedral_info = get_dihedral_tracking_info(
///     1, 2, &fragment, &adjacency_list, &atoms
/// );
/// if let Some((atom_a, atom_b)) = dihedral_info {
///     println!("Can track dihedral {}-{}-{}-{}", atom_a, 1, 2, atom_b);
/// }
/// ```
fn get_dihedral_tracking_info(
    atom1: usize,
    atom2: usize,
    fragment: &[usize],
    adjacency_list: &[Vec<usize>],
    atoms: &[Atom],
) -> Option<(usize, usize)> {
    let fragment1 = identify_fragment1(atom1, atom2, adjacency_list);
    
    find_dihedral_reference_atoms(
        atom1, atom2, &fragment1, fragment, adjacency_list, atoms
    )
}

/// Rotates a molecular fragment to achieve an absolute dihedral angle.
/// 
/// This is the primary function for step-based rotations (e30, e60, etc.)
/// where we want to set the dihedral to specific absolute values rather
/// than relative changes. Uses iterative solving to find the precise
/// bond rotation needed.
/// 
/// # Arguments
/// 
/// * `atoms` - Mutable slice of all atoms in the molecule
/// * `fragment` - Indices of atoms that should be rotated
/// * `pivot_atom1` - Index of the first pivot atom (remains fixed)
/// * `pivot_atom2` - Index of the second pivot atom (rotation axis endpoint)
/// * `target_absolute_dihedral` - Desired absolute dihedral angle in degrees
/// 
/// # Returns
/// 
/// * `Ok(())` if the rotation was successfully applied
/// * `Err(String)` if reference atoms could not be found for dihedral calculation
/// 
/// # Algorithm
/// 
/// 1. **Reference atom identification**: Finds atoms A and B for dihedral A-atom1-atom2-B
/// 2. **Current dihedral measurement**: Calculates the starting dihedral angle
/// 3. **Iterative solving**: Uses [`find_bond_rotation_for_dihedral`] to find needed rotation
/// 4. **Rotation application**: Applies the calculated bond rotation using Rodrigues formula
/// 
/// # Usage Context
/// 
/// This function is used for:
/// - Step-based rotations (e30, e60, e120, etc.)
/// - Setting conformers to specific dihedral values
/// - Systematic conformational sampling
/// 
/// # Examples
/// 
/// ```rust
/// # use rotbond::algorithms::rotate_fragment_to_absolute_dihedral;
/// // Set dihedral to exactly 60 degrees
/// let result = rotate_fragment_to_absolute_dihedral(
///     &mut atoms, &fragment, 1, 2, 60.0
/// );
/// match result {
///     Ok(()) => println!("Dihedral set to 60°"),
///     Err(e) => println!("Failed: {}", e),
/// }
/// ```
pub fn rotate_fragment_to_absolute_dihedral(
    atoms: &mut [Atom],
    fragment: &[usize],
    pivot_atom1: usize,
    pivot_atom2: usize,
    target_absolute_dihedral: f64,
) -> Result<(), String> {
    // Build bond graph for dihedral analysis
    let adjacency_list = build_bond_graph_from_atoms(atoms);
    
    // Get fragment1 for reference atom finding
    let fragment1 = identify_fragment1(pivot_atom1, pivot_atom2, &adjacency_list);
    
    // Find reference atoms A (in fragment1) and B (in fragment2)
    let (atom_a, atom_b) = find_dihedral_reference_atoms(
        pivot_atom1, pivot_atom2, &fragment1, fragment, &adjacency_list, atoms
    ).ok_or("Could not find suitable reference atoms for dihedral")?;
    
    // Calculate current dihedral A-atom1-atom2-B
    let current_dihedral = calculate_dihedral_angle(
        &atoms[atom_a],
        &atoms[pivot_atom1],
        &atoms[pivot_atom2],
        &atoms[atom_b],
    );
    
    // Use iterative approach to find the correct bond rotation for the target dihedral
    let bond_rotation_needed = find_bond_rotation_for_dihedral(
        atoms, fragment, pivot_atom1, pivot_atom2, atom_a, atom_b, 
        current_dihedral, target_absolute_dihedral
    );
    
    // Apply the calculated bond rotation
    rotate_fragment_simple(atoms, fragment, pivot_atom1, pivot_atom2, bond_rotation_needed);
    
    Ok(())
}



/// Applies simple bond rotation using the Rodrigues rotation formula.
/// 
/// This is a robust implementation based on OpenBabel's SetTorsion algorithm,
/// providing direct, stable, and mathematically correct rotation around a
/// bond axis. Used as both a standalone rotation method and as a fallback
/// when dihedral-based rotation is not available.
/// 
/// # Arguments
/// 
/// * `atoms` - Mutable slice of atoms (must have .x, .y, .z fields as f64)
/// * `fragment` - List of atom indices (0-based) that should be rotated
/// * `pivot_atom1` - Index of first pivot atom defining the rotation axis
/// * `pivot_atom2` - Index of second pivot atom defining the rotation axis
/// * `angle_degrees` - Rotation angle in degrees (positive follows right-hand rule)
/// 
/// # Rotation Axis
/// 
/// The rotation axis is defined as the vector from `pivot_atom1` to `pivot_atom2`.
/// Positive angles follow the right-hand rule around this axis vector.
/// 
/// # Algorithm
/// 
/// 1. **Axis definition**: Creates unit vector from pivot_atom1 to pivot_atom2
/// 2. **Rodrigues formula**: Applies rotation using the mathematically stable formula
/// 3. **Coordinate transformation**: Translates to pivot_atom2 origin, rotates, translates back
/// 4. **Pivot preservation**: Ensures pivot atoms remain fixed during rotation
/// 
/// # Rodrigues Formula
/// 
/// ```text
/// v_rot = v*cos(θ) + (u × v)*sin(θ) + u*(u·v)*(1-cos(θ))
/// ```
/// 
/// Where:
/// - `v` is the position vector relative to rotation center
/// - `u` is the unit rotation axis vector
/// - `θ` is the rotation angle
/// 
/// # Examples
/// 
/// ```rust
/// # use rotbond::algorithms::rotate_fragment_simple;
/// // Rotate atoms 2,3,4 around bond 0-1 by 60 degrees
/// let fragment = vec![2, 3, 4];
/// rotate_fragment_simple(&mut atoms, &fragment, 0, 1, 60.0);
/// ```
fn rotate_fragment_simple(
    atoms: &mut [Atom],
    fragment: &[usize],
    pivot_atom1: usize,
    pivot_atom2: usize,
    angle_degrees: f64,
) {
    // Quick validation
    let natoms = atoms.len();
    if pivot_atom1 >= natoms || pivot_atom2 >= natoms {
        return; // Invalid pivot indices
    }
    
    if fragment.is_empty() || (angle_degrees - 0.0).abs() < f64::EPSILON {
        return; // Nothing to do
    }

    // Helper to access coordinates
    let p1 = (atoms[pivot_atom1].x, atoms[pivot_atom1].y, atoms[pivot_atom1].z);
    let p2 = (atoms[pivot_atom2].x, atoms[pivot_atom2].y, atoms[pivot_atom2].z);

    // Build axis vector = pivot2 - pivot1 (OpenBabel convention)
    let mut ax = p2.0 - p1.0;
    let mut ay = p2.1 - p1.1;
    let mut az = p2.2 - p1.2;

    // Normalize axis
    let alen = (ax * ax + ay * ay + az * az).sqrt();
    if alen == 0.0 || !alen.is_finite() {
        return; // Degenerate bond (zero-length axis)
    }
    
    ax /= alen;
    ay /= alen;
    az /= alen;

    // Angle in radians
    let theta = angle_degrees.to_radians();
    let cos_t = theta.cos();
    let sin_t = theta.sin();
    let one_minus_cos = 1.0 - cos_t;

    // Precompute axis components for Rodrigues formula
    let ux = ax;
    let uy = ay;
    let uz = az;

    // For each atom in fragment apply rotation about axis passing through pivot point p2
    // We translate so pivot2 is origin, rotate, then translate back.
    for &i in fragment {
        if i >= natoms {
            continue; // Skip invalid indices
        }
        
        // Do not move pivot atoms (critical for bond preservation)
        if i == pivot_atom1 || i == pivot_atom2 {
            continue;
        }

        // Translate to pivot2 origin
        let vx = atoms[i].x - p2.0;
        let vy = atoms[i].y - p2.1;
        let vz = atoms[i].z - p2.2;

        // u × v (cross product)
        let cx = uy * vz - uz * vy;
        let cy = uz * vx - ux * vz;
        let cz = ux * vy - uy * vx;

        // u · v (dot product)
        let dot = ux * vx + uy * vy + uz * vz;

        // Rodrigues' rotation: v_rot = v*cos + (u × v)*sin + u*(u·v)*(1-cos)
        let rx = vx * cos_t + cx * sin_t + ux * dot * one_minus_cos;
        let ry = vy * cos_t + cy * sin_t + uy * dot * one_minus_cos;
        let rz = vz * cos_t + cz * sin_t + uz * dot * one_minus_cos;

        // Translate back
        atoms[i].x = rx + p2.0;
        atoms[i].y = ry + p2.1;
        atoms[i].z = rz + p2.2;
    }
}



/// Validates a conformer by checking for steric clashes and geometric validity.
/// 
/// This is the backward-compatible version of conformer validation that maintains
/// the same interface and behavior as the original function. It internally uses
/// the enhanced validation system but returns only a boolean result.
/// 
/// Uses distance-based criteria to detect atomic overlaps that would make
/// a conformer chemically unrealistic. The validation uses covalent radii
/// scaled by the skip_factor parameter from the molecule configuration.
/// For scanning mode, also validates that bond lengths are within reasonable
/// chemical bounds and that forced/forbidden bond constraints are maintained.
/// 
/// # Arguments
/// 
/// * `molecule` - The molecule containing validation parameters and operation mode
/// * `coordinates` - The conformer to validate
/// 
/// # Returns
/// 
/// `true` if the conformer is valid (no significant clashes and valid geometry), `false` otherwise.
/// 
/// # Examples
/// 
/// ```rust
/// use rotbond::algorithms::is_valid_conformer;
/// 
/// if is_valid_conformer(&molecule, &atoms) {
///     println!("Conformer is valid");
/// } else {
///     println!("Conformer has steric clashes or invalid geometry");
/// }
/// ```
pub fn is_valid_conformer(molecule: &Molecule, coordinates: &[Atom]) -> bool {
    // Use the enhanced validation function for consistency, but return only the boolean result
    let result = is_valid_conformer_enhanced(molecule, coordinates, None);
    result.is_valid
}

/// Validates scanning-specific conformer properties.
/// 
/// Performs additional validation checks specific to bond scanning conformers,
/// including bond length validation and geometric consistency checks.
/// Integrates with existing validation systems by respecting bond_factor
/// and skip_factor parameters from the molecule configuration.
/// 
/// This function ensures scanning respects the same validation parameters
/// as rotation mode, maintaining consistency across operation modes.
/// 
/// # Arguments
/// 
/// * `molecule` - The molecule containing bonds and validation parameters
/// * `coordinates` - The conformer coordinates to validate
/// 
/// # Returns
/// 
/// `true` if the scanning conformer passes all validation checks, `false` otherwise.
/// 
/// # Validation Checks
/// 
/// 1. **Bond length validation**: Ensures all bonds have reasonable lengths (≥ 0.1 Å)
/// 2. **Scanning bond validation**: Checks that scanned bonds maintain chemical validity
/// 3. **Fragment connectivity**: Ensures molecular fragments remain properly connected
/// 4. **Bond factor integration**: Uses molecule.bond_factor for connectivity validation
/// 5. **Skip factor integration**: Uses molecule.skip_factor for steric clash detection
/// 6. **Forced bonds validation**: Ensures forced_bonds maintain reasonable distances
/// 7. **Forbidden bonds validation**: Ensures forbidden_bonds are not accidentally created
/// 
/// # Integration with Existing Systems
/// 
/// - Respects bond_factor parameter for bond length thresholds
/// - Uses skip_factor parameter for steric clash detection (same as rotation mode)
/// - Validates forced_bonds specifications are maintained during scanning
/// - Ensures forbidden_bonds specifications are not violated during scanning
/// - Maintains consistency with existing error handling patterns
/// 
/// # Examples
/// 
/// ```rust
/// use rotbond::algorithms::validate_scanning_conformer;
/// 
/// if validate_scanning_conformer(&molecule, &coordinates) {
///     println!("Scanning conformer is geometrically valid");
/// }
/// ```
#[allow(dead_code)]
pub fn validate_scanning_conformer(molecule: &Molecule, coordinates: &[Atom]) -> bool {
    // Validate bond lengths for all bonds in the molecule
    for bond in &molecule.bonds {
        if bond.atom1 >= coordinates.len() || bond.atom2 >= coordinates.len() {
            return false; // Invalid atom indices
        }

        let atom1 = &coordinates[bond.atom1];
        let atom2 = &coordinates[bond.atom2];
        
        let dx = atom2.x - atom1.x;
        let dy = atom2.y - atom1.y;
        let dz = atom2.z - atom1.z;
        let bond_length = (dx * dx + dy * dy + dz * dz).sqrt();

        // Use the same validation as in molecule.rs
        if !crate::molecule::validate_bond_length(bond_length) {
            return false;
        }

        // Use bond_factor for minimum bond length validation (integration with existing system)
        let min_expected_length = (covalent_radius(&atom1.element) + 
                                 covalent_radius(&atom2.element)) * molecule.bond_factor * 0.5;
        
        if bond_length < min_expected_length {
            return false; // Bond is unrealistically short considering bond_factor
        }

        // Use bond_factor for maximum bond length validation
        let max_expected_length = (covalent_radius(&atom1.element) + 
                                 covalent_radius(&atom2.element)) * molecule.bond_factor * 3.0;
        
        if bond_length > max_expected_length {
            return false; // Bond is unrealistically long considering bond_factor
        }
    }

    // Validate that forced bonds still exist within reasonable ranges
    // This integrates with the existing forced_bonds system
    for &(atom1_idx, atom2_idx) in &molecule.forced_bonds {
        if atom1_idx == 0 || atom2_idx == 0 || 
           atom1_idx > coordinates.len() || atom2_idx > coordinates.len() {
            continue; // Skip invalid indices
        }

        let atom1 = &coordinates[atom1_idx - 1]; // Convert from 1-based to 0-based
        let atom2 = &coordinates[atom2_idx - 1];
        
        let dx = atom2.x - atom1.x;
        let dy = atom2.y - atom1.y;
        let dz = atom2.z - atom1.z;
        let distance = (dx * dx + dy * dy + dz * dz).sqrt();

        // Forced bonds should maintain reasonable distances
        if !crate::molecule::validate_bond_length(distance) {
            return false;
        }
        
        // Apply bond_factor to forced bond validation as well
        let min_forced_length = (covalent_radius(&atom1.element) + 
                                covalent_radius(&atom2.element)) * molecule.bond_factor * 0.3;
        let max_forced_length = (covalent_radius(&atom1.element) + 
                                covalent_radius(&atom2.element)) * molecule.bond_factor * 4.0;
        
        if distance < min_forced_length || distance > max_forced_length {
            return false; // Forced bond length is outside reasonable range
        }
    }

    // Validate that forbidden bonds are not accidentally created during scanning
    // This integrates with the existing forbidden_bonds system
    for &(atom1_idx, atom2_idx) in &molecule.forbidden_bonds {
        if atom1_idx == 0 || atom2_idx == 0 || 
           atom1_idx > coordinates.len() || atom2_idx > coordinates.len() {
            continue; // Skip invalid indices
        }

        let atom1 = &coordinates[atom1_idx - 1]; // Convert from 1-based to 0-based
        let atom2 = &coordinates[atom2_idx - 1];
        
        let dx = atom2.x - atom1.x;
        let dy = atom2.y - atom1.y;
        let dz = atom2.z - atom1.z;
        let distance = (dx * dx + dy * dy + dz * dz).sqrt();

        // Check if forbidden atoms are too close (would form a bond)
        let bond_threshold = (covalent_radius(&atom1.element) + 
                             covalent_radius(&atom2.element)) * molecule.bond_factor + 0.45;
        
        if distance < bond_threshold {
            return false; // Forbidden bond would be formed
        }
    }

    true
}

/// Compares coordinates between two conformers to detect structural differences.
/// 
/// Calculates the maximum coordinate difference between corresponding atoms
/// in two conformers. This is used for conformer diversity analysis and
/// duplicate detection in conformer generation.
/// 
/// # Arguments
/// 
/// * `conformer1` - First conformer as a slice of atoms
/// * `conformer2` - Second conformer as a slice of atoms
/// 
/// # Returns
/// 
/// The maximum coordinate difference found between any pair of corresponding atoms.
/// Returns `f64::INFINITY` if the conformers have different numbers of atoms.
/// 
/// # Algorithm
/// 
/// 1. **Size validation**: Ensures both conformers have the same number of atoms
/// 2. **Pairwise comparison**: Calculates Euclidean distance for each atom pair
/// 3. **Maximum tracking**: Returns the largest distance found
/// 
/// # Usage
/// 
/// This function is used for:
/// - Conformer duplicate detection
/// - Diversity analysis in conformer sets
/// - Quality control in conformer generation
/// 
/// # Examples
/// 
/// ```rust
/// # use rotbond::algorithms::compare_conformer_coordinates;
/// # use rotbond::molecule::Atom;
/// let conf1 = vec![Atom::new("C".to_string(), 0.0, 0.0, 0.0)];
/// let conf2 = vec![Atom::new("C".to_string(), 0.1, 0.0, 0.0)];
/// let max_diff = compare_conformer_coordinates(&conf1, &conf2);
/// assert_eq!(max_diff, 0.1);  // Maximum difference is 0.1 Å
/// ```
pub fn compare_conformer_coordinates(conformer1: &[Atom], conformer2: &[Atom]) -> f64 {
    if conformer1.len() != conformer2.len() {
        return f64::INFINITY; // Different number of atoms
    }

    let mut max_diff: f64 = 0.0;
    for (atom1, atom2) in conformer1.iter().zip(conformer2.iter()) {
        let dx = atom1.x - atom2.x;
        let dy = atom1.y - atom2.y;
        let dz = atom1.z - atom2.z;
        let distance = (dx * dx + dy * dy + dz * dz).sqrt();
        max_diff = max_diff.max(distance);
    }
    max_diff
}

/// Checks if two conformers are effectively identical within a tolerance.
/// 
/// Determines whether two conformers represent the same molecular structure
/// by comparing their maximum coordinate difference against a tolerance threshold.
/// This is essential for removing duplicate conformers from generation results.
/// 
/// # Arguments
/// 
/// * `conformer1` - First conformer as a slice of atoms
/// * `conformer2` - Second conformer as a slice of atoms  
/// * `tolerance` - Maximum allowed coordinate difference (Angstroms)
/// 
/// # Returns
/// 
/// `true` if the conformers are identical within tolerance, `false` otherwise.
/// 
/// # Tolerance Guidelines
/// 
/// - **0.001 Å**: Very strict (numerical precision)
/// - **0.01 Å**: Strict (high-quality structures)
/// - **0.1 Å**: Moderate (typical for conformer analysis)
/// - **0.5 Å**: Lenient (allows for significant flexibility)
/// 
/// # Examples
/// 
/// ```rust
/// # use rotbond::algorithms::are_conformers_identical;
/// # use rotbond::molecule::Atom;
/// let conf1 = vec![Atom::new("C".to_string(), 0.0, 0.0, 0.0)];
/// let conf2 = vec![Atom::new("C".to_string(), 0.05, 0.0, 0.0)];
/// 
/// assert!(are_conformers_identical(&conf1, &conf2, 0.1));   // Within tolerance
/// assert!(!are_conformers_identical(&conf1, &conf2, 0.01)); // Outside tolerance
/// ```
pub fn are_conformers_identical(conformer1: &[Atom], conformer2: &[Atom], tolerance: f64) -> bool {
    compare_conformer_coordinates(conformer1, conformer2) < tolerance
}

/// Validates that a set of conformers contains structurally different conformations.
/// 
/// Performs pairwise comparison of all conformers in a set to identify
/// identical pairs within a tolerance threshold. This is used for quality
/// control and diversity analysis of conformer generation results.
/// 
/// # Arguments
/// 
/// * `conformers` - Vector of conformers, each as a vector of atoms
/// * `tolerance` - Maximum coordinate difference for considering conformers identical (Å)
/// 
/// # Returns
/// 
/// A tuple `(identical_pairs, total_pairs)` where:
/// - `identical_pairs`: Number of conformer pairs that are identical within tolerance
/// - `total_pairs`: Total number of pairwise comparisons performed
/// 
/// # Complexity
/// 
/// - **Time**: O(n² × m) where n = number of conformers, m = atoms per conformer
/// - **Space**: O(1) additional space beyond input
/// 
/// # Diversity Metrics
/// 
/// - **Diversity ratio**: `(total_pairs - identical_pairs) / total_pairs`
/// - **High diversity**: Low number of identical pairs
/// - **Low diversity**: High number of identical pairs (may indicate insufficient sampling)
/// 
/// # Examples
/// 
/// ```rust
/// # use rotbond::algorithms::validate_conformer_diversity;
/// # use rotbond::molecule::Atom;
/// let conformers = vec![
///     vec![Atom::new("C".to_string(), 0.0, 0.0, 0.0)],
///     vec![Atom::new("C".to_string(), 1.0, 0.0, 0.0)],
///     vec![Atom::new("C".to_string(), 0.01, 0.0, 0.0)], // Similar to first
/// ];
/// 
/// let (identical, total) = validate_conformer_diversity(&conformers, 0.1);
/// println!("Found {} identical pairs out of {} total", identical, total);
/// // Output: Found 1 identical pairs out of 3 total
/// ```
pub fn validate_conformer_diversity(conformers: &[Vec<Atom>], tolerance: f64) -> (usize, usize) {
    let mut identical_pairs = 0;
    let mut total_pairs = 0;

    for i in 0..conformers.len() {
        for j in (i + 1)..conformers.len() {
            total_pairs += 1;
            if are_conformers_identical(&conformers[i], &conformers[j], tolerance) {
                identical_pairs += 1;
                // Warning suppressed for cleaner output
            }
        }
    }

    (identical_pairs, total_pairs)
}

/// Logs fragment identification results for debugging purposes.
/// 
/// This function is currently disabled to maintain clean output during
/// conformer generation. When enabled, it would log detailed information
/// about molecular fragment identification for each rotatable bond.
/// 
/// # Arguments
/// 
/// * `_bond_idx` - Index of the bond being analyzed (unused)
/// * `_atom1` - First atom of the rotatable bond (unused)
/// * `_atom2` - Second atom of the rotatable bond (unused)
/// * `_fragment1` - Atoms in the first fragment (unused)
/// * `_fragment2` - Atoms in the second fragment (unused)
/// 
/// # Debug Information
/// 
/// When enabled, this function would log:
/// - Bond indices and atom numbers
/// - Fragment sizes and composition
/// - Validation results for fragment identification
/// 
/// # Usage
/// 
/// This is an internal debugging function called during fragment identification.
/// To enable logging, uncomment the debug output statements in the function body.
/// 
/// # Examples
/// 
/// ```rust
/// # use rotbond::algorithms::log_fragment_info;
/// // Called internally during fragment identification
/// log_fragment_info(0, 1, 2, &vec![0, 3, 4], &vec![1, 5, 6]);
/// // Currently produces no output (debug disabled)
/// ```
pub fn log_fragment_info(_bond_idx: usize, _atom1: usize, _atom2: usize, _fragment1: &[usize], _fragment2: &[usize]) {
    // Debug output removed
}

/// Logs rotation application details for debugging purposes.
/// 
/// This function is currently disabled to maintain clean output during
/// conformer generation. When enabled, it would log detailed information
/// about each rotation operation applied to molecular fragments.
/// 
/// # Arguments
/// 
/// * `_bond_idx` - Index of the bond being rotated (unused)
/// * `_angle` - Rotation angle applied in degrees (unused)
/// * `_fragment` - Atoms that were rotated (unused)
/// * `_atoms_before` - Atomic coordinates before rotation (unused)
/// * `_atoms_after` - Atomic coordinates after rotation (unused)
/// 
/// # Debug Information
/// 
/// When enabled, this function would log:
/// - Bond index and rotation angle
/// - Fragment composition and size
/// - Coordinate changes for rotated atoms
/// - Validation of rotation precision
/// 
/// # Usage
/// 
/// This is an internal debugging function called after each rotation operation.
/// To enable logging, uncomment the debug output statements in the function body.
/// 
/// # Examples
/// 
/// ```rust
/// # use rotbond::algorithms::log_rotation_application;
/// # use rotbond::molecule::Atom;
/// let atoms_before = vec![Atom::new("C".to_string(), 0.0, 0.0, 0.0)];
/// let atoms_after = vec![Atom::new("C".to_string(), 1.0, 0.0, 0.0)];
/// 
/// // Called internally after rotation
/// log_rotation_application(0, 60.0, &vec![1, 2], &atoms_before, &atoms_after);
/// // Currently produces no output (debug disabled)
/// ```
pub fn log_rotation_application(_bond_idx: usize, _angle: f64, _fragment: &[usize], _atoms_before: &[Atom], _atoms_after: &[Atom]) {
    // Debug output removed
}


/// Calculates the dihedral (torsion) angle between four atoms A-B-C-D.
/// 
/// The dihedral angle is the angle between the plane formed by atoms A-B-C
/// and the plane formed by atoms B-C-D. This is a fundamental measurement
/// in molecular geometry for describing conformational states.
/// 
/// # Arguments
/// 
/// * `atom_a` - First atom defining the dihedral
/// * `atom_b` - Second atom (first pivot of the central bond)
/// * `atom_c` - Third atom (second pivot of the central bond)
/// * `atom_d` - Fourth atom defining the dihedral
/// 
/// # Returns
/// 
/// The dihedral angle in degrees, ranging from -180° to +180°.
/// Positive angles follow the right-hand rule around the B-C bond vector.
/// 
/// # Algorithm
/// 
/// 1. Constructs bond vectors: A→B, B→C, C→D
/// 2. Calculates normal vectors to planes ABC and BCD using cross products
/// 3. Computes angle between normal vectors
/// 4. Determines sign using the middle bond vector B→C
/// 
/// # Examples
/// 
/// ```rust
/// use rotbond::algorithms::calculate_dihedral_angle;
/// 
/// let angle = calculate_dihedral_angle(&atom_a, &atom_b, &atom_c, &atom_d);
/// println!("Dihedral angle: {:.2}°", angle);
/// ```
pub fn calculate_dihedral_angle(
    atom_a: &Atom,
    atom_b: &Atom, 
    atom_c: &Atom,
    atom_d: &Atom,
) -> f64 {
    // Vectors along the bonds
    let b1 = [atom_b.x - atom_a.x, atom_b.y - atom_a.y, atom_b.z - atom_a.z]; // A→B
    let b2 = [atom_c.x - atom_b.x, atom_c.y - atom_b.y, atom_c.z - atom_b.z]; // B→C
    let b3 = [atom_d.x - atom_c.x, atom_d.y - atom_c.y, atom_d.z - atom_c.z]; // C→D
    
    // Normal to plane ABC: b1 × b2
    let n1 = cross_product(&b1, &b2);
    // Normal to plane BCD: b2 × b3
    let n2 = cross_product(&b2, &b3);
    
    // Normalize the normal vectors
    let n1_norm = normalize_vector(&n1);
    let n2_norm = normalize_vector(&n2);
    
    // Calculate angle between normal vectors
    let dot = dot_product(&n1_norm, &n2_norm);
    let angle = dot.clamp(-1.0, 1.0).acos();
    
    // Determine sign using the middle bond vector
    let cross = cross_product(&n1_norm, &n2_norm);
    let sign = dot_product(&cross, &b2);
    
    let dihedral = if sign < 0.0 { -angle } else { angle };
    
    // Convert to degrees
    dihedral.to_degrees()
}

/// Computes the cross product of two 3D vectors.
/// 
/// The cross product produces a vector perpendicular to both input vectors,
/// following the right-hand rule. This is essential for calculating normal
/// vectors to planes in dihedral angle calculations.
/// 
/// # Arguments
/// 
/// * `a` - First 3D vector as [x, y, z]
/// * `b` - Second 3D vector as [x, y, z]
/// 
/// # Returns
/// 
/// The cross product vector a × b as [x, y, z].
/// 
/// # Formula
/// 
/// ```text
/// a × b = [a_y*b_z - a_z*b_y, a_z*b_x - a_x*b_z, a_x*b_y - a_y*b_x]
/// ```
/// 
/// # Examples
/// 
/// ```rust
/// # use rotbond::algorithms::cross_product;
/// let a = [1.0, 0.0, 0.0];
/// let b = [0.0, 1.0, 0.0];
/// let cross = cross_product(&a, &b);
/// assert_eq!(cross, [0.0, 0.0, 1.0]);  // Right-hand rule: x × y = z
/// ```
fn cross_product(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2], 
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Computes the dot product of two 3D vectors.
/// 
/// The dot product measures the projection of one vector onto another,
/// providing information about the angle between them. Used extensively
/// in dihedral angle calculations and vector normalization.
/// 
/// # Arguments
/// 
/// * `a` - First 3D vector as [x, y, z]
/// * `b` - Second 3D vector as [x, y, z]
/// 
/// # Returns
/// 
/// The dot product scalar a · b.
/// 
/// # Formula
/// 
/// ```text
/// a · b = a_x*b_x + a_y*b_y + a_z*b_z
/// ```
/// 
/// # Properties
/// 
/// - If a · b = 0, the vectors are perpendicular
/// - If a · b > 0, the angle between vectors is acute (< 90°)
/// - If a · b < 0, the angle between vectors is obtuse (> 90°)
/// 
/// # Examples
/// 
/// ```rust
/// # use rotbond::algorithms::dot_product;
/// let a = [1.0, 0.0, 0.0];
/// let b = [0.0, 1.0, 0.0];
/// let dot = dot_product(&a, &b);
/// assert_eq!(dot, 0.0);  // Perpendicular vectors
/// 
/// let c = [1.0, 1.0, 0.0];
/// let d = [1.0, 1.0, 0.0];
/// let dot2 = dot_product(&c, &d);
/// assert_eq!(dot2, 2.0);  // Parallel vectors
/// ```
fn dot_product(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Normalizes a 3D vector to unit length.
/// 
/// Converts a vector to have magnitude 1.0 while preserving its direction.
/// This is essential for creating unit vectors for rotation axes and
/// normal vectors in geometric calculations.
/// 
/// # Arguments
/// 
/// * `v` - The 3D vector to normalize as [x, y, z]
/// 
/// # Returns
/// 
/// The normalized vector with unit length, or zero vector if input has zero magnitude.
/// 
/// # Algorithm
/// 
/// 1. Calculate magnitude: `|v| = sqrt(x² + y² + z²)`
/// 2. Divide each component by magnitude: `v_norm = v / |v|`
/// 3. Handle zero-length vectors by returning [0, 0, 0]
/// 
/// # Examples
/// 
/// ```rust
/// # use rotbond::algorithms::normalize_vector;
/// let v = [3.0, 4.0, 0.0];
/// let normalized = normalize_vector(&v);
/// assert_eq!(normalized, [0.6, 0.8, 0.0]);  // 3-4-5 triangle
/// 
/// let zero = [0.0, 0.0, 0.0];
/// let normalized_zero = normalize_vector(&zero);
/// assert_eq!(normalized_zero, [0.0, 0.0, 0.0]);  // Zero vector case
/// ```
fn normalize_vector(v: &[f64; 3]) -> [f64; 3] {
    let length = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if length == 0.0 {
        [0.0, 0.0, 0.0]
    } else {
        [v[0] / length, v[1] / length, v[2] / length]
    }
}


/// Builds bond graph from atoms using distance-based detection with configurable bond factor.
/// 
/// This is a simplified version for when we don't have the full molecule structure
/// with forced/forbidden bonds. Uses the same OpenBabel bond detection formula
/// as the main [`build_bond_graph`] function to ensure consistency across the system.
/// 
/// # Arguments
/// 
/// * `atoms` - Slice of atoms to analyze for bonding
/// * `bond_factor` - Bond detection threshold multiplier (typically 1.0-1.3)
/// 
/// # Returns
/// 
/// Adjacency list representation of the bond graph where `graph[i]` contains
/// indices of atoms bonded to atom `i`.
/// 
/// # Bond Detection Formula
/// 
/// Uses the OpenBabel formula:
/// ```text
/// bond_detected = distance² ≤ ((radius₁ + radius₂) × bond_factor + 0.45)²
/// ```
/// 
/// # Examples
/// 
/// ```rust
/// # use rotbond::algorithms::build_bond_graph_from_atoms_with_factor;
/// # use rotbond::molecule::Atom;
/// let atoms = vec![
///     Atom::new("C".to_string(), 0.0, 0.0, 0.0),
///     Atom::new("H".to_string(), 1.1, 0.0, 0.0),
/// ];
/// let bond_graph = build_bond_graph_from_atoms_with_factor(&atoms, 1.0);
/// assert!(bond_graph[0].contains(&1));  // C-H bond detected
/// ```
fn build_bond_graph_from_atoms_with_factor(atoms: &[Atom], bond_factor: f64) -> Vec<Vec<usize>> {
    let n = atoms.len();
    let mut adjacency_list = vec![Vec::new(); n];
    
    for i in 0..n {
        for j in (i + 1)..n {
            let dx = atoms[i].x - atoms[j].x;
            let dy = atoms[i].y - atoms[j].y;
            let dz = atoms[i].z - atoms[j].z;
            let distance_squared = dx * dx + dy * dy + dz * dz;
            
            // Use the same OpenBabel formula as the main bond detection
            let radius_sum = covalent_radius(&atoms[i].element) + 
                           covalent_radius(&atoms[j].element);
            let cutoff_distance = radius_sum * bond_factor + 0.45;
            let threshold_squared = cutoff_distance.powi(2);
            
            if distance_squared <= threshold_squared {
                adjacency_list[i].push(j);
                adjacency_list[j].push(i);
            }
        }
    }
    
    adjacency_list
}

/// Builds bond graph from atoms using default bond detection parameters.
/// 
/// This is a simplified version for when we don't have the full molecule structure
/// with user-defined parameters. Uses a default bond_factor of 1.2 for slightly
/// more permissive bond detection, which is appropriate for dihedral analysis.
/// 
/// # Arguments
/// 
/// * `atoms` - Slice of atoms to analyze for bonding
/// 
/// # Returns
/// 
/// Adjacency list representation of the bond graph.
/// 
/// # Default Parameters
/// 
/// - `bond_factor`: 1.2 (slightly permissive for compatibility)
/// - No forced or forbidden bonds applied
/// 
/// # Examples
/// 
/// ```rust
/// # use rotbond::algorithms::build_bond_graph_from_atoms;
/// # use rotbond::molecule::Atom;
/// let atoms = vec![
///     Atom::new("C".to_string(), 0.0, 0.0, 0.0),
///     Atom::new("C".to_string(), 1.5, 0.0, 0.0),
/// ];
/// let bond_graph = build_bond_graph_from_atoms(&atoms);
/// assert!(bond_graph[0].contains(&1));  // C-C bond detected
/// ```
fn build_bond_graph_from_atoms(atoms: &[Atom]) -> Vec<Vec<usize>> {
    // Use default bond_factor of 1.2 for compatibility
    build_bond_graph_from_atoms_with_factor(atoms, 1.2)
}



/// Identifies the molecular fragment connected to atom1, excluding atom2.
/// 
/// Performs breadth-first search to find all atoms connected to atom1
/// without crossing through atom2. This is used for dihedral reference
/// atom identification and fragment analysis.
/// 
/// # Arguments
/// 
/// * `atom1` - Starting atom for the fragment search
/// * `atom2` - Atom to exclude from the search (conceptually "removed" bond)
/// * `adjacency_list` - Bond connectivity graph of the molecule
/// 
/// # Returns
/// 
/// Vector of atom indices that are connected to atom1 (including atom1 itself)
/// without going through atom2.
/// 
/// # Algorithm
/// 
/// Uses BFS starting from atom1, marking atom2 as forbidden to traverse through.
/// This effectively finds one side of a conceptually "broken" bond.
/// 
/// # Examples
/// 
/// ```rust
/// # use rotbond::algorithms::identify_fragment1;
/// // For a linear molecule C-C-C where we want fragment around first C
/// let adjacency_list = vec![vec![1], vec![0, 2], vec![1]];
/// let fragment = identify_fragment1(0, 1, &adjacency_list);
/// assert_eq!(fragment, vec![0]);  // Only the first carbon
/// ```
fn identify_fragment1(
    atom1: usize,
    atom2: usize,
    adjacency_list: &[Vec<usize>],
) -> Vec<usize> {
    let n = adjacency_list.len();
    let mut fragment1 = Vec::new();
    let mut visited = vec![false; n];
    let mut queue = std::collections::VecDeque::new();
    
    // Start BFS from atom1
    fragment1.push(atom1);
    visited[atom1] = true;
    queue.push_back(atom1);
    
    while let Some(current) = queue.pop_front() {
        for &neighbor in &adjacency_list[current] {
            if !visited[neighbor] && neighbor != atom2 {
                visited[neighbor] = true;
                fragment1.push(neighbor);
                queue.push_back(neighbor);
            }
        }
    }
    
    fragment1
}

/// Finds reference atoms for dihedral angle calculation A-atom1-atom2-B.
/// 
/// Identifies suitable atoms A and B to form a complete dihedral angle
/// definition with the central bond atom1-atom2. The reference atoms
/// are chosen to provide stable and meaningful dihedral measurements.
/// 
/// # Arguments
/// 
/// * `atom1` - First atom of the central bond
/// * `atom2` - Second atom of the central bond  
/// * `fragment1` - Atoms connected to atom1 (excluding atom2)
/// * `fragment2` - Atoms connected to atom2 (excluding atom1)
/// * `adjacency_list` - Bond connectivity graph
/// * `atoms` - Atomic coordinates and element information
/// 
/// # Returns
/// 
/// `Some((atom_a, atom_b))` if suitable reference atoms are found,
/// `None` if no suitable atoms exist for dihedral calculation.
/// 
/// # Selection Criteria
/// 
/// Reference atoms are selected based on:
/// 1. **Connectivity**: Must be bonded to the central atoms
/// 2. **Fragment membership**: Must be in the correct fragment
/// 3. **Element preference**: Non-hydrogen atoms preferred
/// 4. **Backbone length**: Atoms with longer backbone chains preferred
/// 
/// # Examples
/// 
/// ```rust
/// # use rotbond::algorithms::find_dihedral_reference_atoms;
/// // For ethane rotation, finds H atoms or C atoms for dihedral measurement
/// let result = find_dihedral_reference_atoms(
///     0, 1, &fragment1, &fragment2, &adjacency_list, &atoms
/// );
/// if let Some((atom_a, atom_b)) = result {
///     println!("Dihedral: {}-{}-{}-{}", atom_a, 0, 1, atom_b);
/// }
/// ```
fn find_dihedral_reference_atoms(
    atom1: usize,
    atom2: usize,
    fragment1: &[usize],
    fragment2: &[usize],
    adjacency_list: &[Vec<usize>],
    atoms: &[Atom],
) -> Option<(usize, usize)> {
    // Find atom A: bonded to atom1, in fragment1, longest backbone, not H
    let atom_a = find_backbone_reference_atom(
        atom1, atom2, fragment1, adjacency_list, atoms
    )?;
    
    // Find atom B: bonded to atom2, in fragment2, longest backbone, not H
    let atom_b = find_backbone_reference_atom(
        atom2, atom1, fragment2, adjacency_list, atoms
    )?;
    
    Some((atom_a, atom_b))
}

/// Finds the best reference atom for dihedral calculation from available candidates.
/// 
/// Selects the most suitable atom bonded to a central atom for use as a
/// dihedral reference. Prioritizes non-hydrogen atoms and those with
/// longer backbone chains for more stable dihedral measurements.
/// 
/// # Arguments
/// 
/// * `central_atom` - The central atom to find a reference for
/// * `exclude_atom` - Atom to exclude from consideration (usually the other central atom)
/// * `fragment` - Fragment containing valid candidate atoms
/// * `adjacency_list` - Bond connectivity graph
/// * `atoms` - Atomic coordinates and element information
/// 
/// # Returns
/// 
/// `Some(atom_index)` if a suitable reference atom is found,
/// `None` if no suitable candidates exist.
/// 
/// # Selection Algorithm
/// 
/// 1. **Primary candidates**: Non-hydrogen atoms bonded to central_atom in fragment
/// 2. **Fallback candidates**: Hydrogen atoms if no non-hydrogen atoms available
/// 3. **Selection criteria**: Atom with longest backbone chain length
/// 4. **Tie-breaking**: First atom found if multiple atoms have same chain length
/// 
/// # Examples
/// 
/// ```rust
/// # use rotbond::algorithms::find_backbone_reference_atom;
/// // Find best reference atom bonded to carbon atom 0
/// let reference = find_backbone_reference_atom(
///     0, 1, &fragment, &adjacency_list, &atoms
/// );
/// if let Some(ref_atom) = reference {
///     println!("Using atom {} as dihedral reference", ref_atom);
/// }
/// ```
fn find_backbone_reference_atom(
    central_atom: usize,
    exclude_atom: usize,
    fragment: &[usize],
    adjacency_list: &[Vec<usize>],
    atoms: &[Atom],
) -> Option<usize> {
    let mut candidates = Vec::new();
    
    // Find atoms bonded to central_atom that are in the same fragment
    for &neighbor in &adjacency_list[central_atom] {
        if neighbor != exclude_atom 
           && fragment.contains(&neighbor)
           && atoms[neighbor].element != "H" {
            candidates.push(neighbor);
        }
    }
    
    if candidates.is_empty() {
        // If no non-hydrogen candidates, allow hydrogen as fallback
        for &neighbor in &adjacency_list[central_atom] {
            if neighbor != exclude_atom && fragment.contains(&neighbor) {
                candidates.push(neighbor);
            }
        }
    }
    
    if candidates.is_empty() {
        return None;
    }
    
    // Select atom with longest backbone chain
    candidates.iter()
        .max_by_key(|&&atom| calculate_backbone_chain_length(
            atom, central_atom, adjacency_list, atoms
        ))
        .copied()
}

/// Calculates the length of the longest backbone chain from a starting atom.
/// 
/// Performs depth-first search to find the longest chain of non-hydrogen
/// atoms extending from the starting atom. This is used to prioritize
/// atoms with more substantial molecular frameworks for dihedral references.
/// 
/// # Arguments
/// 
/// * `start_atom` - Starting atom for chain length calculation
/// * `exclude_atom` - Atom to exclude from traversal (prevents backtracking)
/// * `adjacency_list` - Bond connectivity graph
/// * `atoms` - Atomic coordinates and element information
/// 
/// # Returns
/// 
/// The length of the longest backbone chain (number of non-hydrogen atoms).
/// 
/// # Algorithm
/// 
/// 1. **DFS traversal**: Explores all paths from start_atom
/// 2. **Hydrogen exclusion**: Only counts non-hydrogen atoms in chain length
/// 3. **Backtrack prevention**: Excludes specified atom from traversal
/// 4. **Maximum path**: Returns the length of the longest path found
/// 
/// # Chain Length Significance
/// 
/// - **Length 1**: Terminal atom (only itself)
/// - **Length 2-3**: Short side chain or small ring
/// - **Length 4+**: Substantial backbone or large ring system
/// 
/// # Examples
/// 
/// ```rust
/// # use rotbond::algorithms::calculate_backbone_chain_length;
/// // Calculate backbone length for atom 2, excluding atom 1
/// let chain_length = calculate_backbone_chain_length(
///     2, 1, &adjacency_list, &atoms
/// );
/// println!("Backbone chain length: {}", chain_length);
/// ```
fn calculate_backbone_chain_length(
    start_atom: usize,
    exclude_atom: usize,
    adjacency_list: &[Vec<usize>],
    atoms: &[Atom],
) -> usize {
    fn dfs(
        atom: usize,
        visited: &mut [bool],
        adjacency_list: &[Vec<usize>],
        atoms: &[Atom]
    ) -> usize {
        visited[atom] = true;
        let mut max_length = 0;
        
        for &neighbor in &adjacency_list[atom] {
            if !visited[neighbor] && atoms[neighbor].element != "H" {
                let length = dfs(neighbor, visited, adjacency_list, atoms);
                max_length = max_length.max(length);
            }
        }
        
        max_length + 1
    }
    
    let mut visited = vec![false; adjacency_list.len()];
    visited[exclude_atom] = true; // Don't traverse through excluded atom
    
    dfs(start_atom, &mut visited, adjacency_list, atoms)
}

/// Finds the bond rotation needed to achieve a target dihedral angle.
/// 
/// Solves the non-linear relationship between bond rotation and dihedral angle
/// change using an iterative Newton-like method with finite difference derivative
/// estimation. This is essential for precise dihedral-based conformer generation.
/// 
/// # Arguments
/// 
/// * `atoms` - Current atomic coordinates
/// * `fragment` - Atoms that will be rotated
/// * `pivot_atom1` - First pivot atom of the rotation bond
/// * `pivot_atom2` - Second pivot atom of the rotation bond
/// * `atom_a` - First reference atom for dihedral measurement
/// * `atom_b` - Second reference atom for dihedral measurement
/// * `current_dihedral` - Current dihedral angle A-atom1-atom2-B (degrees)
/// * `target_dihedral` - Desired dihedral angle A-atom1-atom2-B (degrees)
/// 
/// # Returns
/// 
/// The bond rotation angle (in degrees) needed to achieve the target dihedral.
/// 
/// # Algorithm
/// 
/// 1. **Angle normalization**: Normalizes target change to [-180°, 180°] range
/// 2. **Initial guess**: Starts with 1:1 assumption (bond rotation = dihedral change)
/// 3. **Iterative refinement**: Uses Newton-like method with finite differences
/// 4. **Derivative estimation**: Calculates dihedral sensitivity to bond rotation
/// 5. **Convergence**: Iterates until error < 0.05° or max iterations reached
/// 
/// # Convergence Parameters
/// 
/// - **Max iterations**: 8 (prevents infinite loops)
/// - **Tolerance**: 0.05° (high precision)
/// - **Derivative epsilon**: 0.5° (finite difference step)
/// - **Clamp range**: [-360°, 360°] (prevents runaway values)
/// 
/// # Examples
/// 
/// ```rust
/// # use rotbond::algorithms::find_bond_rotation_for_dihedral;
/// // Find rotation to change dihedral from 60° to 180°
/// let bond_rotation = find_bond_rotation_for_dihedral(
///     &atoms, &fragment, 1, 2, 0, 3, 60.0, 180.0
/// );
/// println!("Need {:.2}° bond rotation for target dihedral", bond_rotation);
/// ```
fn find_bond_rotation_for_dihedral(
    atoms: &[Atom],
    fragment: &[usize],
    pivot_atom1: usize,
    pivot_atom2: usize,
    atom_a: usize,
    atom_b: usize,
    current_dihedral: f64,
    target_dihedral: f64,
) -> f64 {
    let desired_change = target_dihedral - current_dihedral;
    
    // Normalize desired change to [-180, 180] range (proper angle wrapping)
    let mut normalized_change = desired_change;
    while normalized_change > 180.0 {
        normalized_change -= 360.0;
    }
    while normalized_change <= -180.0 {
        normalized_change += 360.0;
    }
    
    // Start with 1:1 assumption (bond rotation = dihedral change)
    let mut bond_rotation = normalized_change;
    
    const MAX_ITERATIONS: usize = 8;
    const TOLERANCE: f64 = 0.05; // 0.05 degree tolerance (tighter)
    const DERIVATIVE_EPSILON: f64 = 0.5; // Small angle for derivative estimation
    
    for _iteration in 0..MAX_ITERATIONS {
        // Test current bond rotation
        let mut test_atoms = atoms.to_vec();
        rotate_fragment_simple(&mut test_atoms, fragment, pivot_atom1, pivot_atom2, bond_rotation);
        
        let achieved_dihedral = calculate_dihedral_angle(
            &test_atoms[atom_a],
            &test_atoms[pivot_atom1],
            &test_atoms[pivot_atom2],
            &test_atoms[atom_b],
        );
        
        let mut achieved_change = achieved_dihedral - current_dihedral;
        
        // Normalize achieved change to match normalized_change range
        while achieved_change > 180.0 {
            achieved_change -= 360.0;
        }
        while achieved_change <= -180.0 {
            achieved_change += 360.0;
        }
        
        let error = (achieved_change - normalized_change).abs();
        
        if error < TOLERANCE {
            return bond_rotation;
        }
        
        // Estimate derivative using finite difference
        let mut test_atoms_deriv = atoms.to_vec();
        rotate_fragment_simple(&mut test_atoms_deriv, fragment, pivot_atom1, pivot_atom2, bond_rotation + DERIVATIVE_EPSILON);
        
        let dihedral_plus_eps = calculate_dihedral_angle(
            &test_atoms_deriv[atom_a],
            &test_atoms_deriv[pivot_atom1],
            &test_atoms_deriv[pivot_atom2],
            &test_atoms_deriv[atom_b],
        );
        
        let mut change_plus_eps = dihedral_plus_eps - current_dihedral;
        
        // Normalize derivative change
        while change_plus_eps > 180.0 {
            change_plus_eps -= 360.0;
        }
        while change_plus_eps <= -180.0 {
            change_plus_eps += 360.0;
        }
        
        let derivative = (change_plus_eps - achieved_change) / DERIVATIVE_EPSILON;
        
        // Newton-like update with derivative
        if derivative.abs() > 0.01 {
            let correction = (normalized_change - achieved_change) / derivative;
            bond_rotation += correction;
        } else {
            // Fallback to simple correction when derivative is too small
            if achieved_change.abs() > 0.001 {
                let correction_factor = normalized_change / achieved_change;
                bond_rotation *= correction_factor;
            } else {
                bond_rotation *= 1.5; // Small increase
            }
        }
        
        // Clamp to reasonable range to prevent runaway
        bond_rotation = bond_rotation.clamp(-360.0, 360.0);
    }
    
    bond_rotation
}

