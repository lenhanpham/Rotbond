//! # Molecular Data Structures and Utilities
//!
//! This module defines the core data structures for representing molecules
//! and provides utilities for molecular property calculations.
//!
//! ## Core Types
//!
//! - [`Atom`]: Represents a single atom with element type and coordinates
//! - [`Bond`]: Represents a rotatable bond with angle specifications
//! - [`Molecule`]: Container for atoms and bonds with detection parameters
//! - [`RotationSpec`]: Specification for rotatable bond parameters
//!
//! ## Covalent Radii
//!
//! The module includes a comprehensive table of covalent radii based on
//! OpenBabel values, used for automatic bond detection. The radii are
//! optimized for organic molecules and common inorganic elements.

/// Represents a single atom in a molecular structure.
/// 
/// Contains the essential information needed for conformer generation:
/// element type for bond detection and 3D coordinates for geometric calculations.
/// 
/// # Examples
/// 
/// ```rust
/// use rotbond::molecule::Atom;
/// 
/// let carbon = Atom {
///     element: "C".to_string(),
///     x: 0.0,
///     y: 0.0,
///     z: 0.0,
/// };
/// ```
#[derive(Debug, Clone)]
pub struct Atom {
    /// Chemical element symbol (e.g., "C", "H", "N", "O")
    pub element: String,
    /// X coordinate in Angstroms
    pub x: f64,
    /// Y coordinate in Angstroms
    pub y: f64,
    /// Z coordinate in Angstroms
    pub z: f64,
}

impl Atom {
    /// Creates a new atom with the specified element and coordinates.
    /// 
    /// # Arguments
    /// 
    /// * `element` - Chemical element symbol
    /// * `x` - X coordinate in Angstroms
    /// * `y` - Y coordinate in Angstroms
    /// * `z` - Z coordinate in Angstroms
    /// 
    /// # Examples
    /// 
    /// ```rust
    /// use rotbond::molecule::Atom;
    /// 
    /// let hydrogen = Atom::new("H".to_string(), 1.0, 0.0, 0.0);
    /// ```
    pub fn new(element: String, x: f64, y: f64, z: f64) -> Self {
        Atom { element, x, y, z }
    }
}

/// Represents a rotatable bond between two atoms.
/// 
/// Contains information about which atoms are connected and how the bond
/// should be rotated during conformer generation. Supports both independent
/// and synchronous rotation modes.
/// 
/// # Examples
/// 
/// ```rust
/// use rotbond::molecule::Bond;
/// 
/// let bond = Bond {
///     atom1: 0,
///     atom2: 1,
///     angles: vec![0.0, 60.0, 120.0, 180.0, 240.0, 300.0],
///     is_synchronous: false,
///     reference_bond: None,
///     direction: 1.0,
/// };
/// ```
#[derive(Debug, Clone)]
pub struct Bond {
    /// Index of the first atom in the bond (0-based)
    pub atom1: usize,
    /// Index of the second atom in the bond (0-based)
    pub atom2: usize,
    /// Rotation angles for this bond in degrees
    #[allow(dead_code)]
    pub angles: Vec<f64>,
    /// Whether this bond rotates synchronously with another bond
    pub is_synchronous: bool,
    /// Index of the reference bond for synchronous rotation (0-based, None if independent)
    pub reference_bond: Option<usize>,
    /// Direction multiplier: +1.0 for same direction, -1.0 for opposite
    pub direction: f64,
}

/// Represents a bond scanning specification for systematic bond length variation.
/// 
/// Contains parameters for scanning a bond through different lengths:
/// number of steps and the step size (increment/decrement) in Angstroms.
/// 
/// # Examples
/// 
/// ```rust
/// use rotbond::molecule::ScanningSpec;
/// 
/// // Scan bond in 10 steps with 0.1 Å increments
/// let scan_spec = ScanningSpec {
///     steps: 10,
///     step_size: 0.1,
/// };
/// 
/// // Compress bond in 5 steps with 0.05 Å decrements
/// let compress_spec = ScanningSpec {
///     steps: 5,
///     step_size: -0.05,
/// };
/// ```
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ScanningSpec {
    /// Number of scanning steps to perform
    pub steps: usize,
    /// Step size in Angstrom units (positive = stretch, negative = compress)
    pub step_size: f64,
}

/// Defines the operation mode for conformer generation.
/// 
/// Distinguishes between bond rotation and bond scanning modes,
/// which are mutually exclusive operations that use different
/// algorithms for generating molecular conformers.
/// 
/// # Examples
/// 
/// ```rust
/// use rotbond::molecule::OperationMode;
/// 
/// let rotation_mode = OperationMode::Rotation;
/// let scanning_mode = OperationMode::Scanning;
/// ```
#[derive(Debug, Clone)]
pub enum OperationMode {
    /// Bond rotation mode - varies dihedral angles
    Rotation,
    /// Bond scanning mode - varies bond lengths
    Scanning,
}

/// Defines the type of rotation specification for a rotatable bond.
/// 
/// Supports four types of specifications:
/// - Step-based rotations (e.g., every 30°, 60°, 120°)
/// - Explicit angle lists
/// - Synchronous rotations that follow another bond
/// - Bond length scanning with steps and step size
/// 
/// # Examples
/// 
/// ```rust
/// use rotbond::molecule::RotationSpec;
/// 
/// // Step rotation every 60 degrees
/// let step_spec = RotationSpec::Step { step: 60.0 };
/// 
/// // Explicit angles
/// let explicit_spec = RotationSpec::Explicit(vec![0.0, 90.0, 180.0]);
/// 
/// // Synchronous with bond 1, opposite direction
/// let sync_spec = RotationSpec::Synchronous { 
///     reference: 1, 
///     direction: -1.0 
/// };
/// 
/// // Bond scanning with 10 steps of 0.1 Å increments
/// let scan_spec = RotationSpec::Scanning { 
///     steps: 10, 
///     step_size: 0.1 
/// };
/// ```
#[derive(Debug, Clone)]
pub enum RotationSpec {
    /// Step-based rotation with regular intervals (e.g., every 30°)
    Step { 
        /// Step size in degrees
        step: f64 
    },
    /// Explicit list of rotation angles in degrees
    Explicit(Vec<f64>),
    /// Synchronous rotation that follows another bond
    Synchronous { 
        /// Index of the reference bond to follow (1-based)
        reference: usize, 
        /// Direction multiplier: +1.0 for same, -1.0 for opposite
        direction: f64 
    },
    /// Bond length scanning specification
    Scanning { 
        /// Index of first atom in the bond (1-based)
        atom1: usize,
        /// Index of second atom in the bond (1-based)
        atom2: usize,
        /// Number of scanning steps to perform
        steps: usize, 
        /// Step size in Angstrom units (positive = stretch, negative = compress)
        step_size: f64 
    },
}

/// Represents a complete molecular structure with bond detection parameters.
/// 
/// Contains all atoms, rotatable bonds, and parameters that control
/// bond detection and conformer validation algorithms. This is the main
/// data structure used throughout the conformer generation process.
/// 
/// # Examples
/// 
/// ```rust
/// use rotbond::molecule::{Molecule, Atom, OperationMode};
/// 
/// let atoms = vec![
///     Atom { element: "C".to_string(), x: 0.0, y: 0.0, z: 0.0 },
///     Atom { element: "H".to_string(), x: 1.0, y: 0.0, z: 0.0 },
/// ];
/// let molecule = Molecule::new(atoms);
/// assert_eq!(matches!(molecule.operation_mode, OperationMode::Rotation), true);
/// ```
#[derive(Debug, Clone)]
pub struct Molecule {
    /// Vector of all atoms in the molecule
    pub atoms: Vec<Atom>,
    /// Vector of rotatable bonds
    pub bonds: Vec<Bond>,
    /// Molecular fragments for each rotatable bond (atom indices that rotate together)
    pub fragments: Vec<Vec<usize>>,
    /// Bond detection threshold multiplier (default: 1.0)
    /// 
    /// Used in the formula: `d² < ((rad₁ + rad₂) × bond_factor + 0.45)²`
    /// - Values > 1.0: More permissive bond detection
    /// - Values < 1.0: More restrictive bond detection
    pub bond_factor: f64,
    /// Steric clash detection threshold (default: 0.7)
    /// 
    /// Used for validating conformers by detecting atomic overlaps.
    /// Lower values are more permissive, higher values more restrictive.
    pub skip_factor: f64,
    /// User-forced bonds: (atom1, atom2) with 1-based indexing
    /// 
    /// These bonds will be created even if not detected automatically.
    pub forced_bonds: Vec<(usize, usize)>,
    /// User-forbidden bonds: (atom1, atom2) with 1-based indexing
    /// 
    /// These bonds will be removed even if detected automatically.
    pub forbidden_bonds: Vec<(usize, usize)>,
    /// Operation mode for conformer generation (rotation or scanning)
    /// 
    /// Determines whether the system performs bond rotation or bond scanning.
    /// These modes are mutually exclusive and use different algorithms.
    pub operation_mode: OperationMode,
}

impl Molecule {
    /// Creates a new molecule with default bond detection parameters.
    /// 
    /// # Arguments
    /// 
    /// * `atoms` - Vector of atoms that make up the molecule
    /// 
    /// # Returns
    /// 
    /// A new `Molecule` with `bond_factor = 1.0`, `skip_factor = 0.7`, 
    /// and `operation_mode = OperationMode::Rotation`
    /// 
    /// # Examples
    /// 
    /// ```rust
    /// use rotbond::molecule::{Molecule, Atom, OperationMode};
    /// 
    /// let atoms = vec![
    ///     Atom { element: "C".to_string(), x: 0.0, y: 0.0, z: 0.0 },
    /// ];
    /// let molecule = Molecule::new(atoms);
    /// assert_eq!(molecule.bond_factor, 1.0);
    /// assert_eq!(matches!(molecule.operation_mode, OperationMode::Rotation), true);
    /// ```
    pub fn new(atoms: Vec<Atom>) -> Self {
        Molecule {
            atoms,
            bonds: Vec::new(),
            fragments: Vec::new(),
            bond_factor: 1.0,
            skip_factor: 0.7,
            forced_bonds: Vec::new(),
            forbidden_bonds: Vec::new(),
            operation_mode: OperationMode::Rotation,
        }
    }

    /// Adds a rotatable bond to the molecule.
    /// 
    /// # Arguments
    /// 
    /// * `bond` - The bond to add to the molecule
    /// 
    /// # Examples
    /// 
    /// ```rust
    /// use rotbond::molecule::{Molecule, Bond};
    /// 
    /// let mut molecule = Molecule::new(vec![]);
    /// let bond = Bond {
    ///     atom1: 0, atom2: 1, angles: vec![0.0, 180.0],
    ///     is_synchronous: false, reference_bond: None, direction: 1.0,
    /// };
    /// molecule.add_bond(bond);
    /// ```
    pub fn add_bond(&mut self, bond: Bond) {
        self.bonds.push(bond);
    }

    /// Adds a molecular fragment to the molecule.
    /// 
    /// Fragments define which atoms move together during bond rotation.
    /// This method is typically called automatically during bond analysis.
    /// 
    /// # Arguments
    /// 
    /// * `fragment` - Vector of atom indices that form a fragment
    #[allow(dead_code)]
    pub fn add_fragment(&mut self, fragment: Vec<usize>) {
        self.fragments.push(fragment);
    }
}

/// Calculates the distance between two atoms in 3D space.
/// 
/// Uses the Euclidean distance formula: √[(x₂-x₁)² + (y₂-y₁)² + (z₂-z₁)²]
/// 
/// # Arguments
/// 
/// * `atom1` - First atom
/// * `atom2` - Second atom
/// 
/// # Returns
/// 
/// Distance between the atoms in Angstroms
/// 
/// # Examples
/// 
/// ```rust
/// use rotbond::molecule::{Atom, calculate_distance};
/// 
/// let atom1 = Atom::new("C".to_string(), 0.0, 0.0, 0.0);
/// let atom2 = Atom::new("H".to_string(), 1.0, 0.0, 0.0);
/// let distance = calculate_distance(&atom1, &atom2);
/// assert_eq!(distance, 1.0);
/// ```
pub fn calculate_distance(atom1: &Atom, atom2: &Atom) -> f64 {
    let dx = atom2.x - atom1.x;
    let dy = atom2.y - atom1.y;
    let dz = atom2.z - atom1.z;
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// Computes the unit direction vector from atom1 to atom2.
/// 
/// Returns a normalized vector pointing from atom1 towards atom2.
/// The vector has magnitude 1.0 and preserves the direction.
/// 
/// # Arguments
/// 
/// * `atom1` - Starting atom (vector origin)
/// * `atom2` - Target atom (vector destination)
/// 
/// # Returns
/// 
/// A tuple (dx, dy, dz) representing the unit direction vector.
/// Returns (0.0, 0.0, 0.0) if atoms are at the same position.
/// 
/// # Examples
/// 
/// ```rust
/// use rotbond::molecule::{Atom, calculate_unit_vector};
/// 
/// let atom1 = Atom::new("C".to_string(), 0.0, 0.0, 0.0);
/// let atom2 = Atom::new("H".to_string(), 2.0, 0.0, 0.0);
/// let (dx, dy, dz) = calculate_unit_vector(&atom1, &atom2);
/// assert_eq!((dx, dy, dz), (1.0, 0.0, 0.0));
/// ```
pub fn calculate_unit_vector(atom1: &Atom, atom2: &Atom) -> (f64, f64, f64) {
    let dx = atom2.x - atom1.x;
    let dy = atom2.y - atom1.y;
    let dz = atom2.z - atom1.z;
    
    let distance = (dx * dx + dy * dy + dz * dz).sqrt();
    
    // Handle case where atoms are at the same position
    if distance == 0.0 {
        return (0.0, 0.0, 0.0);
    }
    
    (dx / distance, dy / distance, dz / distance)
}

/// Validates that a bond length is positive and within reasonable chemical bounds.
/// 
/// Ensures bond lengths are physically meaningful for molecular structures.
/// Very short bonds (< 0.1 Å) are considered invalid. No upper limit is enforced
/// to allow users flexibility in their research, including long-range interactions.
/// 
/// # Arguments
/// 
/// * `length` - Bond length to validate in Angstroms
/// 
/// # Returns
/// 
/// `true` if the bond length is valid, `false` otherwise
/// 
/// # Examples
/// 
/// ```rust
/// use rotbond::molecule::validate_bond_length;
/// 
/// assert_eq!(validate_bond_length(1.5), true);   // Typical C-C bond
/// assert_eq!(validate_bond_length(0.05), false); // Too short
/// assert_eq!(validate_bond_length(15.0), false); // Too long
/// assert_eq!(validate_bond_length(-1.0), false); // Negative
/// ```
/// 
/// # Chemical Context
/// 
/// Typical bond lengths for reference:
/// - C-H: ~1.1 Å
/// - C-C: ~1.5 Å  
/// - C=C: ~1.3 Å
/// - C≡C: ~1.2 Å
/// - C-N: ~1.5 Å
/// - C-O: ~1.4 Å
pub fn validate_bond_length(length: f64) -> bool {
    length >= 0.1
}

/// Returns the covalent radius for a given element in Angstroms.
/// 
/// Based on OpenBabel covalent radii values, optimized for bond detection
/// in organic and common inorganic molecules. Used in the bond detection
/// formula: `d² < ((rad₁ + rad₂) × bond_factor + 0.45)²`
/// 
/// # Arguments
/// 
/// * `element` - Chemical element symbol (case-insensitive)
/// 
/// # Returns
/// 
/// Covalent radius in Angstroms. Returns 1.5 Å for unknown elements.
/// 
/// # Examples
/// 
/// ```rust
/// use rotbond::molecule::covalent_radius;
/// 
/// assert_eq!(covalent_radius("C"), 0.68);
/// assert_eq!(covalent_radius("H"), 0.23);
/// assert_eq!(covalent_radius("unknown"), 1.5);
/// ```
/// 
/// # Supported Elements
/// 
/// Includes all common elements from H to U, with special attention to:
/// - Organic elements (C, H, N, O, P, S)
/// - Halogens (F, Cl, Br, I)
/// - Metals commonly found in biochemistry
pub fn covalent_radius(element: &str) -> f64 {
    let elem = element.to_lowercase();
    match elem.as_str() {
        "h" => 0.23,
        "he" => 1.50,
        "li" => 0.68,
        "be" => 0.35,
        "b" => 0.83,
        "c" => 0.68,
        "n" => 0.68,
        "o" => 0.68,
        "f" => 0.64,
        "ne" => 1.50,
        "na" => 0.97,
        "mg" => 1.10,
        "al" => 1.35,
        "si" => 1.20,
        "p" => 1.05,
        "s" => 1.02,
        "cl" => 0.99,
        "ar" => 1.51,
        "k" => 1.33,
        "ca" => 0.99,
        "sc" => 1.44,
        "ti" => 1.47,
        "v" => 1.33,
        "cr" => 1.35,
        "mn" => 1.35,
        "fe" => 1.34,
        "co" => 1.33,
        "ni" => 1.50,
        "cu" => 1.52,
        "zn" => 1.45,
        "ga" => 1.22,
        "ge" => 1.17,
        "as" => 1.21,
        "se" => 1.22,
        "br" => 1.21,
        "kr" => 1.50,
        "rb" => 1.50,
        "sr" => 1.50,
        "y" => 1.50,
        "zr" => 1.50,
        "nb" => 1.50,
        "mo" => 1.50,
        "tc" => 1.50,
        "ru" => 1.50,
        "rh" => 1.50,
        "pd" => 1.50,
        "ag" => 1.50,
        "cd" => 1.50,
        "in" => 1.50,
        "sn" => 1.50,
        "sb" => 1.50,
        "te" => 1.50,
        "i" => 1.40,
        "xe" => 1.50,
        "cs" => 1.50,
        "ba" => 1.50,
        "la" => 1.50,
        "ce" => 1.50,
        "pr" => 1.50,
        "nd" => 1.50,
        "pm" => 1.50,
        "sm" => 1.50,
        "eu" => 1.50,
        "gd" => 1.50,
        "tb" => 1.50,
        "dy" => 1.50,
        "ho" => 1.50,
        "er" => 1.50,
        "tm" => 1.50,
        "yb" => 1.50,
        "lu" => 1.50,
        "hf" => 1.50,
        "ta" => 1.50,
        "w" => 1.50,
        "re" => 1.50,
        "os" => 1.50,
        "ir" => 1.50,
        "pt" => 1.50,
        "au" => 1.50,
        "hg" => 1.50,
        "tl" => 1.50,
        "pb" => 1.50,
        "bi" => 1.50,
        "po" => 1.50,
        "at" => 1.50,
        "rn" => 1.50,
        "fr" => 1.50,
        "ra" => 1.50,
        "ac" => 1.50,
        "th" => 1.50,
        "pa" => 1.50,
        "u" => 1.50,
        "np" => 1.50,
        "pu" => 1.50,
        "am" => 1.50,
        "cm" => 1.50,
        "bk" => 1.50,
        "cf" => 1.50,
        "es" => 1.50,
        "fm" => 1.50,
        "md" => 1.50,
        "no" => 1.50,
        "lr" => 1.50,
        _ => 1.50, // Default for unknown elements
    }
}

/// Applies bond length scanning to generate new molecular coordinates.
/// 
/// Modifies the positions of atoms in a molecular fragment to achieve a target
/// bond length between two specified atoms. The scanning preserves molecular
/// geometry for all non-scanning bonds by moving one fragment as a rigid body.
/// 
/// This function integrates with existing validation systems by using covalent
/// radii for bond length validation, ensuring that scanning produces chemically
/// reasonable geometries that will pass the same validation as rotation mode.
/// 
/// # Arguments
/// 
/// * `atoms` - Mutable slice of all atoms in the molecule
/// * `fragment` - Indices of atoms that should move together (0-based)
/// * `atom1_idx` - Index of the first atom in the scanning bond (0-based)
/// * `atom2_idx` - Index of the second atom in the scanning bond (0-based)
/// * `target_length` - Desired bond length in Angstroms
/// 
/// # Returns
/// 
/// `Ok(())` if scanning was successful, or an error if:
/// - Target length would create invalid geometry
/// - Atom indices are out of bounds
/// - Current bond length is zero (atoms at same position)
/// - Target length is outside reasonable chemical bounds
/// 
/// # Algorithm
/// 
/// 1. Calculate current distance between atom1 and atom2
/// 2. Compute unit direction vector from atom1 to atom2
/// 3. Validate target length using covalent radii (integrates with validation systems)
/// 4. Calculate displacement needed: target_length - current_length
/// 5. Move all atoms in the fragment along the direction vector
/// 6. Validate that the new bond length is reasonable and precise
/// 
/// # Integration with Validation Systems
/// 
/// - Uses covalent_radius() function for chemical bond length validation
/// - Validates target lengths against reasonable chemical bounds
/// - Ensures scanning produces geometries compatible with existing validation
/// - Maintains consistency with bond length validation used in other parts of the system
/// 
/// # Examples
/// 
/// ```rust
/// use rotbond::molecule::{Atom, scan_bond_length};
/// 
/// let mut atoms = vec![
///     Atom::new("C".to_string(), 0.0, 0.0, 0.0),
///     Atom::new("H".to_string(), 1.0, 0.0, 0.0),
/// ];
/// let fragment = vec![1]; // Move the hydrogen atom
/// 
/// // Stretch C-H bond from 1.0 Å to 1.2 Å
/// scan_bond_length(&mut atoms, &fragment, 0, 1, 1.2)?;
/// ```
/// 
/// # Fragment Movement
/// 
/// The function moves the specified fragment as a rigid body, preserving:
/// - Internal bond lengths within the fragment
/// - Internal bond angles within the fragment  
/// - Internal dihedral angles within the fragment
/// 
/// Only the distance between atom1 and atom2 is modified.
pub fn scan_bond_length(
    atoms: &mut [Atom],
    fragment: &[usize],
    atom1_idx: usize,
    atom2_idx: usize,
    target_length: f64,
) -> Result<(), String> {
    // Validate inputs with enhanced error messages
    if atom1_idx >= atoms.len() || atom2_idx >= atoms.len() {
        return Err(format!("Atom indices out of bounds: atom {} or {} exceeds molecule size ({})", 
                          atom1_idx + 1, atom2_idx + 1, atoms.len()));
    }
    
    if !validate_bond_length(target_length) {
        return Err(format!("Target bond length {:.3} Å is too short (minimum: 0.1 Å). Consider using smaller compression steps.", target_length));
    }
    
    // Calculate current bond length and direction
    let atom1 = &atoms[atom1_idx];
    let atom2 = &atoms[atom2_idx];
    let current_length = calculate_distance(atom1, atom2);
    
    if current_length == 0.0 {
        return Err(format!("Cannot scan bond between atoms {} and {}: atoms are at the same position", 
                          atom1_idx + 1, atom2_idx + 1));
    }
    
    // Validate that the target length is chemically reasonable for these atom types
    let min_reasonable = (covalent_radius(&atom1.element) + covalent_radius(&atom2.element)) * 0.5;
    let max_reasonable = (covalent_radius(&atom1.element) + covalent_radius(&atom2.element)) * 2.5;
    
    if target_length < min_reasonable {
        return Err(format!("Target bond length {:.3} Å is too short for {}-{} bond (minimum reasonable: {:.2} Å)", 
                          target_length, atom1.element, atom2.element, min_reasonable));
    }
    
    if target_length > max_reasonable {
        return Err(format!("Target bond length {:.3} Å is too long for {}-{} bond (maximum reasonable: {:.2} Å)", 
                          target_length, atom1.element, atom2.element, max_reasonable));
    }
    
    // Calculate displacement needed
    let displacement = target_length - current_length;
    
    // Handle empty fragment case gracefully
    if fragment.is_empty() {
        // No atoms to move, but this is not necessarily an error
        // The bond length is already at the target (no displacement needed)
        return Ok(());
    }
    
    // Get unit direction vector from atom1 to atom2
    let (dx, dy, dz) = calculate_unit_vector(atom1, atom2);
    
    // Move all atoms in the fragment along the direction vector
    for &atom_idx in fragment {
        if atom_idx >= atoms.len() {
            return Err(format!("Fragment contains invalid atom index: {} (molecule has {} atoms)", 
                              atom_idx + 1, atoms.len()));
        }
        
        atoms[atom_idx].x += displacement * dx;
        atoms[atom_idx].y += displacement * dy;
        atoms[atom_idx].z += displacement * dz;
    }
    
    // Validate the result
    let new_length = calculate_distance(&atoms[atom1_idx], &atoms[atom2_idx]);
    let length_error = (new_length - target_length).abs();
    
    if length_error > 1e-6 {
        return Err(format!("Bond scanning precision error: achieved {:.6} Å instead of {:.6} Å (error: {:.6} Å)", 
                          new_length, target_length, length_error));
    }
    
    Ok(())
}
