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

/// Defines the type of rotation specification for a rotatable bond.
/// 
/// Supports three types of rotation specifications:
/// - Step-based rotations (e.g., every 30°, 60°, 120°)
/// - Explicit angle lists
/// - Synchronous rotations that follow another bond
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
/// use rotbond::molecule::{Molecule, Atom};
/// 
/// let atoms = vec![
///     Atom { element: "C".to_string(), x: 0.0, y: 0.0, z: 0.0 },
///     Atom { element: "H".to_string(), x: 1.0, y: 0.0, z: 0.0 },
/// ];
/// let molecule = Molecule::new(atoms);
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
    /// A new `Molecule` with `bond_factor = 1.0` and `skip_factor = 0.7`
    /// 
    /// # Examples
    /// 
    /// ```rust
    /// use rotbond::molecule::{Molecule, Atom};
    /// 
    /// let atoms = vec![
    ///     Atom { element: "C".to_string(), x: 0.0, y: 0.0, z: 0.0 },
    /// ];
    /// let molecule = Molecule::new(atoms);
    /// assert_eq!(molecule.bond_factor, 1.0);
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
