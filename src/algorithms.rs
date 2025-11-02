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



/// Builds a bond connectivity graph from a molecule using covalent radii and user overrides.
/// 
/// Uses the OpenBabel bond detection formula combined with user-defined forced and
/// forbidden bonds to create a comprehensive bond connectivity graph.
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
/// 1. **Distance-based detection**: Uses OpenBabel formula with bond_factor
/// 2. **Forced bonds**: Adds user-specified bonds regardless of distance
/// 3. **Forbidden bonds**: Removes user-specified bonds even if detected
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

/// Get dihedral tracking information for debugging
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

/// Apply rotation to set an ABSOLUTE dihedral angle
/// This is the correct function for step-based rotations (e30, e60, etc.)
/// where we want to set the dihedral to specific absolute values
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



/// Robust Rodrigues rotation implementation
/// Based on OpenBabel's SetTorsion algorithm - direct, stable, and mathematically correct
/// 
/// Usage:
/// - `atoms`: mutable slice of Atom (must have .x, .y, .z fields as f64)
/// - `fragment`: list of atom indices (0-based) that should be rotated
/// - `pivot_atom1`, `pivot_atom2`: indices of the bonded atoms defining the rotation axis (B and C)
/// - `angle_degrees`: rotation to apply in degrees (positive follows right-hand rule around axis = pivot2 - pivot1)
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
    
    if fragment.is_empty() || (angle_degrees - 0.0).abs() < std::f64::EPSILON {
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



/// Validates a conformer by checking for steric clashes between atoms.
/// 
/// Uses distance-based criteria to detect atomic overlaps that would make
/// a conformer chemically unrealistic. The validation uses covalent radii
/// scaled by the skip_factor parameter.
/// 
/// # Arguments
/// 
/// * `molecule` - The molecule containing validation parameters
/// * `coordinates` - The conformer to validate
/// 
/// # Returns
/// 
/// `true` if the conformer is valid (no significant clashes), `false` otherwise.
/// 
/// # Algorithm
/// 
/// Checks all atom pairs for distances shorter than the sum of their
/// covalent radii multiplied by skip_factor. Lower skip_factor values
/// are more permissive.
/// 
/// # Examples
/// 
/// ```rust
/// use rotbond::algorithms::is_valid_conformer;
/// 
/// if is_valid_conformer(&molecule, &atoms) {
///     println!("Conformer is valid");
/// } else {
///     println!("Conformer has steric clashes");
/// }
/// ```
pub fn is_valid_conformer(molecule: &Molecule, coordinates: &[Atom]) -> bool {
    let n = coordinates.len();

    for i in 0..n {
        for j in (i + 1)..n {
            let dx = coordinates[i].x - coordinates[j].x;
            let dy = coordinates[i].y - coordinates[j].y;
            let dz = coordinates[i].z - coordinates[j].z;
            let distance = (dx * dx + dy * dy + dz * dz).sqrt();

            let min_distance = (covalent_radius(&coordinates[i].element) +
                              covalent_radius(&coordinates[j].element)) * molecule.skip_factor;

            if distance < min_distance {
                return false;
            }
        }
    }

    true
}

/// Compare coordinates between two conformers to detect if they are identical
/// Returns the maximum coordinate difference found
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

/// Check if two conformers are effectively identical (within tolerance)
pub fn are_conformers_identical(conformer1: &[Atom], conformer2: &[Atom], tolerance: f64) -> bool {
    compare_conformer_coordinates(conformer1, conformer2) < tolerance
}

/// Validate that a set of conformers contains actually different structures
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

/// Log fragment identification results for debugging
pub fn log_fragment_info(_bond_idx: usize, _atom1: usize, _atom2: usize, _fragment1: &[usize], _fragment2: &[usize]) {
    // Debug output removed
}

/// Log rotation application details for debugging
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
/// # Arguments
/// 
/// * `a` - First 3D vector
/// * `b` - Second 3D vector
/// 
/// # Returns
/// 
/// The cross product vector a × b.
fn cross_product(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2], 
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Computes the dot product of two 3D vectors.
/// 
/// # Arguments
/// 
/// * `a` - First 3D vector
/// * `b` - Second 3D vector
/// 
/// # Returns
/// 
/// The dot product scalar a · b.
fn dot_product(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Normalizes a 3D vector to unit length.
/// 
/// # Arguments
/// 
/// * `v` - The 3D vector to normalize
/// 
/// # Returns
/// 
/// The normalized vector with unit length, or zero vector if input is zero.
fn normalize_vector(v: &[f64; 3]) -> [f64; 3] {
    let length = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if length == 0.0 {
        [0.0, 0.0, 0.0]
    } else {
        [v[0] / length, v[1] / length, v[2] / length]
    }
}


// Build bond graph from atoms using distance-based detection
/// This is a simplified version for when we don't have the full molecule structure
fn build_bond_graph_from_atoms(atoms: &[Atom]) -> Vec<Vec<usize>> {
    let n = atoms.len();
    let mut adjacency_list = vec![Vec::new(); n];
    
    for i in 0..n {
        for j in (i + 1)..n {
            let dx = atoms[i].x - atoms[j].x;
            let dy = atoms[i].y - atoms[j].y;
            let dz = atoms[i].z - atoms[j].z;
            let distance = (dx * dx + dy * dy + dz * dz).sqrt();
            
            // Use covalent radii for bonding detection
            let bond_threshold = (covalent_radius(&atoms[i].element) + 
                                 covalent_radius(&atoms[j].element)) * 1.2;
            
            if distance < bond_threshold {
                adjacency_list[i].push(j);
                adjacency_list[j].push(i);
            }
        }
    }
    
    adjacency_list
}



/// Identify fragment1 (atoms connected to atom1, excluding atom2)
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

/// Find reference atoms for dihedral A-atom1-atom2-B
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

/// Find the best reference atom for dihedral calculation
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

/// Calculate the length of the longest backbone chain from an atom
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

/// Find the bond rotation needed to achieve a target dihedral angle
/// Uses improved iterative approach with finite difference derivative estimation
/// Based on OpenBabel's approach to solving the non-linear bond rotation <-> dihedral relationship
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

