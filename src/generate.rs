//! # Conformer Generation Engine
//!
//! This module contains the core conformer generation algorithms that systematically
//! explore conformational space by rotating around rotatable bonds.
//!
//! ## Generation Process
//!
//! 1. **Bond Analysis**: Identify rotatable bonds and molecular fragments
//! 2. **Combinatorial Generation**: Create all possible angle combinations
//! 3. **Rotation Application**: Apply dihedral rotations using Rodrigues formula
//! 4. **Validation**: Filter conformers using steric clash detection
//! 5. **Output**: Write valid conformers to XYZ files
//!
//! ## Key Features
//!
//! - **Systematic exploration**: Generates all theoretically possible conformers
//! - **Synchronous rotations**: Supports coordinated bond movements
//! - **Steric filtering**: Removes conformers with atomic overlaps
//! - **Progress tracking**: Reports generation statistics and success rates
//! - **Multiple outputs**: Individual files and trajectory format
//!
//! ## Performance
//!
//! The generation scales as O(n^k) where n is the average number of angles
//! per bond and k is the number of rotatable bonds. For molecules with many
//! rotatable bonds, consider using larger step sizes (e.g., e120 instead of e30).

use crate::molecule::{Molecule, Bond, OperationMode, RotationSpec, scan_bond_length};
use crate::algorithms::{build_bond_graph, identify_fragments, rotate_fragment, is_valid_conformer, 
                        log_fragment_info, log_rotation_application, validate_conformer_diversity};
use crate::io::{write_trajectory_file, write_individual_xyz_file, format_conformer_filename};

/// Generates all possible conformers based on rotation or scanning specifications.
/// 
/// This is the main conformer generation function that systematically explores
/// conformational space by applying all combinations of specified rotations or
/// bond length variations. Supports both independent and synchronous bond rotations,
/// as well as single and multi-dimensional bond scanning.
/// 
/// # Arguments
/// 
/// * `molecule` - The molecular structure with bond detection parameters and operation mode
/// * `angle_sets` - Pre-computed angle sets for rotatable bonds (empty for scanning bonds)
/// * `rotation_specs` - Original rotation specifications for extracting scanning parameters
/// * `base_name` - Base name for output files (e.g., "ethane" → "ethane_1.xyz")
/// * `max_conformers` - Optional limit on number of conformers to generate
/// 
/// # Returns
/// 
/// A tuple containing:
/// - Total number of theoretical combinations attempted
/// - Number of valid conformers successfully generated
/// 
/// # Errors
/// 
/// Returns an error if:
/// - Bond graph construction fails
/// - Fragment identification fails
/// - File writing operations fail
/// - Scanning parameters are invalid
/// 
/// # Process
/// 
/// ## Rotation Mode
/// 1. **Preparation**: Build bond graph and identify molecular fragments
/// 2. **Bond Classification**: Separate independent from synchronous bonds
/// 3. **Combination Generation**: Create all possible angle combinations
/// 4. **Rotation Application**: Apply rotations using dihedral algorithms
/// 5. **Validation**: Check for steric clashes and geometric validity
/// 6. **Output**: Write valid conformers to individual XYZ files
/// 7. **Trajectory**: Combine all conformers into a single trajectory file
/// 
/// ## Scanning Mode
/// 1. **Preparation**: Build bond graph and identify molecular fragments
/// 2. **Parameter Extraction**: Extract scanning parameters from rotation specs
/// 3. **Combination Generation**: Create all possible bond length combinations
/// 4. **Scanning Application**: Apply bond length changes using fragment movement
/// 5. **Validation**: Check for steric clashes and geometric validity
/// 6. **Output**: Write valid conformers to individual XYZ files
/// 7. **Trajectory**: Combine all conformers into a single trajectory file
/// 
/// # Examples
/// 
/// ```rust
/// use rotbond::generate::generate_conformers;
/// 
/// // Rotation mode
/// let (total, valid) = generate_conformers(&mut molecule, &angle_sets, &rotation_specs, "ethane", None)?;
/// println!("Generated {}/{} conformers ({:.1}% success)", 
///          valid, total, 100.0 * valid as f64 / total as f64);
/// 
/// // Scanning mode (angle_sets will be empty, scanning params extracted from rotation_specs)
/// let (total, valid) = generate_conformers(&mut scanning_molecule, &[], &scanning_specs, "ethane_scan", Some(100))?;
/// ```
/// 
/// # Performance Notes
/// 
/// - **Rotation**: Complexity O(n^k) where n = avg angles per bond, k = number of bonds
/// - **Scanning**: Complexity O(s₁ × s₂ × ... × sₖ) where sᵢ = steps for bond i
/// - Memory usage scales with the number of valid conformers
/// - For large molecules, consider using larger step sizes or fewer scanning steps
pub fn generate_conformers(
    molecule: &mut Molecule,
    angle_sets: &[Vec<f64>],
    rotation_specs: &[RotationSpec],
    base_name: &str,
    max_conformers: Option<usize>,
) -> Result<(usize, usize), Box<dyn std::error::Error>> {
    // Detect operation mode and delegate to appropriate generation function
    match molecule.operation_mode {
        OperationMode::Rotation => {
            generate_rotation_conformers(molecule, angle_sets, base_name, max_conformers)
        }
        OperationMode::Scanning => {
            generate_scanning_conformers(molecule, rotation_specs, base_name, max_conformers)
        }
    }
}

/// Generates conformers using bond rotation mode.
/// 
/// This function handles the traditional dihedral rotation approach to conformer
/// generation. It builds bond graphs, identifies fragments, and systematically
/// explores all combinations of rotation angles.
/// 
/// # Arguments
/// 
/// * `molecule` - The molecular structure with bond detection parameters
/// * `angle_sets` - Pre-computed angle sets for each rotatable bond
/// * `base_name` - Base name for output files
/// * `max_conformers` - Optional limit on number of conformers to generate
/// 
/// # Returns
/// 
/// A tuple containing (total_combinations_attempted, valid_conformers_generated).
/// 
/// # Process
/// 
/// 1. Build bond graph and identify molecular fragments
/// 2. Separate independent from synchronous bonds
/// 3. Generate all combinations of rotation angles
/// 4. Apply rotations using dihedral algorithms
/// 5. Validate conformers and write output files
fn generate_rotation_conformers(
    molecule: &mut Molecule,
    angle_sets: &[Vec<f64>],
    base_name: &str,
    max_conformers: Option<usize>,
) -> Result<(usize, usize), Box<dyn std::error::Error>> {
    // Build bond graph
    let adjacency_list = build_bond_graph(molecule);

    // Identify fragments for each bond
    molecule.fragments.clear();
    for (i, bond) in molecule.bonds.iter().enumerate() {
        let (frag1, frag2) = identify_fragments(bond.atom1, bond.atom2, &adjacency_list);
        
        // Log fragment information for debugging
        log_fragment_info(i, bond.atom1, bond.atom2, &frag1, &frag2);
        
        // Use fragment2 (atoms attached to atom2) for rotation
        molecule.fragments.push(frag2);
    }

    // Separate independent bonds from synchronous bonds
    let mut independent_indices = Vec::new();
    let mut synchronous_bonds = Vec::new();

    for (i, bond) in molecule.bonds.iter().enumerate() {
        if bond.is_synchronous {
            synchronous_bonds.push(i);
        } else {
            independent_indices.push(i);
        }
    }

    // Estimate total conformers
    let mut total_estimated = 1;
    for &idx in &independent_indices {
        total_estimated *= angle_sets[idx].len().max(1);
    }

    // Generate cartesian product of independent bonds
    let (total_combinations, valid_conformers) = generate_conformer_combinations(
        molecule,
        &angle_sets,
        &independent_indices,
        &synchronous_bonds,
        base_name,
        total_estimated,
        max_conformers,
    )?;

    Ok((total_combinations, valid_conformers))
}

/// Generates conformers using bond length scanning mode.
/// 
/// This function handles the bond scanning approach to conformer generation.
/// It extracts scanning parameters from rotation specifications, builds bond
/// graphs, identifies fragments, and systematically explores all combinations
/// of bond length variations.
/// 
/// Integrates with existing validation systems by:
/// - Using bond_factor and skip_factor from molecule configuration
/// - Respecting forced_bonds and forbidden_bonds specifications
/// - Applying maxgen limits for conformer generation
/// - Using same steric clash validation as rotation mode
/// 
/// # Arguments
/// 
/// * `molecule` - The molecular structure with scanning specifications in bonds
/// * `rotation_specs` - Rotation specifications containing scanning parameters
/// * `base_name` - Base name for output files
/// * `max_conformers` - Optional limit on number of conformers to generate (from maxgen)
/// 
/// # Returns
/// 
/// A tuple containing (total_combinations_attempted, valid_conformers_generated).
/// 
/// # Process
/// 
/// 1. Extract scanning parameters from bond rotation specifications
/// 2. Build bond graph using molecule's bond_factor and forced/forbidden bonds
/// 3. Identify molecular fragments for each scanning bond
/// 4. Generate all combinations of bond length variations
/// 5. Apply bond length changes using fragment movement
/// 6. Validate conformers using molecule's skip_factor and existing validation
/// 7. Write output files with proper conformer limits
/// 
/// # Errors
/// 
/// Returns an error if:
/// - No scanning specifications are found
/// - Bond graph construction fails (respecting bond_factor)
/// - Fragment identification fails
/// - Scanning parameters are invalid
/// - Forced/forbidden bond constraints are violated
fn generate_scanning_conformers(
    molecule: &mut Molecule,
    rotation_specs: &[RotationSpec],
    base_name: &str,
    max_conformers: Option<usize>,
) -> Result<(usize, usize), Box<dyn std::error::Error>> {
    // Build bond graph using molecule's bond_factor and forced/forbidden bonds
    // This ensures scanning respects the same bond detection parameters as rotation
    let adjacency_list = build_bond_graph(molecule);

    // Extract scanning parameters from rotation specifications and identify fragments
    molecule.fragments.clear();
    let mut scanning_bonds = Vec::new();
    let mut scanning_params = Vec::new();
    
    for (_i, spec) in rotation_specs.iter().enumerate() {
        if let RotationSpec::Scanning { atom1, atom2, steps, step_size } = spec {
            // This is a scanning bond - create the bond with correct atom indices
            let bond = Bond {
                atom1: *atom1 - 1, // Convert from 1-based to 0-based indexing
                atom2: *atom2 - 1,
                angles: Vec::new(),
                is_synchronous: false,
                reference_bond: None,
                direction: 1.0,
            };
            
            // Identify fragments for this bond
            let (frag1, frag2) = identify_fragments(bond.atom1, bond.atom2, &adjacency_list);
            
            // Log fragment information for debugging
            log_fragment_info(molecule.bonds.len(), bond.atom1, bond.atom2, &frag1, &frag2);
            
            // Add the bond to the molecule
            molecule.bonds.push(bond);
            
            // Use fragment2 (atoms attached to atom2) for scanning
            molecule.fragments.push(frag2);
            
            scanning_bonds.push(molecule.bonds.len() - 1); // Use the actual bond index
            scanning_params.push((*steps, *step_size));
        } else {
            // Non-scanning bond, add empty fragment for consistency
            molecule.fragments.push(Vec::new());
        }
    }
    
    if scanning_bonds.is_empty() {
        return Err("No scanning bonds found for scanning mode".into());
    }

    // Calculate total combinations for multi-dimensional scanning
    let mut total_estimated = 1;
    for &(steps, _) in &scanning_params {
        total_estimated *= steps;
    }

    // Apply same warning system as rotation for large conformer counts
    // This integrates with the maxgen parameter from conformer configuration
    // Inherits maxgen and autoconfirm parameters for scanning conformer limits
    let effective_limit = max_conformers.unwrap_or(500);
    
    if total_estimated > effective_limit {
        println!("WARNING: Large scanning job detected - {} total combinations", total_estimated);
        println!("This exceeds the current limit of {} conformers.", effective_limit);
        println!("This may take significant time and memory.");
        println!("Consider reducing the number of steps or bonds to scan.");
        println!("Or increase maxgen in your .rp file to allow more conformers.");
        
        // Note: autoconfirm parameter would be handled at the main application level
        // The max_conformers parameter passed here already reflects the maxgen setting
    } else if total_estimated > 500 {
        println!("INFO: Large scanning job detected - {} total combinations", total_estimated);
        println!("This may take significant time and memory.");
        println!("Consider reducing the number of steps or bonds to scan.");
    }

    // Generate all combinations of scanning parameters using recursive pattern similar to rotation
    let (total_combinations, valid_conformers) = generate_scanning_combinations(
        molecule,
        &scanning_bonds,
        &scanning_params,
        base_name,
        total_estimated,
        max_conformers,
    )?;

    Ok((total_combinations, valid_conformers))


}

/// Generates all combinations of conformers using cartesian product approach.
/// 
/// This function handles the combinatorial explosion of conformer generation
/// by systematically exploring all possible combinations of rotation angles
/// for independent bonds, while properly handling synchronous bond relationships.
/// 
/// # Arguments
/// 
/// * `molecule` - The molecular structure with atoms and bonds
/// * `angle_sets` - Pre-computed rotation angles for each bond
/// * `independent_indices` - Indices of bonds that rotate independently
/// * `synchronous_bond_indices` - Indices of bonds that follow other bonds
/// * `base_name` - Base name for output files
/// * `total_estimated` - Estimated total number of conformers for progress tracking
/// 
/// # Returns
/// 
/// A tuple containing (total_combinations_attempted, valid_conformers_generated).
/// 
/// # Algorithm
/// 
/// Uses recursive generation to explore the full combinatorial space while
/// maintaining memory efficiency. Each combination is validated immediately
/// and discarded if invalid, preventing memory overflow for large search spaces.
fn generate_conformer_combinations(
    molecule: &Molecule,
    angle_sets: &[Vec<f64>],
    independent_indices: &[usize],
    synchronous_bond_indices: &[usize],
    base_name: &str,
    total_estimated: usize,
    max_conformers: Option<usize>,
) -> Result<(usize, usize), Box<dyn std::error::Error>> {
    let num_independent = independent_indices.len();

    if num_independent == 0 {
        return Ok((0, 0));
    }

    // For each independent bond index combination
    let mut valid_conformers = Vec::new();

    // Generate all combinations using recursion
    let mut current_combination = vec![0.0; num_independent];
    let mut current_bond_indices = vec![0; num_independent];

    let total_combinations = calculate_total_combinations(angle_sets, independent_indices);

    // Generating combinations...

    let mut processed = 0;

    generate_recursive(
        0,
        molecule,
        angle_sets,
        independent_indices,
        synchronous_bond_indices,
        &mut current_combination,
        &mut current_bond_indices,
        &mut processed,
        total_combinations,
        &mut valid_conformers,
        base_name,
        total_estimated,
        max_conformers,
    )?;

    // Generated conformers

    // Validate conformer diversity
    if valid_conformers.len() > 1 {
        let (_identical_pairs, _total_pairs) = validate_conformer_diversity(&valid_conformers, 0.001);
        // Conformer diversity validation completed
    }

    // Write all conformers to files
    write_all_conformers(base_name, &valid_conformers, total_estimated)?;

    Ok((total_combinations, valid_conformers.len()))
}

/// Recursively generates all angle combinations for conformer exploration.
/// 
/// This is the core recursive function that explores the full combinatorial
/// space of rotation angles. It applies rotations depth-first, validates
/// each complete conformer, and collects valid results.
/// 
/// # Arguments
/// 
/// * `depth` - Current recursion depth (bond index being processed)
/// * `molecule` - The molecular structure
/// * `angle_sets` - Rotation angles for each bond
/// * `independent_indices` - Indices of independently rotating bonds
/// * `synchronous_bond_indices` - Indices of synchronously rotating bonds
/// * `current_combination` - Current angle combination being built
/// * `current_bond_indices` - Current bond indices for tracking
/// * `processed` - Counter for processed combinations
/// * `total_combinations` - Total expected combinations
/// * `valid_conformers` - Collection of valid conformers found
/// * `base_name` - Base name for output files
/// * `total_estimated` - Total estimated conformers for progress
/// 
/// # Algorithm
/// 
/// 1. **Base case**: When depth equals number of bonds, apply all rotations
/// 2. **Recursive case**: Try each angle for current bond, recurse deeper
/// 3. **Validation**: Check each complete conformer for steric clashes
/// 4. **Collection**: Store valid conformers for output
fn generate_recursive(
    depth: usize,
    molecule: &Molecule,
    angle_sets: &[Vec<f64>],
    independent_indices: &[usize],
    synchronous_bond_indices: &[usize],
    current_combination: &mut [f64],
    current_bond_indices: &mut [usize],
    processed: &mut usize,
    total_combinations: usize,
    valid_conformers: &mut Vec<Vec< crate::molecule::Atom>>,
    base_name: &str,
    total_estimated: usize,
    max_conformers: Option<usize>,
) -> Result<(), Box<dyn std::error::Error>> {
    if depth == independent_indices.len() {
        *processed += 1;

        // Progress tracking removed for cleaner output

        // Clone original coordinates
        let mut coordinates = molecule.atoms.clone();

        // Apply rotations for independent bonds
        for i in 0..independent_indices.len() {
            let bond_idx = independent_indices[i];
            let angle = current_combination[i];
            let bond = &molecule.bonds[bond_idx];
            let fragment = &molecule.fragments[bond_idx];

            if angle != 0.0 {
                // Store coordinates before rotation for debugging
                let coords_before = coordinates.clone();
                
                rotate_fragment(&mut coordinates, fragment, bond.atom1, bond.atom2, angle);
                
                // Log rotation application for debugging
                log_rotation_application(bond_idx, angle, fragment, &coords_before, &coordinates);
            }
        }

        // Apply synchronous bond rotations
        for &sync_idx in synchronous_bond_indices {
            let sync_bond = &molecule.bonds[sync_idx];

            if let Some(ref reference) = sync_bond.reference_bond {
                let ref_idx = *reference - 1; // Convert to 0-indexed

                // Find the reference bond in independent_indices
                let mut found = false;
                for (i, &indep_idx) in independent_indices.iter().enumerate() {
                    if indep_idx == ref_idx {
                        let angle = current_combination[i] * sync_bond.direction;
                        let fragment = &molecule.fragments[sync_idx];

                        if angle != 0.0 {
                            rotate_fragment(&mut coordinates, fragment, sync_bond.atom1, sync_bond.atom2, angle);
                        }
                        found = true;
                        break;
                    }
                }

                // If reference bond is not independent, it's another synchronous bond
                if !found {
                    for &other_sync_idx in synchronous_bond_indices {
                        if other_sync_idx == ref_idx {
                            // This is a chain of synchronous bonds
                            // For simplicity, skip (should not happen in well-formed input)
                            break;
                        }
                    }
                }
            }
        }

        // Validate conformer
        if is_valid_conformer(molecule, &coordinates) {
            // Check conformer limit before adding
            if let Some(max_limit) = max_conformers {
                if valid_conformers.len() >= max_limit {
                    // Limit reached, stop generation
                    return Ok(());
                }
            }
            valid_conformers.push(coordinates);
        }

        return Ok(());
    }

    let bond_idx = independent_indices[depth];
    let angles = &angle_sets[bond_idx];

    if angles.is_empty() {
        // No angles for this bond (e.g., it's synchronous)
        current_combination[depth] = 0.0;
        current_bond_indices[depth] = 0;

        generate_recursive(
            depth + 1,
            molecule,
            angle_sets,
            independent_indices,
            synchronous_bond_indices,
            current_combination,
            current_bond_indices,
            processed,
            total_combinations,
            valid_conformers,
            base_name,
            total_estimated,
            max_conformers,
        )?;
    } else {
        for (angle_idx, &angle) in angles.iter().enumerate() {
            current_combination[depth] = angle;
            current_bond_indices[depth] = angle_idx;

            generate_recursive(
                depth + 1,
                molecule,
                angle_sets,
                independent_indices,
                synchronous_bond_indices,
                current_combination,
                current_bond_indices,
                processed,
                total_combinations,
                valid_conformers,
                base_name,
                total_estimated,
                max_conformers,
            )?;
        }
    }

    Ok(())
}

/// Calculates the total number of conformer combinations.
/// 
/// Computes the cartesian product size for all independent rotatable bonds
/// to determine the total search space size. This is used for progress
/// tracking and memory estimation.
/// 
/// # Arguments
/// 
/// * `angle_sets` - Rotation angles for each bond
/// * `independent_indices` - Indices of independently rotating bonds
/// 
/// # Returns
/// 
/// Total number of theoretical conformer combinations.
/// 
/// # Examples
/// 
/// ```rust
/// // Two bonds with 6 and 4 angles respectively = 24 combinations
/// let total = calculate_total_combinations(&angle_sets, &[0, 1]);
/// ```
fn calculate_total_combinations(angle_sets: &[Vec<f64>], independent_indices: &[usize]) -> usize {
    let mut total = 1;
    for &idx in independent_indices {
        let len = angle_sets[idx].len().max(1);
        total *= len;
    }
    total
}

/// Writes all valid conformers to individual files and a trajectory file.
/// 
/// Creates properly formatted XYZ files for each conformer and combines
/// them into a single trajectory file for visualization and analysis.
/// Uses smart filename padding for consistent directory ordering.
/// 
/// # Arguments
/// 
/// * `base_name` - Base name for output files
/// * `conformers` - Vector of valid conformer structures
/// * `_total_estimated` - Total estimated conformers (unused, for future progress tracking)
/// 
/// # Returns
/// 
/// `Ok(())` on success, or an error if file operations fail.
/// 
/// # Output Files
/// 
/// - Individual files: `base_name_001.xyz`, `base_name_002.xyz`, etc.
/// - Trajectory file: `base_name_traj.xyz`
/// 
/// # Examples
/// 
/// ```rust
/// write_all_conformers("ethane", &conformers, 24)?;
/// // Creates: ethane_01.xyz, ethane_02.xyz, ..., ethane_traj.xyz
/// ```
fn write_all_conformers(
    base_name: &str,
    conformers: &[Vec< crate::molecule::Atom>],
    _total_estimated: usize,
) -> Result<(), Box<dyn std::error::Error>> {
    if conformers.is_empty() {
        return Ok(());
    }

    // Write trajectory file
    let traj_filename = format!("{}_traj.xyz", base_name);
    write_trajectory_file(&traj_filename, base_name, conformers, conformers.len())?;
    // Trajectory file created

    // Write individual files
    for (idx, conformer) in conformers.iter().enumerate() {
        let conformer_num = idx + 1;
        let filename = format_conformer_filename(base_name, conformer_num, conformers.len());
        write_individual_xyz_file(&filename, conformer)?;
    }

    // Individual conformer files created

    Ok(())
}

/// Generates all combinations of scanning conformers using cartesian product approach.
/// 
/// This function handles the combinatorial explosion of scanning conformer generation
/// by systematically exploring all possible combinations of bond length variations
/// for multiple scanning bonds.
/// 
/// # Arguments
/// 
/// * `molecule` - The molecular structure with atoms and bonds
/// * `scanning_bonds` - Indices of bonds that will be scanned
/// * `scanning_params` - Scanning parameters (steps, step_size) for each bond
/// * `base_name` - Base name for output files
/// * `total_estimated` - Estimated total number of conformers for progress tracking
/// * `max_conformers` - Optional limit on number of conformers to generate
/// 
/// # Returns
/// 
/// A tuple containing (total_combinations_attempted, valid_conformers_generated).
/// 
/// # Algorithm
/// 
/// Uses recursive generation to explore the full combinatorial space while
/// maintaining memory efficiency. Each combination is validated immediately
/// and discarded if invalid, preventing memory overflow for large search spaces.
fn generate_scanning_combinations(
    molecule: &Molecule,
    scanning_bonds: &[usize],
    scanning_params: &[(usize, f64)], // (steps, step_size)
    base_name: &str,
    total_estimated: usize,
    max_conformers: Option<usize>,
) -> Result<(usize, usize), Box<dyn std::error::Error>> {
    if scanning_bonds.is_empty() {
        return Ok((0, 0));
    }

    let mut valid_conformers = Vec::new();
    let total_combinations = calculate_scanning_combinations(scanning_params);

    // Generate all combinations using recursion
    let mut current_steps = vec![0; scanning_bonds.len()];
    let mut processed = 0;

    generate_scanning_recursive(
        0,
        molecule,
        scanning_bonds,
        scanning_params,
        &mut current_steps,
        &mut processed,
        total_combinations,
        &mut valid_conformers,
        base_name,
        total_estimated,
        max_conformers,
    )?;

    // Validate conformer diversity
    if valid_conformers.len() > 1 {
        let (_identical_pairs, _total_pairs) = validate_conformer_diversity(&valid_conformers, 0.001);
    }

    // Write all conformers to files
    write_all_conformers(base_name, &valid_conformers, total_estimated)?;

    Ok((total_combinations, valid_conformers.len()))
}

/// Recursively generates all scanning combinations for conformer exploration.
/// 
/// This is the core recursive function that explores the full combinatorial
/// space of bond length variations. It applies scanning depth-first, validates
/// each complete conformer, and collects valid results.
/// 
/// # Arguments
/// 
/// * `depth` - Current recursion depth (bond index being processed)
/// * `molecule` - The molecular structure
/// * `scanning_bonds` - Indices of bonds being scanned
/// * `scanning_params` - Scanning parameters for each bond
/// * `current_steps` - Current step combination being built
/// * `processed` - Counter for processed combinations
/// * `total_combinations` - Total expected combinations
/// * `valid_conformers` - Collection of valid conformers found
/// * `base_name` - Base name for output files
/// * `total_estimated` - Total estimated conformers for progress
/// * `max_conformers` - Optional limit on conformers to generate
/// 
/// # Algorithm
/// 
/// 1. **Base case**: When depth equals number of bonds, apply all scannings
/// 2. **Recursive case**: Try each step for current bond, recurse deeper
/// 3. **Validation**: Check each complete conformer for steric clashes
/// 4. **Collection**: Store valid conformers for output
fn generate_scanning_recursive(
    depth: usize,
    molecule: &Molecule,
    scanning_bonds: &[usize],
    scanning_params: &[(usize, f64)],
    current_steps: &mut [usize],
    processed: &mut usize,
    total_combinations: usize,
    valid_conformers: &mut Vec<Vec<crate::molecule::Atom>>,
    base_name: &str,
    total_estimated: usize,
    max_conformers: Option<usize>,
) -> Result<(), Box<dyn std::error::Error>> {
    if depth == scanning_bonds.len() {
        *processed += 1;

        // Clone original coordinates
        let mut coordinates = molecule.atoms.clone();

        // Apply scanning for each bond
        for i in 0..scanning_bonds.len() {
            let bond_idx = scanning_bonds[i];
            let bond = &molecule.bonds[bond_idx];
            let fragment = &molecule.fragments[i]; // Use scanning-specific fragment index
            let (_steps, step_size) = scanning_params[i];
            let current_step = current_steps[i];

            // Calculate target bond length for this step
            let original_length = crate::molecule::calculate_distance(
                &molecule.atoms[bond.atom1], 
                &molecule.atoms[bond.atom2]
            );
            let target_length = original_length + (current_step as f64 * step_size);

            // Apply bond length scanning with enhanced error handling
            if let Err(_e) = scan_bond_length(
                &mut coordinates,
                fragment,
                bond.atom1,
                bond.atom2,
                target_length,
            ) {
                // Log the scanning failure for debugging (optional)
                // eprintln!("Scanning failed for bond {} at step {}: {}", bond_idx + 1, current_step, e);
                
                // Skip this conformer if scanning fails - this is expected behavior
                // for some parameter combinations that create invalid geometry
                return Ok(());
            }
        }

        // Validate conformer using integrated validation systems
        // This uses molecule's skip_factor for steric clash detection
        // and applies scanning-specific validation including forced/forbidden bonds
        if is_valid_conformer(molecule, &coordinates) {
            // Check conformer limit before adding (integrates with maxgen parameter)
            if let Some(max_limit) = max_conformers {
                if valid_conformers.len() >= max_limit {
                    // Limit reached, stop generation (respects maxgen setting)
                    return Ok(());
                }
            }
            valid_conformers.push(coordinates);
        }

        return Ok(());
    }

    let (steps, _step_size) = scanning_params[depth];

    for step in 0..steps {
        current_steps[depth] = step;

        generate_scanning_recursive(
            depth + 1,
            molecule,
            scanning_bonds,
            scanning_params,
            current_steps,
            processed,
            total_combinations,
            valid_conformers,
            base_name,
            total_estimated,
            max_conformers,
        )?;
    }

    Ok(())
}

/// Calculates the total number of scanning combinations.
/// 
/// Computes the cartesian product size for all scanning bonds
/// to determine the total search space size.
/// 
/// # Arguments
/// 
/// * `scanning_params` - Scanning parameters (steps, step_size) for each bond
/// 
/// # Returns
/// 
/// Total number of theoretical scanning combinations.
fn calculate_scanning_combinations(scanning_params: &[(usize, f64)]) -> usize {
    let mut total = 1;
    for &(steps, _) in scanning_params {
        total *= steps;
    }
    total
}
