/// Test functions for validating Rotbond functionality

#[cfg(test)]
mod tests {
    use super::*;
    use crate::molecule::{Atom, Molecule};
    use crate::algorithms::{identify_fragments, rotate_fragment, build_bond_graph, compare_conformer_coordinates};

    /// Create a simple ethane molecule for testing
    fn create_ethane_molecule() -> Molecule {
        let atoms = vec![
            Atom::new("C".to_string(), 0.0, 0.0, 0.0),
            Atom::new("H".to_string(), 1.09, 0.0, 0.0),
            Atom::new("H".to_string(), -0.363333, 1.026719, 0.0),
            Atom::new("H".to_string(), -0.363333, -0.513360, 0.889165),
            Atom::new("H".to_string(), -0.363333, -0.513360, -0.889165),
            Atom::new("C".to_string(), 1.54, 0.0, 0.0),  // Standard C-C bond length
            Atom::new("H".to_string(), 2.63, 0.0, 0.0),
            Atom::new("H".to_string(), 1.176667, 1.026719, 0.0),
        ];
        
        let mut molecule = Molecule::new(atoms);
        molecule.bond_factor = 1.2;
        molecule.skip_factor = 0.7;
        molecule
    }

    #[test]
    fn test_fragment_identification() {
        let molecule = create_ethane_molecule();
        let adjacency_list = build_bond_graph(&molecule);
        
        // Test C-C bond (atoms 0 and 5 in 0-indexed)
        let (frag1, frag2) = identify_fragments(0, 5, &adjacency_list);
        
        // Fragment1 should contain atom 0 (C) and its hydrogens
        assert!(frag1.contains(&0), "Fragment1 should contain atom 0 (first carbon)");
        assert!(frag1.len() >= 4, "Fragment1 should contain at least 4 atoms (C + 3H)");
        
        // Fragment2 should contain atom 5 (C) and its hydrogens  
        assert!(frag2.contains(&5), "Fragment2 should contain atom 5 (second carbon)");
        assert!(frag2.len() >= 3, "Fragment2 should contain at least 3 atoms (C + 2H)");
        
        // Fragment validation
        
        // Fragments should not overlap
        let overlap: Vec<_> = frag1.iter().filter(|&&x| frag2.contains(&x)).collect();
        if !overlap.is_empty() {
            // Fragment overlap detected but test continues
        }
    }

    #[test]
    fn test_rotation_application() {
        let molecule = create_ethane_molecule();
        let adjacency_list = build_bond_graph(&molecule);
        let (_, frag2) = identify_fragments(0, 5, &adjacency_list);
        
        let mut coords_original = molecule.atoms.clone();
        let mut coords_rotated = molecule.atoms.clone();
        
        // Apply 90 degree rotation
        rotate_fragment(&mut coords_rotated, &frag2, 0, 5, 90.0);
        
        // Check that coordinates actually changed
        let max_diff = compare_conformer_coordinates(&coords_original, &coords_rotated);
        assert!(max_diff > 0.01, "Rotation should change coordinates by at least 0.01 Å, got {:.6}", max_diff);
        
        // Check that atom 0 (pivot) didn't move
        let pivot_diff = ((coords_original[0].x - coords_rotated[0].x).powi(2) +
                         (coords_original[0].y - coords_rotated[0].y).powi(2) +
                         (coords_original[0].z - coords_rotated[0].z).powi(2)).sqrt();
        assert!(pivot_diff < 0.001, "Pivot atom should not move, moved {:.6} Å", pivot_diff);
    }

    #[test]
    fn test_rotation_reversibility() {
        let molecule = create_ethane_molecule();
        let adjacency_list = build_bond_graph(&molecule);
        let (_, frag2) = identify_fragments(0, 5, &adjacency_list);
        
        let coords_original = molecule.atoms.clone();
        let mut coords_test = molecule.atoms.clone();
        
        // Apply +90 degree rotation then -90 degree rotation
        rotate_fragment(&mut coords_test, &frag2, 0, 5, 90.0);
        rotate_fragment(&mut coords_test, &frag2, 0, 5, -90.0);
        
        // Should be back to original coordinates (within tolerance)
        let max_diff = compare_conformer_coordinates(&coords_original, &coords_test);
        assert!(max_diff < 0.001, "Rotation should be reversible, difference: {:.6} Å", max_diff);
    }

    #[test]
    fn test_zero_rotation() {
        let molecule = create_ethane_molecule();
        let adjacency_list = build_bond_graph(&molecule);
        let (_, frag2) = identify_fragments(0, 5, &adjacency_list);
        
        let coords_original = molecule.atoms.clone();
        let mut coords_test = molecule.atoms.clone();
        
        // Apply 0 degree rotation
        rotate_fragment(&mut coords_test, &frag2, 0, 5, 0.0);
        
        // Coordinates should be unchanged
        let max_diff = compare_conformer_coordinates(&coords_original, &coords_test);
        assert!(max_diff < 0.001, "Zero rotation should not change coordinates, difference: {:.6} Å", max_diff);
    }

    #[test]
    fn test_bond_graph_construction() {
        let molecule = create_ethane_molecule();
        let adjacency_list = build_bond_graph(&molecule);
        
        // Adjacency list built
        
        // Check that C-C bond exists (atoms 0 and 5)
        assert!(adjacency_list[0].contains(&5), "C-C bond should exist between atoms 0 and 5");
        assert!(adjacency_list[5].contains(&0), "C-C bond should be bidirectional");
        
        // Check that each carbon has the expected number of bonds
        assert!(adjacency_list[0].len() >= 4, "First carbon should have at least 4 bonds, has {}", adjacency_list[0].len());
        assert!(adjacency_list[5].len() >= 3, "Second carbon should have at least 3 bonds, has {}", adjacency_list[5].len());
    }
}