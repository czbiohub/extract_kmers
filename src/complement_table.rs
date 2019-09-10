use std::collections::HashMap;

#[inline]
pub fn complement(n: &u8) -> u8 {
    //! Returns the complementary base for a given IUPAC base code.
    //!
    //! Does not work for RNA sequences (maybe we should raise an error or something?)
    match *n as char {
        'a' => 't' as u8, 'A' => 'T' as u8,
        'c' => 'g' as u8, 'C' => 'G' as u8,
        'g' => 'c' as u8, 'G' => 'C' as u8,
        't' => 'a' as u8, 'T' => 'A' as u8,

        // IUPAC codes
        'r' => 'y' as u8, 'y' => 'r' as u8,
        'k' => 'm' as u8, 'm' => 'k' as u8,
        'b' => 'v' as u8, 'v' => 'b' as u8,
        'd' => 'h' as u8, 'h' => 'd' as u8,
        's' => 's' as u8, 'w' => 'w' as u8,
        'R' => 'Y' as u8, 'Y' => 'R' as u8,
        'K' => 'M' as u8, 'M' => 'K' as u8,
        'B' => 'V' as u8, 'V' => 'B' as u8,
        'D' => 'H' as u8, 'H' => 'D' as u8,
        'S' => 'S' as u8, 'W' => 'W' as u8,

        // anything else just pass through
        // 'u' | 'U' => panic!("Does not support complements of U"),
        x => x as u8,
    }
}
