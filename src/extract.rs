use std::str;
use std::collections::HashSet;
extern crate lazy_static;
extern crate needletail;
use needletail::{fastx};


// Local code
//mod codon_table;
//mod complement_table;

pub fn extract_kmers(sequence_files: Vec<&str>, ksize: u8) {

    let mut all_kmers = HashSet::new();


    for sequence_file in sequence_files {
        eprintln!("Counting k-mers in file: {}", sequence_file);

        let mut n_bases = 0;
        fastx::fastx_cli(&sequence_file[..], |_| {}, |seq| {
// seq.id is the name of the record
// seq.seq is the base sequence
// seq.qual is an optional quality score

// keep track of the total number of bases
            n_bases += seq.seq.len();

// keep track of the number of AAAA (or TTTT via canonicalization) in the
// file (normalize makes sure every base is capitalized for comparison)
            for (_, kmer, _) in seq.normalize(false).kmers(ksize, false) {
                let kmer = str::from_utf8(&kmer).unwrap();
                all_kmers.insert(kmer.to_owned());
//            println!("{}", kmer);
            }
        }).expect(&format!("Could not read {}", sequence_file));


        eprintln!("There are {} bases in your file.", n_bases);
    }
    eprintln!("There are {} k-mers in your file.", all_kmers.len());

    for kmer in all_kmers {
        println!("{}", kmer);
    }
}
