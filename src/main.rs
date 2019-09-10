use std::str;
use std::collections::HashSet;
extern crate lazy_static;
extern crate needletail;
use needletail::{fastx};
extern crate clap;
use clap::{Arg, App, value_t};

// Local code
mod codon_table;

fn main() {
    let matches = App::new("My Super Program")
        .version("1.0")
        .author("Olga Botvinnik <olga.botvinnik@czbiohub.org>")
        .about("Does awesome things")
        .arg(Arg::with_name("files")
            .short("f")
            .long("files")
            .value_name("FILES")
            .help("Input sequence file")
            .required(true)
            .multiple(true)
            .takes_value(true))
        .arg(Arg::with_name("ksize")
            .short("k")
            .long("ksize")
            .value_name("KSIZE")
            .help("k-mer size")
            // Clap currently only handles strings. Have to convert to int later
            .default_value("21")
            .takes_value(true))
        .get_matches();

    let files = matches.values_of("files").unwrap();

    // Convert ksize string argument to integer
    let ksize = value_t!(matches, "ksize", u8).unwrap_or_else(|e| e.exit());
    let mut all_kmers = HashSet::new();


    for file in files {
        eprintln!("Counting k-mers in file: {}", file);

        let mut n_bases = 0;
        fastx::fastx_cli(&file[..], |_| {}, |seq| {
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
        }).expect(&format!("Could not read {}", file));


        eprintln!("There are {} bases in your file.", n_bases);
    }
    eprintln!("There are {} k-mers in your file.", all_kmers.len());

    for kmer in all_kmers {
        println!("{}", kmer);
    }

}
