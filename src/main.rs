use std::path::Path;
extern crate lazy_static;
extern crate needletail;
use exitfailure::ExitFailure;
extern crate clap;
use clap::{App, load_yaml, value_t};

use std::path::Path;

mod extract;
mod classify_reads;

fn main() -> Result<(), ExitFailure> {
    // let kmer_path = Path::new("../test-data/primers_R1.fasta");
    // classify_reads::classify(
        
    //     kmer_path);

    let yml = load_yaml!("must.yml");
    let m = App::from_yaml(yml).get_matches();

    // Vary the output based on how many times the user used the "verbose" flag
    // (i.e. 'myprog -v -v -v' or 'myprog -vvv' vs 'myprog -v'
    match m.occurrences_of("v") {
        0 => println!("No verbose info"),
        1 => println!("Some verbose info"),
        2 => println!("Tons of verbose info"),
        3 | _ => println!("Don't be crazy"),
    }
    let verbosity = m.occurrences_of("v") as usize;

    match m.subcommand_name() {
        Some("extract") => {
            let cmd = m.subcommand_matches("compute").unwrap();
            let sequence_files = cmd
                .values_of("sequence_files")
                .map(|vals| vals.collect::<Vec<_>>())
                .unwrap();

            // Convert ksize string argument to integer
            let ksize: u8 = value_t!(cmd, "ksize", u8).unwrap_or_else(|e| e.exit());
            extract::extract_kmers(sequence_files, ksize);
        }
        Some("classify") => {
            let cmd = m.subcommand_matches("compute").unwrap();
            let sequence_files = cmd
                .values_of("sequence_files")
                .map(|vals| vals.collect::<Vec<_>>())
                .unwrap();

            // Convert ksize string argument to integer
            let ksize: u8 = value_t!(cmd, "ksize", u8).unwrap_or_else(|e| e.exit());
            println!("{}", ksize);

            let coding_kmer_file = Path::new(cmd.value_of("coding_kmers").unwrap());
            let non_coding_kmer_file = Path::new(cmd.value_of("non_coding_kmers").unwrap());

            // Convert ksize string argument to integer
            let ksize: u8 = value_t!(cmd, "ksize", u8).unwrap_or_else(|e| e.exit());

            classify_reads::classify(sequence_files,
                                     coding_kmer_file,
                                     non_coding_kmer_file,
                                     ksize, verbosity);
        }
        _ => {
            println!("{:?}", m);
        }
    }
    Ok(())
}
