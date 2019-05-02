extern crate bio;
extern crate debruijn;
extern crate debruijn_mapping;
extern crate failure;
extern crate hamming;

use debruijn::dna_string::DnaString;
use std::env;
use failure::Error;
use bio::io::{fasta, fastq};
use debruijn_mapping::{config, utils};
use debruijn_mapping::{build_index::build_index,
                       pseudoaligner::process_reads,
                       mappability::analyze_graph};

use std::{
   fs::File,
   io::{prelude::*, BufReader},
   path::Path,
};

use std::collections::HashSet;

fn read_barcodes(filename: impl AsRef<Path>) -> HashSet<String> {
   let file = File::open(filename).expect("no such file");
   let buffer = BufReader::new(file);
   let barcodes: Vec<String> = buffer.lines().map(|l| l.expect("Could not parse line")).collect();
   let mut barcode_set = HashSet::new();
   for barcode in barcodes {
       barcode_set.insert(barcode);
   }
   barcode_set
}

fn main() -> Result<(), Error>  {
    println!("\n**** scrna differential expression tss shit ****\n");
    let args: Vec<String> = env::args().collect();

    let fasta =  &args[1];
    let fastq =  &args[2];
    let output = String::from("tss_coverage.tsv");

    println!("Reading Reference...");
    let reference = fasta::Reader::from_file(fasta)?;
    println!("Finished Reading Reference.\n");

    println!("Reading R1 Fastq...");
    let reads = fastq::Reader::from_file(fastq)?;
    println!("Finished Reading R1.\n");

    println!("Hashing Barcode Set...");
    let barcode_set = read_barcodes("test/737K-august-2016.txt");
    println!("Finished Hashing Barcode Set.\n");

    for result in reads.records() {
        let record = result?;
        let sequence = DnaString::from_acgt_bytes_hashn(record.seq(), record.id().as_bytes());
        let barcode = sequence.slice(0, 16);
        let barcode_str = barcode.to_string();
        let umi = sequence.slice(16, 26);
        let umi_str = umi.to_string();
        if barcode_set.contains(&barcode_str) {
            println!("{:?} - 0", barcode);
        } else {
            for known_barcode in barcode_set.iter() {
                let mut mismatch = 0;
                for (known_nt, nt) in known_barcode.chars().zip(barcode_str.chars()) {
                    if known_nt != nt {
                        mismatch = mismatch + 1;
                    }
                }
                if mismatch == 1 {
                    println!("{:?} - {:?}", barcode, mismatch);
                    break;
                }
            }
        }
    }

    let (seqs, tx_names, tx_gene_map) = utils::read_transcripts(reference)?;
    let index = build_index::<config::KmerType>(&seqs, &tx_names, &tx_gene_map)?;

    let mapped = process_reads::<config::KmerType, _>(reads, &index, output)?;
    for read_data in mapped.iter() {
        if read_data.3 != 0 {
            println!("Read: {:?}", read_data);
        }
    }

    Ok(())
}
