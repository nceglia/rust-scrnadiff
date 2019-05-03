extern crate bio;
extern crate debruijn;
extern crate debruijn_mapping;
extern crate failure;
extern crate rayon;

use debruijn::dna_string::{DnaString,DnaStringSlice};
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

use rayon::prelude::*;
use std::collections::HashSet;
use std::collections::HashMap;
use std::iter::FromIterator;


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

fn mismatches(sequence : &String, reference : &String) -> u32 {
    let mut mismatch : u32 = 0;
    for (known_nt, nt) in reference.chars().zip(sequence.chars()) {
        if known_nt != nt {
            mismatch = mismatch + 1;
        }
    }
    mismatch
}

fn parse_barcode(record: &bio::io::fastq::Record) -> String {
    let seqid = record.id().to_owned();
    let sequence = DnaString::from_acgt_bytes_hashn(record.seq(), record.id().as_bytes());
    let barcode = sequence.slice(0, 16);
    barcode.to_string()
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
    let records = fastq::Reader::from_file(fastq)?;
    println!("Finished Reading R1.\n");

    println!("Hashing Barcode Set...");
    let barcode_set = read_barcodes("test/3M-february-2018.txt");
    println!("Finished Hashing Barcode Set.\n");

    // let mut valid_barcodes = HashMap::new();
    let mut read_records = Vec::new();

    for record_ref in records.records() {
        let record = record_ref?;
        read_records.push(record.to_owned());
    }
    println!("\nParsing Sequence Barcodes..");
    let barcodes = read_records.par_iter().map(|x| parse_barcode(x));
    println!("Finished.");

    println!("\nFinding Valid Barcode Set...");
    let mut valid_barcodes : Vec<String> = barcodes.clone().filter(|x| barcode_set.contains(x)).collect();
    valid_barcodes.par_sort_unstable();
    valid_barcodes.dedup();
    println!("Perfectly Matched Barcodes: {:?}", valid_barcodes.len());

    println!("\nFinding Invalid Barcode Set...");
    let mut invalid_barcodes : Vec<String> = barcodes.clone().filter(|x| !barcode_set.contains(x)).collect();
    invalid_barcodes.par_sort_unstable();
    invalid_barcodes.dedup();
    println!("Invalid Barcodes: {:?}", invalid_barcodes.len());

    let barcode_to_mismatches = invalid_barcodes.par_iter().map(|x| barcode_set.par_iter().map(move |y| mismatches(x,y))).collect();

    // for record_ref in records.records() {
    //     let record = record_ref?;
    //     let seqid = record.id().to_owned();
    //     let sequence = DnaString::from_acgt_bytes_hashn(record.seq(), record.id().as_bytes());
    //     let barcode = sequence.slice(0, 16);
    //     let barcode_str = barcode.to_string();
    //     if barcode_set.contains(&barcode_str) {
    //         valid_barcodes.insert(seqid,barcode_str);
    //     } else {
    //         let one_dist_corrected = barcode_set.par_iter().map(|x| mismatches(barcode_str.clone(), x.to_string())).filter(|&x| x == 1).collect()?;
            // for bcs in one_dist_corrected {
            //     println!("{:?} Records with 1 mismatch.", bcs);
            // }

            // for known_barcode in  {
            //     let mut mismatch = mismatches(barcode_str, known_barcode)
            // if mismatch == 1 {
            //
            //     break;
            // }
            // }
    //     }
    // }
    // println!("\n{:?} - Valid Barcodes", valid_barcodes.keys().len());


    // let (seqs, tx_names, tx_gene_map) = utils::read_transcripts(reference)?;
    // let index = build_index::<config::KmerType>(&seqs, &tx_names, &tx_gene_map)?;
    //
    // let mapped = process_reads::<config::KmerType, _>(reads, &index, output)?;
    // for read_data in mapped.iter() {
    //     if read_data.3 != 0 {
    //         println!("Read: {:?}", read_data);
    //     }
    // }

    Ok(())
}
