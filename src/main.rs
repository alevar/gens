extern crate clap;

use clap::{Command,Arg, ArgAction};

use std::collections::HashMap;

/*
========================================
    Module for computing gene densities
========================================
*/
mod compute_densities {
    use std::collections::HashMap;
    use std::fs;
    use std::io::{self, BufRead};


    pub fn load_faidx(faidx_fname: &str) -> HashMap<String, usize> {
        let mut seqid_lengths = HashMap::new();

        if let Ok(file) = fs::File::open(faidx_fname) {
            let reader = io::BufReader::new(file);

            for line in reader.lines() {
                if let Ok(line) = line {
                    let lcs: Vec<&str> = line.split('\t').collect();
                    if lcs.len() == 2 {
                        let seqid = lcs[0].to_string();
                        let length = lcs[1].trim().parse().unwrap_or(0);
                        seqid_lengths.insert(seqid, length);
                    }
                }
            }
        }

        seqid_lengths
    }

    pub fn load_gtf(gtf_file: &str) -> HashMap<String, Vec<(usize, usize)>> {
        // load specified features (TODO: add new argument for this)
        // into an interval tree
        unimplemented!()
    }

    pub fn compute_density(seqid_lengths: &HashMap<String, usize>, genes: &HashMap<String, Vec<(usize, usize)>>) -> HashMap<String, Vec<u32>> {
        // iterate over the contents of the interval tree in sorted order
        // and compute the densities

        // iterate over faidx
        // for every chromosome iterate in intervals of specified resolution
        // for every interval pull slice of the interval tree
        // compute total number of genes in the slice

        // after everything has been loaded - normalize to [0-1] based on max observed value
        unimplemented!()
    }

    pub fn report_results(densities: &HashMap<String, Vec<u32>>) {
        for (seqid, density) in densities.iter() {
            println!("Chromosome {}: {:?}", seqid, density);
        }
    }
}

/*
========================================
    Module for manipulating densities
========================================
*/
mod zoom_densities {
    use std::collections::HashMap;

    pub fn zoom_in(densities: &HashMap<String, Vec<u32>>, zoom_factor: &f64) -> HashMap<String, Vec<u32>> {
        unimplemented!()
    }
}

fn main() {
    let matches = Command::new("SAM 2 GTF")
        .version("0.1.0")
        .author("Ales Varabyou")
        .about("Compute gene densities over the genome")
        .subcommands( [
            Command::new("compute")
            .about("Controls configuration functionality")
            .arg(
                Arg::new("gtf")
                .short('i')
                .long("gtf")
                .help("Gene annotation in GTF format.")
                .required(true)
                .value_name("FILE")
            )
            .arg(
                Arg::new("output")
                .short('o')
                .long("output")
                .help("Output file containing densities in the format: chromosome_name\tstart_position\tend_position\tdensity")
                .required(true)
                .value_name("FILE")
            )
            .arg(
                Arg::new("faidx")
                .short('f')
                .long("faidx")
                .help("FASTA index file")
                .required(true)
                .value_name("FILE")
            )
            .arg(
                Arg::new("resolution")
                .short('r')
                .long("resolution")
                .help("Resolution for which densities will be computed. If 1000000 is specified, densities will be computed for 1Mb windows.")
                .default_value("1000000")
                .value_name("INT")
                .value_parser(clap::value_parser!(usize))
            ),
            Command::new("zoom")
            .about("Changes resolution of the density map")
            .arg(
                Arg::new("zoom_factor")
                .short('x')
                .long("zoom_factor")
                .help("Magnification factor")
                .required(true)
                .value_name("FLOAT")
                .value_parser(clap::value_parser!(f64))
            )
            .arg(
                Arg::new("input")
                .short('i')
                .long("input")
                .help("Input file containing densities in the format: chromosome_name\tstart_position\tend_position\tdensity")
                .required(true)
                .value_name("FILE")
            )
            .arg(
                Arg::new("output")
                .short('o')
                .long("output")
                .help("Output file containing densities in the format: chromosome_name\tstart_position\tend_position\tdensity")
                .required(true)
                .value_name("FILE")
            )]
        )
        .after_help("--help or -h")
        .get_matches();

    // // Extract and use the value of the "input" argument
    // let gtf: &String = matches.get_one("gtf").unwrap();
    // let faidx: &String = matches.get_one("faidx").unwrap();
    // let output_fname: &String = matches.get_one("output").unwrap();

    match matches.subcommand() {
        Some(("compute",  sub_matches)) => {
            let faidx_fname: &String = sub_matches.get_one("faidx").unwrap();
            let gtf_fname: &String = sub_matches.get_one("gtf").unwrap();
            let out_fname: &String = sub_matches.get_one("output").unwrap();
            let resolution: &usize = sub_matches.get_one::<usize>("resolution").unwrap();

            println!("Input GTF file: {}", gtf_fname);
            println!("Input FASTA index file: {}", faidx_fname);
            println!("Output file: {}", out_fname);
            println!("Resolution: {}", resolution);


            let seqid_lengths = compute_densities::load_faidx(&faidx_fname);
            let gene_positions = compute_densities::load_gtf(&gtf_fname);
            let densities = compute_densities::compute_density(&seqid_lengths, &gene_positions);
            compute_densities::report_results(&densities);
        },
        Some(("zoom", sub_matches)) => {
            let zoom_factor: &f64 = sub_matches.get_one::<f64>("zoom_factor")
                                                .unwrap();
                                            
        }
        _ => println!("Please use one of the available subcommands: compute, zoom"),
    }
}
