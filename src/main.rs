extern crate clap;

use clap::{Command,Arg};
/*
========================================
    Module for computing gene densities
========================================
*/
mod compute_densities {
    use std::collections::{HashMap, HashSet};
    use std::fs;
    use std::io::{self, BufRead, Write};
    use bio::utils::Interval;
    use bio::data_structures::interval_tree::ArrayBackedIntervalTree;

    pub fn load_faidx(faidx_fname: &str) -> HashMap<String, usize> {
        let mut faidx = HashMap::new();

        if let Ok(file) = fs::File::open(faidx_fname) {
            let reader = io::BufReader::new(file);

            for line in reader.lines() {
                if let Ok(line) = line {
                    let lcs: Vec<&str> = line.split('\t').collect();
                    let seqid = lcs[0].to_string();
                    let length = lcs[1].trim().parse().unwrap();
                    faidx.insert(seqid, length);
                }
            }
        }

        faidx
    }

    pub fn load_gtf(gtf_file: &str, feature: &str) -> HashMap<String, ArrayBackedIntervalTree<u64, String>> {
        // into an interval tree
        let mut exon_map: HashMap<String, ArrayBackedIntervalTree<u64, String>> = HashMap::new();

        if let Ok(file) = fs::File::open(gtf_file) {
            let reader = io::BufReader::new(file);

            for line in reader.lines() {
                if let Ok(line) = line {
                    let lcs: Vec<&str> = line.split('\t').collect();
                    if lcs.len() == 9 && lcs[2] == feature {
                        let seqid = lcs[0].to_string();
                        let start: u64 = lcs[3].trim().parse().unwrap_or(0);
                        let end: u64 = lcs[4].trim().parse().unwrap_or(0);

                        // initialize interval tree for the chromosome
                        let map_entry = exon_map.entry(seqid.clone()).or_insert_with(|| ArrayBackedIntervalTree::new());

                        // find gene_id in attributes. if not found - skip
                        let attributes: Vec<&str> = lcs[8].split(';').collect();
                        match attributes.iter().find(|&&x| x.trim().starts_with("gene_id ")) {
                            Some(gene_id) => {
                                let gene_id: Vec<&str> = gene_id.trim().split(' ').collect();
                                let gene_id = gene_id[1].trim_matches('"').to_string();
                                map_entry.insert(Interval::new(start..end).unwrap(), gene_id.clone());
                            },
                            None => continue,
                        }
                    }
                }
            }
        }

        // perform indexing of the interval trees
        for (_seqid, tree) in exon_map.iter_mut() {
            tree.index();
        }
        exon_map
    }

    pub fn compute_density(faidx: &HashMap<String, usize>, genes: &HashMap<String, ArrayBackedIntervalTree<u64,String>>, interval: &u32) -> HashMap<String, ArrayBackedIntervalTree<u64,u32>> {
        // iterate over faidx
        // for every chromosome iterate in intervals of specified resolution
        // for every interval pull slice of the interval tree
        // compute total number of genes in the slice

        // after everything has been loaded - normalize to [0-1] based on max observed value
        
        let mut densities = HashMap::new();
        
        // iterate over the faidx
        for (seqid, length) in faidx.iter() {
            let seqid_densities = densities.entry(seqid.clone()).or_insert_with(|| ArrayBackedIntervalTree::new());
            let mut start = 0;
            let end = *length as u64;

            // get reference to the interval tree for the chromosome
            let gene_tree = genes.get(seqid).unwrap();

            // iterate over the intervals of specified resolution/interval
            while start < end {
                // get all intervals in the gene_tree that overlap current interval
                let overlapping_intervals = gene_tree.find(start..(start + interval.clone() as u64));

                // count number of unique gene_ids
                let mut unique_genes: HashSet<String> = HashSet::new();
                
                for overlap in overlapping_intervals {
                    // get gene_id
                    let gene_id = overlap.data();
                    unique_genes.insert(gene_id.to_string());
                }

                // update gene_count
                seqid_densities.insert(Interval::new(start..(start + interval.clone() as u64)).unwrap(), unique_genes.len() as u32);

                // update start
                start += interval.clone() as u64;
            }
        }
        
        densities
    }

    pub fn report(densities: &HashMap<String, ArrayBackedIntervalTree<u64,u32>>, outfname: Option<&String>) {
        match outfname {
            Some(outfname) => {
                let mut out_fp = fs::File::create(outfname).unwrap();
                for (seqid, seqid_density) in densities.iter() {
                    for interval in seqid_density {
                        let start = interval.interval().start;
                        let end = interval.interval().end;
                        let seqid_density = interval.data();
                        writeln!(out_fp, "{}\t{}\t{}\t{}", seqid, start, end, seqid_density).unwrap();
                    }
                }
            }, 
            None => {
                for (seqid, seqid_density) in densities.iter() {
                    for interval in seqid_density {
                        let start = interval.interval().start;
                        let end = interval.interval().end;
                        let seqid_density = interval.data();
                        println!("{}\t{}\t{}\t{}", seqid, start, end, seqid_density);
                    }
                }
            }
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
                .short('g')
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
                .required(false)
                .value_name("FILE")
            )
            .arg(
                Arg::new("faidx")
                .short('i')
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
                .value_parser(clap::value_parser!(u32))
            )
            .arg(
                Arg::new("feature")
                .short('f')
                .long("feature")
                .help("Feature to use for computing densities. Default: exon")
                .default_value("exon")
                .value_name("STRING")
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

    match matches.subcommand() {
        Some(("compute",  sub_matches)) => {
            let faidx_fname: &String = sub_matches.get_one("faidx").unwrap();
            let gtf_fname: &String = sub_matches.get_one("gtf").unwrap();
            let out_fname: Option<&String> = sub_matches.get_one("output");
            let resolution: &u32 = sub_matches.get_one::<u32>("resolution").unwrap();
            let feature: &String = sub_matches.get_one("feature").unwrap();

            let seqid_lengths = compute_densities::load_faidx(&faidx_fname);
            for(seqid, &length) in seqid_lengths.iter() {
                println!("{}: {}", seqid, length);
            }
            let gene_positions = compute_densities::load_gtf(&gtf_fname, feature);
            let densities = compute_densities::compute_density(&seqid_lengths, &gene_positions, &resolution);
            compute_densities::report(&densities,out_fname);
        },
        Some(("zoom", sub_matches)) => {
            let zoom_factor: &f64 = sub_matches.get_one::<f64>("zoom_factor")
                                                .unwrap();
                                            
        }
        _ => println!("Please use one of the available subcommands: compute, zoom"),
    }
}
