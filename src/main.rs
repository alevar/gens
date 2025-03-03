extern crate clap;

use clap::{Command,Arg, ArgAction};
/*
========================================
    Module for computing gene densities
========================================
*/

mod compute_densities {
    use std::collections::{HashMap, HashSet};
    use std::fs;
    use std::io::{self, BufRead, Write};
    use std::cmp::Ordering;
    use bio::utils::Interval;
    use bio::data_structures::interval_tree::{ArrayBackedIntervalTree, EntryT};

    #[derive(Debug, Clone)]
    pub struct IntervalEntry<T> {
        pub data: T,
        pub interval: Interval<u64>,
    }

    impl<T> IntervalEntry<T> {
        pub fn new(data: T, interval: Interval<u64>) -> Self {
            IntervalEntry { data, interval }
        }

        pub fn data(&self) -> &T {
            &self.data
        }
    }

    // Implement EntryT for IntervalEntry
    impl<T> EntryT for IntervalEntry<T> {
        type N = u64;
        fn interval(&self) -> &Interval<Self::N> {
            &self.interval
        }
    }

    // Implement EntryT for references to IntervalEntry
    impl<'a, T> EntryT for &'a IntervalEntry<T> {
        type N = u64;

        fn interval(&self) -> &Interval<Self::N> {
            &self.interval
        }
    }

    // Implement PartialEq for IntervalEntry
    impl<T: PartialEq> PartialEq for IntervalEntry<T> {
        fn eq(&self, other: &Self) -> bool {
            self.data == other.data && self.interval == other.interval
        }
    }

    impl<T: Eq> Eq for IntervalEntry<T> {}

    // Implement Ord for IntervalEntry
    impl<T: Ord> Ord for IntervalEntry<T> {
        fn cmp(&self, other: &Self) -> Ordering {
            match self.interval.start.cmp(&other.interval.start) {
                Ordering::Equal => match self.interval.end.cmp(&other.interval.end) {
                    Ordering::Equal => self.data.cmp(&other.data),
                    other => other,
                },
                other => other,
            }
        }
    }

    // Implement PartialOrd for IntervalEntry
    impl<T: Ord> PartialOrd for IntervalEntry<T> {
        fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
            Some(self.cmp(other))
        }
    }

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

    pub fn load_gtf(gtf_file: &str, feature: &str) -> HashMap<String, ArrayBackedIntervalTree<IntervalEntry<String>>> {
        // into an interval tree
        let mut exon_map: HashMap<String, ArrayBackedIntervalTree<IntervalEntry<String>>> = HashMap::new();

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
                                let gene_entry = IntervalEntry::new(gene_id.clone(), Interval::new(start..end).unwrap());
                                map_entry.insert(gene_entry.clone());
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

    pub fn compute_density(faidx: &HashMap<String, usize>, genes: &HashMap<String, ArrayBackedIntervalTree<IntervalEntry<String>>>, interval: &f64) -> HashMap<String, ArrayBackedIntervalTree<IntervalEntry<f64>>> {
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
            let mut dummy_empty_tree: ArrayBackedIntervalTree<IntervalEntry<String>> = ArrayBackedIntervalTree::new();
            dummy_empty_tree.index();
            let gene_tree: &ArrayBackedIntervalTree<IntervalEntry<String>> = match genes.get(seqid) {
                Some(value) => value,
                None => &dummy_empty_tree,
            };

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
                let density_entry = IntervalEntry::new(unique_genes.len() as f64, Interval::new(start..(start + interval.clone() as u64)).unwrap());
                seqid_densities.insert(density_entry.clone());
                // update start
                start += interval.clone() as u64;
            }
        }
        
        densities
    }

    pub fn normalize(densities: &HashMap<String, ArrayBackedIntervalTree<IntervalEntry<f64>>>) -> HashMap<String, ArrayBackedIntervalTree<IntervalEntry<f64>>> {
        // find maximum value
        let mut max_value = 0.0;
        for (_seqid, seqid_density) in densities.iter() {
            for interval in seqid_density {
                let seqid_density = interval.data();
                if seqid_density > &max_value {
                    max_value = *seqid_density as f64;
                }
            }
        }

        // normalize
        let mut normalized_densities = HashMap::new();

        for (seqid, seqid_density) in densities.iter() {
            let seqid_normalized_density = normalized_densities.entry(seqid.clone()).or_insert_with(|| ArrayBackedIntervalTree::new());
            for interval in seqid_density {
                let start = interval.interval().start;
                let end = interval.interval().end;
                let seqid_density = interval.data();
                let density_entry = IntervalEntry::new(seqid_density / max_value, Interval::new(start..end).unwrap());
                seqid_normalized_density.insert(density_entry.clone());
            }
        }

        normalized_densities
    }

    pub fn report(densities: &HashMap<String, ArrayBackedIntervalTree<IntervalEntry<f64>>>, outfname: Option<&String>) {
        match outfname {
            Some(outfname) => {
                let mut out_fp = fs::File::create(outfname).unwrap();
                for (seqid, seqid_density) in densities.iter() {
                    for interval in seqid_density {
                        let start = interval.interval().start;
                        let end = interval.interval().end;
                        let seqid_density = interval.data();
                        writeln!(out_fp, "{}\t{}\t{}\t{}\t{}\t{}", seqid, start, end, '-', seqid_density, '.').unwrap();
                    }
                }
            }, 
            None => {
                for (seqid, seqid_density) in densities.iter() {
                    for interval in seqid_density {
                        let start = interval.interval().start;
                        let end = interval.interval().end;
                        let seqid_density = interval.data();
                        println!("{}\t{}\t{}\t{}\t{}\t{}", seqid, start, end, '-', seqid_density, '.');
                    }
                }
            }
        }
    }
}

fn main() {

    let matches = Command::new("SAM 2 GTF")
        .version("0.1.0")
        .author("Ales Varabyou")
        .about("Compute gene densities over the genome")
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
            .help("Output file containing densities in the bed format: chromosome_name\tstart_position\tend_position\tblank_name\tdensity\tblank_strand")
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
            .default_value("10000000")
            .value_name("FLOAT")
            .value_parser(clap::value_parser!(f64))
        )
        .arg(
            Arg::new("feature")
            .short('f')
            .long("feature")
            .help("Feature to use for computing densities. Default: exon")
            .default_value("exon")
            .value_name("STRING")
        )
        .arg(Arg::new("normalize")
            .short('n')
            .long("normalize")
            .action(ArgAction::SetTrue)
            .help("Normalize densities to [0-1] range")
        )
        .after_help("--help or -h")
        .get_matches();

    let faidx_fname: &String = matches.get_one("faidx").unwrap();
    let gtf_fname: &String = matches.get_one("gtf").unwrap();
    let out_fname: Option<&String> = matches.get_one("output");
    let resolution: &f64 = matches.get_one::<f64>("resolution").unwrap();
    let feature: &String = matches.get_one("feature").unwrap();
    let normalize: bool = matches.get_flag("normalize");

    let seqid_lengths = compute_densities::load_faidx(&faidx_fname);
    let gene_positions = compute_densities::load_gtf(&gtf_fname, feature);
    let mut densities = compute_densities::compute_density(&seqid_lengths, &gene_positions, &resolution);
    
    if normalize {
        densities = compute_densities::normalize(&densities);
    }
    
    compute_densities::report(&densities,out_fname);
}
