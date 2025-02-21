import subprocess
import os
import glob
import sys

def get_cardinality(sketch):
    """
    Get the cardinality of a given HLL sketch using Dashing.
    
    :param sketch: Path to the HLL sketch file.
    :return: Estimated cardinality (integer).
    """
    result = subprocess.run(["dashing", "card", "--presketched", sketch], 
                            capture_output=True, text=True, check=True)
    result = result.stdout.strip().split('\t')[-1]
    return float(result)

def union_sketch(bin_sketch, genome_sketch, output_sketch):
    """
    Compute the union of two HLL sketches and save it as a new sketch.
    """
    subprocess.run(["dashing", "union", "-o", output_sketch, bin_sketch, genome_sketch],
                    check=True, stdout=subprocess.DEVNULL)

def extract_genome_name(filepath):
    """
    Extract the genome name from a given file path.

    Example:
    path/to/file/SAMEA897824.fa.gz.w.31.spacing.10.hll → SAMEA897824
    """
    filename = os.path.basename(filepath)  # Get the file name
    parts = filename.split(".")  # Split by '.'
    # Find the first part that does not contain known extensions
    for part in parts:
        if part not in ["fa", "fq", "gz", "w", "31", "spacing", "10", "hll"]:
            return part  # Return the first valid part as genome name

    return filename  # If nothing found, return original filename

def greedy_partitioning(numbins : int, sketches_dir="tmp/sketches", bin_dir="tmp/bins"):
    
    os.makedirs(bin_dir, exist_ok=True)
    bin_result = []  # A list of lists, each sub-list is a bin containing genome names
    bin_paths = []   # A list storing paths of the bin sketches
    bin_cardinalities = []
    
    # Get all .hll files in tmp/sketches
    genome_sketches = sorted(glob.glob(os.path.join(sketches_dir, "*.hll")))
    
    if not genome_sketches:
        return bin_result,bin_cardinalities     
    
    # Initialize bins
    bin_result = [[] for _ in range(numbins)]
    bin_cardinalities = [0] * numbins
    bin_paths = [os.path.join(bin_dir, f"bin_{i}.hll") for i in range(numbins)]
    
    # using an algorithm called greedy number partitioning
    # the algorithm works as follows: it loops over the items from large to small, 
    # and puts the next item into a bin that currently contains the smallest total size.    

    # pseudocode:
    # calculate the distinct kmers (cardinality/card) of all genomes
    # Compute cardinalities and create a list of genome data [(name,sketch,card)]
    genome_data = []
    for sketch in genome_sketches:
        genome_name = extract_genome_name(sketch)
        cardinality = get_cardinality(sketch)
        genome_data.append((genome_name, sketch, cardinality))   

    # sort genomes based on the nb of distinct kmers from large to small 
    # Sort genomes by decreasing cardinality
    genome_data.sort(key=lambda x: x[2], reverse=True)

    # Assign first `numbins` genomes directly to bins
    for i in range(min(numbins, len(genome_data))):
        genome_name, genome_sketch, cardinality = genome_data[i]
        bin_result[i].append(genome_name)
        bin_cardinalities[i] = cardinality
        subprocess.run(["cp", genome_sketch, bin_paths[i]], check=True)    

    # get the current bin with the current smallest card
    # put the genome into the bin
    # update the bin sketch with the genome sketch using dashing union
    # calculate the card of the new merged bin
    # loop over and the genomes and do the same
    # Assign remaining genomes using greedy binning

    tmp_union_sketch = os.path.join(bin_dir, "tmp_union_sketch.hll")  # Temporary file

    for genome_name, genome_sketch, cardinality in genome_data[numbins:]:
        min_bin_index = bin_cardinalities.index(min(bin_cardinalities))
        # Compute union sketch
        union_sketch(bin_paths[min_bin_index], genome_sketch, tmp_union_sketch)
        subprocess.run(["mv", tmp_union_sketch, bin_paths[min_bin_index]], check=True)
        bin_result[min_bin_index].append(genome_name)
        bin_cardinalities[min_bin_index] = get_cardinality(bin_paths[min_bin_index])

    return bin_result,bin_cardinalities


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python greedy_partitioning_hyperloglog.py <num_bins> <file_name>")
        sys.exit(1)

    num_bins = int(sys.argv[1])
    file_name = sys.argv[2]

    os.makedirs("output", exist_ok=True)
    os.makedirs("tmp/completion", exist_ok=True)

    bins, cardinalities = greedy_partitioning(num_bins)

    # Save bin assignments
    with open(f'output/{file_name}_bin_assignment.txt', "w") as f:
        for i, (bin_content, bin_card) in enumerate(zip(bins, cardinalities)):
            f.write(f"Bin {i}: {', '.join(bin_content)}; Cardinality: {bin_card}\n")

    # Create completion marker
    open(f'tmp/completion/{file_name}_binned.done', "w").close()

    print(f"Binning completed. Results saved to output/{file_name}_bin_assignment.txt")