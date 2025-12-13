def windows2bedgraph(args):
    """
    Convert windows output file to bedGraph format.
    
    Reads the input file line by line and writes to output without loading
    entire file into memory. Assumes the input file is already correctly ordered by genomic coordinates.
    
    The input file should contain columns: chrom, wstart, wend, name, wEI, strand, q_value
    The output bedGraph will have columns: chrom, start, end, value (wEI)
    
    Args:
        input_file (str): Path to input windows file
        output_file (str): Path to output bedGraph file
    """
    input_file = args["input"]
    output_file = args["output"]
    try:
        
        line_count = 0
        with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
            for line in fin:
                line = line.strip()
                if not line:
                    continue
                
                # Skip header line
                if line.startswith("chrom"):
                    continue
                
                # Parse tab-separated fields
                fields = line.split('\t')
                
                if len(fields) < 5:
                    continue
                
                try:
                    chrom = fields[0]
                    start = int(fields[1])
                    end = int(fields[2])
                    value = float(fields[4])
                    # Write bedGraph line
                    fout.write(f"{chrom}\t{start}\t{end}\t{value}\n")
                    line_count += 1
                    
                except (ValueError, IndexError) as e:
                    continue
        
    except FileNotFoundError:
        print(f"Input file not found: {input_file}")
        raise
    except IOError as e:
        print(f"I/O error during conversion: {str(e)}")
        raise
    except Exception as e:
        print(f"Error during conversion: {str(e)}")
        raise
