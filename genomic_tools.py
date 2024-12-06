"""
This program includes a set of tools for DNA sequence analyses. 
To use the toolset, provide the path and filename to the class:
    dna_tools = dna_tool_sets("../data/example.fasta")
Then call class methods, e.g., to check the DNA sequence length:
    length = dna_tools.check_length()
"""

class dna_tool_sets:
    """
    A class for DNA sequence analyses, including:
        a. count_records
        b. check_length
        c. orf_identifier
        d. repeats_identifier

    Parameter: filename (including path)
    Note: The file should be in FASTA format.
    """
    def __init__(self, file_name):
        """Initialize the class by loading and parsing the FASTA file."""
        self.file_name = file_name
        self.dict = {}
        with open(self.file_name) as f_reader:
            for line in f_reader:
                line = line.strip("\n")
                if line.startswith(">"):
                    header = line
                    self.dict[header] = ""
                else:
                    self.dict[header] += line

    def count_records(self):
        """Count the number of records in the file."""
        number_of_records = len(self.dict)
        print(f"Number of records in the multi-FASTA file: {number_of_records}")

    def check_length(self):
        """Check the length of each record and identify the longest/shortest."""
        length_dict = {key: len(value) for key, value in self.dict.items()}
        lengths = list(length_dict.values())
        max_length = max(lengths)
        min_length = min(lengths)
        record_max_length = [item for item in length_dict if length_dict[item] == max_length]
        record_min_length = [item for item in length_dict if length_dict[item] == min_length]
        
        print(f"The longest sequence length: {max_length}")
        print(f"Number of longest sequences: {len(record_max_length)}")
        print(f"The shortest sequence length: {min_length}")
        print(f"Number of shortest sequences: {len(record_min_length)}")

    def find_pos(self, dna):
        """Identify start positions and lengths of ORFs in each reading frame."""
        start_code = "ATG"
        stop_codes = ["TAA", "TAG", "TGA"]
        pos_dict = {}

        for i in range(3):
            pos = []
            frame = [dna[j:j+3] for j in range(i, len(dna), 3)]
            start_pos = [m for m, y in enumerate(frame) if y == start_code]
            stop_pos = [n for n, x in enumerate(frame) if x in stop_codes and n > min(start_pos, default=-1)]

            while start_pos:
                start = min(start_pos)
                try:
                    end = min([stop for stop in stop_pos if stop > start])
                    s_pos = len("".join(frame[:start])) + 1
                    pos.append((s_pos, (end - start + 1) * 3))
                    start_pos.remove(start)
                except ValueError:
                    break
            pos_dict[f"frame{i+1}"] = pos if pos else [(-1, 0)]
        
        return pos_dict

    def revs_complement(self, dna):
        """Generate the reverse complement of a DNA sequence."""
        pairs = {"A": "T", "C": "G", "G": "C", "T": "A"}
        return "".join(pairs[base] for base in dna)[::-1]

    def orf_identifier(self):
        """Identify ORFs and report details like longest ORF length and start position."""
        orf = {header: self.find_pos(seq) for header, seq in self.dict.items()}
        all_frames = [frame for dict_value in orf.values() for frames in dict_value.values() for frame in frames]
        frame2 = [frame for frames in all_frames for frame in frames if "frame2" in frames]
        
        max_frame2 = max(frame2, key=lambda x: x[1], default=(-1, 0))
        print(f"Longest ORF in frame 2: Length={max_frame2[1]}")

        # Extend logic for other ORF-related queries as needed.

    def find_repeats(self, dna, n):
        """Find and count repeats of length n in a DNA sequence."""
        repeats = {}
        for i in range(len(dna) - n + 1):
            repeat = dna[i:i+n]
            repeats[repeat] = repeats.get(repeat, 0) + 1
        return repeats


if __name__ == "__main__":
