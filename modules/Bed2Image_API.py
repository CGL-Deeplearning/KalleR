from modules.ImageCreator import ImageCreator
from modules.BamHandler import BamHandler
from modules.FastaHandler import FastaHandler


class Bed2ImageAPI:
    """
    Works as a main class and handles user interaction with different modules.
    """
    def __init__(self, bam_file_path, reference_file_path):
        # --- initialize handlers ---
        self.bam_handler = BamHandler(bam_file_path)
        self.fasta_handler = FastaHandler(reference_file_path)

    @staticmethod
    def create_image(bam_handler, fasta_handler, bed_record):
        """
        Create an image from a bed record
        :param bam_handler: Handles bam file
        :param fasta_handler: Handles fasta file
        :param bed_record: Bed record
        :return: Imagearray, label
        """
        chromosome_name, start_position, end_position, ref, alts, genotype = tuple(bed_record.rstrip().split('\t'))
        start_position = int(start_position)
        end_position = int(end_position)
        genotype = int(genotype)

        reads = bam_handler.get_reads(chromosome_name=chromosome_name, start=start_position, stop=end_position+1)
        image_creator = ImageCreator(fasta_handler, chromosome_name, start_position, end_position)

        image_creator.process_reads(reads)

        image_array = image_creator.generate_image(start_position, alts)
        return image_array, genotype