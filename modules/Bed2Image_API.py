from modules.ImageCreator import ImageCreator
from modules.BamHandler import BamHandler
from modules.FastaHandler import FastaHandler
import os

class Bed2ImageAPI:
    """
    Works as a main class and handles user interaction with different modules.
    """
    def __init__(self, bam_file_path, reference_file_path):
        # --- initialize handlers ---
        self.bam_handler = BamHandler(bam_file_path)
        self.fasta_handler = FastaHandler(reference_file_path)

    @staticmethod
    def create_image(bam_handler, fasta_handler, bed_record, output_dir, file_name):
        """
        Create an image from a bed record
        :param bam_handler: Handles bam file
        :param fasta_handler: Handles fasta file
        :param bed_record: Bed record
        :return: Imagearray, label
        """
        chromosome_name, start_position, end_position, ref, alts, genotype, qual, g_filter, in_conf = \
            tuple(bed_record.rstrip().split('\t'))
        start_position = int(start_position)
        genotype = int(genotype)

        pileups = bam_handler.get_pileupcolumns_aligned_to_a_site(chromosome_name, start_position)
        image_creator = ImageCreator(fasta_handler, pileups, chromosome_name, start_position, genotype, alts)

        image_array, image_shape = image_creator.create_image(start_position, ref, alts)
        image_creator.save_image_as_png(image_array, output_dir, file_name)

        return image_array, genotype, image_shape