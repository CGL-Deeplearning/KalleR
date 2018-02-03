import argparse

from modules.ImageCreator import ImageCreator
from modules.BamHandler import BamHandler
from modules.FastaHandler import FastaHandler
from modules.BedHandler import BedHandler


class Bed2ImageAPI:
    """
    Works as a main class and handles user interaction with different modules.
    """
    def __init__(self, bam_file_path, reference_file_path, bed_file_path):
        # --- initialize handlers ---
        self.bam_handler = BamHandler(bam_file_path)
        self.fasta_handler = FastaHandler(reference_file_path)
        self.bed_handler = BedHandler(bed_file_path)

    @staticmethod
    def create_image(bam_handler, fasta_handler, bed_record):
        """
        Iterate through all the reads that fall in a region, find candidates, label candidates and output a bed file.
        :param start_position: Start position of the region
        :param end_position: End position of the region
        :return:
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

    def test(self):
        """
        Test the API
        :return:
        """
        for record in self.bed_handler.all_bed_records:
            self.create_image(self.bam_handler, self.fasta_handler, record)


if __name__ == '__main__':
    '''
    Processes arguments and performs tasks to generate the pileup.
    '''
    parser = argparse.ArgumentParser()
    parser.register("type", "bool", lambda v: v.lower() == "true")
    parser.add_argument(
        "--ref",
        type=str,
        required=True,
        help="Reference corresponding to the BAM file."
    )
    parser.add_argument(
        "--bam",
        type=str,
        required=True,
        help="BAM file containing reads of interest."
    )
    parser.add_argument(
        "--bed",
        type=str,
        required=True,
        help="bed file path."
    )
    FLAGS, unparsed = parser.parse_known_args()

    view = Bed2ImageAPI(FLAGS.bam, FLAGS.ref, FLAGS.bed)
    view.test()

