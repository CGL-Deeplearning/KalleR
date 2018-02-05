import argparse
import time
import os
from scipy import misc
import sys
sys.path.insert(0,'..')

from modules.ImageCreatorRGB import ImageCreatorRGB
from modules.BamHandler import BamHandler
from modules.FastaHandler import FastaHandler
from modules.BedHandler import BedHandler


class Bed2ImageAPI:
    """
    Works as a main class and handles user interaction with different modules.
    """
    def __init__(self, bam_file_path, reference_file_path, bed_file_path, output_dir):
        # --- initialize handlers ---
        self.bam_handler = BamHandler(bam_file_path)
        self.fasta_handler = FastaHandler(reference_file_path)
        self.bed_handler = BedHandler(bed_file_path)
        self.output_dir = output_dir

    @staticmethod
    def save_image_rgb(image_array, chr_name, start_position, alts, genotype, output_dir):
        output_file_name = output_dir + chr_name + "_" + str(start_position) + "_" + '_'.join(alts) + "_" + str(genotype)
        misc.imsave(output_file_name + ".png", image_array, format="PNG")

    @staticmethod
    def create_image_rgb(bam_handler, fasta_handler, bed_record, output_dir):
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
        image_creator = ImageCreatorRGB(fasta_handler, chromosome_name, start_position, end_position)

        image_creator.process_reads(reads)

        image_array = image_creator.generate_image(start_position, alts)
        Bed2ImageAPI.save_image_rgb(image_array, chromosome_name, start_position, alts, genotype, output_dir)

    def test(self):
        """
        Test the API
        :return:
        """
        for record in self.bed_handler.all_bed_records:
            print(record)
            self.create_image_rgb(self.bam_handler, self.fasta_handler, record, self.output_dir)


def handle_output_directory(output_dir):
    """
    Process the output directory and return a valid directory where we save the output
    :param output_dir: Output directory path
    :return:
    """
    # process the output directory
    if output_dir[-1] != "/":
        output_dir += "/"
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # create an internal directory so we don't overwrite previous runs
    timestr = time.strftime("%m%d%Y_%H%M%S")
    internal_directory = "run_" + timestr + "/"
    output_dir = output_dir + internal_directory

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    return output_dir


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
    parser.add_argument(
        "--save_img",
        type=bool,
        help="If true the images will be saved in output directory"
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        default="output/",
        help="Path to output directory."
    )

    FLAGS, unparsed = parser.parse_known_args()
    if FLAGS.save_img is True:
        FLAGS.output_dir = handle_output_directory(FLAGS.output_dir)

    view = Bed2ImageAPI(FLAGS.bam, FLAGS.ref, FLAGS.bed, FLAGS.output_dir)
    view.test()

"""
python3 bed2RGBimgAPI.py --bam ~/Kishwar/Whole_chr3_data/illumina/vcf_whole_chr/chr3.bam \
--ref ~/Kishwar/Whole_chr3_data/illumina/vcf_whole_chr/chr3.fa \
--bed ~/Kishwar/software/KalleR/test_bed_file/sampled_train.bed \
--output_dir ../output/ \
--save_img True 

"""
