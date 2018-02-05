import argparse
import random

from modules.BedHandler import BedHandler

"""
Subsample the bed file for homozygous cases for training faster
"""


class BedSubsampler:
    """
    Downsamples a bed file
    """
    def __init__(self, bed_file_path):
        """
        Initialize a bed subsampler
        :param bed_file_path:
        """
        # --- initialize handlers ---
        self.bed_handler = BedHandler(bed_file_path)

    @staticmethod
    def select_or_not(bed_record, downsample_rate):
        """
        Determines if a bed record should be selected given a downsampling rate
        :param bed_record: A bed record
        :param downsample_rate: A downsampling probability
        :return: Boolean
        """
        chromosome_name, start_position, end_position, ref, alts, genotype = tuple(bed_record.rstrip().split('\t'))
        genotype = int(genotype)
        # if not homozygous, always pick
        if genotype != 0:
            return True

        # else do a sampling based on probability
        random_chance = random.uniform(0, 1)
        if random_chance <= downsample_rate:
            return True
        return False

    def downsample_bed_file(self):
        """
        Downsample the bed file
        :return:
        """
        # calculate the downsample rate based on distribution of three classes
        downsample_rate = max(self.bed_handler.total_het, self.bed_handler.total_hom_alt) / self.bed_handler.total_hom

        for i, record in enumerate(self.bed_handler.all_bed_records):
            if BedSubsampler.select_or_not(record, downsample_rate) is True:
                print(record, end='')

    def print_bed_stats(self):
        """
        Print the distribution of the bed file
        :return:
        """
        print("Total records: ")
        print(self.bed_handler.total_hom, self.bed_handler.total_het, self.bed_handler.total_hom_alt)
        total_cases = self.bed_handler.total_hom + self.bed_handler.total_het + self.bed_handler.total_hom_alt
        print("Percent homozygous records:\t", int(self.bed_handler.total_hom * 100 / total_cases))
        print("Percent Heterozygous records:\t", int(self.bed_handler.total_het * 100 / total_cases))
        print("Percent Hom-alt records:\t", int(self.bed_handler.total_hom_alt * 100 / total_cases))


if __name__ == '__main__':
    '''
    Processes arguments and performs tasks to generate the pileup.
    '''
    parser = argparse.ArgumentParser()
    parser.register("type", "bool", lambda v: v.lower() == "true")
    parser.add_argument(
        "--bed",
        type=str,
        required=True,
        help="bed file path."
    )
    parser.add_argument(
        "--stats",
        type=bool,
        help="If true, will print bed file stats."
    )
    FLAGS, unparsed = parser.parse_known_args()
    down_sampler = BedSubsampler(FLAGS.bed)

    if FLAGS.stats is True:
        down_sampler.print_bed_stats()
    else:
        down_sampler.downsample_bed_file()

