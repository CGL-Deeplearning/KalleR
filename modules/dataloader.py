import torch
from modules.BedHandler import BedHandler
from modules.Bed2Image_API import Bed2ImageAPI

from torch.utils.data import Dataset
import os
from PIL import Image, ImageOps
import numpy as np


class DataSetLoader(Dataset):
    def __init__(self, bam_file_path, fasta_file_path, bed_file_path, img_output_dir, transform, img_w=300, img_h=300, img_c=7):
        self.bam_file_path = bam_file_path
        self.fasta_file_path = fasta_file_path
        self.bed_handler = BedHandler(bed_file_path)

        self.all_bed_records = self.bed_handler.all_bed_records
        self.img_output_dir = img_output_dir
        self.transform = transform
        self.shape = (img_w, img_h, img_c)

        self.generated_files = {}

    def __getitem__(self, index):
        bed_record = self.all_bed_records[index]
        contig, pos_s, pos_e, ref, alt, genotype, qual, gen_filter, in_confident = bed_record.rstrip().split('\t')
        file_name = contig + "_" + pos_s + "_" + alt + "_" + genotype
        summary_string = ''
        if os.path.isfile(os.path.abspath(self.img_output_dir + file_name) + ".png"):
            # read the file
            file = os.path.abspath(self.img_output_dir + file_name) + ".png"
            label = int(genotype)
            shape = self.shape

            img = Image.open(file)
            np_array_of_img = np.array(img.getdata())

            img = np.reshape(np_array_of_img, shape)
            img = np.transpose(img, (0, 1, 2))
            label = torch.LongTensor([label])
        else:
            # save the file
            api_object = Bed2ImageAPI(self.bam_file_path, self.fasta_file_path)
            img, label, img_shape = api_object.create_image(api_object.bam_handler, api_object.fasta_handler,
                                                 bed_record, self.img_output_dir, file_name)
            label = torch.LongTensor([label])
            summary_string += os.path.abspath(self.img_output_dir + file_name) + ".png," + str(genotype) + ',' \
                              + ','.join(map(str, img_shape)) + "," + str(qual) + "," + str(gen_filter) + "," \
                              + str(in_confident) + '\n'

        if self.transform is not None:
            img = self.transform(img)
        return img, label, bed_record, summary_string

    def __len__(self):
        return len(self.all_bed_records)
