import torch
from modules.BedHandler import BedHandler
from modules.Bed2Image_API import Bed2ImageAPI

from torch.utils.data import Dataset


class DataSetLoader(Dataset):
    def __init__(self, bam_file_path, fasta_file_path, bed_file_path, transform=None):
        self.bam_file_path = bam_file_path
        self.fasta_file_path = fasta_file_path
        self.bed_handler = BedHandler(bed_file_path)

        self.all_bed_records = self.bed_handler.all_bed_records

        self.transform = transform

    def __getitem__(self, index):
        bed_record = self.all_bed_records[index]
        api_object = Bed2ImageAPI(self.bam_file_path, self.fasta_file_path)
        img, label = api_object.create_image(api_object.bam_handler, api_object.fasta_handler, bed_record)
        label = torch.LongTensor([label])
        if self.transform is not None:
            img = self.transform(img)
        return img, label, bed_record

    def __len__(self):
        return len(self.all_bed_records)
