#! /usr/bin/env python2

import vcf.utils
from cyvcf2 import VCF
import hgvs.parser
import os
import subprocess

__author__ = "Josef Moser"


class Assignment3:
    def __init__(self):
        print("PyVCF version: %s" % vcf.VERSION)
        print("HGVS version: %s" % hgvs.__version__)
        self.son = self.get_file("ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/"
                                 "analysis/IonTorrent_TVC_03162015/AmpliseqExome.20141120.NA24385.vcf")
        self.mother = self.get_file("ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/"
                                    "analysis/IonTorrent_TVC_03162015/AmpliseqExome.20141120.NA24143.vcf")
        self.father = self.get_file("ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/"
                                    "analysis/IonTorrent_TVC_03162015/AmpliseqExome.20141120.NA24149.vcf")
        # check dependencies:
        with open(os.devnull, "w") as fnull:
            x = subprocess.call("vcf-merge -h", stdout=fnull, stderr=subprocess.PIPE, shell=True)
            y = subprocess.call("tabix -h", stdout=fnull, stderr=subprocess.PIPE, shell=True)
            if x == 127:
                print "vcftools is not installed - run: sudo install_vcftools.py"
                exit()
            if y == 127:
                print "tabix is not installed: try sudo apt-get install tabix"
                exit()

    @staticmethod
    def get_file(ftp_link):
        bn = os.path.basename(ftp_link)
        if not os.path.isfile(bn):
            subprocess.call(["wget", ftp_link])
        if not os.path.isfile(bn + ".gz"):
            subprocess.call(" ".join(["bgzip", "-c", bn, ">", bn + ".gz"]), shell=True)
        if not os.path.isfile(bn + ".gz.tbi"):
            subprocess.call(["tabix", "-p", "vcf", bn + ".gz"])

        return bn + ".gz"

    def get_total_number_of_variants_mother(self):
        """
        Return the total number of identified variants in the mother
        :return:
        """
        print "Variants in mother: %d" % len(list(VCF(self.mother)))

    def get_total_number_of_variants_father(self):
        """
        Return the total number of identified variants in the father
        :return:
        """
        print "Variants in father: %d" % len(list(VCF(self.father)))

    @staticmethod
    def compare_parent_child(parent, child):
        both = []
        for p, c in vcf.utils.walk_together(parent, child):
            if p and c \
                    and p.CHROM == c.CHROM \
                    and p.ALT == c.ALT \
                    and p.REF == c.REF \
                    and p.POS == c.POS \
                    and p.INFO.get("GT") == c.INFO.get("GT"):
                both.append(c)

        return both

    def get_variants_shared_by_father_and_son(self):
        """
        Return the number of identified variants shared by father and son
        :return:
        """
        self.father_and_son = self.compare_parent_child(VCF(self.father), VCF(self.son))
        print "Variants shared by father and son: %d" % len(self.father_and_son)

    def get_variants_shared_by_mother_and_son(self):
        """
        Return the number of identified variants shared by mother and son
        :return:
        """
        self.mother_and_son = self.compare_parent_child(VCF(self.mother), VCF(self.son))
        print "Variants shared by mother and son: %d" % len(self.mother_and_son)

    def get_variants_shared_by_trio(self):
        """
        Return the number of identified variants shared by father, mother and son
        :return:
        """
        m = VCF(self.mother)
        f = VCF(self.father)
        s = VCF(self.son)
        all = []
        for x in vcf.utils.walk_together(m, f, s):
            m, f, s = x[0], x[1], x[2]
            if m and f and s \
                    and m.CHROM == f.CHROM == s.CHROM \
                    and m.ALT == f.ALT == s.ALT \
                    and m.REF == f.REF == s.REF \
                    and m.POS == f.POS == s.POS \
                    and m.INFO.get("GT") == f.INFO.get("GT") == s.INFO.get("GT"):
                all.append(m)

        print "Variants shared by father, mother and son: %d" % len(all)

    def merge_mother_father_son_into_one_vcf(self):
        """
        Creates one VCF containing all variants of the trio (merge VCFs)
        :return:
        """

        subprocess.call(" ".join(["vcf-merge", self.father, self.mother, self.son, ">", "test.vcf"]), shell=True)

    def convert_first_variants_of_son_into_HGVS(self):
        """
        Convert the first 100 variants identified in the son into the corresponding transcript HGVS.
        Each variant should be mapped to all corresponding transcripts. Pointer:
        - https://hgvs.readthedocs.io/en/master/examples/manuscript-example.html#project-genomic-variant-to-a-new-transcript
        :return:
        """


    def print_summary(self):
        print(__author__)
        self.get_total_number_of_variants_mother()
        self.get_total_number_of_variants_father()
        self.get_variants_shared_by_father_and_son()
        self.get_variants_shared_by_mother_and_son()
        self.get_variants_shared_by_trio()
        self.merge_mother_father_son_into_one_vcf()


if __name__ == '__main__':
    print("Assignment 3")
    assignment1 = Assignment3()
    assignment1.print_summary()
