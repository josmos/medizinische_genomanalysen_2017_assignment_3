#! /usr/bin/env python2

import os
import subprocess
import hgvs.assemblymapper
import hgvs.dataproviders.uta
import hgvs.exceptions
import hgvs.parser
import vcf.utils
from bioutils.assemblies import make_name_ac_map
from cyvcf2 import VCF

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
        hdp = hgvs.dataproviders.uta.connect()  # Connect to UTA
        vm = hgvs.assemblymapper.AssemblyMapper(hdp)  # Used to get the transcripts
        hp = hgvs.parser.Parser()  # Used for parsing
        ## Now parse the variant
        ## http://hgvs.readthedocs.io/en/master/modules/io.html?highlight=parser_hgvs
        genome_hgvs = []
        cds_transcipt_hgvs = []
        noncoding_transcript_hgvs = []
        unmappable = []
        outfile = "Transcipt_mapping_results.txt"
        out_fh = open(outfile, "w")

        def map_transcript(g):
            for tr in vm.relevant_transcripts(g):
                try:
                    c = vm.g_to_c(g, tr)  # coding transcript
                    out_fh.writelines("\t%s\n" % c)
                    cds_transcipt_hgvs.append(c)
                except hgvs.exceptions.HGVSError:
                    try:
                        n = vm.g_to_n(g, tr)  # noncoding transcript
                        noncoding_transcript_hgvs.append(n)
                        out_fh.writelines("\t %s Noncoding transcript!\n" % n)
                    except hgvs.exceptions.HGVSError:
                        unmappable.append((g, tr))
                        out_fh.writelines("Variant %s can't be mapped to Transcript %s\n" % (g, tr))

        def parse_variant(g_hgvs):
            g = hp.parse_hgvs_variant(g_hgvs)
            out_fh.writelines("%s\n" %g)
            genome_hgvs.append(g)
            map_transcript(g)

        for i, v in enumerate(VCF(self.son)):
            if i < 100:
                refseq_nc_number = make_name_ac_map("GRCh37.p13")[v.CHROM[3:]]
                for alt in v.ALT:
                    # http://varnomen.hgvs.org/
                    if len(v.REF) == 1 and len(alt) == 1:  # substitution
                        string = "%s:g.%s%s>%s" % (refseq_nc_number, str(v.POS), str(v.REF), str(alt))
                        parse_variant(string)

                    elif len(v.REF) == 1 and len(alt) > 1:  # insertion
                        start = str(v.POS)
                        end = str(v.POS + 1)
                        ins = str(alt[1:])
                        string = "%s:g.%s_%sins%s" % (refseq_nc_number, start, end, ins)
                        parse_variant(string)

                    elif len(v.REF) > 1 and len(alt) == 1:  # deletion
                        start = str(v.POS + 1)
                        end = str(v.POS + len(v.REF[1:]))
                        delet = str(v.REF[1:])
                        string = "%s:g.%s_%sdel%s" % (refseq_nc_number, start, end, delet)
                        parse_variant(string)

                    elif len(v.REF) > 1 and len(alt) > 1 and len(alt) == len(v.REF):  # complex
                        start = str(v.POS)
                        end = str(v.POS + len(v.REF[1:]))
                        string = "%s:g.%s_%sdelins%s" % (refseq_nc_number, start, end, alt)
                        parse_variant(string)

                    else:
                        raise Exception("Case not implemented!")
            else:
                break

        out_fh.close()

        print "The first %i Variants contain" % len(genome_hgvs)
        print "\t%i CDS Transcripts" % len(cds_transcipt_hgvs)
        print "\t%i Noncoding Transcripts" % len(noncoding_transcript_hgvs)
        print "\t%i Unmappable Transcripts (NCBI Refseq outdated!)" % len(unmappable)
        print "Transkript mapping results are written to file %s!" % (outfile)

    def print_summary(self):
        print(__author__)
        # self.get_total_number_of_variants_mother()
        # self.get_total_number_of_variants_father()
        # self.get_variants_shared_by_father_and_son()
        # self.get_variants_shared_by_mother_and_son()
        # self.get_variants_shared_by_trio()
        # self.merge_mother_father_son_into_one_vcf()
        self.convert_first_variants_of_son_into_HGVS()


if __name__ == '__main__':
    print("Assignment 3")
    assignment1 = Assignment3()
    assignment1.print_summary()
