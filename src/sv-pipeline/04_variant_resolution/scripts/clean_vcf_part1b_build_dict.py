"""
Remove CNVs that are improperly genotyped by depth because they are nested
within a real CNV
"""

import logging
import pybedtools
import pysam
import sys
import json

SVTYPE = "SVTYPE"
BLANK_SAMPLES = "blanksample"


class SVType:
    DUP = "DUP"
    DEL = "DEL"


class VariantFormatTypes:
    # Predicted copy state
    RD_CN = "RD_CN"
    # Classes of evidence supporting final genotype
    EV = "EV"


class VCFReviser:
    def __init__(self):
        self.rd_cn = {}
        self.ev = {}
        self.sample_indices_dict = {}
        self.sample_list = []

    def _update_rd_cn_ev(self, variant):
        self.rd_cn[variant.id] = [v[VariantFormatTypes.RD_CN] for _, v in variant.samples.items()]
        self.ev[variant.id] = [v[VariantFormatTypes.EV] for _, v in variant.samples.items()]

    @staticmethod
    def get_wider(f):
        if int(f[2]) - int(f[1]) >= int(f[8]) - int(f[7]):
            return f[0:6], f[6:12]
        else:
            return f[6:12], f[0:6]

    @staticmethod
    def get_coverage(wider, narrower):
        n_start = int(narrower[1])
        n_stop = int(narrower[2])
        w_start = int(wider[1])
        w_stop = int(wider[2])

        coverage = 0
        if w_start <= n_stop and n_start <= w_stop:
            intersection_size = min(n_stop, w_stop) - max(n_start, w_start)
            coverage = intersection_size / (n_stop - n_start)
        return coverage

    def get_geno_normal_revise(self, vcf_file, bed_file):
        overlap_test_text = {}
        with pysam.VariantFile(vcf_file, "r") as f:
            header = f.header
            i = -1
            for sample in header.samples:
                i += 1
                self.sample_indices_dict[sample] = i
                self.sample_list.append(sample)

            logging.info('Getting RD/CN/EV')
            for variant in f:
                self._update_rd_cn_ev(variant)

            logging.info("Filtering intersect results")
            bed = pybedtools.BedTool(bed_file)
            for interval in bed.intervals:
                wider, narrower = self.get_wider(interval.fields)
                if wider[5] == BLANK_SAMPLES:
                    continue

                coverage = self.get_coverage(wider, narrower)
                if coverage >= 0.5:
                    wider_samples = set(wider[5].split(","))
                    narrower_samples = set(narrower[5].split(","))
                    non_common_samples = wider_samples - narrower_samples
                    for x in non_common_samples:
                        overlap_test_text[f"{narrower[3]}@{x}"] = \
                            [f"{narrower[3]}@{x}", wider[3], wider[4]]



        logging.info('Generating geno_normal_revise_dict')
        geno_normal_revise_dict = {}
        for k, v in overlap_test_text.items():
            var_id = k.split("@")[0]
            sample_index = self.sample_indices_dict[k.split("@")[1]]
            new_val = None
            if v[2] == SVType.DUP and \
                    self.rd_cn[var_id][sample_index] == 2 and \
                    self.rd_cn[v[1]][sample_index] == 3:
                new_val = 1
            elif v[2] == SVType.DEL and \
                    self.rd_cn[var_id][sample_index] == 2 \
                    and self.rd_cn[v[1]][sample_index] == 1:
                new_val = 3

            if new_val:
                if var_id not in geno_normal_revise_dict:
                    geno_normal_revise_dict[var_id] = {}
                sample_id = self.sample_list[sample_index]
                geno_normal_revise_dict[var_id][sample_id] = new_val

        return geno_normal_revise_dict


def main(args):
    logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)
    logging.info('Starting script')
    reviser = VCFReviser()
    filtered_vcf = args[1]
    intersected_bed = args[2]
    geno_normal_revise_dict = reviser.get_geno_normal_revise(filtered_vcf, intersected_bed)
    logging.info('Dumping dictionary')
    sys.stdout.write(json.dumps(geno_normal_revise_dict))
    logging.info('Done')


if __name__ == '__main__':
    main(sys.argv)
