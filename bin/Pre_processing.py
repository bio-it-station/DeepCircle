import os    
import argparse

def mainf():
    ### store every boundary of each chromosome into a dictionary
    output_dir_DNABERT = './DNABERT/eccdna'
    bound = {}
    with open(f'{output_dir_DNABERT}/genome/human/{args.boundary}', 'r') as input_bound:
        for line in input_bound:
            bound[line.split('\t')[0]] = int(line.split('\t')[1])

    output_dir_CNN = './output'
    if not os.path.exists(output_dir_CNN):
        os.mkdir(output_dir_CNN)

    output_dir_CNN += '/{}'.format(args.datatype)
    if not os.path.exists(output_dir_CNN):
        os.mkdir(output_dir_CNN)

    with open(f'{output_dir_DNABERT}/genome/human/{args.gap}', 'r') as gap_f:
        with open(f'{output_dir_DNABERT}/db/human/{args.eccdna}', 'r') as input_f:
            with open(f'{output_dir_DNABERT}/db/human/{args.datatype}_circleseq_eccdna_filt_uniq_seq_{args.extend * 2}.bed', 'w') as output_f:
                with open(f'{output_dir_DNABERT}/db/human/{args.datatype}_circleseq_eccdna_filt_uniq_excl_{args.extend * 2}.bed', 'w') as output_ex:
                    ### write genome gap into exclude .bed file (for generating negative label .bed)
                    for line in gap_f:
                        name, start, end = line.split('\t')
                        output_ex.writelines('%s\t%s\t%s' % (name, start, end))
                    
                    ### select sequence shorter than constraint you set, and extend bidirectionally from their center 
                    for line in input_f:
                        name, start, end = line.split('\t')
                        mid = (int(end) + int(start)) // 2
                        if int(end) - int(start) <= args.limit:
                            ### if sequence extend out of boundary, just put this sequence into excl .bed
                            if mid - args.extend < 0 or mid + args.extend > bound[name]:
                                output_ex.writelines('%s\t%s\t%s' % (name, start, end))
                            else:
                                output_f.writelines('%s\t%d\t%d\n' % (name, mid-args.extend, mid+args.extend))
                                output_ex.writelines('%s\t%d\t%d\n' % (name, mid-args.extend, mid+args.extend))
                        else:
                            output_ex.writelines('%s\t%s\t%s' % (name, start, end))

    ### Generate negative label .bed file
    os.system(f"bedtools shuffle \
        -i {output_dir_DNABERT}/db/human/{args.datatype}_circleseq_eccdna_filt_uniq_seq_{args.extend * 2}.bed \
        -g {output_dir_DNABERT}/genome/human/{args.boundary} \
        -excl {output_dir_DNABERT}/db/human/{args.datatype}_circleseq_eccdna_filt_uniq_seq_{args.extend * 2}.bed \
        > {output_dir_DNABERT}/db/human/{args.datatype}_circleseq_eccdna_filt_uniq_comp_{args.extend * 2}.bed"
    )    
    print("Finish generating negative label data")
    ### Check if negative label .bed is intersect with positive label .bed
    os.system(f"bedtools intersect \
        -a {output_dir_DNABERT}/db/human/{args.datatype}_circleseq_eccdna_filt_uniq_seq_{args.extend * 2}.bed \
        -b {output_dir_DNABERT}/db/human/{args.datatype}_circleseq_eccdna_filt_uniq_comp_{args.extend * 2}.bed"
    )

    with open(f"{output_dir_DNABERT}/db/human/{args.datatype}_circleseq_eccdna_filt_uniq_seq_{args.extend * 2}.bed", 'r') as DNABERT_positive:
        with open(f"{output_dir_DNABERT}/db/human/{args.datatype}_circleseq_eccdna_filt_uniq_comp_{args.extend * 2}.bed", 'r') as DNABERT_negative:
            with open(f"{output_dir_CNN}/{args.datatype}_circleseq_eccdna_filt_uniq_{args.window}_window.bed", 'w') as CNN_positive:
                with open(f"{output_dir_CNN}/rand_{args.datatype}_circleseq_eccdna_filt_uniq_{args.window}_window.bed", 'w') as CNN_negative:
                    difference = args.extend - args.window // 2
                    for line in DNABERT_positive:
                        name, start, end = line.split('\t')
                        CNN_positive.writelines('%s\t%d\t%d\n' % (name, int(start)+difference, int(end)-difference))
                    for line in DNABERT_negative:
                        name, start, end = line.split('\t')
                        CNN_negative.writelines('%s\t%d\t%d\n' % (name, int(start)+difference, int(end)-difference))

    ### Acquire its corresponding FASTA file (DNABERT)
    os.system(f"bedtools getfasta -fi {output_dir_DNABERT}/genome/human/{args.genome} -bed {output_dir_DNABERT}/db/human/{args.datatype}_circleseq_eccdna_filt_uniq_seq_{args.extend * 2}.bed -fo {output_dir_DNABERT}/output/human/human_{args.datatype}_positive_{args.extend * 2}_limit.fa.out")
    os.system(f"bedtools getfasta -fi {output_dir_DNABERT}/genome/human/{args.genome} -bed {output_dir_DNABERT}/db/human/{args.datatype}_circleseq_eccdna_filt_uniq_comp_{args.extend * 2}.bed -fo {output_dir_DNABERT}/output/human/human_{args.datatype}_negative_{args.extend * 2}_limit.fa.out")

    ### Acquire its corresponding FASTA file (CNN)
    os.system(f"bedtools getfasta -fi {output_dir_DNABERT}/genome/human/{args.genome} -bed {output_dir_CNN}/{args.datatype}_circleseq_eccdna_filt_uniq_{args.window}_window.bed -fo {output_dir_CNN}/{args.datatype}_circleseq_eccdna_filt_uniq_{args.window}_window.fa")
    os.system(f"bedtools getfasta -fi {output_dir_DNABERT}/genome/human/{args.genome} -bed {output_dir_CNN}/rand_{args.datatype}_circleseq_eccdna_filt_uniq_{args.window}_window.bed -fo {output_dir_CNN}/rand_{args.datatype}_circleseq_eccdna_filt_uniq_{args.window}_window.fa")

    print("Finish converting bed files to fasta files")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--extend", type=int, default = 512)                    # sequence length you want to extend from center
    parser.add_argument("--window", type=int, default = 1000)                   # window size of sequence in CNN
    parser.add_argument("--limit", type=int, default = 1000)                    # sequence limit
    parser.add_argument("--datatype", type=str)                                 # species name of eccdna
    parser.add_argument("--boundary", type=str)                                 # genome boundary file path
    parser.add_argument("--gap", type=str)                                      # gap bedfile path
    parser.add_argument("--genome", type=str)                                   # genome reference file path
    parser.add_argument("--eccdna", type=str)                                   # eccdna bedfile path
    args = parser.parse_args()
    mainf()   