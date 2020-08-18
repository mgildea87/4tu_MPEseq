#!/usr/bin/env python2.7
import argparse, subprocess, os, gzip, logging, time, collections, HTSeq, pysam, sys
from itertools import izip_longest
from collections import defaultdict

logging.basicConfig(filename='spliceSeq.log', filemode='w', format='%(asctime)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=logging.INFO)

def main(args):

	logging.info(' '.join(sys.argv))
	global_start_time = time.time()
	samples = []

	logging.info('Reading intron ranges from %s' % (args.intron_intervals))
	targets, intron_set, fiveSS, threeSS, Branches, Branchto3ss = readTargetFeatures(args.intron_intervals, args.Branch_windows, args.Branch_to3ss)

	for n, i in enumerate(args.input_files):

		sample_start_time = time.time()

		infile = i.split('/')[-1]
		logging.info(' Processing %s' % (infile))
		
		samples.append(infile)

		# Create folder for each sample
		if not os.path.exists(infile):
			os.makedirs(infile)
		
		# Trim adapter sequences
		if not args.skipTrim:
			logging.info('\tTrimming Adapter Sequences')
			trimAdapters(i, infile)
		else:
			logging.info('\tSkipping Trimming of Adapter Sequences')

		# Alignment
		if not args.skipAlignment:
			logging.info('\tAligning Reads to Genome')
			align(infile, args.STARIndex, args.read_1_length)
		else:
			logging.info('\tSkipping Alignment of Reads to Genome')

		# Pool Reads
		if not args.skipPool:
			logging.info('\tPooling Reads')
			pool(infile, targets, intron_set, fiveSS, threeSS, Branches, Branchto3ss)
		else:
			logging.info('\tSkipping Pooling of Reads')

	#combine count files from each sample
	if not args.skipCombine:
		logging.info('\tCombining count files')
		suffix = '_Splicing_counts.txt'
		sumCol = False
		unspliced = False
		spliced = False
		SI = False
		branch = False
		
		output = 'Combined_total_counts'
		sumCol = True
		combine(suffix, samples, output, sumCol, unspliced, spliced, SI, branch)
		
		sumCol = False
		unspliced = True
		output = 'Combined_unspliced_counts'
		combine(suffix, samples, output, sumCol, unspliced, spliced, SI, branch)

		unspliced = False
		SI = True
		output = 'Combined_SI'
		combine(suffix, samples, output, sumCol, unspliced, spliced, SI, branch)	

		SI = False
		spliced = True
		output = 'Combined_spliced_counts'
		combine(suffix, samples, output, sumCol, unspliced, spliced, SI, branch)

		suffix = '_Concordant_splicing_counts.txt'
		spliced = True
		unspliced = False
		output = 'Combined_concordant_spliced_counts'
		combine(suffix, samples, output, sumCol, unspliced, spliced, SI, branch)
		
		suffix = '_Concordant_splicing_counts.txt'
		spliced = False
		unspliced = True
		output = 'Combined_concordant_unspliced_counts'
		combine(suffix, samples, output, sumCol, unspliced, spliced, SI, branch)

		suffix = '_Branch_to3ss_counts.txt'
		unspliced = False
		branch = True
		output = 'Combined_branch_to3ss_counts'
		combine(suffix, samples, output, sumCol, unspliced, spliced, SI, branch)

		suffix = '_Branched_counts.txt'
		output = 'Combined_branched_counts'
		combine(suffix, samples, output, sumCol, unspliced, spliced, SI, branch)

	else:
		logging.info('\tSkipping count file combine')
	
	logging.info('Total Time: %ds' % (time.time()-global_start_time))

def readTargetFeatures(interval, Branch_windows, Branch_to3ss):
	intron_set = set()
	fiveSS = {}
	threeSS = {}
	targets = HTSeq.GenomicArrayOfSets('auto', stranded=True)
#strands are switched here. MPE-seq reads are on the opposite strand
	for line in open(interval):
		fields = line.rstrip().split('\t')
		if fields[5] == '+':
			fields[5] = '-'
		else:
			fields[5] = '+'
#2 is subtracted and 1 is added such that a read must go atleast 3bp into the intron to be called unspliced		
		if fields[5] == '-':
			iv = HTSeq.GenomicInterval(fields[0], int(fields[1])-1, int(fields[2])-2, fields[5])
		else:
			iv = HTSeq.GenomicInterval(fields[0], int(fields[1])+1, int(fields[2]), fields[5])
		targets[iv] += fields[3]
		intron_set.add(fields[3])
#5'SS and 3'SS are swithed here because MPE-seq reads are on the opposite strand and the orientation was switched above 
		if fields[5] == '+':
			fiveSS[(fields[0], int(fields[2]), fields[5])] = tuple(fields[3].split(';'))
			threeSS[(fields[0], int(fields[1]), fields[5])] = tuple(fields[3].split(';'))
		else:
			fiveSS[(fields[0], int(fields[1]), fields[5])] = tuple(fields[3].split(';'))
			threeSS[(fields[0], int(fields[2]), fields[5])] = tuple(fields[3].split(';'))

	Branches = HTSeq.GenomicArrayOfSets("auto", stranded=False)
	for line in open(Branch_windows):
		fields = line.rstrip().split('\t')
		iv = HTSeq.GenomicInterval(fields[1], int(fields[2]), int(fields[3]))
		Branches[iv] += fields[0]

	
	Branchto3ss = HTSeq.GenomicArrayOfSets("auto", stranded=False)
	for line in open(Branch_to3ss):
		fields = line.rstrip().split('\t')
		iv = HTSeq.GenomicInterval(fields[1], int(fields[2]), int(fields[3]))
		Branchto3ss[iv] += fields[0]


	return targets, intron_set, fiveSS, threeSS, Branches, Branchto3ss

def trimAdapters(full_path, infile):
	cmd = 'fastp -i %s_R1.fastq.gz -I %s_R2.fastq.gz -o %s/%s_Trimmed_R1_P.fastq.gz -O %s/%s_Trimmed_R2_P.fastq.gz --disable_trim_poly_g --adapter_sequence CTGTCTCTTATACACATCT --adapter_sequence_r2 CTGTCTCTTATACACATCT --html %s/%s.html --json %s/%s.json --thread 2' % (full_path, full_path, infile, infile, infile, infile, infile, infile, infile, infile)
	proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	(out, err) = proc.communicate()


def align(infile, index, read_1_length):
	cmd = 'STAR --genomeDir %s --readFilesIn %s/%s_Trimmed_R1_P.fastq.gz %s/%s_Trimmed_R2_P.fastq.gz --readFilesCommand zcat --clip5pNbases 7 0 --peOverlapNbasesMin 5 --peOverlapMMp 0.1 --outFilterMultimapNmax 1 --alignIntronMin 10 --alignIntronMax 1100 --outSAMattributes All --runThreadN 8 --outSAMunmapped Within KeepPairs --alignSJoverhangMin 3 --alignSplicedMateMapLmin 3 --alignMatesGapMax 3000 --outFilterMismatchNmax 999 --alignEndsType EndToEnd --outFileNamePrefix %s/%s_ --outSAMtype BAM Unsorted' % (index, infile, infile, infile, infile, infile, infile)
	proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	(out, err) = proc.communicate()
	with open('%s/%s_alignment.log' % (infile, infile), 'w') as alignOut:
		alignOut.write(err)
	
#removes read pairs with insert sizes <33bp
	bam = pysam.AlignmentFile('%s/%s_Aligned.out.bam' % (infile, infile) , 'rb')
	newbam = pysam.AlignmentFile('%s/%s_temp.bam' % (infile, infile), "wb", template=bam)
	newbam2 = pysam.AlignmentFile('%s/%s_tooshort.bam' % (infile, infile), "wb", template=bam)	
	reads = set()
	short_reads = set()
	for read in bam:
		if read.is_read1 == True and len(read.query_sequence) > read_1_length:
			newbam.write(read)
			reads.add(read.query_name)
		elif read.is_read1 == True:
			newbam2.write(read)		
			short_reads.add(read.query_name)
	logging.info('Number of read 1 alignments longer than 33nt = %s' % (len(reads)))
	bam.close()
	bam = pysam.AlignmentFile('%s/%s_Aligned.out.bam' % (infile, infile) , 'rb')
	for read in bam:
		if read.is_read2 == True and read.query_name in reads:
			newbam.write(read)
		elif read.is_read2 == True and read.query_name in short_reads:
			newbam2.write(read)
	bam.close()
	newbam.close()
	newbam2.close()
	cmd = 'samtools sort -n %s/%s_temp.bam  -o %s/%s.bam' % (infile, infile, infile, infile)
	proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	(out, err) = proc.communicate()
	os.remove('%s/%s_temp.bam' % (infile, infile))
	os.remove('%s/%s_Aligned.out.bam' % (infile, infile))

def pool(infile, targets, intron_set, fiveSS, threeSS, Branches, Branchto3ss):
	
	SI_counts = defaultdict(int)
	junction_counts = defaultdict(int)
	
	for f, s in HTSeq.pair_SAM_alignments_with_buffer(HTSeq.BAM_Reader('%s/%s.bam' % (infile, infile))):

		if f != None and f.aligned == True and f.aQual > 5:
			chrome = f.iv.chrom
			start = f.iv.start
			end = f.iv.end
			strand = f.iv.strand
			if strand == '+':
				geneint = HTSeq.GenomicPosition(chrome, start, strand)
			else:	
				geneint = HTSeq.GenomicPosition(chrome, end, strand)
			if len(targets[geneint]) == 0: 
				introns = set()
				junctions = set()

				for i, cigop in enumerate(f.cigar):
					if cigop.type == 'M':
						for iv, val in targets[cigop.ref_iv].steps():
							introns |= val

					elif cigop.type == 'N':
						if f.cigar[i-1].type =='M' and f.cigar[i-1].size > 3 and f.cigar[i+1].type =='M' and f.cigar[i+1].size > 3:
							for iv, val in targets[cigop.ref_iv].steps():
								junctions |= val

							chrom = cigop.ref_iv.chrom
							if cigop.ref_iv.strand == '+':
								first = cigop.ref_iv.end
								second = cigop.ref_iv.start + 1
								strand  = "+"
							else:
								first = cigop.ref_iv.start + 1
								second = cigop.ref_iv.end
								strand = '-'

							if (chrom, first, strand) in fiveSS and (chrom, second, strand) in threeSS:
								up = fiveSS[chrom, first, strand]
								down = threeSS[chrom, second, strand]
								if up[0] == down[0]:
									if up[1] == down[1]:
										junction_counts[(infile, up[0], int(up[1]), int(down[1])+1, "Constituitive")] += 1
									else:
										junction_counts[(infile, up[0], int(up[1]), int(down[1])+1, "Exon Skipping")] += 1
							elif (chrom, first, strand) in fiveSS:
								junction_counts[(infile, up[0], int(up[1]), int(down[1])+1, "Alternative 3'")] += 1
							elif (chrom, second, strand) in threeSS:
								junction_counts[(infile, up[0], int(up[1]), int(down[1])+1, "Alternative 5'")] += 1

				
				intron_num_mat = {}
				intron_num_pre = {}
				intron = ''
				junction = ''

				if len(introns) > 0:
					for i in introns:
						a = i.split(';')
						intron_num_pre[i] = a[1]
					intron = max(intron_num_pre.items(), key=lambda x: x[1])
					intron = intron[0]
				
				if len(junctions) > 0:
					for i in junctions:
						a = i.split(';')
						intron_num_mat[i] = a[1]
					junction = max(intron_num_mat.items(), key=lambda x: x[1])
					junction = junction[0]
				
				if junction == intron:
					intron = ''
					junction = ''

				if junction and intron:
					if junction.split(';')[1] > intron.split(';')[1]:
						intron = ''
					else:
						junction = ''

				candidate_genes = set()
				for i in introns:
					candidate_genes.add(i.split(';')[0])
				for i in junctions:
					candidate_genes.add(i.split(';')[0])

				if len(candidate_genes) == 1:
					if junction:
						SI_counts[('mature', junction)] += 1
					if intron:
						SI_counts[('premature', intron)] +=1
					if f.proper_pair == True and s.proper_pair == True and s.aligned == True and s.aQual > 5:
						if junction:
							SI_counts[('concordant_mature', junction)] += 1
						if intron:
							SI_counts[('concordant_premature', intron)] +=1

#read 2 counting
				if intron > 0 and s.aligned == True and s.proper_pair == True and s.aQual > 5:		
					chrome = s.iv.chrom
					start = s.iv.start
					end = s.iv.end
					strand = s.iv.strand
					if strand == '+':
						geneint = HTSeq.GenomicPosition(chrome, start, strand)
					else:	
						geneint = HTSeq.GenomicPosition(chrome, end, strand)
					if intron in Branches[ geneint ] and len(Branches[ geneint ]) == 1:
						SI_counts[('Branched', intron)] += 1
					if intron in Branchto3ss[ geneint ] and len(Branchto3ss[ geneint ]) == 1:
						SI_counts[('Branch_to3ss', intron)] += 1

	with open('%s/%s_Splicing_counts.txt' % (infile, infile), 'w') as out:
		for intron in sorted(intron_set):
			out.write('%s\t%d\t%d\n' % (intron, SI_counts[('mature', intron)], SI_counts[('premature', intron)]))

	with open('%s/%s_Concordant_splicing_counts.txt' % (infile, infile), 'w') as out:
		for intron in sorted(intron_set):
			out.write('%s\t%d\t%d\n' % (intron, SI_counts[('concordant_mature', intron)], SI_counts[('concordant_premature', intron)]))

	with open('%s/%s_Branched_counts.txt' % (infile, infile), 'w') as out:
		for intron in sorted(intron_set):
			out.write('%s\t%d\n' % (intron, SI_counts[('Branched', intron)]))

	with open('%s/%s_Branch_to3ss_counts.txt' % (infile, infile), 'w') as out:
		for intron in sorted(intron_set):
			out.write('%s\t%d\n' % (intron, SI_counts[('Branch_to3ss', intron)]))			

	with open('%s/%s_junctionCounts.txt' % (infile, infile), 'w') as out:
		out.write('Gene\tUpstream\tDownstream\tType\tCount\n')
		for junc in sorted(junction_counts):
			out.write('%s\t%d\t%d\t%s\t%d\n' % (junc[1], junc[2], junc[3], junc[4], junction_counts[junc]))


def combine (suffix, prefix_list, output_prefix, sumCol, unspliced, spliced, SI, branch):
	repos = {}
	gene = 'gene'
	length = len(prefix_list)
	a = 0
	for i in prefix_list:
		if gene not in repos:
			repos[gene] = [i]
		else:	
			repos[gene].append(i)
		for line in open('%s/%s%s' % (i,i, suffix)):
			cur = line.rstrip().split('\t')
			if len(cur) > 1:
				if sumCol == True:	
					if cur[0] not in repos:
						repos[cur[0]] = [0] * length
					repos[cur[0]][a] = repos[cur[0]][a] + float(cur[2]) + float(cur[1])
		 		elif SI == True:
		 			if cur[0] not in repos:
						repos[cur[0]] = [0] * length
					if float(cur[2]) == 0 or float(cur[1]) == 0:
						repos[cur[0]][a] = 0
					else:
						repos[cur[0]][a] = repos[cur[0]][a] + float(cur[2]) / float(cur[1])
		 		elif unspliced == True:
					if cur[0] not in repos:
						repos[cur[0]] = [0] * length
					repos[cur[0]][a] = repos[cur[0]][a] + float(cur[2])
		 		elif spliced == True or branch == True:
					if cur[0] not in repos:
						repos[cur[0]] = [0] * length
					repos[cur[0]][a] = repos[cur[0]][a] + float(cur[1])								
		a +=1	
	with open('%s.txt' % (output_prefix), 'w') as Out:
		for i in sorted(repos.keys(), reverse=True):
			output = ''
			for col in range(len(repos[i])):
				output = output + '\t' + str(repos[i][col])
			Out.write('%s%s\n' % (i, output.rstrip()))

def parseArguments():
	
	parser = argparse.ArgumentParser(prog="MPE_PipeLine_PE", description='', usage='%(prog)s [options]')
	required = parser.add_argument_group('required arguments')
	trim = parser.add_argument_group('Trimming Options')
	alignment = parser.add_argument_group('Alignment Options')
	pool = parser.add_argument_group('Pooling Options')
	combine = parser.add_argument_group('File Combine Options')
	
	required.add_argument('-i', '--input_files', nargs='+', required=True, help=' Basename of files to run. (fastq.gz)', metavar='', dest='input_files')
	trim.add_argument('--skip_Trim', action='store_true', help=' Skip the trimming of adapter sequences. Assumes files already exist.', dest='skipTrim')
	alignment.add_argument('--skip_Align', action='store_true', help=' Skip the alignment of reads to genome. Assumes files already exist.', dest='skipAlignment')
	alignment.add_argument('--STAR_index', default='/home/mag456/genomes/concat_cere_pombe/STAR_genome_2_7/', help=' Location of STAR index files.', dest='STARIndex')
	alignment.add_argument('--Read_1_length_min', type = int, default=40, help='Minimum alignment length of read 1. The purpose of which is to remove unextended primer reads', dest='read_1_length')
	pool.add_argument('--skip_Pool', action='store_true', help=' Skip the pooling of reads. Assumes files already exist.', dest='skipPool')
	pool.add_argument('--intron_intervals', default='/home/mag456/genomes/concat_cere_pombe/concat_intron_ranges.bed', help=' File containing intervals for which to count reads.', dest='intron_intervals')
	pool.add_argument('--Branch_windows', default='/home/mag456/genomes/cere/cere_annotations/Branch_window.txt', help=' File containing branch intervals for which to count reads.', dest='Branch_windows')
	pool.add_argument('--Branch_to3ss', default='/home/mag456/genomes/cere/cere_annotations/3SSto6bpDownStreamofBranch.txt', help=' File containing branch intervals for which to count reads.', dest='Branch_to3ss')
	pool.add_argument('--junction_anchor', default=1, type=int, help=' Required length of anchor on each side of junction.', dest='junction_anchor')
	combine.add_argument('--skip_Combine', action='store_true', help='Skip combining count files', dest='skipCombine')
	
	return parser.parse_args()

args = parseArguments()
main(args)
