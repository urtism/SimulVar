import argparse
import subprocess
import json
import os
import textwrap
import random as rm
import datetime
import pysam
import random as r
import set_tools_config as set
import tools_functions as tool

def prRed(prt): print("\033[91m {}\033[00m" .format(prt))
def prGreen(prt): print("\033[92m {}\033[00m" .format(prt))

def Increase_target(bed,incr_len,opts):
	for line in open(bed,'r'):
		line=line.rstrip()
		if line.startswith('chr'):
			chr=line.split('\t')[0]
			start=int(line.split('\t')[1])-int(incr_len)
			stop=int(line.split('\t')[2])+int(incr_len)
			out.write('\t'.join([chr,str(start),str(stop)]) +'\n')
		else:
			try:
				chr=line.split('\t')[0]
				start=int(line.split('\t')[1])-int(incr_len)
				stop=int(line.split('\t')[2])+int(incr_len)
				out.write('\t'.join([chr,str(start),str(stop)]) +'\n')
			except:
				continue

		


def calc_gt(maf_ref,maf_alt):
	gt='0/1'
	if [maf_ref,maf_alt] == ['.','.']:
		return gt
	else:
		r=rm.randrange(1,10000000)
		p_omo=float(maf_alt)*10000000
		if p_omo >= r:
			gt='1/1'
	return gt

def vars_from_db(db,num_snv,num_indel,out):

	database=open(db,'r').readlines()
	vcf=open(out,'w')
	start = [database.index(x) for x in database if x.startswith("#CHROM")][0]
	varianti=database[start:]
	j=0
	simulate=[['chr1','0']]
	while j < int(num_snv):
		rand=r.randrange(len(varianti))
		var=varianti[rand].rstrip().split('\t')
		del(varianti[rand])
		alt=var[4].split(',')[0]
		chr=var[0]
		pos=var[1]
		info=var[7].split(';')
		maf_ref,maf_alt=['.','.']
		if len(var[3])==1 and len(alt)==1:
			for i in info:
				if i.startswith('CAF='):
					caf=i.split('=')[1].split('[')[1].split(']')[0].split(',')
					if '.' in caf:
						del caf[caf.index('.')]
					try:
						maf_ref,maf_alt=caf
					except:
						print caf
			for v in simulate:
				if chr == v[0] and int(pos) < int(v[1]) + 150 and chr == v[0] and int(pos) > int(v[1]) - 150:
					#print v,var
					pass
				else:
					j+=1
					gt=calc_gt(maf_ref,maf_alt)
					vcf.write('\t'.join(var[:4]+[alt,'.','.','.','GT',gt])+'\n')
					#print var
					simulate+=[[chr,pos]]
					break
		else:	
			continue 
	varianti=database[start:]
	j=0
	while j < int(num_indel):
		rand=r.randrange(len(varianti))
		var=varianti[rand].rstrip().split('\t')
		del(varianti[rand])
		chr=var[0]
		pos=var[1]
		alt=var[4].split(',')[0]
		info=var[7].split(';')
		maf_ref,maf_alt=['.','.']
		if len(var[3])>1 or len(alt)>1:
			for v in simulate:
				for i in info:
					if i.startswith('CAF='):
						caf=i.split('=')[1].split('[')[1].split(']')[0].split(',')
						if '.' in caf:
							del caf[caf.index('.')]
						try:
							maf_ref,maf_alt=caf
						except:
							print caf
				if chr == v[0] and int(pos) < int(v[1]) + 150 and int(pos) > int(v[1]) - 150:
					#print v,var
					pass
				else:
					j+=1
					gt=calc_gt(maf_ref,maf_alt)
					vcf.write('\t'.join(var[:4]+[alt,'.','.','.','GT',gt])+'\n')
					#print var
					simulate+=[[chr,pos]]
					break
		else:
			continue
	vcf.close()
	return out

def check_simul_vars(log_bs_dir,bam,vars,outpath):

	for dir in os.listdir(log_bs_dir):
		if bam.split('/')[-1] in dir:
			log_dir=dir
		else:
			continue
		bam_log=outpath+'/'+bam.split('/')[-1]+'.log'
		varianti= open(vars,'r')
		log= open(bam_log,'w')
		for line in varianti:
			line=line.rstrip()
			id_var='_'.join(line.split('\t')[:3])
			for file in os.listdir(log_bs_dir+'/'+log_dir):
				if id_var in file:
					try:
						line_to_check=(open(log_bs_dir+'/'+log_dir+'/'+file,'r').readlines())[-1]
					except:
						var=line.split('\t')[:3]
						log.write('\t'.join(var+['.','NOT SIMULATED']) +'\n')
						print ':'.join(line.split('\t')[:2]),'--> variant not simulated'
						continue
					if line_to_check.startswith('indel'):
						vaf=str(round(float(line_to_check.split('\t')[-2]),3))
						type=line.split('\t')[4]
						var=line.split('\t')[:3]
						log.write('\t'.join(var+[vaf,type,'SIMULATED']) +'\n')
					elif line_to_check.startswith('snv'):
						vaf=str(round(float(line_to_check.split('\t')[-2]),3))
						var=line.split('\t')[:3]
						log.write('\t'.join(var+[vaf,'SIMULATED']) +'\n')
					else:
						var=line.split('\t')[:3]
						log.write('\t'.join(var+['.','NOT SIMULATED']) +'\n')
						print ':'.join(line.split('\t')[:2]),'--> variant not simulated'
		varianti.close()
		log.close()

def print_header(vcf):
	header=True
	for line in vcf:
		if line.startswith('#'):
			header=False
			break
		else:
			header=True
	if header:
		vcf.write(textwrap.dedent("""\
		##fileformat=VCFv4.1
		##phasing=none
		##INDIVIDUAL=TRUTH
		##SAMPLE=<ID=TRUTH,Individual="TRUTH",Description="bamsurgeon spike-in">
		##INFO=<ID=Somatic,Number=0,Type=Flag,Description="Somatic mutation in primary">
		##INFO=<ID=Germline,Number=0,Type=Flag,Description="Germline mutation in primary">
		##INFO=<ID=VAF,Number=1,Type=Float,Description="Variant Allele Frequency">
		##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
		#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSPIKEIN""")+'\n')

def print_vcf(analisi,log_bs_dir,outpath,reference):
	vcf=open(outpath+'/'+analisi+'.vcf','a+')
	print_header(vcf)

	for pathname in os.listdir(log_bs_dir):
		if os.path.isdir(log_bs_dir+'/'+pathname):
			for filename in os.listdir(log_bs_dir+'/'+pathname):
				#print filename
				if filename.endswith('.log'):
					with open(log_bs_dir +'/'+pathname + '/' + filename, 'r') as infile:
						for line in infile:
							if line.startswith('snv'):
								#chrom, pos, mut = line.strip().split()
								c = line.strip().split()
								chrom = c[1].split(':')[0]
								pos = c[3]
								mut = c[4]
								dpr = c[6]
								vaf = str(round(float(line.split('\t')[-2]),3))
								if float(vaf) >= 0.75:
									gt='1/1'
								elif float(vaf) > 0.00 and float(vaf) < 0.75:
									gt='0/1'

								ref,alt = mut.split('-->')
								try:
									ref=ref.upper()
								except:
									pass
								try:
									alt=alt.upper()
								except:
									pass
								vcf.write( "\t".join((chrom,pos,'.',ref,alt,'.','PASS', analisi+';VAF=' + vaf,'GT',gt))+'\n')

							elif line.startswith('indel'):
								vaf=str(round(float(line.split('\t')[-2]),3))
								if float(vaf) >= 0.75:
									gt='1/1'
								elif float(vaf) > 0.00 and float(vaf) < 0.75:
									gt='0/1'
								fa = pysam.Fastafile(reference)
								indelinfo = line.strip().split()[1].split(':')
								if indelinfo[0] == 'DEL':
									chrom, start, end = indelinfo[1:4]
									ref = fa.fetch(chrom, int(start)-1, int(end)).upper()
									alt = ref[0].upper()

								if indelinfo[0] == 'INS':
									chrom, start, seq = indelinfo[1:4]
									ref = fa.fetch(chrom, int(start)-1, int(start)).upper()
									alt = ref.upper() + seq.upper()
									
								assert ref != '' and alt != '' and start != ''

								vcf.write('\t'.join((chrom, start, '.', ref, alt, '.', 'PASS', analisi+';VAF=' + vaf , 'GT', gt))+'\n')
					subprocess.call("mv " + log_bs_dir+'/'+pathname+ '/'+filename + ' ' + log_bs_dir+'/'+pathname+ '/'+filename +'.checked', shell=True)

def Vcf_to_bamsurgeon(vars,min,max,err):
	print vars
	snp=open('/'.join(vars.split('.')[:-1])+'.snp','w')
	indel=open('/'.join(vars.split('.')[:-1])+'.indel','w')
	vcf=open(vars,'r')
	for line in vcf:
		line=(line.rstrip()).split('\t')
		if line[0].startswith('#'):
			continue
		else:
			chr=line[0]
			pos=line[1]
			id=line[2]
			ref=line[3]
			alt=line[4]
			format=line[-2]
			samp_format=line[-1]
			if min != None and max != None:
				freq=Freq_calc(min,max,err)
			else:
				try:
					freq=samp_format.split(':')[format.split(':').index('VAF')]
				except:
					try:
						gt=samp_format.split(':')[format.split(':').index('GT')]
					except:
						print "ERROR: Set VAF or/and GT in format"
						exit(1)
						
					if gt == '0/1':
						freq=Freq_calc(0.5,0.51,err)
					elif gt == '1/1':
						freq=Freq_calc(1.0,1.1,err)

			if len(ref)==1 and len(alt)==1:
				if alt =='.':
					snp.write('\t'.join([chr,pos,pos,freq]) + '\n')
				else:
					snp.write('\t'.join([chr,pos,pos,freq,alt]) + '\n')
			elif len(ref)>1 and len(alt)==1:
				del_pos=str(int(pos) + len(ref) - len(alt))
				indel.write('\t'.join([chr,pos,del_pos,freq,'DEL']) + '\n')
			elif len(ref)==1 and len(alt)>1:
				INS=alt[1:]
				indel.write('\t'.join([chr,pos,str(int(pos)+1),freq,'INS',INS]) + '\n')
	snp.close()
	indel.close()
	vcf.close()
	return '/'.join(vars.split('.')[:-1])+'.snp','/'.join(vars.split('.')[:-1])+'.indel'

def Freq_calc(min,max,err):
	freq=float(rm.randrange(int(float(min)*100),int(float(max)*100)))/100.0
	if err != None:
		delta=float(err)*freq
		freq_rand=rm.randrange(int((freq-delta)*10000),int((freq+delta)*10000))
		if freq_rand/10000.0 > 1.0:
			freq_rand=10000.0
		freq=freq_rand/10000.0
	return str(freq)


###################################################################################################################################################################################################

	
if __name__ == '__main__':

	parser = argparse.ArgumentParser('A Pipeline to simulate Germline and Somatic variants')
	parser.add_argument('-s','--sample', help="Name of sample",default='simulated')
	parser.add_argument('--bam', help="Bam to simulate the variants",default=None)
	parser.add_argument('-fq1', '--fastq1', help="Fastq 1 of paired sample",default=None)
	parser.add_argument('-fq2', '--fastq2', help="Fastq 2 of paired sample",default=None)
	parser.add_argument('-fa', '--fasta', help="Fasta to generate fastq using pirs or ART",default=None)
	parser.add_argument('-r', '--ref', help="Reference.fasta",default=None)
	parser.add_argument('-b', '--bed', help="Bed file to generate fastq files",default=None)
	parser.add_argument('--vars', help="Germline variant list to simulate in vcf format",default=None)
	parser.add_argument('--vars_som', help="Somatic variant list to simulate in vcf format",default=None)
	parser.add_argument('--dbsnp', help="dbsnp path for variant random extraction",default=None)
	parser.add_argument('--num_snv_dbsnp', help="number of snv to extract from db_snp",default=None)
	parser.add_argument('--num_indel_dbsnp', help="number of indel to extract from db_snp",default=None)
	parser.add_argument('--cosmic', help="cosmic path for variant random extraction",default=None)
	parser.add_argument('--num_snv_cosmic', help="number of snv to extract from cosmic",default=None)
	parser.add_argument('--num_indel_cosmic', help="number of indel to extract from cosmic",default=None)
	parser.add_argument('--minfreq', help="Min threshold of simulating frequence for Germline variants",default=None)
	parser.add_argument('--maxfreq', help="Max threshold of simulating frequence for Germline variants",default=None)
	parser.add_argument('--minfreq_som', help="Min threshold of simulating frequence for Somatic variants",default='0.05')
	parser.add_argument('--maxfreq_som', help="Max threshold of simulating frequence for Somatic variants",default='0.25')
	parser.add_argument('--err', help="Simulating frequence error. Variants will be simulated with freq between freq +- err*freq. [float]",default=None)
	parser.add_argument('-a', '--analysis',choices=['Somatic','Germline'],help="Somatic or Germline",default='Germline')
	parser.add_argument('-o', '--out_path',help="output path")
	parser.add_argument('-c', '--cfg',help="configuration file")
	parser.add_argument('--amplicon',help="Amplicon design for fastq simulation",action='store_true')

	global opts
	opts = parser.parse_args()

	if opts.cfg != None:
		conf = json.loads((open(opts.cfg).read()).encode('utf8'))
	else:
		conf = json.loads((open(os.path.dirname(os.path.abspath(__file__)) + '/CFG/SimulVar.config.json').read()).encode('utf8'))
 

	fastq_path = opts.out_path + '/FASTQ'
	bam_path = opts.out_path + '/BAM/'
	log_bs_dir = opts.out_path + '/BS_log/'
	sim_log_dir = opts.out_path + '/Simulation_LOGS/'
	log_path = opts.out_path + '/' + str(datetime.datetime.now())
	log = open(log_path+'.log','w')
	print '\nLog file will be generated: '+ log_path+'\n'

	tool.makedirs([fastq_path,bam_path])
	
	## Setting path e arguments for each software #####
	path_art,art_args = set.Set_art(conf,opts)
	path_pirs,pirs_args = set.Set_pirs_args(conf,opts)
	path_bwa,bwa_args = set.Set_bwa(conf,opts)
	path_picard,picard_args = set.Set_picard(conf,opts)
	addsnv,addindel,bs_args = set.Set_bamsurgeon(conf,opts)
	samtools_path = set.Set_samtools(conf,opts)
	fasta = opts.fasta
	if opts.bed != None:
		# Cutting reference.fasta with bed to generate intervals.fasta
		fasta = tool.getfasta(opts,log)
		tool.Index_fasta(fasta,samtools_path,log)


	if fasta != None:

		# Generating paired ends fastq from given sequence.fasta
		
		if opts.amplicon:
			print "- Reads simulation using ART"
			tool.Simulate_fastq_art(path_art,art_args,path_picard,picard_args,fasta,log)
		else:
			print "- Reads simulation using Pirs"
			tool.Simulate_fastq_pirs(path_pirs,pirs_args,fasta,log)
		
		for file in os.listdir(opts.out_path + '/FASTQ'):
			if '.sam' in file:
				continue
				#status = subprocess.call("rm " + fastq_path + '/' + file, shell=True)
			if '1.f' in file:
				fastq1 = fastq_path + '/' + file
			elif '2.f' in file:
				fastq2 = fastq_path + '/' + file

		# Alignment of generated fastq in bam file, sorting and indexing
		out_bwa = tool.Alignment_bwa(path_bwa,bwa_args,fastq1,fastq2,log,opts)
		outbam = tool.SamFormatConverter(path_picard,picard_args,out_bwa,log)
		status = subprocess.call("rm " + out_bwa, shell=True)
		bam = tool.SortSam(path_picard,picard_args,outbam,log)
		status = subprocess.call("rm " + outbam , shell=True)
		tool.Index_bam(path_picard,picard_args,bam,log)
		print "Sample simulated: " + bam + '\n'
		
	elif opts.fastq1 != None or opts.fastq2 != None :
		# if fastqs are given
		fastq1 = opts.fastq1
		fastq2 = opts.fastq2
		# Alignment of given fastq in bam file, sorting and indexing
		out_bwa = tool.Alignment_bwa(path_bwa,bwa_args,fastq1,fastq2,log,opts)
		outbam = tool.SamFormatConverter(path_picard,picard_args,out_bwa,log)
		status = subprocess.call("rm "+out_bwa, shell=True)
		bam = tool.SortSam(path_picard,picard_args,outbam,log)
		status = subprocess.call("rm "+outbam , shell=True)
		tool.Index_bam(path_picard,picard_args,bam,log)
		print "Sample simulated: " + bam + '\n'

	elif opts.bam != None:
		# if bam file is given
		bam = opts.bam

	vars = opts.vars
	if opts.dbsnp != None:
		vars = vars_from_db(opts.dbsnp,opts.num_snv_dbsnp,opts.num_indel_dbsnp,opts.out_path+'/Germline_dbsnp.vcf')
	if vars != None:
		tool.makedirs([log_bs_dir,sim_log_dir])
		print "\nStarting simulation of Germline variants:"
		# if a file containing variant to simulate is given
		#write variant from vcf format to bam surgeon format (bed format)
		snp,indel=Vcf_to_bamsurgeon(vars,opts.minfreq,opts.maxfreq,opts.err)
		#if list contains snps
		if os.stat(snp).st_size != 0:
			outbam = bam_path + bam.split('/')[-1].split('.')[0] + '.snp.bam'
			os.chdir(log_bs_dir)
			tool.Bam_surgeon(addsnv,bs_args,bam,outbam,snp,log,opts)
			check_simul_vars(log_bs_dir,outbam,snp,sim_log_dir)
			print_vcf('Germline',log_bs_dir,sim_log_dir,opts.ref)
			#status = subprocess.call("rm "+bam , shell=True)
			#status = subprocess.call("rm "+bam + '.bai' , shell=True)
			bam = tool.SortSam(path_picard,picard_args,outbam,log)
			status = subprocess.call("rm " + outbam, shell=True)
			tool.Index_bam(path_picard,picard_args,bam,log)
		#if list contains indels
		if os.stat(indel).st_size != 0:

			outbam = bam_path + bam.split('/')[-1].split('.')[0] + '.indel.bam'
			os.chdir(log_bs_dir)
			tool.Bam_surgeon(addindel,bs_args,bam,outbam,indel,log,opts)
			check_simul_vars(log_bs_dir,outbam,indel,sim_log_dir)
			print_vcf('Germline',log_bs_dir,sim_log_dir,opts.ref)
			status = subprocess.call("rm "+ bam , shell=True)
			status = subprocess.call("rm "+ bam + '.bai' , shell=True)
			bam = tool.SortSam(path_picard,picard_args,outbam,log)
			status = subprocess.call("rm "+ outbam, shell=True)
			tool.Index_bam(path_picard,picard_args,bam,log)

		germlinebam = bam_path + bam.split('/')[-1].split('.')[0] + '.Germline.bam'
		status = subprocess.call("mv " + bam + ' ' + germlinebam , shell=True)
		status = subprocess.call("mv " + bam + '.bai '+ germlinebam + '.bai', shell=True)
		bam = germlinebam
		print "Sample Germline simulated: " + germlinebam + '\n'

	#if analysis is somatic
	if opts.analysis == 'Somatic':
		tool.makedirs([log_bs_dir,sim_log_dir])
		vars_som = opts.vars_som
		if opts.cosmic != None:
			vars_som = vars_from_db(opts.cosmic,opts.num_snv_cosmic,opts.num_indel_cosmic,opts.out_path+'/Somatic_cosmic.vcf')
		print "Starting simulation of Somatic variants:"
		#write variant from vcf format to bam surgeon format (bed format)
		snp,indel = Vcf_to_bamsurgeon(vars_som,opts.minfreq_som,opts.maxfreq_som,opts.err)
		somaticbam = bam.split('.')[0] + '.Somatic.bam'
		if os.stat(snp).st_size != 0:
			outbam = bam_path + bam.split('/')[-1].split('.')[0] + '.snp.Somatic.bam'
			os.chdir(log_bs_dir)
			tool.Bam_surgeon(addsnv,bs_args,bam,outbam,snp,log,opts)
			check_simul_vars(log_bs_dir,outbam,snp,sim_log_dir)
			print_vcf('Somatic',log_bs_dir,sim_log_dir,opts.ref)
			bam = tool.SortSam(path_picard,picard_args,outbam,log)
			status = subprocess.call("rm " + outbam , shell=True)
			tool.Index_bam(path_picard,picard_args,bam,log)
		
		if os.stat(indel).st_size != 0:
			outbam = bam_path + bam.split('/')[-1].split('.')[0] + '.indel.Somatic.bam'
			os.chdir(log_bs_dir)
			tool.Bam_surgeon(addindel,bs_args,bam,outbam,indel,log,opts)
			check_simul_vars(log_bs_dir,outbam,indel,sim_log_dir)
			print_vcf('Somatic',log_bs_dir,sim_log_dir,opts.ref)
			if '.sort' in bam:
				status = subprocess.call("rm " + bam , shell=True)
				status = subprocess.call("rm " + bam + '.bai' , shell=True)
			bam = tool.SortSam(path_picard,picard_args,outbam,log)
			status = subprocess.call("rm " + outbam, shell=True)
			tool.Index_bam(path_picard,picard_args,bam,log)

		status = subprocess.call("mv " + bam + ' ' + somaticbam , shell=True)
		status = subprocess.call("mv " + bam + '.bai ' + somaticbam + '.bai', shell=True)
		print "Somatic sample: " + somaticbam
		print "Control sample: " + germlinebam

	#subprocess.call("rm -r " + log_bs_dir, shell=True)





