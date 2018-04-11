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
import regex as re

def prRed(prt): print("\033[91m {}\033[00m" .format(prt))
def prGreen(prt): print("\033[92m {}\033[00m" .format(prt))

def Increase_target(bed,incr_len,opts):
	outpath = opts.out_path+ '/FILES/'+'.'.join((bed.split('/')[-1]).split('.')[:-1]+['Incr'+str(incr_len)]+['bed'])
	out = open( outpath,'w')
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
	return outpath
		
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


def vars_from_db(vars,db,num_snv,num_indel,out):

	database=open(db,'r').readlines()
	vcf=open(out,'w')
	start = [database.index(x) for x in database if x.startswith("#CHROM")][0]
	vcf.write("##fileformat=VCFv4.2\n")
	vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
	simulate=[]
	if vars != None:
		vars_file = open(vars,'r')
		for line in vars_file:
			if line.startswith('#'):
				continue
			else:
				chr,pos,id,ref,alt=line.rstrip().split('\t')[:5]
				vcf.write(line)
				simulate+=[[chr,pos,ref,alt]]
		vars_file.close()

	varianti=database[start:]
	j=0
	while j < int(num_indel):
		rand=r.randrange(len(varianti))
		var=varianti[rand].rstrip().split('\t')
		del(varianti[rand])
		chr=var[0]
		pos=var[1]
		ref=var[3]
		alt=var[4].split(',')
		info=var[7].split(';')
		maf_ref,maf_alt=['.','.']
		
		for i in info:
			if i.startswith('CAF='):
				caf=i.split('=')[1].split('[')[1].split(']')[0].split(',')

				if '.' in caf:
					try:
						del alt[caf.index('.')]
					except:
						pass
					del caf[caf.index('.')]
				try:
					maf_ref,maf_alt=caf
				except:
					pass
		alt=alt[0]
		try:
			if len(ref)>1 and len(alt)==1 or len(alt)>1 and len(ref)==1:
				for schr,spos,sref,salt in simulate:
					len_svar= abs(len(sref)-len(salt))
					len_var= abs(len(ref)-len(alt))
					if schr == chr and int(pos) >= int(spos) and int(pos) <= (int(spos)+int(len_svar)):
						raise Exception()
					elif schr == chr and int(spos) >= int(pos) and int(spos) <= (int(pos)+int(len_var)):
						raise Exception()
				j+=1
				gt=calc_gt(maf_ref,maf_alt)
				vcf.write('\t'.join(var[:4]+[alt,'.','.','.','GT',gt])+'\n')
				simulate+=[[chr,pos,ref,alt]]
			else:
				continue
		except Exception:
			continue


	varianti=database[start:]
	
	j=0		
	while j < int(num_snv):
		rand=r.randrange(len(varianti))
		var=varianti[rand].rstrip().split('\t')
		del(varianti[rand])
		chr=var[0]
		pos=var[1]
		alt=var[4].split(',')
		info=var[7].split(';')
		maf_ref,maf_alt=['.','.']
		
		for i in info:
			if i.startswith('CAF='):
				caf=i.split('=')[1].split('[')[1].split(']')[0].split(',')

				if '.' in caf:
					try:
						del alt[caf.index('.')]
					except:
						pass
					del caf[caf.index('.')]
				try:
					maf_ref,maf_alt=caf
				except:
					pass
		alt=alt[0]
		#print alt,maf_alt
		try:
			if len(var[3])==1 and len(alt)==1:
				for schr,spos,sref,salt in simulate:
					len_svar= abs(len(sref)-len(salt))
					if schr == chr and int(pos) >= int(spos) and int(pos) <= (int(spos)+int(len_svar)):
						raise Exception()
				j+=1
				gt=calc_gt(maf_ref,maf_alt)
				vcf.write('\t'.join(var[:4]+[alt,'.','.','.','GT',gt])+'\n')
			else:	
				continue
		except Exception:
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
						print ':'.join(line.split('\t')[:2]),'--> ' + type +' simulated'
					elif line_to_check.startswith('snv'):
						vaf=str(round(float(line_to_check.split('\t')[-2]),3))
						var=line.split('\t')[:3]
						log.write('\t'.join(var+[vaf,'SIMULATED']) +'\n')
						print ':'.join(line.split('\t')[:2]),'--> SNV simulated'
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
		##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
		##FORMAT=<ID=VAF,Number=1,Type=Float,Description="Variant Allele Frequency">
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
								vcf.write( "\t".join((chrom,pos,'.',ref,alt,'.','PASS', analisi,'GT:VAF',gt+':'+vaf))+'\n')

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

								vcf.write('\t'.join((chrom, start, '.', ref, alt, '.', 'PASS',  analisi,'GT:VAF',gt+':'+vaf))+'\n')
					subprocess.call("mv " + log_bs_dir+'/'+pathname+ '/'+filename + ' ' + log_bs_dir+'/'+pathname+ '/'+filename +'.checked', shell=True)

def Vcf_to_bamsurgeon(vars,min,max,err):
	snpfile = opts.out_path+ '/FILES/'+'.'.join((vars.split('/')[-1]).split('.')[:-1]+['snp'])
	indelfile = opts.out_path+ '/FILES/'+'.'.join((vars.split('/')[-1]).split('.')[:-1]+['indel'])
	svfile = opts.out_path+ '/FILES/'+'.'.join((vars.split('/')[-1]).split('.')[:-1]+['sv'])
	snp=open(snpfile,'w')
	indel=open(indelfile,'w')
	#sv=open(svfile,'w')
	vcf=open(vars,'r')
	for line in vcf:
		line=(line.rstrip()).split('\t')
		if line[0].startswith('#'):
			continue
		else:
			chr,pos,id,ref,alt,qual,filter,info,format,samp_format=line
			# chr=line[0]
			# pos=line[1]
			# id=line[2]
			# ref=line[3]
			# alt=line[4]
			# format=line[-2]
			# samp_format=line[-1]
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
	return snpfile,indelfile

def SV_vcf_to_bamsurgeon(vars):
	svfile = opts.out_path+ '/FILES/'+'.'.join((vars.split('/')[-1]).split('.')[:-1]+['sv'])
	sv=open(svfile,'w')
	vcf=open(vars,'r')
	alufile='/home/jarvis/NGS_TOOLS/bamsurgeon-master/scripts/ins_seqs/alu.fa'

	for line in vcf:
		line=(line.rstrip()).split('\t')
		if line[0].startswith('#'):
			continue
		else:
			chr,pos,id,ref,alt,qual,filter,info,format,samp_format=line
			alt=re.sub('<|>','',alt)
			print line
			try:
				gt=samp_format.split(':')[format.split(':').index('GT')]
				gt=re.sub('\|','/',gt)
			except:
				gt='0/1'
			if gt == '0/1':
				af='0.5'
			elif gt == '1/1':
				af='1.0'
			else:
				continue
			svtype=''
			svend=''
			svlen=''
			ciend=''
			cipos=''
			mstart=''
			mend=''
			for i in info.split(';'):
				if i.startswith('TSD='):
					tsd=i.split('=')[1]
				if i.startswith('SVTYPE='):
					svtype=i.split('=')[1]
				if i.startswith('END='):
					svend=i.split('=')[1]
				if i.startswith('CIEND'):
					ciend=i.split('=')[1].split(',')
				if i.startswith('CIPOS'):
					cipos=i.split('=')[1].split(',')
				if i.startswith('SVLEN'):
					svlen=i.split('=')[1]
			

			print svtype,alt

			if 'ALU' in svtype and 'INS' in alt:
				svend = str(int(pos)+1)
				for i in info.split(';'):
					if i.startswith('MEINFO'):
						aluname,alustart,aluend,polarity=i.split('=')[1].split(',')
				
				seq=find_seq_in_file(aluname,alustart,aluend,alufile,None)

				if tsd != 'null' and seq.startswith(tsd) and seq.endswith(tsd):
					sv.write('\t'.join([chr,pos,svend,'INS',seq,str(len(tsd))]) + '\n')
				else:
					sv.write('\t'.join([chr,pos,svend,'INS',seq]) + '\n')


			elif 'SVA' in svtype and 'INS' in alt:

				svafasta=os.path.isfile(opts.out_path+ '/FILES/'+'.'.join((vars.split('/')[-1]).split('.')[:-1]+['SVA',chr+':'+pos,fa]))
				svend = str(int(pos)+1)
				svaname,svastart,svaend=['SVA',0,'all']
				for i in info.split(';'):
					if i.startswith('MEINFO'):
						svaname,svastart,svaend,polarity=i.split('=')[1].split(',')
				
				seq=find_seq_in_file(svaname,svastart,svaend,svafile,svafasta)

				if tsd != 'null' and seq.startswith(tsd) and seq.endswith(tsd):
					sv.write('\t'.join([chr,pos,svend,'INS',svafasta,str(len(tsd))]) + '\n')
				else:
					sv.write('\t'.join([chr,pos,svend,'INS',svafasta]) + '\n')

			elif 'LINE1' in svtype and 'INS' in alt:

				l1fasta=os.path.isfile(opts.out_path+ '/FILES/'+'.'.join((vars.split('/')[-1]).split('.')[:-1]+['LINE1',chr+':'+pos,fa]))
				svend = str(int(pos)+1)
				l1name,l1start,l1end=['LINE1',0,'all']
				for i in info.split(';'):
					if i.startswith('MEINFO'):
						l1name,l1start,l1end,polarity=i.split('=')[1].split(',')
				
				seq=find_seq_in_file(l1name,l1start,l1end,l1file,l1fasta)

				if tsd != 'null' and seq.startswith(tsd) and seq.endswith(tsd):
					sv.write('\t'.join([chr,pos,svend,'INS',l1fasta,str(len(tsd))]) + '\n')
				else:
					sv.write('\t'.join([chr,pos,svend,'INS',l1fasta]) + '\n')
				

			elif 'INS' in svtype and 'MT' in alt:
				mtfasta=os.path.isfile(opts.out_path+ '/FILES/'+'.'.join((vars.split('/')[-1]).split('.')[:-1]+['MT',chr+':'+pos,fa]))
				svend = str(int(pos)+1)
				mtname,l1start,l1end=['MT',0,'all']
				for i in info.split(';'):
					if i.startswith('MSTART'):
						mtstart=i.split('=')[1]
					if i.startswith('MEND'):
						mtend=i.split('=')[1]
				
				seq=find_seq_in_file(mtname,mtstart,mtend,mtfile,mtfasta)
				sv.write('\t'.join([chr,pos,svend,'INS',mtfasta]) + '\n')

			elif 'INS' in svtype and 'MT' in alt:

			if ciend != '':
				ciend=rm.randrange(int(ciend[0]),int(ciend[1]))
				svend=str(int(svend)+ciend)

			if cipos != '':
				cipos=rm.randrange(int(cipos[0]),int(cipos[1]))
				pos=str(int(pos)+cipos)

			if 'CNV' in svtype:
				alt=re.sub('<|>','',alt)
				alt=alt.split(',')[rm.randrange(1,len(alt.split(',')))]
				if alt == 'INV':
					svtype='INV'
				elif alt == 'CN0':
					svtype='DEL'
				else:
					svtype='DUP'
					
			if 'INV' in svtype:
				sv.write('\t'.join([chr,pos,svend,'INV']) + '\n')

			if 'DEL'in svtype:
				sv.write('\t'.join([chr,pos,svend,'DEL','1.0',af]) + '\n')
			
			if 'DUP' in svtype:
				alt=re.sub('<|>','',alt)
				nalt=str(int(float(af)*2)*(int(re.sub('CN','',alt))-1))
				print nalt
				sv.write('\t'.join([chr,pos,svend,'DUP',nalt]) + '\n')
	sv.close()
	vcf.close()
	return svfile

def find_seq_in_file(name,start,end,seqfile,svfasta):
	seqf=open(seqfile,'r')
	seqline=[]
	for line in seqf:
		line=line.rstrip()
		if line == '>'+name:
			ok=1
			continue
		elif line.startswith('>'):
			ok=0
			continue
		if ok==1:
			seqline+=[line]

	if end == 'all':
		seq= ''.join(seqline)
	else:
		try:
			seq= ''.join(seqline)[int(start):int(end)]
		except:
			seq= ''.join(seqline)[int(start):]
	if svtype =='ALU':
		return seq.upper()
	elif svtype == 'SVA':
		svfasta = open(svfasta,'w')
		svfa.write('>SVA\n')
		svfa.write(seq.upper())
		svfa.close()
		return seq.upper()
	elif svtype == 'LINE1':
		svfasta = open(svfasta,'w')
		svfa.write('>LINE1\n')
		svfa.write(seq.upper())
		svfa.close()
		return seq.upper()
		

def split_neighbors(v_file,distance,out_path):
	vfile=open(v_file,'r')
	left_variants=vfile.readlines()
	all_variants=[x for x in left_variants]
	i=0
	splitted_variants=[]
	while len(left_variants)!=0:
		added=[]
		if i>10:
			quit(1)

		newpath= out_path+ '/'+'.'.join((v_file.split('/')[-1]).split('.')[:-1]+['split'+str(i+1)] +[v_file.split('.')[-1]])
		splitted_variants+= [newpath]
		out=open(splitted_variants[i],'w')
		for line in all_variants:
			var=(line.rstrip()).split('\t')
			if line in left_variants:
				chr=var[0]
				start=var[1]
				stop=var[2]
				add=1
				if len(added) == 0:
					pass
				else:
					for vline in added:
						v=(vline.rstrip()).split('\t')
						vchr=v[0]
						vstart=v[1]
						vstop=v[2]
						if vchr == chr:
							if int(stop) <= int(vstart)-int(distance) and int(start) <= int(vstart)-int(distance) or int(stop) >= int(vstop)+int(distance) and int(start) >= int(vstop)+int(distance):
								pass
							else:
								add=0
								break

				if add == 1:
					added+=[line]
					out.write(line)
					left_variants.remove(line)
				
		i+=1
	return splitted_variants



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
	parser.add_argument('-v', '--version', action='version', version='SimulVar v1.2.0')

	parser.add_argument('-r', '--ref', help="Reference.fasta", required=True)
	parser.add_argument('-o', '--out_path',help="output path", default='~')
	parser.add_argument('-c', '--cfg',help="configuration file", required=True)
	parser.add_argument('-s','--sample', help="Name of sample", default='simulated')

	parser.add_argument('--Generate_fastq',help="It allows fastq generation using: a)bedfile + reference, b) custom fasta file; pIRS (Enrichment design) or ART (Amplicon design) will be used for fastq generation ",action='store_true')
	parser.add_argument('-b', '--bed', help="Bed file used to filter reference using targets",default=None)
	parser.add_argument('-fa', '--fasta', help="Custom fasta file used to generate fastq",default=None)
	parser.add_argument('-d', '--design',choices=['Enrichment','Amplicon'],help="The fastq design: Amplicon or Enrichment",default='Enrichment')

	parser.add_argument('--Germline_simulation',help='''It allows simulation of germline variants starting from: a) paired fastq files, b) bam file. Variants will be taken from a vcf (--vcfvariants) or dbSNP (--dbsnp). 
						Allele frequency can be set: a) using freq thresholds (--minfreq and --maxfreq), b) using VAF in the FORMAT field in vcf file, c) using GT in the FORMAT field in vcf file.''',action='store_true')
	parser.add_argument('--vcfvariants', help="Germline variant list to simulate in vcf format", default=None)
	parser.add_argument('--sv_vcf', help="Germline Structural variation list to simulate in vcf format", default=None)
	parser.add_argument('--dbsnp', help="dbsnp path for variant random extraction",default=None)
	parser.add_argument('--num_snv_dbsnp', help="number of snv to extract from db_snp", type=int, default=None)
	parser.add_argument('--num_indel_dbsnp', help="number of indel to extract from db_snp", type=int, default=None)
	parser.add_argument('--minfreq', help="Min threshold of simulating frequence for Germline variants", type=float, default=None)
	parser.add_argument('--maxfreq', help="Max threshold of simulating frequence for Germline variants", type=float, default=None)

	parser.add_argument('--Somatic_simulation',help='''It allows simulation of somatic variants starting from: a) paired fastq files, b) bam file. Variants will be taken from a vcf (--vcfvariants_som) or COSMIC db (--cosmic). 
						Allele frequency can be set: a) using freq thresholds (--minfreq_som and --maxfreq_som), b) using VAF in the FORMAT field in vcf file.''',action='store_true')
	parser.add_argument('--vcfvariants_som', help="Somatic variant list to simulate in vcf format", default=None)
	parser.add_argument('--cosmic', help="cosmic path for variant random extraction", default=None)
	parser.add_argument('--num_snv_cosmic', help="number of snv to extract from cosmic", type=int, default=None)
	parser.add_argument('--num_indel_cosmic', help="number of indel to extract from cosmic", type=int, default=None)
	parser.add_argument('--minfreq_som', help="Min threshold of simulating frequence for Somatic variants", type=float, default='0.05')
	parser.add_argument('--maxfreq_som', help="Max threshold of simulating frequence for Somatic variants", type=float, default='0.25')

	parser.add_argument('-fq1', '--fastq1', help="Fastq 1 of paired sample",default=None)
	parser.add_argument('-fq2', '--fastq2', help="Fastq 2 of paired sample",default=None)
	parser.add_argument('--bam', help="Bam to simulate the variants",default=None)
	parser.add_argument('--err', help="Simulating frequence error. Variants will be simulated with freq between freq +- err*freq.",type=float, default=None)
	parser.add_argument('--fastq_out', help="Indicate output fastq file of simulated sample after variant spikein", action='store_true')
	
	start_time = datetime.datetime.now()
	
	global opts
	opts = parser.parse_args()

	if opts.cfg != None:
		conf = json.loads((open(opts.cfg).read()).encode('utf8'))
	else:
		conf = json.loads((open(os.path.dirname(os.path.abspath(__file__)) + '/CFG/SimulVar.config.json').read()).encode('utf8'))
 
	tool.makedirs([opts.out_path])

	fastq_path = opts.out_path + '/FASTQ'
	bam_path = opts.out_path + '/BAM/'
	files_path = opts.out_path + '/FILES/'
	log_bs_dir = opts.out_path + '/BAMSURGEON_LOGS/'
	sim_log_dir = opts.out_path + '/SIMULATION_LOGS/'
	log_path = files_path + '/' + str(start_time)
	
	tool.makedirs([fastq_path,bam_path,files_path])
	log = open(log_path+'.log','w')
	print '\nLog file will be generated: '+ log_path+'\n'
	
	## Setting path e arguments for each software #####
	path_art,art_args = set.Set_art(conf,opts)
	path_pirs,pirs_args = set.Set_pirs_args(conf,opts)
	path_bwa,bwa_args = set.Set_bwa(conf,opts)
	path_picard,picard_args = set.Set_picard(conf,opts)
	addsnv,addindel,bs_args = set.Set_bamsurgeon(conf,opts)
	samtools_path = set.Set_samtools(conf,opts)
	bedtools_path = set.Set_bedtools(conf,opts)
	

	if opts.Generate_fastq:
		fasta = opts.fasta
		if opts.bed != None:
			bed=opts.bed
			# Cutting reference.fasta with bed to generate intervals.fasta
			if opts.design == "Enrichment":
				bed = Increase_target(opts.bed,150,opts)
			fasta = tool.getfasta(bed,opts,bedtools_path,log)
			tool.Index_fasta(fasta,samtools_path,log)

		if fasta != None:
			# Generating paired ends fastq from given sequence.fasta
			
			if opts.design == "Amplicon":
				print "- Reads simulation using ART"
				tool.Simulate_fastq_art(path_art,art_args,path_picard,picard_args,fasta,log)

			elif opts.design == "Enrichment":
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


	if opts.Germline_simulation or opts.Somatic_simulation:

		if opts.fastq1 != None or opts.fastq2 != None :
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

		elif opts.bam != None:
			# if bam file is given
			bam = opts.bam
			tool.Index_bam(path_picard,picard_args,bam,log)

	if opts.Germline_simulation:

		vars = None
		if opts.vcfvariants != None:
			vars = opts.vcfvariants
		if opts.dbsnp != None:
			dbsnp=opts.dbsnp
			if opts.bed != None:
				dbsnp=tool.Intersect_DB_bed(opts.bed,opts.dbsnp,bedtools_path,opts,log)
			vars = vars_from_db(vars,dbsnp,opts.num_snv_dbsnp,opts.num_indel_dbsnp,opts.files_path+'/Spikein.Germline.vcf')
		
		if vars == None:
			print "If you want to simulate germline variants you have to indicate a vcf file or dbSNP path."
			sv=SV_vcf_to_bamsurgeon(opts.sv_vcf)
			exit(1)

		tool.makedirs([log_bs_dir,sim_log_dir])
		print "\nStarting Germline spike-in..."
		# if a file containing variant to simulate is given
		#write variant from vcf format to bam surgeon format (bed format)
		snp,indel=Vcf_to_bamsurgeon(vars,opts.minfreq,opts.maxfreq,opts.err)
		sv=SV_vcf_to_bamsurgeon(opts.sv_vcf)
		#if list contains snps
		if os.stat(snp).st_size != 0:
			for snpfile in split_neighbors(snp,150,files_path):
				outbam = bam_path + bam.split('/')[-1].split('.')[0] + '.snp.bam'
				os.chdir(log_bs_dir)
				tool.Bam_surgeon(addsnv,bs_args,bam,outbam,snpfile,log,opts)
				bam = tool.SortSam(path_picard,picard_args,outbam,log)
				status = subprocess.call("rm " + outbam, shell=True)
				status = subprocess.call("rm "+ snpfile, shell=True)
				tool.Index_bam(path_picard,picard_args,bam,log)
			check_simul_vars(log_bs_dir,outbam,snp,sim_log_dir)
			print_vcf('Germline',log_bs_dir,sim_log_dir,opts.ref)			
			
		
		#if list contains indels
		if os.stat(indel).st_size != 0:
			for indelfile in split_neighbors(indel,150,files_path):
				outbam = bam_path + bam.split('/')[-1].split('.')[0] + '.indel.bam'
				os.chdir(log_bs_dir)
				tool.Bam_surgeon(addindel,bs_args,bam,outbam,indelfile,log,opts)
				status = subprocess.call("rm "+ bam , shell=True)
				status = subprocess.call("rm "+ bam + '.bai' , shell=True)
				bam = tool.SortSam(path_picard,picard_args,outbam,log)
				status = subprocess.call("rm "+ outbam, shell=True)
				status = subprocess.call("rm "+ indelfile, shell=True)
				tool.Index_bam(path_picard,picard_args,bam,log)
			check_simul_vars(log_bs_dir,outbam,indel,sim_log_dir)
			print_vcf('Germline',log_bs_dir,sim_log_dir,opts.ref)

		germlinebam = bam_path + outbam.split('/')[-1].split('.')[0] + '.Germline.bam'
		status = subprocess.call("mv " + bam + ' ' + germlinebam , shell=True)
		status = subprocess.call("mv " + bam + '.bai '+ germlinebam + '.bai', shell=True)
		bam = germlinebam
		if opts.fastq_out:
			tool.BamtoFastq(germlinebam,path_picard,picard_args,fastq_path,log)
		print "Germline spike-in: DONE\n"
		print "Sample Germline simulated: " + germlinebam + '\n'


	#if somatic simulation is enabled
	if opts.Somatic_simulation:
		tool.makedirs([log_bs_dir,sim_log_dir])
		germlinebam=bam
		
		vars_som = None
		if opts.vcfvariants_som != None:
			vars_som = opts.vcfvariants_som
		if opts.cosmic != None:
			cosmic=opts.cosmic
			if opts.bed != None:
				cosmic=tool.Intersect_DB_bed(opts.bed,opts.cosmic,bedtools_path,opts,log)
			vars_som = vars_from_db(vars_som,cosmic,opts.num_snv_cosmic,opts.num_indel_cosmic,files_path+'/Spikein.Somatic.vcf')

		if vars_som == None:
			print "If you want to simulate somatic variants you have to indicate a vcf file or COSMICdb path."
			exit(1)

		print "Starting Somatic spike-in..."
		#write variant from vcf format to bam surgeon format (bed format)
		snp,indel = Vcf_to_bamsurgeon(vars_som,opts.minfreq_som,opts.maxfreq_som,opts.err)
		if os.stat(snp).st_size != 0:
			for snpfile in split_neighbors(snp,150,files_path):
				outbam = bam_path + bam.split('/')[-1].split('.')[0] + '.snp.Somatic.bam'
				os.chdir(log_bs_dir)
				tool.Bam_surgeon(addsnv,bs_args,bam,outbam,snpfile,log,opts)

				bam = tool.SortSam(path_picard,picard_args,outbam,log)
				status = subprocess.call("rm " + outbam , shell=True)
				status = subprocess.call("rm "+ snpfile, shell=True)
				tool.Index_bam(path_picard,picard_args,bam,log)
			check_simul_vars(log_bs_dir,outbam,snp,sim_log_dir)
			print_vcf('Somatic',log_bs_dir,sim_log_dir,opts.ref)

		
		if os.stat(indel).st_size != 0:
			for indelfile in split_neighbors(indel,150,files_path):
				outbam = bam_path + bam.split('/')[-1].split('.')[0] + '.indel.Somatic.bam'
				os.chdir(log_bs_dir)
				tool.Bam_surgeon(addindel,bs_args,bam,outbam,indelfile,log,opts)
				status = subprocess.call("rm "+ bam , shell=True)
				status = subprocess.call("rm "+ bam + '.bai' , shell=True)
				bam = tool.SortSam(path_picard,picard_args,outbam,log)
				status = subprocess.call("rm " + outbam, shell=True)
				status = subprocess.call("rm "+ indelfile, shell=True)
				tool.Index_bam(path_picard,picard_args,bam,log)
			check_simul_vars(log_bs_dir,outbam,indel,sim_log_dir)
			print_vcf('Somatic',log_bs_dir,sim_log_dir,opts.ref)

		somaticbam = bam_path + outbam.split('/')[-1].split('.')[0] + '.Somatic.bam'
		status = subprocess.call("mv " + bam + ' ' + somaticbam , shell=True)
		status = subprocess.call("mv " + bam + '.bai ' + somaticbam + '.bai', shell=True)
		tool.BamtoFastq(somaticbam,path_picard,picard_args,fastq_path,log)
		print "Somatic spike-in: DONE\n"
		print "Somatic sample: " + somaticbam
		print "Germline sample: " + germlinebam

	#subprocess.call("rm -rf " + log_bs_dir, shell=True)
	elapsed_time = divmod((datetime.datetime.now() - start_time).total_seconds(),60)

	print "\nTime elapsed: %d min, %d sec" % elapsed_time




