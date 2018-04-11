import subprocess
import os

def prRed(prt): print("\033[91m {}\033[00m" .format(prt))
def prGreen(prt): print("\033[92m {}\033[00m" .format(prt))

def makedirs(dirs):
    for d in dirs:
        if not os.path.exists(d):
			os.makedirs(d)


def Bam_surgeon(path_bs,opts_bs,inbam,outbam,vars,log,opts):
	
	tmpdir='/'.join(outbam.split('/')[:-1]+['tmp'])
	makedirs([tmpdir])
	args = ['python',path_bs,'-v',vars,'-f', inbam,'-r', opts.ref,'-o',outbam,'--tmpdir',tmpdir] + opts_bs
	success = subprocess.call(args,stdout=log,stderr=log)
	try:
		os.rmdir(tmpdir)
	except:
		pass
	
	if success:
		prRed('Error in Variant spikein. Check log file.')
		exit(1)


def Simulate_fastq_pirs(path_pirs,opts_pirs,fasta,log):
	print '- Fastq simulation from: '+fasta
	if '--diploid' in opts_pirs:
		opts_pirs+=[fasta,fasta]
	else:
		opts_pirs+=[fasta]
	#print opts_pirs
	args = [path_pirs,'simulate'] + opts_pirs
	success = subprocess.call(args,stdout=log,stderr=log)
	if not success:
		print 'Fastq simulation: DONE'
	else:
		prRed('Error in Fastq simulation. Check log file.')
		exit(1)


def Simulate_fastq_art(path_art,opts_art,path_picard,opts_picard,fasta,log):
	print '- Fastq simulation from: '+fasta
	
	args = [path_art,'-sam','-i',fasta] + opts_art
	success = subprocess.call(args,stdout=log,stderr=log)

	if "--errfree" in opts_art:
		status = subprocess.call("rm " + opts_art[-1]+'1.fq' , shell=True)
		status = subprocess.call("rm " + opts_art[-1]+'2.fq' , shell=True)
		status = subprocess.call("rm " + opts_art[-1]+'.sam' , shell=True)
		args = ['java',opts_picard[0],'-jar',path_picard,'SamToFastq','I='+opts_art[-1]+'_errFree.sam','FASTQ='+opts_art[-1]+'_R1.fastq','SECOND_END_FASTQ='+opts_art[-1]+'_R2.fastq']
		success = subprocess.call(args,stdout=log,stderr=log)
		status = subprocess.call("rm " + opts_art[-1]+'_errFree.sam' , shell=True)
	else:
		status = subprocess.call("mv " + opts_art[-1]+'1.fq ' + opts_art[-1]+'_R1.fastq', shell=True)
		status = subprocess.call("mv " + opts_art[-1]+'2.fq ' + opts_art[-1]+'_R2.fastq' , shell=True)
		status = subprocess.call("rm " + opts_art[-1]+'.sam' , shell=True)

	if not success:
		print 'Fastq simulation: DONE'
	else:
		prRed('Error in Fastq simulation. Check log file.')
		exit(1)


def Alignment_bwa(path_bwa,opts_bwa,fastq1,fastq2,log,opts):
	
	sam=opts_bwa[-1]+'.sam'
	sam_file=open(sam,'w+')
	if fastq1 == None:
		print '- Alignment using BWA:'+fastq2
		args = [path_bwa,opts_bwa[0],opts.ref,fastq2,opts_bwa[1]]
	elif fastq2 == None:
		print '- Alignment using BWA:'+fastq1
		args = [path_bwa,opts_bwa[0],opts.ref,fastq1,opts_bwa[1]]
	else:
		print '- Alignment using BWA:'+fastq1 + ' '+ fastq2
		args = [path_bwa,opts_bwa[0],opts.ref,fastq1,fastq2,opts_bwa[1]]
	
	success = subprocess.call(args, stdout=sam_file,stderr=log)

	if not success:
		print 'Alignment: DONE ---> Sam: '+sam
		return sam
	else:
		prRed('Error in Alignment. Check log file.')
		exit(1)
	

def Index_bam(path_picard,opts_picard,bam,log):
	bai=bam+'.bai'
	if not os.path.exists(bai):
		args = ['java',opts_picard[0],'-jar',path_picard,'BuildBamIndex','I='+bam,'O='+bai,'VALIDATION_STRINGENCY=LENIENT']
		success = subprocess.call(args,stdout=log,stderr=log)
		if not success:
			pass
		else:
			prRed('Error in Indexing. Check log file.')
			exit(1)


def SamFormatConverter(path_picard,opts_picard,sam,log):

	print '- From Sam to Bam'
	bam='.'.join(sam.split('.')[:-1])+'.bam'
	args = ['java',opts_picard[0],'-jar',path_picard,'SamFormatConverter','I='+sam,'O='+bam]
	success = subprocess.call(args,stdout=log,stderr=log)
	if not success:
		return bam
	else:
		prRed('Error in Conversion. Check log file.')
		exit(1)
	

def SortSam(path_picard,opts_picard,bam,log):

	sort='.'.join(bam.split('.')[:-1])+'.sort.bam'
	args = ['java',opts_picard[0],'-jar',path_picard,'SortSam','I='+bam,'O='+sort,'SORT_ORDER=coordinate']
	success = subprocess.call(args,stdout=log,stderr=log)
	if not success:
		return sort
	else:
		prRed('Error in Sorting. Check log file.')
		exit(1)


def getfasta(bed,opts,bedtools_path,log):
	print '- Getting Fasta from bed'

	outfasta=opts.out_path+ '/FILES/'+'.'.join((opts.ref.split('/')[-1]).split('.')[:-1])+'.bed.fasta'
	args = [bedtools_path, 'getfasta', '-fo', outfasta, '-fi', opts.ref, '-bed' ,bed]
	success = subprocess.call(args,stdout=log,stderr=log)
	
	if not success:
		print 'New Fasta generated from reference and bed: '+ outfasta+'\n'
		return outfasta
	else:
		prRed('Error generating new fasta. Check log file.')
		exit(1)


def Index_fasta(fasta,st_path,log):
	args = [st_path, 'faidx',fasta]
	success = subprocess.call(args,stdout=log,stderr=log)


def Intersect_DB_bed(bed,db,bedtools_path,opts,log):
	outbd=opts.out_path+ '/FILES/'+'.'.join((db.split('/')[-1]).split('.')[:-1])+'.bedfilter.vcf'
	out=open(outbd,'w+')
	args = [bedtools_path, 'intersect','-header', '-a',db,'-b',bed]
	success = subprocess.call(args,stdout=out,stderr=log)
	out.close()
	return outbd


def BamtoFastq(bam,path_picard,opts_picard,out_path,log):

	fastq1=out_path+ '/'+'.'.join((bam.split('/')[-1]).split('.')[:-1])+'_R1.fastq'
	fastq2=out_path+ '/'+'.'.join((bam.split('/')[-1]).split('.')[:-1])+'_R2.fastq'
	args = ['java',opts_picard[0],'-jar',path_picard,'SamToFastq','I='+bam,'FASTQ='+fastq1,'SECOND_END_FASTQ='+fastq2]
	success = subprocess.call(args,stdout=log,stderr=log)
	args = ['gzip',fastq1]
	success = subprocess.call(args,stdout=log,stderr=log)
	args = ['gzip',fastq2]
	success = subprocess.call(args,stdout=log,stderr=log)