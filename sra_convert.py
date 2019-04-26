#!/usr/bin/env python

from argparse import ArgumentParser
import os
import pypiper
import sys


CONVERT = "convert"
DELETE_SRA = "delete_sra"
DELETE_BAM = "delete_bam"
REMOVAL_OPTNAME = "--remove"


def _parse_cmdl(cmdl):
	parser = ArgumentParser(description="Automatic GEO SRA data downloader")
	
	parser.add_argument(
			"-b", "--bamfolder", 
			default=safe_echo("SRABAM"),
			help="Optional: Specify a location to store bam files "
			"[Default: $SRABAM:" + safe_echo("SRABAM") + "]")
	
	parser.add_argument(
			"-s", "--srafolder", default=safe_echo("SRARAW"),
			help="Optional: Specify a location to store pipeline output "
			"[Default: $SRARAW:" + safe_echo("SRARAW") + "]")

	parser.add_argument(
			"-m", "--mode", choices=[CONVERT, DELETE_SRA, DELETE_BAM], default=CONVERT,
			help="Main action to be performed for each input file; e.g. SRA-to-BAM conversion, SRA removal, or BAM removal")

	parser.add_argument(REMOVAL_OPTNAME, action="store_true", help="Once converted, remove original file.")

	# parser.add_argument(
	# 		"--picard", dest="picard_path", default=safe_echo("PICARD"),
	# 		help="Specify a path to the picard jar, if you want to convert "
	# 		"fastq to bam [Default: $PICARD:" + safe_echo("PICARD") + "]")
	
	parser.add_argument(
			"-r", "--srr", required=True, nargs="+", help="SRR files")

	parser = pypiper.add_pypiper_args(
		parser, groups=["pypiper", "config"], args=["output-parent"])
	return parser.parse_args(cmdl)


def safe_echo(var):
	""" Returns an environment variable if it exists, or an empty string if not"""
	return os.getenv(var, "")


if __name__ == "__main__":
	args = _parse_cmdl(sys.argv[1:])

	def do_rm():
		return getattr(args, REMOVAL_OPTNAME, None)

	if do_rm() and args.mode != CONVERT:
		print("Warning: {rm} applies only for mode {conv} (not {curr})".format(rm=REMOVAL_OPTNAME, conv=CONVERT, curr=args.mode))

	key = args.srr[0]
	outfolder = os.path.join(args.srafolder, "sra_convert_pipeline")
	pm = pypiper.PipelineManager(name="sra_convert", outfolder=outfolder, args=args)

	nfiles = len(args.srr)
	for i in range(nfiles):
		print("Processing " + str(i+1) + " of " + str(nfiles))
		infile = args.srr[i]
		srr_acc = os.path.splitext(os.path.basename(args.srr[i]))[0]
		outfile = os.path.join(args.bamfolder, srr_acc + ".bam")
		if not os.path.isfile(infile):
			infile = os.path.join(args.srafolder, args.srr[i] + ".sra")
			outfile = os.path.join(args.bamfolder, args.srr[i] + ".bam")

		if args.mode == CONVERT:
			cmd = "sam-dump -u {data_source} | samtools view -bS - > {outfile}".format(
				data_source=infile, outfile=outfile)
			target = outfile
			pm.run(cmd, target=target, follow=(lambda: os.unlink(infile)) if do_rm() else None)
		elif args.mode == DELETE_SRA:
			cmd = "rm {data_source}".format(data_source=infile)
			target = "rm_" + infile
			pm.run(cmd, target=target, nofail=True)
		elif args.mode == DELETE_BAM:
			cmd = "rm {data_source}".format(data_source=outfile)
			target = "rm_" + outfile
			pm.run(cmd, target=target, nofail=True)

	pm.stop_pipeline()
