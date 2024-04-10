#!/usr/bin/env python

import argparse
import os
import numpy

parser = argparse.ArgumentParser()
parser.add_argument('--input', '-i', help='Input file')
parser.add_argument('--cutoff', '-c', type=int, help='Cutoff value')
parser.add_argument('--percentile', '-p', type=float, help='Percentile value to calculate cutoffs')
parser.add_argument('what_shall_i_write_here', nargs=argparse.REMAINDER)
args = parser.parse_args()

def main():
	areas = parse_bed(args.input)
	potentialHotspots = []
	if args.percentile:
		cutoff = calc_percentile_cutoff(areas,args.percentile)
	elif args.cutoff:
		cutoff = int(args.cutoff)
	else:
		cutoff = 0

	for area in areas:
		# Check if the current area is the first or the last on the reference. If so then it is skipped
		if not areas.index(area) == 0 and not areas.index(area) == areas.index(areas[-1]) and area["coverage"] >= cutoff:
			currentIndex = areas.index(area)
			previousArea = areas[currentIndex - 1]
			followingArea = areas[currentIndex + 1]
			startConnected = True if area["start"] == previousArea["end"] else False
			endConnected = True if area["end"] == followingArea["start"] else False
			# Filter out all areas that are not connected to any other area -> might need to reconsider 
			if not startConnected and not endConnected:
				pass
			elif startConnected and not endConnected:
				if area["coverage"] >= previousArea["coverage"]:
					area["local_high"] = True
			elif not startConnected and endConnected:
				if area["coverage"] >= followingArea["coverage"]:
					area["local_high"] = True
			elif startConnected and endConnected:
				if area["coverage"] >= previousArea["coverage"] and area["coverage"] >= followingArea["coverage"] : 
					area["local_high"] = True			
			# Prolong the current area if other areas are surrounding it
			if area["local_high"]:
				tmpStart = area
				tmpEnd = area
				# Expand the start of the area
				if startConnected:
					while tmpStart["start"] == areas[areas.index(tmpStart) - 1]["end"] and areas[areas.index(tmpStart) - 1]["coverage"] >= cutoff:
						tmpStart = areas[areas.index(tmpStart) - 1]
						if areas.index(tmpStart) == 0:
							break
				# Expand the end of the area
				if endConnected:
					while tmpEnd["end"] == areas[areas.index(tmpEnd) + 1]["start"] and areas[areas.index(tmpEnd) + 1]["coverage"] >= cutoff:
						tmpEnd = areas[areas.index(tmpEnd) + 1]
						if areas.index(tmpEnd) == areas.index(areas[-1]):
							break

				start = int(tmpStart["start"])
				end = int(tmpEnd["end"])
				potentialHotspots.append({"start": start, "end": end, "coverage": area["coverage"], "length": end - start})
		else:
			pass
	outputString = 'start\tend\tcoverage\tlength'

	for i in potentialHotspots:
		outputString += '\n%s\t%s\t%s\t%s' % (i["start"],i["end"],i['coverage'],i["length"])
	print(outputString)

def parse_bed( bedFile ):
	areas = []
	with open ( bedFile, 'r' ) as bed:
		for line in bed:
			line = line.strip()
			(origin, start, end, cov) = line.split("\t")
			areas.append({"start": int(start), "end": int(end), "coverage": int(cov), "local_high": False})
	return areas

def calc_percentile_cutoff( areas,percentile ):
	values = []
	for i in areas:
		if(i["coverage"] > 0):
			values.append([i["coverage"]]*(i["end"] - i["start"]))
	values = sum(values,[])
	percentileCutoff = numpy.percentile(values, percentile)

	return percentileCutoff

main()