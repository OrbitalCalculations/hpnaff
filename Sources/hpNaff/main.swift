//
//  main.swift
//  hpnaff
//
//  Created by Heiko Pälike on 29/04/2015.
//  Copyright (c) 2015 Heiko Pälike. All rights reserved.
//

import Foundation
import Accelerate
//import Surge
//import CommandLineKit
import Rainbow
//import Upsurge
//import Numerics


let version = Bundle.main.infoDictionary?["CFBundleShortVersionString"] as? String ?? "unknown"
let version2 = Bundle.main.infoDictionary?["CFBundleVersion"] as? String ?? "unknown"
hpNAFFLog(("Welcome to hpNAFF Analysis".green.bold + "(version \(version.red), build \(version2.red)).").blue.bold + " ©2015-2019 Heiko Pälike. 3rd party licenses: see --licenses option")

// MARK: Defaults
let kDecimate = 1
let kmaxFrequencies = 30
let kfreqCycleLimit = 100
let kfracRMSChangeLimit = 0.000002
let kfracFreqAccuracyLimit = 0.00000001
let klowerFreqLimit = 0.0
let kupperFreqLimit = 0.5
let kcolumnSeparator = "\t"
let kskipHeaderLines = 0
let kabscissaColumn = 1
let kordinateColumn = 2
let kordinateComplexColumn = 3

let kdeltaT = 1.0
//let kchunks = 1

// MARK: Command Line setup
let cli = CommandLine()
cli.formatOutput = { s, type in
	var str: String
	switch(type) {
		
	case .error:
		str = s.red.bold
	case .optionFlag:
		str = s.magenta.underline
	case .optionHelp:
		str = s.blue
	default:
		str = s
	}
	return cli.defaultFormat(s: str, type: type)
}
let help_option = BoolOption(shortFlag: "h", longFlag: "help", required: false, helpMessage: "Displays help message.")
let fracRMSChangeLimit_option = DoubleOption(shortFlag: "r", longFlag: "fracRMSChangeLimit", required: false, helpMessage: "fracRMSChangeLimit. Default: 0.000002")
let fracFreqAccuracyLimit_option = DoubleOption(shortFlag: "f", longFlag: "fracFreqAccuracyLimit", required: false, helpMessage: "fracFreqAccuracyLimit. Default: 0.00000001")
let maxFrequencies_option = IntOption(shortFlag: "m", longFlag: "maxFrequencies", required: false, helpMessage: "maxFrequencies. Default: \(kmaxFrequencies)")
let freqCycleLimit_option = IntOption(shortFlag: "c", longFlag: "freqCycleLimit", required: false, helpMessage: "freqCycleLimit. Default: \(kfreqCycleLimit)")
let lowerFreqLimit_option = DoubleOption(shortFlag: "l", longFlag: "lowerFreqLimit", required: false, helpMessage: "lowerFreqLimit. Default: 0.0")
let upperFreqLimit_option = DoubleOption(shortFlag: "u", longFlag: "upperFreqLimit", required: false, helpMessage: "upperFreqLimit. Default: 0.5 (NYQUIST)")
let detrending_option = IntOption(shortFlag: "d", longFlag: "detrendingOrder", required: false, helpMessage: "Order for polynomial detrending. 0 removes DC only. >0 removes a polynomial fit. Default: not set")
let detrending_option2 = BoolOption(shortFlag: "dC", longFlag: "detrendChunks", required: false, helpMessage: "Option to apply detrending (-d) to individual chunks. Default: false")
let filterAbscissaLower = DoubleOption(shortFlag: "fAL", longFlag: "filterAbscissaLower", required: false, helpMessage: "filter on Abscissa (x-axis) data greater or equal than value specified.")
let filterAbscissaUpper = DoubleOption(shortFlag: "fAU", longFlag: "filterAbscissaUpper", required: false, helpMessage: "filter on Abscissa (x-axis) data lesser or equal than value specified.")

let columnSeparator_option = StringOption(shortFlag: "z", longFlag: "columnSeparator", required: false, helpMessage: "Use this as columns separtor for input file instead of default tab \\t. Default: \\t")
let skipHeaderLines_option = IntOption(shortFlag: "k", longFlag: "skipHeaderLines", required: false, helpMessage: "Specify number of (header)lines to be skipped. Default: 0")
//abscissa
let abscissa_option = IntOption(shortFlag: "a", longFlag: "abscissaColumn", required: false, helpMessage: "Specify column number for abscissa (time or depth). Default: 1")
let ordinate_option = IntOption(shortFlag: "o", longFlag: "ordinateColumn", required: false, helpMessage: "Specify column number for ordinate (data). Default: 2")
let complex_ordinate_option = IntOption(shortFlag: "Co", longFlag: "ComplexOrdinateColumn", required: false, helpMessage: "Specify column number for ordinate (data), complex part. Default: 3 (if C option used)")
let complex_option = BoolOption(shortFlag: "C", longFlag: "complex", required: false, helpMessage: "Compute spectra for complex input. Default: false. If true, by default takes values from abscissa column + 1.")
let deltat_option = DoubleOption(shortFlag: "t", longFlag: "deltaT", required: false, helpMessage: "Specify time/depth offset between data. Default: 1.0")
let chunking_option = IntOption(shortFlag: "s", longFlag: "chunkSize", required: false, helpMessage: "Split the input file into chunks of this size. Default: not set, no chunking")
let windowoffset_option = IntOption(shortFlag: "w", longFlag: "window_offsets", required: false, helpMessage: "Specify windowing offset for evolutive analysis. Default: chunkSize if set, or None.")
let ratio_option = BoolOption(shortFlag: "x", longFlag: "ratios", required: false, helpMessage: "Compute and output frequency ratio matrix. Default: false.")
let ratio_list = StringOption(longFlag: "ratio_list", required: false, helpMessage: "If ratio_option is set, provides comma separated list of frequency ratios to highlight in output.")
let output_detrended = StringOption(longFlag: "output_detrended_data", required: false, helpMessage: "Output detrended data to this filename. Not implemented yet.")
let display_licenses = BoolOption(longFlag: "licenses", required: false, helpMessage: "Shows licenses of code components used in this software.")
let display_references = BoolOption(longFlag: "references", required: false, helpMessage: "Shows references to papers and sources.")
let version_option = BoolOption(longFlag: "version", required: false, helpMessage: "Shows current version and build of this software")
let decimate_option = IntOption(shortFlag: "n", longFlag: "decimate", required: false, helpMessage: "Specify decimation of input data. Default: 1 (take every value)")


let verbosity = CounterOption(shortFlag: "v", longFlag: "verbose", required: false,
							  helpMessage: "Print verbose messages. Specify multiple times to increase verbosity.")

cli.addOptions(help_option, fracRMSChangeLimit_option, fracFreqAccuracyLimit_option, maxFrequencies_option, freqCycleLimit_option, lowerFreqLimit_option, upperFreqLimit_option, detrending_option, detrending_option2, filterAbscissaLower, filterAbscissaUpper, chunking_option, windowoffset_option, columnSeparator_option, skipHeaderLines_option, decimate_option, abscissa_option, ordinate_option, complex_option, complex_ordinate_option, deltat_option,ratio_option, ratio_list, output_detrended, version_option, verbosity, display_references, display_licenses)

var fileNames = [String]()
do {
	try cli.parse()
	fileNames = cli.unparsedArguments
} catch {
	//cli.printUsage(error)
	cli.printCustomUsage("Usage: hpnaff [options] file\nfile format: tab separated columns, data must be equidistantly spaced.  Instead of file can specify '-' for stdin.")
	exit(EX_USAGE)
}

if help_option.wasSet {
	cli.printCustomUsage("Usage: hpnaff [options] file\nfile format: tab separated columns, data must be equidistantly spaced. Instead of file can specify '-' for stdin.")
	exit(0)
}

if version_option.wasSet {
	print("hpNAFF (version \(version.red), build \(version2.red)) ©2015-2019 Heiko Pälike. 3rd party licenses: see --licenses option")
	exit(0)
}

if display_licenses.wasSet {
	print("hpNAFF Licenses:".bold)
	print(licenses.blue)
	exit(0)
}

if display_references.wasSet {
	print("hpNAFF Notes and References:".bold)
	print(references.blue)
	exit(0)
}


let fracRMSChangeLimit = fracRMSChangeLimit_option.value ?? kfracRMSChangeLimit // 0.000002
let fracFreqAccuracyLimit = fracFreqAccuracyLimit_option.value ?? kfracFreqAccuracyLimit //0.00000001
let maxFrequencies = maxFrequencies_option.value ?? kmaxFrequencies // 30
let freqCycleLimit = maxFrequencies_option.value ?? kfreqCycleLimit // 100
let lowerFreqLimit = lowerFreqLimit_option.value ?? klowerFreqLimit // 0.0
let upperFreqLimit = upperFreqLimit_option.value ?? kupperFreqLimit // 0.5
let detrendOrder : Int? = detrending_option.value
let detrendChunks : Bool = detrending_option2.value
let complex : Bool = complex_option.value
let decimate = decimate_option.value ?? kDecimate


let columnSeparator = columnSeparator_option.value ?? kcolumnSeparator
let skipHeaderLines = skipHeaderLines_option.value ?? kskipHeaderLines

let abscissaColumn = (abscissa_option.value ?? kabscissaColumn) - 1 //transform to 0 counted columns
let ordinateColumn = (ordinate_option.value ?? kordinateColumn) - 1
let ordinateComplexColumn = (complex_ordinate_option.value ?? kordinateComplexColumn) - 1

var deltaT = deltat_option.value ?? kdeltaT


var ratios = [Double]()
if let ratio_list_string = ratio_list.value {
	if ratio_option.wasSet == false {
		hpNAFFLog("WARNING: ratio_list option was set but ratio_option was not set", messageKind: .error)
	} else {
		ratios = ratio_list_string.split(separator: ",").compactMap{Double($0)}
	}
}

if let _ = windowoffset_option.value {
	if chunking_option.wasSet == false {
		hpNAFFLog("WARNING: windowing option was set but chunking_option was not set", messageKind: .error)
	}
}


guard (Set([abscissaColumn, ordinateColumn, ordinateComplexColumn]).count > 1) else {
	hpNAFFLog("ERROR: Specified abscissa(\(abscissaColumn+1), ordinate(\(ordinateColumn+1), and complex ordinate (\(ordinateComplexColumn+1)) columns cannot be the same ", messageKind: .error)
	exit(-1)
}

//print(ExecutableUUID())

guard fileNames.count == 1,
	let fileName = fileNames.first else {
		hpNAFFLog("Please provide exactly one file name argument", messageKind: .error)
		cli.printCustomUsage("Usage: hpnaff [options] file\nfile format: tab separated columns, data must be equidistantly spaced. Instead of file can specify '-' for stdin.")
		exit(-1)
}
hpNAFFLog( "Chosen File: \(fileName == "-" ? "STDIN" : fileName)", messageKind: .info, verbosity: verbosity.value)

enum hpNAFFError: Error {
	case unknownError
	case missingInputFileError(String)
}


// MARK: Data reading
let inputurl = URL(fileURLWithPath: fileName)
//let data = NSData(contentsOfURL: inputurl)
var errorptr: NSErrorPointer = nil

var time = [Double]()
var data = [Double]()
var cplxdata = [Double]()
// TODO
do {
	var lineCounter = 0
	
	var mystring: String = ""
	
	if (fileName == "-") {
		while let newline = readLine(strippingNewline: false) {
			//print(newline)
			mystring.append(newline)
		}
	} else {
		mystring = try String(contentsOf: inputurl, encoding: .utf8)
	}
	
	let dataString = mystring.components(separatedBy: CharacterSet.newlines)
	let dataLines = Array(dataString[skipHeaderLines...].lazy.filter({!$0.isEmpty }))
	
	let totalDataCount = dataLines.count
	let finalData = stride(from: 0, to: (totalDataCount), by: decimate).map{dataLines[$0]}
	
	
	if (decimate != 1) {
		hpNAFFLog("Decimate option set to \(decimate). Adjusting deltaT from \(deltaT) to \(deltaT*Double(decimate))", messageKind: .info, verbosity: verbosity.value)
		deltaT *= Double(decimate)
	}
	
	for line in finalData {
		lineCounter += 1
		let columns = line.trimmingCharacters(in: .whitespaces).split(separator: Character(columnSeparator), omittingEmptySubsequences: true).map{$0.trimmingCharacters(in: .whitespaces)}
		
		if lineCounter <= skipHeaderLines {
			continue //skips this line
		}
		guard columns.count > 1  else {
			hpNAFFLog("ERROR: file \(inputurl.standardizedFileURL), line #\(lineCounter) has less than 2 columns. Aborting.\n Original line = '\(mystring)'", messageKind: .error)
			exit(-1)
		}
		
		let validRange = 0..<columns.count
		
		
		
		if (complex == true) {
			guard (validRange.contains(abscissaColumn) && validRange.contains(ordinateColumn) && validRange.contains(ordinateComplexColumn)) else {
				hpNAFFLog("ERROR: file \(inputurl.standardizedFileURL), line #\(lineCounter) does not contain columns for requested abscissa(\(abscissaColumn+1)), ordinate(\(ordinateColumn+1)) and complex ordinate(\(ordinateComplexColumn+1)). Line has \(columns.count) columns. \n Original line = '\(mystring)'",messageKind: .error)
				exit(-1)
			}
		} else {
			guard (validRange.contains(abscissaColumn) && validRange.contains(ordinateColumn)) else {
				hpNAFFLog("ERROR: file \(inputurl.standardizedFileURL), line #\(lineCounter) does not contain columns for requested abscissa(\(abscissaColumn+1)), and ordinate(\(ordinateColumn+1)). Line has \(columns.count) columns. \n Original line = '\(mystring)'",messageKind: .error)
				exit(-1)
			}
		}
		
		
		if (complex == true) {
			guard columns.count > 2  else {
				hpNAFFLog("ERROR: file \(inputurl.standardizedFileURL), line #\(lineCounter) has less than 3 columns. Aborting", messageKind: .error)
				exit(-1)
			}
		} else {
			guard columns.count > 1  else {
				hpNAFFLog("ERROR: file \(inputurl.standardizedFileURL), line #\(lineCounter) has less than 2 columns. Aborting", messageKind: .error)
				exit(-1)
			}
		}
		
		if (complex == true) {
			guard let timeval = Double(columns[abscissaColumn].trimmingCharacters(in: CharacterSet.whitespaces)),
				let dataval = Double(columns[ordinateColumn].trimmingCharacters(in: CharacterSet.whitespaces)),
				let cplxval = Double(columns[ordinateComplexColumn].trimmingCharacters(in: CharacterSet.whitespaces)) else {
					hpNAFFLog("ERROR: file \(inputurl.standardizedFileURL), line #\(lineCounter) has column data that cannot be parsed as a number [\(columns[abscissaColumn]), \(columns[ordinateColumn]), \(columns[ordinateComplexColumn])]. Aborting", messageKind: .error)
					exit(-1)
			}
			time.append(timeval)
			data.append(dataval)
			cplxdata.append(cplxval)
		} else {
			guard let timeval = Double(columns[abscissaColumn].trimmingCharacters(in: CharacterSet.whitespaces)),
				let dataval = Double(columns[ordinateColumn].trimmingCharacters(in: CharacterSet.whitespaces)) else {
					hpNAFFLog("ERROR: file \(inputurl.standardizedFileURL), line #\(lineCounter) has column data that cannot be parsed as a number [\(columns[abscissaColumn]), \(columns[ordinateColumn])]. Aborting", messageKind: .error)
					exit(-1)
			}
			time.append(timeval)
			data.append(dataval)
		}
		
		
	}
	
	let count = data.count
	hpNAFFLog("read \(count) lines", messageKind: .info, verbosity: verbosity.value)
	
} catch let error as NSError {
	errorptr?.pointee = error
}

let filterLower : Double = filterAbscissaLower.value ?? time.min() ?? 0.0
let filterUpper : Double = filterAbscissaUpper.value ?? time.max() ?? 0.0 //max(time)

let startIndex1 = time.firstIndex(where: {$0 >= filterLower})
let startIndex2 = time.lastIndex(where: {$0 <= filterLower})
let endIndex1 = time.lastIndex(where: {$0 <= filterUpper})
let endIndex2 = time.firstIndex(where: {$0 >= filterUpper})

var startIndex: Int
var endIndex: Int

let filterDiff = fabs(filterLower-filterUpper)

guard let startIndex1 = startIndex1, let endIndex1 = endIndex1,
	let startIndex2 = startIndex2, let endIndex2 = endIndex2
	else {
		fatalError("cannot obtain indices")
//		exit(-1)
}

if fabs(time[endIndex1]-time[startIndex1]) <= filterDiff {
	startIndex = min(startIndex1,endIndex1)
	endIndex = max(startIndex1,endIndex1)
} else {
	startIndex = min(startIndex2,endIndex2)
	endIndex = max(startIndex2,endIndex2)
}

let lowerIndex = startIndex
let upperIndex = endIndex

hpNAFFLog("filter lower,upper = \(filterLower), \(filterUpper), \(lowerIndex), \(upperIndex)", messageKind: .info, verbosity: verbosity.value)


// MARK: Detrending

let origData = data
let origCplx = cplxdata
//print("detrendChunks = \(detrendChunks)")

time = Array(time[lowerIndex...upperIndex])
data = Array(data[lowerIndex...upperIndex])
if (complex == true) { cplxdata =  Array(cplxdata[lowerIndex...upperIndex]) }

if detrendChunks == false {
	if let detrendOrder = detrending_option.value {
		if (complex == true) {
      let (realresult,realcoefficients) =  (data, [1.0]) //TODO: implement detrend(data, degree: detrendOrder)//detrend(real: data, imag: cplxdata, degree: detrendOrder)
      let (imagresult,imagcoefficients) =  (data, [1.0]) //TODO: implement detrend(cplxdata, degree: detrendOrder)//detrend(real: data, imag: cplxdata, degree: detrendOrder)
			
			data = realresult
			cplxdata = imagresult
			hpNAFFLog("Detrending coefficients (REAL) order (high to low) \(detrendOrder) : \(realcoefficients)", messageKind: .info, verbosity: verbosity.value)
			hpNAFFLog("Detrending coefficients (IMAG) order (high to low) \(detrendOrder) : \(imagcoefficients)", messageKind: .info, verbosity: verbosity.value)
			
		}
		else {
			let (detrendedSignal, coefficients) = (data, [1.0]) //TODO: implement detrend(data, degree: detrendOrder)
			hpNAFFLog("Detrending coefficients order (high to low) \(detrendOrder) : \(coefficients)", messageKind: .info, verbosity: verbosity.value)
			data = detrendedSignal
		}
		//        let signal_to_detrend = complex ? zip(data,cplxdata).map{hypot($0.0,$0.1)} : data
		//let
		
		/*if let detrendedFilePath = output_detrended.value {
		//write_detrended(detrendedFilePath, origData, data)
		}*/
	} else {
		if output_detrended.wasSet {
			hpNAFFLog("Option to output detrended value set, but no detrending requested.", messageKind: .warn, verbosity: verbosity.value)
		}
	}
}

// MARK: Chunking

let totalDataCount = data.count


var chunkingSize = totalDataCount //Default
var windowOffset = totalDataCount //Default

if let requestedChunkingSize = chunking_option.value,
	requestedChunkingSize <= totalDataCount/2
{
	chunkingSize = requestedChunkingSize
	windowOffset = windowoffset_option.value ?? chunkingSize
}
else {
	if chunking_option.wasSet {
		hpNAFFLog("WARNING: requested chunk size \(String(describing: chunking_option.value)) > data.count/2 (\(data.count/2)); Ignoring and using whole data set", messageKind: .error)
		
	}
	
}



let dataArrays = data.chunked(by: chunkingSize, windowOffset: windowOffset)
let imagArrays = cplxdata.chunked(by: chunkingSize, windowOffset: windowOffset)

let timeArrays = time.chunked(by: chunkingSize, windowOffset: windowOffset)
let cplxDataArrays : [[Double]]

if (complex == true) {
	cplxDataArrays = cplxdata.chunked(by: chunkingSize, windowOffset: windowOffset)
} else {
	cplxDataArrays = [[]]
}

let numberChunks = dataArrays.count

let maxSize = dataArrays.map{ $0.count }.max() ?? 0
var weights: [Double] //vDSP_DFT_SetupD?


let realweights = [Double]() //TODO: vDSP_DFT_zrop_CreateSetupD(nil, vDSP_Length(maxSize), .FORWARD)
let cplxWeights = [Double]() //TODO: vDSP_DFT_zop_CreateSetupD(weights, vDSP_Length(maxSize), .FORWARD)
if (complex == true) {
	weights = cplxWeights
} else {
	weights = realweights
}

for chunk in 0..<numberChunks {
	var localData : [Double] = []
	var localCplxData : [Double] = []
	
	if detrendChunks == true {
		if (detrendOrder == nil) {
			hpNAFFLog("Detrending coefficients order not set, using default of 0", messageKind: .warn, verbosity: verbosity.value)
		}
    let (detrendedSignal, coefficients) = (dataArrays[chunk],[1.0]) //TODO: implement detrend(dataArrays[chunk], degree: detrendOrder ?? 0)
		hpNAFFLog("Chunk #\(chunk), detrending coefficients order \(detrendOrder ?? 0) : \(coefficients)", messageKind: .info, verbosity: verbosity.value)
		localData = detrendedSignal
		localCplxData = detrendedSignal
		if (complex == true) {
      let (detrendedImagSignal, imagCoefficients) = (cplxDataArrays[chunk], [1.0]) //TODO: detrend(cplxDataArrays[chunk], degree: detrendOrder ?? 0)
			hpNAFFLog("Chunk #\(chunk), detrending coefficients (IMAG) order \(detrendOrder ?? 0) : \(imagCoefficients)", messageKind: .info, verbosity: verbosity.value)
			localCplxData = detrendedImagSignal
		}
	} else {
		localData = dataArrays[chunk]
		if (complex == true) {
			localCplxData = cplxDataArrays[chunk]
		}
	}
	
	let time = timeArrays[chunk]
	let points = data.count
	
	if points < 2 {
		continue
	}
	guard let t0 = time.first else {
		hpNAFFLog("ERROR: Chunk #\(chunk) does not have a valid time value at index #0. Aborting", messageKind: .error)
		exit(-1)
	}
	let dt = deltaT
	
	let result = performNAFF(data: localData, cplxData: localCplxData, dt: dt, nfreqs: maxFrequencies, t0: t0, maxFrequencies: maxFrequencies, fracRMSChangeLimit: fracRMSChangeLimit, freqCycleLimit: freqCycleLimit, fracFreqAccuracyLimit: fracFreqAccuracyLimit, lowerFreqLimit: lowerFreqLimit, upperFreqLimit: upperFreqLimit, weights: weights, complex: complex)
	let (freq,ampl,phase,sig) = result
	
	let count = freq.count
	//let headerstr = ["    #","chunk","      t0      ","points","   frequency","    amplitude","       phase","significance", "period"].joined(separator: " | ").blue
	let headerstr1 = ["    #","chunk","      t0      ","points","   frequency","    amplitude","       phase","significance", "      period","   arcsec/yr"," phase (deg)"].joined(separator: " | ")
	let headerstr = headerstr1.blue
	
	
	let separatorline = [String](repeating:"=",count:headerstr1.count).joined(separator: "").blue
	print(separatorline)
	print(headerstr)
    let indices = 0..<count
        /*zip(0..<count, freq).sorted { (lhs: (Int, Double), rhs: (Int, Double)) -> Bool in
        return rhs.1 > lhs.1
    }.map{$0.0}*/
    
	for i in 0..<count {
		let chunkstr = String(format:"% 5d",chunk)
		let timestr = String(format:"% 14.4f",t0)
		let pointsstr = String(format:"% 6d",localData.count) //points
		let index = String(format:"% 5d",indices[i])
		let f = String(format:"% 011.9f",freq[indices[i]])
		let a = String(format:"%  13.8f",ampl[indices[i]])
		let p = String(format:"% 011.9f",phase[indices[i]])
		let s = String(format:"% 011.9f",sig[indices[i]])
		let per = String(format:"% 12.7f",1.0/freq[indices[i]])
		let arcsec = String(format:"% 12.7f",360.0*3.600*freq[indices[i]])
		let pdeg = String(format:"% 12.7f",phase[indices[i]]*180.0/Double.pi)
		
		print(index,chunkstr,timestr,pointsstr,f, a, p, s,per, arcsec,pdeg, separator:" | ")
	}
	
	if ratio_option.value == true {
		print(separatorline)
		let oldindices = freq.indices.map{Int($0)}
		let sortedfreqs  = zip(freq[1...],oldindices[1...]).sorted(by: {$0.0<$1.0})
		//print(sortedfreqs)
		let sortedindices = sortedfreqs.map{$0.1}
		let sortedfreqsonly = sortedfreqs.map{$0.0}
		
		let lineseparation = ("=====" + [String](repeating:"===========",count: count-1).joined(separator: "")).green
		//let indices = (1..<count).map{String(format:"% 8d",$0).red}
		let indices = sortedindices.map{String(format:"% 8d",$0).red}
		
		print("    #", indices.joined(separator: " | "), separator: " | ")
		print(lineseparation)
		//for i in 1..<count {
		
		for i in 0..<sortedindices.count {
			let index = String(format:"% 5d",sortedindices[i]).red
			var ratioterms = [String]()
			for j in 0..<sortedindices.count {
				let ratio = sortedfreqsonly[j]/sortedfreqsonly[i]
				var str = String(format:"%8.3f",ratio)
				if (i==j) {str = str.lightBlue}
				for r in ratios {
					if abs(ratio-r)<0.05 {
						str = str.red
					}
				}
				
				ratioterms.append(str)
			}
			//let formattedratioterms = [String](repeating:"        ",count:i-1) + ratioterms
			let formattedratioterms = ratioterms
			//let ratiotermstring = ratioterms.joined(separator: " | ")
			print(index, formattedratioterms.joined(separator: " | "), separator: " | ")
		}
		print(lineseparation)
	}
}

//vDSP_DFT_DestroySetupD(weights)
