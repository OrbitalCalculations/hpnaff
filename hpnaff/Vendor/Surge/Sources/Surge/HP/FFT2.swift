//
//  FFT2.swift
//  Surge
//
//  Created by Heiko Pälike on 28/05/2015.
//  Copyright (c) 2015 Mattt Thompson. All rights reserved.
//

import Accelerate

// MARK: Fast Fourier Transform
public func paddedData(_ data: [Double]) -> ContiguousArray<Double>{
	let original_length = data.count
	let centeridx = Int(ceil(Double(original_length)/2.0))
	let dftLength = nextDFTn(original_length)
	let paddingLength = dftLength - original_length
	let padding = Array<Double>(repeating: 0.0, count: paddingLength)
	// zero phase padding  http://dsp.stackexchange.com/questions/18938/merits-of-zero-phase-zero-padding/
	let result = ContiguousArray([Array(data[centeridx..<original_length]), padding, Array(data[0..<centeridx])].joined())
	return result
}

public func matlabrealfft(_ real_in: [Double]) -> [Double] {

	//let dftres : [DSPDoubleComplex] = realfft(real_in)
	
	
	let original_length = real_in.count
	let dftLength = nextDFTn(original_length)
	
	let paddingLength = dftLength - original_length
	let padding = Array<Double>(repeating: 0.0, count: paddingLength)
	let paddingcenteridx = Int(ceil(Double(paddingLength)/2.0))
	
	let length = dftLength
	// zero phase padding  http://dsp.stackexchange.com/questions/18938/merits-of-zero-phase-zero-padding/
	let inputsignal = Array(ContiguousArray([Array(padding[0..<paddingcenteridx]), real_in, Array(padding[paddingcenteridx..<paddingLength])].joined()))

	
	let halflength = length/2
	
	
	
	let forward = vDSP_DFT_zrop_CreateSetupD(nil, vDSP_Length(length), .FORWARD)
	let inverse = vDSP_DFT_zrop_CreateSetupD(forward, vDSP_Length(length), .INVERSE)
	assert((forward != nil) && (inverse != nil),"FFT Setup Failed")

	var even = [Double](repeating: 0.0, count:halflength)
	var odd = [Double](repeating: 0.0, count:halflength)
	var complex = DSPDoubleSplitComplex(realp:&even, imagp: &odd)

	UnsafePointer(inputsignal).withMemoryRebound(to: DSPDoubleComplex.self, capacity: halflength) { inputsignal in
		vDSP_ctozD(inputsignal, vDSP_Stride(2), &complex, vDSP_Stride(1), vDSP_Length(halflength))
	}

	var real = [Double](repeating: 0.0, count:length)
	var imag = [Double](repeating: 0.0, count: length)
	
	vDSP_DFT_ExecuteD(forward!, &even, &odd, &real, &imag)
	
	var scale = 0.5 //0.5 to account for vDSP zrop   impl = 2 * math

	var realout = real
	var imagout = imag
	vDSP_vsmulD(&real, 1, &scale, &realout, 1, vDSP_Length(length))
	vDSP_vsmulD(&imag, 1, &scale, &imagout, 1, vDSP_Length(length))


	for i in (1..<halflength) { //vDSP only outputs first half of mirrored set
		let mirror = length-i
		realout[mirror] = realout[i]
		imagout[mirror] = -imagout[i] //complex conjugate
	}

	realout[halflength] = imagout[0]; imagout[0] = 0.0
	//var split_complex_out = DSPDoubleSplitComplex(realp:&real, imagp: &imag)
	
	let result = zip(realout,imagout).flatMap { [$0, $1] }

	return result
}

public func realfftmag(_ real_in: [Double], weights inputweights: vDSP_DFT_SetupD? = nil) -> [Double] {
    
	let dftres : [DSPDoubleComplex] = realfft(real_in, weights: inputweights)
	let length = dftres.count
	let FFTFreqs = length/2+1
	var real = [Double](repeating: 0.0, count: length)
	var imag = [Double](repeating: 0.0, count: length)

	var complex = DSPDoubleSplitComplex(realp:&real, imagp: &imag)
	vDSP_ctozD(dftres, vDSP_Stride(2), &complex, vDSP_Stride(1), vDSP_Length(length))

	var scale = 1.0/Double(length) //0.5 to account for vDSP zrop   impl = 2 * math
	//vDSP_vsmulD(&out_real, 1, &scale, &out_real, 1, vDSP_Length(length))
	//vDSP_vsmulD(&out_imag, 1, &scale, &out_imag, 1, vDSP_Length(length))
    var realinp = real
    var imaginp = imag
	vDSP_vsmulD(&realinp, 1, &scale, &real, 1, vDSP_Length(length))
	vDSP_vsmulD(&imaginp, 1, &scale, &imag, 1, vDSP_Length(length))

	real = sqr(real)
	imag = sqr(imag)
	
	let magnitudesq = add(real,imag)
	let result = Array(magnitudesq[0..<FFTFreqs])
    return result
}

extension Collection where Index: Comparable {
	subscript(back i: Int) -> Iterator.Element {
		let backBy = i + 1
		return self[self.index(self.endIndex, offsetBy: -backBy)]
	}
}

public func cplxfftmag(_ real_in: [Double], _ cplx_in: [Double], weights inputweights: vDSP_DFT_SetupD? = nil) -> [Double] {
    
    let dftres : [DSPDoubleComplex] = cplxfftnopad(real_in, cplx_in, weights: inputweights)
    let length = dftres.count
   // let FFTFreqs = length/2+1
	
    var real = [Double](repeating: 0.0, count: length)
    var imag = [Double](repeating: 0.0, count: length)
    
    var complex = DSPDoubleSplitComplex(realp:&real, imagp: &imag)
    vDSP_ctozD(dftres, vDSP_Stride(2), &complex, vDSP_Stride(1), vDSP_Length(length))
    
    real = sqr(real)
    imag = sqr(imag)
	
	//var magnitudesq = [Double]()
    //magnitudesq.reserveCapacity(length) //FFTFreqs
	
	/*magnitudesq.append(real[0]+imag[0])
	for i in 1..<length/2 {
		magnitudesq.append(real[i]+real[back: i-1]+imag[i]+imag[back: i-1])
	}
	magnitudesq.append(real[length/2]+imag[length/2])*/
	//magnitudesq = real + imag
    let magnitudesq = add(real,imag)
    //let result = Array(magnitudesq[0..<FFTFreqs])
    return magnitudesq
}

public func realfft(_ real_in: [Double], weights inputweights: vDSP_DFT_SetupD? = nil) -> [DSPDoubleComplex] {
	
	let original_length = real_in.count
	let dftLength = nextDFTn(original_length)
	
	#if DEBUG
		//print("entering realfft, original length of input = \(original_length), using padded length \(dftLength)")
	#endif
	
	let paddingLength = dftLength - original_length
	let padding = Array<Double>(repeating: 0.0, count: paddingLength)
	let paddingcenteridx = Int(ceil(Double(paddingLength)/2.0))
	
	let length = dftLength
	// zero phase padding  http://dsp.stackexchange.com/questions/18938/merits-of-zero-phase-zero-padding/
	let inputsignal = Array(ContiguousArray([Array(padding[0..<paddingcenteridx]), real_in, Array(padding[paddingcenteridx..<paddingLength])].joined()))

	//let inputsignal = [Double](real_in[0..<centeridx])+ padding
	
	let halflength = length/2
	
	
    let weights: vDSP_DFT_SetupD = inputweights ?? fftPlanR(dftLength)
	//let forward = vDSP_DFT_zrop_CreateSetupD(nil, vDSP_Length(length), .FORWARD)
	//let inverse = vDSP_DFT_zrop_CreateSetupD(forward, vDSP_Length(length), .INVERSE)
	//assert((forward != nil) && (inverse != nil),"FFT Setup Failed")

	var even = [Double](repeating: 0.0, count:halflength)
	var odd = [Double](repeating: 0.0, count:halflength)
	var complex = DSPDoubleSplitComplex(realp:&even, imagp: &odd)

	UnsafePointer(inputsignal).withMemoryRebound(to: DSPDoubleComplex.self, capacity: halflength) { inputsignal in
		vDSP_ctozD(inputsignal, vDSP_Stride(2), &complex, vDSP_Stride(1), vDSP_Length(halflength))
	}
	//vDSP_ctozD(UnsafePointer<DSPDoubleComplex>(inputsignal), vDSP_Stride(2), &complex, vDSP_Stride(1), vDSP_Length(halflength))

	var real = [Double](repeating: 0.0, count:length)
	var imag = [Double](repeating: 0.0, count: length)
	
	vDSP_DFT_ExecuteD(weights, &even, &odd, &real, &imag)
	//  At this point, we have the forward real-to-complex DFT, which agrees with
	//  MATLAB up to a factor of two.  Since we want to double all but DC and NY
	//  as part of the Hilbert transform anyway, I'm not going to bother to
	//  unscale the rest of the frequencies -- they're already the values that
	//  we really want.  So we just need to move NY into the "right place",
	//  and scale DC and NY by 0.5.  The reflection frequencies are already
	//  zeroed out because the real-to-complex DFT only writes to the first n/2
	//  elements of real and imag.
	real[0] *= 0.5;	real[length/2] = 0.5 * imag[0];	imag[0] = 0.0
	//real[length/2] = imag[0];	imag[0] = 0.0
	
//	var out_real = [Double](repeating: 0.0, count: length)
//	var out_imag = [Double](repeating: 0.0, count: length)
	//  Now we have the completed hilbert transform up to a scale factor of n.
	//  We can unscale using vDSP_vsmulD.
	// vDSP_DFT_ExecuteD(inverse!, &real, &imag, &out_real, &out_imag)
	//  Now we have the completed hilbert transform up to a scale factor of n.
	//  We can unscale using vDSP_vsmulD.
	//var scale = 1.0/Double(length) //0.5 to account for vDSP zrop   impl = 2 * math
//	vDSP_vsmulD(&real, 1, &scale, &real, 1, vDSP_Length(length))
//	vDSP_vsmulD(&imag, 1, &scale, &imag, 1, vDSP_Length(length))
	

	var split_complex_out = DSPDoubleSplitComplex(realp:&real, imagp: &imag)
    //var result = DSPDoubleSplitComplex(realp:&real, imagp: &imag)

	//vDSP_DFT_ExecuteD(forward, &real, &imag, &real_out, &imag_out)
	
	var result : [DSPDoubleComplex] = [DSPDoubleComplex](repeating: DSPDoubleComplex(real: 0.0, imag: 0.0), count: length)
	vDSP_ztocD(&split_complex_out, 1, &result, 2, vDSP_Length(length))
	
	//if (hasOddLength == true){
	//  result.removeLast()
	//}
    if inputweights == nil {
        vDSP_DFT_DestroySetupD(weights)
    }
	return result
	
}

public func cplxfft(_ real_in: [Double], _ cplx_in: [Double], weights inputweights: vDSP_DFT_SetupD? = nil) -> [DSPDoubleComplex] {
    
    let original_length = real_in.count
    let dftLength = nextDFTn(original_length)
    
    #if DEBUG
    //print("entering realfft, original length of input = \(original_length), using padded length \(dftLength)")
    #endif
    
    let paddingLength = dftLength - original_length
    let padding = Array<Double>(repeating: 0.0, count: paddingLength)
    let paddingcenteridx = Int(ceil(Double(paddingLength)/2.0))
    
    let length = dftLength
    // zero phase padding  http://dsp.stackexchange.com/questions/18938/merits-of-zero-phase-zero-padding/
    var inputsignal = Array(ContiguousArray([Array(padding[0..<paddingcenteridx]), real_in, Array(padding[paddingcenteridx..<paddingLength])].joined()))
    var inputsignalcplx = Array(ContiguousArray([Array(padding[0..<paddingcenteridx]), cplx_in, Array(padding[paddingcenteridx..<paddingLength])].joined()))
    
    //let inputsignal = [Double](real_in[0..<centeridx])+ padding
    
    //let halflength = length/2
    
    
    let weights: vDSP_DFT_SetupD = inputweights ?? fftPlanC(dftLength)
    //let forward = vDSP_DFT_zrop_CreateSetupD(nil, vDSP_Length(length), .FORWARD)
    //let inverse = vDSP_DFT_zrop_CreateSetupD(forward, vDSP_Length(length), .INVERSE)
    //assert((forward != nil) && (inverse != nil),"FFT Setup Failed")
    
    //var even = [Double](repeating: 0.0, count:halflength)
    //var odd = [Double](repeating: 0.0, count:halflength)
    //var complex = DSPDoubleSplitComplex(realp:&even, imagp: &odd)
    
    //UnsafePointer(inputsignal).withMemoryRebound(to: DSPDoubleComplex.self, capacity: halflength) { inputsignal in
    //    vDSP_ctozD(inputsignal, vDSP_Stride(2), &complex, vDSP_Stride(1), vDSP_Length(halflength))
    //}
    let complex = DSPDoubleSplitComplex(realp: &inputsignal, imagp: &inputsignalcplx)
    //vDSP_ctozD(UnsafePointer<DSPDoubleComplex>(inputsignal), vDSP_Stride(2), &complex, vDSP_Stride(1), vDSP_Length(halflength))
    
    var real = [Double](repeating: 0.0, count: length)
    var imag = [Double](repeating: 0.0, count: length)
    
    vDSP_DFT_ExecuteD(weights, complex.realp, complex.imagp, &real, &imag)
    //vDSP_DFT_ExecuteD(<#T##__Setup: OpaquePointer##OpaquePointer#>, <#T##__Ir: UnsafePointer<Double>##UnsafePointer<Double>#>, <#T##__Ii: UnsafePointer<Double>##UnsafePointer<Double>#>, <#T##__Or: UnsafeMutablePointer<Double>##UnsafeMutablePointer<Double>#>, <#T##__Oi: UnsafeMutablePointer<Double>##UnsafeMutablePointer<Double>#>)
    //  At this point, we have the forward real-to-complex DFT, which agrees with
    //  MATLAB up to a factor of two.  Since we want to double all but DC and NY
    //  as part of the Hilbert transform anyway, I'm not going to bother to
    //  unscale the rest of the frequencies -- they're already the values that
    //  we really want.  So we just need to move NY into the "right place",
    //  and scale DC and NY by 0.5.  The reflection frequencies are already
    //  zeroed out because the real-to-complex DFT only writes to the first n/2
    //////real[0] *= 0.5;    real[length/2] = 0.5 * imag[0];    imag[0] = 0.0
    //real[length/2] = imag[0];    imag[0] = 0.0
    
    //    var out_real = [Double](repeating: 0.0, count: length)
    //    var out_imag = [Double](repeating: 0.0, count: length)
    //  Now we have the completed hilbert transform up to a scale factor of n.
    //  We can unscale using vDSP_vsmulD.
    // vDSP_DFT_ExecuteD(inverse!, &real, &imag, &out_real, &out_imag)
    //  Now we have the completed hilbert transform up to a scale factor of n.
    //  We can unscale using vDSP_vsmulD.
    //var scale = 1.0/Double(length) //0.5 to account for vDSP zrop   impl = 2 * math
    //    vDSP_vsmulD(&real, 1, &scale, &real, 1, vDSP_Length(length))
    //    vDSP_vsmulD(&imag, 1, &scale, &imag, 1, vDSP_Length(length))
    
    
    var split_complex_out = DSPDoubleSplitComplex(realp:&real, imagp: &imag)
    //var result = DSPDoubleSplitComplex(realp:&real, imagp: &imag)
    
    //vDSP_DFT_ExecuteD(forward, &real, &imag, &real_out, &imag_out)
    
    var result = [DSPDoubleComplex](repeating: DSPDoubleComplex(real: 0.0, imag: 0.0), count: length)
    vDSP_ztocD(&split_complex_out, 1, &result, 2, vDSP_Length(length))
    
    //if (hasOddLength == true){
    //  result.removeLast()
    //}
    if inputweights == nil {
        vDSP_DFT_DestroySetupD(weights)
    }
    return result
    
}

public func cplxfftnopad(_ real_in: [Double], _ cplx_in: [Double], weights inputweights: vDSP_DFT_SetupD? = nil) -> [DSPDoubleComplex] {
	
	
	let length = real_in.count
	let weights: vDSP_DFT_SetupD = inputweights ?? fftPlanC(real_in.count)
	//let weights: vDSP_DFT_SetupD = fftPlanC(real_in.count)

	var inputsignal = real_in
	var inputsignalcplx = cplx_in
	
	//var realtest=(0...31).map{cos(2.0*Double.pi/32.0*Double($0)*4)}
	//var imagtest=(0...31).map{sin(2.0*Double.pi/32.0*Double($0)*4)}

	let complex = DSPDoubleSplitComplex(realp: &inputsignal, imagp: &inputsignalcplx)
	//let complex = DSPDoubleSplitComplex(realp: &realtest, imagp: &imagtest)

	//print("INPUT")
	//_ = zip(real_in,cplx_in).map{print("\($0.0)\t\($0.1)\t\(sqrt($0.0*$0.0+$0.1*$0.1))")}
	var real = [Double](repeating: 0.0, count: length)
	var imag = [Double](repeating: 0.0, count: length)
	
	vDSP_DFT_ExecuteD(weights, complex.realp, complex.imagp, &real, &imag)
	//print("OUTPUT")
	//_ = zip(real,imag).map{print("\($0.0)\t\($0.1)\t\(sqrt($0.0*$0.0+$0.1*$0.1))")}
	/*var scale = 1.0/Double(length) //0.5 to account for vDSP zrop   impl = 2 * math
	//vDSP_vsmulD(&out_real, 1, &scale, &out_real, 1, vDSP_Length(length))
	//vDSP_vsmulD(&out_imag, 1, &scale, &out_imag, 1, vDSP_Length(length))
	var realinp = real
	var imaginp = imag
	vDSP_vsmulD(&realinp, 1, &scale, &real, 1, vDSP_Length(length))
	vDSP_vsmulD(&imaginp, 1, &scale, &imag, 1, vDSP_Length(length))*/
	
	var split_complex_out = DSPDoubleSplitComplex(realp:&real, imagp: &imag)
	var result = [DSPDoubleComplex](repeating: DSPDoubleComplex(real: 0.0, imag: 0.0), count: length)
	
	
	
	vDSP_ztocD(&split_complex_out, 1, &result, 2, vDSP_Length(length))
	if inputweights == nil {
		vDSP_DFT_DestroySetupD(weights)
	}
	return result
	
}
public func nextDFTn(_ n: Int)-> Int {
	
	// calculate minimum necessary base 2 exponent
	let minbase2 = Int(ceil(log2(Double(n))))
	
	var allowedDFTCases = [Int]()
	
	//if base2 < 4, only 2**base2 allowed (for zrop, 5)
	
	for i in 1..<5 {
		allowedDFTCases.append(1<<i)
	}
	
	//if base2 >= 4, valid cases also include multiples of 1, 3, 5, 15
	if minbase2 >= 5 {
		for i in 5...minbase2 {
			let pow2 = 1<<i
			allowedDFTCases.append(pow2)
			allowedDFTCases.append(3*pow2)
			allowedDFTCases.append(5*pow2)
			allowedDFTCases.append(15*pow2)
		}
	}
	return allowedDFTCases.sorted(by: {$0 < $1}).filter({$0 >= n}).first!
}


public func simplefft2(_ inputarray: [Double]) -> [Double] {
    var inputarray = inputarray
    let length = inputarray.count
    let fft_weights: FFTSetupD = vDSP_create_fftsetupD(vDSP_Length(log2(Float(length))), FFTRadix(kFFTRadix2))!

    var fftMagnitudes = [Double](repeating:0.0, count:inputarray.count)
    var zeroArray = [Double](repeating:0.0, count:inputarray.count)
    
    var splitComplexInput = DSPDoubleSplitComplex(realp: &inputarray, imagp: &zeroArray)
    
    vDSP_fft_zipD(fft_weights, &splitComplexInput, 1, vDSP_Length(log2(CDouble(inputarray.count))), FFTDirection(FFT_FORWARD));
    vDSP_zvmagsD(&splitComplexInput, 1, &fftMagnitudes, 1, vDSP_Length(inputarray.count));
    
    let roots = sqrt(fftMagnitudes) // vDSP_zvmagsD returns squares of the FFT magnitudes, so take the root here
    var normalizedValues = [Double](repeating:0.0, count:inputarray.count)
    
    vDSP_vsmulD(roots, vDSP_Stride(1), [2.0 / Double(inputarray.count)], &normalizedValues, vDSP_Stride(1), vDSP_Length(inputarray.count))
    return normalizedValues
}

public func cfft(_ real: [Double], imag: [Double]) -> [Double] {
    assert(real.count == imag.count, "length of real and imag inputs must be equal")
    var real = [Double](real)
    var imag = [Double](imag)
	let count = real.count
	
   // var imaginary = [Double](repeating: 0.0, count: input.count)
    var splitComplex = DSPDoubleSplitComplex(realp: &real, imagp: &imag)

    let length = vDSP_Length(floor(log2(Float(count))))
    let radix = FFTRadix(kFFTRadix2)
    let weights = vDSP_create_fftsetupD(length, radix)
    withUnsafeMutablePointer(to: &splitComplex) { splitComplex in
        vDSP_fft_zipD(weights!, splitComplex, 1, length, FFTDirection(FFT_FORWARD))
    }

    var magnitudes = [Double](repeating: 0.0, count: real.count)
    withUnsafePointer(to: &splitComplex) { splitComplex in
        magnitudes.withUnsafeMutableBufferPointer { magnitudes in
            vDSP_zvmagsD(splitComplex, 1, magnitudes.baseAddress!, 1, vDSP_Length(count))
        }
    }

    var normalizedMagnitudes = [Double](repeating: 0.0, count: real.count)
    normalizedMagnitudes.withUnsafeMutableBufferPointer { normalizedMagnitudes in
//        vDSP_vsmulD(sqrt(magnitudes), 1, [2.0 / Double(count)], normalizedMagnitudes.baseAddress!, 1, vDSP_Length(count))
        vDSP_vsmulD(sqrt(magnitudes), 1, [1.0], normalizedMagnitudes.baseAddress!, 1, vDSP_Length(count))

    }

    vDSP_destroy_fftsetupD(weights)

    return normalizedMagnitudes
}

private func fftPlanC(_ count: Int, direction: vDSP_DFT_Direction = .FORWARD) -> vDSP_DFT_SetupD {
    let result: vDSP_DFT_SetupD?
    switch direction {
    case .FORWARD:
        result =  vDSP_DFT_zop_CreateSetupD(nil, vDSP_Length(count), .FORWARD)
    case .INVERSE:
        result = vDSP_DFT_zop_CreateSetupD(nil, vDSP_Length(count), .INVERSE)
    }
    assert(result != nil, "FFT Plan Failed")
    return result!
}

private func fftPlanR(_ count: Int, direction: vDSP_DFT_Direction = .FORWARD) -> vDSP_DFT_SetupD {
    let result: vDSP_DFT_SetupD?
    switch direction {
    case .FORWARD:
        result =  vDSP_DFT_zrop_CreateSetupD(nil, vDSP_Length(count), .FORWARD)
    case .INVERSE:
        result = vDSP_DFT_zrop_CreateSetupD(nil, vDSP_Length(count), .INVERSE)
    }
    assert(result != nil, "FFT Plan Failed")
    return result!
}

public func complexfft(_ input: [DSPDoubleComplex], weights inputweights: vDSP_DFT_SetupD? = nil, direction: vDSP_DFT_Direction = .FORWARD) -> [DSPDoubleComplex] {
    
    let inputLength = input.count
    let dftLength = nextDFTn(inputLength)
    
    let weights = inputweights ?? fftPlanC(dftLength, direction: direction)

    let paddingLength = dftLength - inputLength
    let padding = Array<Double>(repeating: 0.0, count: paddingLength)
    let paddingcenteridx = Int(ceil(Double(paddingLength)/2.0))

    // zero phase padding  http://dsp.stackexchange.com/questions/18938/merits-of-zero-phase-zero-padding/
    var paddedReal = Array(ContiguousArray([Array(padding[0..<paddingcenteridx]), input.map{$0.real}, Array(padding[paddingcenteridx..<paddingLength])].joined()))
    var paddedImag = Array(ContiguousArray([Array(padding[0..<paddingcenteridx]), input.map{$0.imag}, Array(padding[paddingcenteridx..<paddingLength])].joined()))
    
    var real = [Double](repeating: 0.0, count: dftLength)
    var imag = [Double](repeating: 0.0, count: dftLength)
    
    
    
    vDSP_DFT_ExecuteD(weights, &paddedReal, &paddedImag, &real, &imag)
    
    var result: [DSPDoubleComplex]
    
    switch direction {
    case .FORWARD:
        result =  zip(real,imag).compactMap{DSPDoubleComplex(real: $0, imag: $1)}
    case .INVERSE:
        result = zip(real,imag).compactMap{DSPDoubleComplex(real: $0/Double(dftLength), imag: $1/Double(dftLength))}
    }
    
    if inputweights == nil {
        vDSP_DFT_DestroySetupD(weights)
    }
    
    return result
}

public func complexfft(real: [Double], imag: [Double], weights inputweights: vDSP_DFT_SetupD? = nil, direction: vDSP_DFT_Direction = .FORWARD) -> (real: [Double], imag:[Double]) {
    let input = zip(real,imag).compactMap({DSPDoubleComplex(real: $0, imag: $1)})
    let fft = complexfft(input, weights: inputweights, direction: direction)
    let real = fft.map{($0.real)}
    let imag = fft.map{($0.imag)}

    return (real, imag)
}
